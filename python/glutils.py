import sys
import os
import random
import re
import glob
from collections import OrderedDict
import csv
from subprocess import check_call, Popen, PIPE
from Bio import SeqIO

import utils

# ----------------------------------------------------------------------------------------
glfo_dir = 'germline-sets'  # always put germline info into a subdir with this name

# single-chain file names
def glfo_csv_fnames():
    return [cdn + '-positions.csv' for cdn in utils.conserved_codons.values()]
def glfo_fasta_fnames(chain):
    return ['ig' + chain + r + algn + '.fasta' for r in utils.regions for algn in ('', '-aligned')]
def glfo_fnames(chain):
    return glfo_csv_fnames() + glfo_fasta_fnames(chain)

# fasta file names including all chains
def all_glfo_fasta_fnames():
    return ['ig' + chain + r + algn + '.fasta' for chain in utils.chains for r in utils.regions for algn in ('', '-aligned')]

#----------------------------------------------------------------------------------------
def read_germline_seqs(datadir, chain):
    seqs, aligned_seqs = {r : OrderedDict() for r in utils.regions}, {r : OrderedDict() for r in utils.regions}
    for fname in glfo_fasta_fnames(chain):
        infodict = aligned_seqs if 'aligned' in fname else seqs
        for seq_record in SeqIO.parse(datadir + '/' + chain + '/' + fname, 'fasta'):
            gene = seq_record.name.split('|')[0]
            seq = str(seq_record.seq).upper()
            if len(seq.strip(''.join(utils.expected_characters))) > 0:  # return the empty string if it only contains expected characters
                raise Exception('unexpected characters in %s (expected %s)' % (seq, ' '.join(utils.expected_characters)))
            infodict[utils.get_region(gene)][gene] = seq
    return seqs, aligned_seqs

#----------------------------------------------------------------------------------------
def clean_up_glfo(glfo, debug=False):
    """ remove any genes in 'aligned-seqs' and codon position info that aren't in 'seqs' """
    desired_genes = set([g for r in utils.regions for g in glfo['seqs'][r]])
    n_aligned_removed = 0
    for gene in [g for r in utils.regions for g in glfo['aligned-seqs'][r]]:
        if gene not in desired_genes:
            n_aligned_removed += 1
            del glfo['aligned-seqs'][utils.get_region(gene)][gene]
    n_codon_removed = {c : 0 for c in utils.conserved_codons.values()}
    for codon in utils.conserved_codons.values():
        for gene in glfo[codon + '-positions'].keys():
            if gene not in desired_genes:
                n_codon_removed[codon] += 1
                del glfo[codon + '-positions'][gene]
    if debug:
        print '    removed %d genes from aligned seqs, %s' % (n_aligned_removed, ' and '.join([('%d from %s positions' % (n_codon_removed[c], c)) for c in utils.conserved_codons.values()]))

# ----------------------------------------------------------------------------------------
def get_new_alignments(glfo, debug=False):
    """
    Pass genes with alignments (from 'aligned-seqs') and those without alignments (anything in 'seqs' that isn't in 'aligned-seqs') into mafft,
    and replace all existing alignments (i.e. 'aligned-seqs') with the resulting msa.
    """
    for region in utils.regions:
        if debug:
            print region
        genes_with_alignments = set(glfo['aligned-seqs'][region])
        genes_without_alignments = set(glfo['seqs'][region]) - set(glfo['aligned-seqs'][region])
        if len(genes_without_alignments) == 0:
            if debug:
                print '  no missing alignments'
            continue

        if debug:
            print '  missing alignments for %d genes: %s' % (len(genes_without_alignments), ' '.join([utils.color_gene(g) for g in genes_without_alignments]))
            print '  existing alignments:'
            for g, seq in glfo['aligned-seqs'][region].items():
                print '    %s   %s' % (seq, utils.color_gene(g))


        # find the longest aligned sequence, so we can pad everybody else with dots on the right out to that length
        biggest_length = None
        for gene in genes_with_alignments:
            if biggest_length is None or len(glfo['aligned-seqs'][region][gene]) > biggest_length:
                biggest_length = len(glfo['aligned-seqs'][region][gene])

        tmpdir = '/tmp'
        already_aligned_fname = tmpdir + '/already-aligned.fasta'
        not_aligned_fname = tmpdir + '/not-aligned.fasta'
        msa_table_fname = tmpdir + '/msa-table.txt'
        aligned_and_not_fnamefname = tmpdir + '/aligned-and-not.fasta'
        mafft_outfname = tmpdir + '/everybody-aligned.fasta'
        with open(already_aligned_fname, 'w') as tmpfile, open(msa_table_fname, 'w') as msafile:
            mysterious_index = 1
            msa_str = ''
            for gene in genes_with_alignments:
                dotstr = '.' * (biggest_length - len(glfo['aligned-seqs'][region][gene]))
                alistr = glfo['aligned-seqs'][region][gene] + dotstr
                tmpfile.write('>%s\n%s\n' % (gene, alistr.replace('.', '-')))
                msa_str += ' ' + str(mysterious_index)
                mysterious_index += 1
            msafile.write('%s # %s\n' % (msa_str, already_aligned_fname))
        with open(not_aligned_fname, 'w') as tmpfile:
            for gene in genes_without_alignments:
                tmpfile.write('>%s\n%s\n' % (gene, glfo['seqs'][region][gene]))

        check_call('cat ' + already_aligned_fname + ' ' + not_aligned_fname + ' >' + aligned_and_not_fnamefname, shell=True)

        # actually run mafft
        cmd = 'mafft --merge ' + msa_table_fname + ' ' + aligned_and_not_fnamefname + ' >' + mafft_outfname  # options=  # "--localpair --maxiterate 1000"
        if debug:
            print '    RUN %s' % cmd
        proc = Popen(cmd, shell=True, stderr=PIPE)
        _, err = proc.communicate()  # debug info goes to err

        if debug:
            # deal with debug info (for err -- out gets redirected to a file)
            err = err.replace('\r', '\n')
            printstrs = []
            for errstr in err.split('\n'):  # remove the stupid progress bar things
                matches = re.findall('[0-9][0-9]* / [0-9][0-9]*', errstr)
                if len(matches) == 1 and errstr.strip() == matches[0]:
                    continue
                if len(errstr) == 0:
                    continue
                printstrs.append(errstr)
            print '        ' + '\n        '.join(printstrs)

        # deal with fasta output
        for seq_record in SeqIO.parse(mafft_outfname, 'fasta'):
            gene = seq_record.name.split('|')[0]
            seq = str(seq_record.seq).upper()
            if gene not in glfo['seqs'][region]:  # only really possible if there's a bug in the preceding fifty lines, but oh well, you can't be too careful
                raise Exception('unexpected gene %s in mafft output' % gene)
            glfo['aligned-seqs'][region][gene] = seq  # overwrite the old alignment with the new one
        if debug:
            print '  new alignments:'
            for g, seq in glfo['aligned-seqs'][region].items():
                print '    %s   %s  %s' % (seq, utils.color_gene(g), '<--- new' if g in genes_without_alignments else '')

        os.remove(already_aligned_fname)
        os.remove(not_aligned_fname)
        os.remove(msa_table_fname)
        os.remove(aligned_and_not_fnamefname)
        os.remove(mafft_outfname)

#----------------------------------------------------------------------------------------
def get_missing_codon_info(glfo, debug=False):
    # ----------------------------------------------------------------------------------------
    def get_n_gaps_up_to_pos(aligned_seq, pos):
        # NOTE I think this duplicates the functionality of count_gaps()
        """ return number of gapped positions in <aligned_seq> before <pos> """
        ipos = 0  # position in unaligned sequence
        n_gaps_passed = 0  # number of gapped positions in the aligned sequence that we pass before getting to <pos> (i.e. while ipos < pos)
        while ipos < pos:
            if aligned_seq[ipos + n_gaps_passed] in utils.gap_chars:
                n_gaps_passed += 1
            else:
                ipos += 1
        return n_gaps_passed

    # ----------------------------------------------------------------------------------------
    def get_pos_in_alignment(codon, aligned_seq, seq, pos):
        """ given <pos> in <seq>, find the codon's position in <aligned_seq> """
        utils.check_codon(codon, seq, pos, debug=True)
        pos_in_alignment = pos + get_n_gaps_up_to_pos(aligned_seq, pos)
        utils.check_codon(codon, aligned_seq, pos_in_alignment, debug=True)
        return pos_in_alignment

    for region, codon in utils.conserved_codons.items():
        missing_genes = set(glfo['seqs'][region]) - set(glfo[codon + '-positions'])
        if len(missing_genes) == 0:
            if debug:
                print '    no missing %s info' % codon
            continue

        if debug:
            print '      missing %d %s positions' % (len(missing_genes), codon)

        if region == 'j':
            raise Exception('missing tryp position for %s, and we can\'t infer it because tryp positions don\'t reliably align to the same position' % ' '.join(missing_genes))

        # existing codon position (this assumes that once aligned, all genes have the same codon position)
        if len(glfo[codon + '-positions']) > 0:
            known_gene, known_pos = glfo[codon + '-positions'].iteritems().next()  # just take the "first" one
            # NOTE for cyst, should be 309 (imgt says 104th codon --> subtract 1 to get zero-indexing, then multiply by three 3 * (104 - 1) = 309
            known_pos_in_alignment = get_pos_in_alignment(codon, glfo['aligned-seqs'][region][known_gene], glfo['seqs'][region][known_gene], known_pos)
            if debug:
                print '      using known position %d (aligned %d) in %s' % (known_pos, known_pos_in_alignment, known_gene)
        elif codon == 'cyst':
            known_pos_in_alignment = 309
            print '      assuming aligned %s position is %d' % (codon, known_pos_in_alignment)
        else:
            raise Exception('no existing %s info, and couldn\'t guess it, either' % codon)

        for gene in missing_genes:
            unaligned_pos = known_pos_in_alignment - utils.count_gaps(glfo['aligned-seqs'][region][gene], istop=known_pos_in_alignment)
            print '      adding %s: %d' % (utils.color_gene(gene), unaligned_pos)
            seq_to_check = glfo['seqs'][region][gene]
            if gene == 'IGHV4-39*03':  # in the imgt database, it ends one base short of finishing the cysteine
                seq_to_check += 'T'
            utils.check_codon(codon, seq_to_check, unaligned_pos, debug=True)
            glfo[codon + '-positions'][gene] = unaligned_pos

#----------------------------------------------------------------------------------------
def add_missing_glfo(glfo, generate_new_alignment=False, debug=False):
    genes_without_alignments = set([g for r in utils.regions for g in set(glfo['seqs'][r]) - set(glfo['aligned-seqs'][r])])
    if len(genes_without_alignments) > 0:
        if generate_new_alignment:  # replace existing alignment, which is missing some genes, with a new one that includes everybody
            get_new_alignments(glfo, debug=debug)
        else:
            raise Exception('missing alignments for %d genes %s (set --generate-new-alignment if you want to attempt to figure out new ones)' % (len(genes_without_alignments), ' '.join(genes_without_alignments)))
    elif debug:
        print '    no missing alignments'

    get_missing_codon_info(glfo, debug=debug)

# ----------------------------------------------------------------------------------------
def read_codon_positions(csvfname):
    positions = {}
    with open(csvfname) as csvfile:
        reader = csv.DictReader(csvfile)
        for line in reader:
            if line['istart'] == '':  # I don't know how these got in the file, but they should probably be removed
                continue
            positions[line['gene']] = int(line['istart'])
    return positions

#----------------------------------------------------------------------------------------
def check_codon_positions(glfo):
    for codon in utils.conserved_codons.values():
        for gene, istart in glfo[codon + '-positions'].items():
            utils.check_codon(codon, glfo['seqs'][utils.get_region(gene)][gene], istart, debug=True)

#----------------------------------------------------------------------------------------
def read_glfo(datadir, chain, only_genes=None, generate_new_alignment=False, debug=False):

    # ----------------------------------------------------------------------------------------
    # warn about backwards compatibility breakage
    if not os.path.exists(datadir + '/' + chain):
        raise Exception("""
        Germline set directory \'%s\' does not exist.
        This is probably because as of v0.6.0 we switched to per-dataset germline sets, which means the germline set fastas have to be *within* --parameter-dir.
        You will, unfortunately, need to delete your existing --parameter-dir and generate a new one (either explicitly with \'cache-parameters\', or just re-run whatever action you were running, and it'll do it automatically when it can\'t find the directory you deleted).
        Sorry! We promise this is better in the long term, though.
        """ % (datadir + '/' + chain))
    # ----------------------------------------------------------------------------------------

    if debug:
        print '  reading %s chain glfo from %s' % (chain, datadir)
    glfo = {'chain' : chain}
    glfo['seqs'], glfo['aligned-seqs'] = read_germline_seqs(datadir, chain)
    for fname in glfo_csv_fnames():
        glfo[utils.get_codon(fname) + '-positions'] = read_codon_positions(datadir + '/' + chain + '/' + fname)
    clean_up_glfo(glfo, debug=debug)  # remove any extra info
    # check_codon_positions(glfo)
    add_missing_glfo(glfo, generate_new_alignment=generate_new_alignment, debug=debug)
    restrict_to_genes(glfo, only_genes, debug=debug)
    if debug:
        print '  read %s' % '  '.join([('%s: %d' % (r, len(glfo['seqs'][r]))) for r in utils.regions])
    return glfo

# ----------------------------------------------------------------------------------------
def stringify_mutfo(mutfo):
    return_str_list = []
    for pos in sorted(mutfo.keys()):
        return_str_list.append(mutfo[pos]['original'] + str(pos) + mutfo[pos]['new'])
    return '.'.join(return_str_list)

# ----------------------------------------------------------------------------------------
def get_mutfo_from_seq(old_seq, new_seq):
    mutfo = {}
    assert len(old_seq) == len(new_seq)
    for inuke in range(len(old_seq)):
        assert old_seq[inuke] in utils.nukes and new_seq[inuke] in utils.nukes
        if new_seq[inuke] != old_seq[inuke]:
            mutfo[inuke] = {'original' : old_seq[inuke], 'new' : new_seq[inuke]}
    return mutfo

# ----------------------------------------------------------------------------------------
def get_mutfo_from_name(gene_name):
    allele_list = utils.allele(gene_name).split('+')
    if len(allele_list) != 2:
        raise Exception('couldn\'t get snp info from gene name %s' % gene_name)
    mutfo = {}
    for mutstr in allele_list[1].split('.'):
        if len(mutstr) < 3:
            raise Exception('couldn\'t extract mutation info from %s' % mutstr)
        original, new = mutstr[0], mutstr[-1]
        if original not in utils.nukes or new not in utils.nukes:
            raise Exception('couldn\'t extract mutation info from %s' % mutstr)
        position = int(mutstr[1:-1])
        if position not in mutfo:
            mutfo[position] = {}
            mutfo[position]['original'] = original  # if it *is* already there, we want to *keep* the old 'original'
        mutfo[position]['new'] = new
        if mutfo[position]['new'] == mutfo[position]['original']:  # reverted back to the original base
            del mutfo[position]

    return mutfo

# ----------------------------------------------------------------------------------------
def get_new_allele_name_and_change_mutfo(template_gene, mutfo):
    if '+' in utils.allele(template_gene):  # template gene was already snp'd
        old_mutfo = get_mutfo_from_name(template_gene)
        for position, info in mutfo.items():
            if position not in old_mutfo:
                old_mutfo[position] = {}
                old_mutfo[position]['original'] = info['original']  # if it *is* already there, we want to *keep* the old 'original'
            old_mutfo[position]['new'] = info['new']
            if old_mutfo[position]['new'] == old_mutfo[position]['original']:  # reverted back to the original base
                del old_mutfo[position]
        final_mutfo = old_mutfo
        assert len(template_gene.split('+')) == 2
        template_gene = template_gene.split('+')[0]  # before we did any snp'ing
    else:
        final_mutfo = mutfo

    final_name = template_gene
    if len(final_mutfo) > 0:
        assert '+' not in final_name
        final_name += '+' + stringify_mutfo(final_mutfo)  # full, but possibly overly verbose
    return final_name, final_mutfo

# ----------------------------------------------------------------------------------------
def generate_snpd_gene(gene, cpos, seq, aligned_seq, positions):
    assert utils.get_region(gene) == 'v'  # others not yet handled
    def choose_position():
        snp_pos = None
        while snp_pos is None or snp_pos in snpd_positions or not utils.check_conserved_cysteine(tmpseq, cpos, debug=True, assert_on_fail=False):
            snp_pos = random.randint(10, len(seq) - 15)  # note that randint() is inclusive
            tmpseq = seq[: snp_pos] + 'X' + seq[snp_pos + 1 :]  # for checking cyst position
        return snp_pos

    snpd_positions = set()  # only used if a position wasn't specified (i.e. was None) in <snps_to_add>
    mutfo = OrderedDict()
    for snp_pos in positions:
        if snp_pos is None:
            snp_pos = choose_position()
        snpd_positions.add(snp_pos)
        new_base = None
        while new_base is None or new_base == seq[snp_pos]:
            new_base = utils.nukes[random.randint(0, len(utils.nukes) - 1)]
        print '        %3d   %s --> %s' % (snp_pos, seq[snp_pos], new_base)
        mutfo[snp_pos] = {'original' : seq[snp_pos], 'new' : new_base}

        seq = seq[: snp_pos] + new_base + seq[snp_pos + 1 :]

        igl = 0  # position in unaligned germline seq
        for ialign in range(len(aligned_seq)):
            if aligned_seq[ialign] in utils.gap_chars:
                continue
            else:
                # assert aligned_seq[ialign] == seq[igl]  # won't necessarily be true after the first mutation
                if igl == snp_pos:
                    aligned_seq = aligned_seq[ : ialign] + new_base + aligned_seq[ialign + 1 :]
                    break
                igl += 1

    utils.check_conserved_cysteine(seq, cpos)
    snpd_name, mutfo = get_new_allele_name_and_change_mutfo(gene, mutfo)
    return {'template-gene' : gene, 'gene' : snpd_name, 'seq' : seq, 'aligned-seq' : aligned_seq}

# ----------------------------------------------------------------------------------------
def restrict_to_genes(glfo, only_genes, debug=False):
    """ remove from <glfo> any genes which are not in <only_genes> """
    if only_genes is None:
        return
    genes_to_remove = set([g for r in utils.regions for g in glfo['seqs'][r]]) - set(only_genes)
    if debug:
        print '    removing %d genes that aren\'t in <only_genes> (%s)' % (len(genes_to_remove), ' '.join(genes_to_remove))
    remove_genes(glfo, genes_to_remove)

# ----------------------------------------------------------------------------------------
def remove_genes(glfo, genes, debug=False):
    """ remove <genes> from <glfo> """
    if debug:
        print '  removing %s from glfo' % ' '.join([utils.color_gene(g) for g in genes])
    for gene in genes:
        remove_gene(glfo, gene)

# ----------------------------------------------------------------------------------------
def remove_gene(glfo, gene, debug=False):
    """ remove <gene> from <glfo> """
    if debug:
        print '  removing %s from glfo' % utils.color_gene(gene)
    region = utils.get_region(gene)
    if region in utils.conserved_codons:
        del glfo[utils.conserved_codons[region] + '-positions'][gene]
    del glfo['seqs'][region][gene]
    del glfo['aligned-seqs'][region][gene]

# ----------------------------------------------------------------------------------------
def add_new_alleles(glfo, newfos, remove_template_genes, debug=False):
    for newfo in newfos:
        add_new_allele(glfo, newfo, remove_template_genes, debug=debug)

# ----------------------------------------------------------------------------------------
def add_new_allele(glfo, newfo, remove_template_genes, debug=False):
    """
    Add a new allele to <glfo>, specified by <newfo> which is of the
    form: {'template-gene' : 'IGHV3-71*01', 'gene' : 'IGHV3-71*01+C35T.T47G', 'seq' : 'ACTG yadda yadda CGGGT', 'aligned-seq' : None}
    If <remove_template_genes>, we also remove 'template-gene' from <glfo>.
    """

    template_gene = newfo['template-gene']
    region = utils.get_region(template_gene)
    if template_gene not in glfo['seqs'][region]:
        raise Exception('unknown template gene %s' % template_gene)

    new_gene = newfo['gene']

    if region == 'v':
        glfo['cyst-positions'][new_gene] = glfo['cyst-positions'][template_gene]
    elif region == 'j':
        glfo['tryp-positions'][new_gene] = glfo['tryp-positions'][template_gene]

    glfo['seqs'][region][new_gene] = newfo['seq']

    if newfo['aligned-seq'] is None:
        # get_new_alignments(glfo, debug=debug)  # don't want to do this -- we'd kind of rather keep the imgt-gapped alignment if we can
        # just copy the template gene's alignment
        template_alignment = glfo['aligned-seqs'][region][template_gene]
        new_alignment = ''
        assert len(glfo['seqs'][region][new_gene]) == len(glfo['seqs'][region][template_gene])  # only works if they're only separated by snps
        n_gaps = 0
        for ipos in range(len(template_alignment)):
            if template_alignment[ipos] in utils.gap_chars:
                new_alignment += template_alignment[ipos]
                n_gaps += 1
            else:
                new_alignment += glfo['seqs'][region][new_gene][ipos - n_gaps]
        glfo['aligned-seqs'][region][new_gene] = new_alignment
    else:
        glfo['aligned-seqs'][region][new_gene] = newfo['aligned-seq']

    if debug:
        print '    adding new allele to glfo:'
        print '      template %s   %s' % (glfo['seqs'][region][template_gene], utils.color_gene(template_gene))
        print '           new %s   %s' % (utils.color_mutants(glfo['seqs'][region][template_gene], newfo['seq']), utils.color_gene(new_gene))

    if remove_template_genes:
        remove_gene(glfo, template_gene, debug=True)

# ----------------------------------------------------------------------------------------
def remove_the_stupid_godamn_template_genes_all_at_once(glfo, templates_to_remove):
    for gene in templates_to_remove:
        remove_gene(glfo, gene, debug=True)

# ----------------------------------------------------------------------------------------
def add_some_snps(snps_to_add, glfo, remove_template_genes=False, debug=False):
    """
    Generate some snp'd genes and add them to glfo, specified with <snps_to_add>.
    e.g. [{'gene' : 'IGHV3-71*01', 'positions' : (35, None)}, ] will add a snp at position 35 and at a random location.
    The resulting snp'd gene will have a name like IGHV3-71*01+C35T.T47G
    """

    templates_to_remove = set()

    for isnp in range(len(snps_to_add)):
        snpinfo = snps_to_add[isnp]
        gene, positions = snpinfo['gene'], snpinfo['positions']
        print '    adding %d %s to %s' % (len(positions), utils.plural_str('snp', len(positions)), gene)
        seq = glfo['seqs'][utils.get_region(gene)][gene]
        aligned_seq = glfo['aligned-seqs'][utils.get_region(gene)][gene]
        assert utils.get_region(gene) == 'v'
        cpos = glfo['cyst-positions'][gene]

        snpfo = None
        itry = 0
        while snpfo is None or snpfo['gene'] in glfo['seqs'][utils.get_region(gene)]:
            if itry > 0:
                print '      already in glfo, try again'
                if itry > 99:
                    raise Exception('too many tries while trying to generate new snps -- did you specify a lot of snps on the same position?')
            snpfo = generate_snpd_gene(gene, cpos, seq, aligned_seq, positions)
            itry += 1

        if remove_template_genes:
            templates_to_remove.add(gene)
        add_new_allele(glfo, snpfo, remove_template_genes=False, debug=debug)  # *don't* remove the templates here, since we don't know if there's another snp later that needs them

    remove_the_stupid_godamn_template_genes_all_at_once(glfo, templates_to_remove)  # works fine with zero-length <templates_to_remove>

# ----------------------------------------------------------------------------------------
def write_glfo(output_dir, glfo, only_genes=None, debug=False):
    if debug:
        print '  writing glfo to %s%s' % (output_dir, '' if only_genes is None else ('  (restricting to %d genes)' % len(only_genes)))
    if os.path.exists(output_dir + '/' + glfo['chain']):
        remove_glfo_files(output_dir, glfo['chain'])  # also removes output_dir
    os.makedirs(output_dir + '/' + glfo['chain'])

    for fname in glfo_fasta_fnames(glfo['chain']):
        glseqfo = glfo['aligned-seqs'] if 'aligned' in fname else glfo['seqs']
        with open(output_dir + '/' + glfo['chain'] + '/' + fname, 'w') as outfile:
            for gene in glseqfo[utils.get_region(fname)]:
                if only_genes is not None and gene not in only_genes:
                    continue
                outfile.write('>' + gene + '\n')
                outfile.write(glseqfo[utils.get_region(fname)][gene] + '\n')
    for fname in glfo_csv_fnames():
        with open(output_dir + '/' + glfo['chain'] + '/' + fname, 'w') as codonfile:
            writer = csv.DictWriter(codonfile, ('gene', 'istart'))
            writer.writeheader()
            for gene, istart in glfo[utils.get_codon(fname) + '-positions'].items():
                if only_genes is not None and gene not in only_genes:
                    continue
                writer.writerow({'gene' : gene, 'istart' : istart})

    # make sure there weren't any files lingering in the output dir when we started
    # NOTE this will ignore the dirs corresponding to any *other* chains (which is what we want now, I think)
    unexpected_files = set(glob.glob(output_dir + '/' + glfo['chain'] + '/*')) - set([output_dir + '/' + glfo['chain'] + '/' + fn for fn in glfo_fnames(glfo['chain'])])
    if len(unexpected_files) > 0:
        raise Exception('unexpected file(s) while writing germline set: %s' % (' '.join(unexpected_files)))

# ----------------------------------------------------------------------------------------
def remove_glfo_files(datadir, chain):
    for fname in glfo_fnames(chain):
        os.remove(datadir + '/' + chain + '/' + fname)
    os.rmdir(datadir + '/' + chain)
    os.rmdir(datadir)  # at the moment, we should only be running on single-chain stuff, so the only dir with info for more than one chain should be data/imgt

