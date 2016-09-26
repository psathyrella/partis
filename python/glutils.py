import tempfile
import copy
import sys
import os
import random
import re
import numpy
import glob
from collections import OrderedDict
import csv
from subprocess import check_call, Popen, PIPE
from Bio import SeqIO

import utils

# ----------------------------------------------------------------------------------------
glfo_dir = 'germline-sets'  # always put germline info into a subdir with this name

dummy_d_genes = {'h' : None, 'k' : 'IGKDx-x*x', 'l' : 'IGLDx-x*x'}  # for light chain

# single-chain file names
extra_fname = 'extras.csv'
def glfo_fasta_fnames(chain):
    if chain == 'h':
        return ['ig' + chain + r + '.fasta' for r in utils.regions]
    else:
        return ['ig' + chain + r + '.fasta' for r in utils.regions if r != 'd']
def glfo_fnames(chain):
    return [extra_fname, ] + glfo_fasta_fnames(chain)

# fasta file names including all chains
def all_glfo_fasta_fnames():
    return ['ig' + chain + r + '.fasta' for chain in utils.chains for r in utils.regions]

csv_headers = ['gene', 'cyst_position', 'tryp_position', 'phen_position', 'aligned_seq']

imgt_info_indices = ('accession-number', 'gene', 'species', 'functionality', '', '', '', '', '', '', '', '', '')  # I think this is the right number of entries, but it doesn't really matter
functionalities = ['F', 'ORF', 'P', '(F)', '[F]', '[P]', '[ORF]']   # not actually sure what the last few mean
pseudogene_funcionalities = ['P', '[P]']

duplicate_names = {
    'v' : [
        set(['IGHV3-23*04', 'IGHV3-23D*02']),
        set(['IGHV3-30*18', 'IGHV3-30-5*01']),
        set(['IGHV1-69*01', 'IGHV1-69D*01']),
        set(['IGHV3-30*02', 'IGHV3-30-5*02']),
        set(['IGHV3-30*04', 'IGHV3-30-3*03']),
        set(['IGHV3-23*01', 'IGHV3-23D*01'])],
    'd' : [
        set(['IGHD1/OR15-1a*01', 'IGHD1/OR15-1b*01']),
        set(['IGHD2/OR15-2a*01', 'IGHD2/OR15-2b*01']),
        set(['IGHD4/OR15-4a*01', 'IGHD4/OR15-4b*01']),
        set(['IGHD3/OR15-3a*01', 'IGHD3/OR15-3b*01']),
        set(['IGHD5/OR15-5a*01', 'IGHD5/OR15-5b*01']),
        set(['IGHD5-18*01', 'IGHD5-5*01']),
        set(['IGHD4-11*01', 'IGHD4-4*01'])],
    'j' : []
}

#----------------------------------------------------------------------------------------
def convert_to_duplicate_name(glfo, gene):
    for equivalence_class in duplicate_names[utils.get_region(gene)]:
        if gene in equivalence_class:
            for alternate_name in equivalence_class:
                if alternate_name != gene and alternate_name in glfo['seqs'][utils.get_region(gene)]:
                    # print 'converting %s --> %s' % (gene, alternate_name)
                    return alternate_name
    raise Exception('couldn\'t find alternate name for %s' % gene)

#----------------------------------------------------------------------------------------
def read_fasta_file(seqs, fname, skip_pseudogenes, aligned=False):
    n_skipped_pseudogenes = 0
    seq_to_gene_map = {}
    for seq_record in SeqIO.parse(fname, 'fasta'):
        linefo = [p.strip() for p in seq_record.description.split('|')]

        # first get gene name
        if linefo[0][:2] != 'IG':  # if it's an imgt file, with a bunch of header info (and the accession number first)
            gene = linefo[imgt_info_indices.index('gene')]
            functionality = linefo[imgt_info_indices.index('functionality')]
            if functionality not in functionalities:
                raise Exception('unexpected functionality %s in %s' % (functionality, fname))
            if skip_pseudogenes and functionality in pseudogene_funcionalities:
                n_skipped_pseudogenes += 1
                continue
        else:  # plain fasta with just the gene name after the '>'
            gene = linefo[0]
        utils.split_gene(gene)  # just to check if it's a valid gene name
        if not aligned and utils.get_region(gene) != utils.get_region(os.path.basename(fname)):  # if <aligned> is True, file name is expected to be whatever
            raise Exception('gene %s from %s has unexpected region %s' % (gene, os.path.basename(fname), utils.get_region(gene)))
        if gene in seqs[utils.get_region(gene)]:
            raise Exception('gene name %s appears twice in %s' % (gene, fname))

        # then the sequence
        seq = str(seq_record.seq).upper()
        if not aligned:
            seq = utils.remove_gaps(seq)
        if 'Y' in seq:
            print '      replacing Y --> N (%d of \'em) in %s' % (seq.count('Y'), utils.color_gene(gene))
            seq = seq.replace('Y', 'N')
        if len(seq.strip(''.join(utils.expected_characters))) > 0:  # return the empty string if it only contains expected characters
            raise Exception('unexpected character %s in %s (expected %s)' % (seq.strip(''.join(utils.expected_characters)), seq, ' '.join(utils.expected_characters)))
        if seq not in seq_to_gene_map:
            seq_to_gene_map[seq] = []
        seq_to_gene_map[seq].append(gene)

        seqs[utils.get_region(gene)][gene] = seq

    tmpcounts = [len(gl) for gl in seq_to_gene_map.values()]  # number of names corresponding to each sequence (should all be ones)
    if tmpcounts.count(1) != len(tmpcounts):
        print '  mutliple names in %s for the following sequences:' % fname
        for seq, genelist in seq_to_gene_map.items():
            if len(genelist) > 1:
                print '    %-50s   %s' % (' '.join(genelist), seq)
        raise Exception('please de-duplicate the fasta and re-run.')

    if n_skipped_pseudogenes > 0:
        print '    skipped %d %s pseudogenes (leaving %d)' % (n_skipped_pseudogenes, utils.get_region(os.path.basename(fname)), len(seqs[utils.get_region(os.path.basename(fname))]))

#----------------------------------------------------------------------------------------
def read_germline_seqs(gldir, chain, skip_pseudogenes):
    seqs = {r : OrderedDict() for r in utils.regions}
    for fname in glfo_fasta_fnames(chain):
        read_fasta_file(seqs, gldir + '/' + chain + '/' + fname, skip_pseudogenes)
    if chain != 'h':
        seqs['d'][dummy_d_genes[chain]] = 'A'  # NOTE this choice ('A') is also set in packages/ham/src/bcrutils.cc
    return seqs

# ----------------------------------------------------------------------------------------
def read_aligned_gl_seqs(fname, glfo):
    aligned_gl_seqs = {r : OrderedDict() for r in utils.regions}
    read_fasta_file(aligned_gl_seqs, fname, skip_pseudogenes=False, aligned=True)

    # pad all the Vs to the same length (imgt fastas just leave them all unequal lengths on the right of the alignment)
    max_aligned_length = max([len(seq) for seq in aligned_gl_seqs['v'].values()])
    for gene in aligned_gl_seqs['v']:
        n_extra_gaps = max_aligned_length - len(aligned_gl_seqs['v'][gene])
        aligned_gl_seqs['v'][gene] += n_extra_gaps * utils.gap_chars[0]

    # check that we got all the genes
    glfo_genes = set([g for r in utils.regions for g in glfo['seqs'][r]]) - set([dummy_d_genes[glfo['chain']], ])
    aligned_genes = set([g for r in utils.regions for g in aligned_gl_seqs[r]])
    if len(glfo_genes - aligned_genes) > 0:
        raise Exception('missing alignments for %s' % ' '.join(glfo_genes - aligned_genes))

    return aligned_gl_seqs

# ----------------------------------------------------------------------------------------
def get_new_alignments(glfo, region, debug=False):
    aligned_seqs = {}

    genes_with_alignments = set(aligned_seqs)  # used to already have some sequences aligned, and may as well keep around the code to handle that case
    genes_without_alignments = set(glfo['seqs'][region]) - set(aligned_seqs)
    if len(genes_without_alignments) == 0:
        if debug:
            print '  no missing %s alignments' % region
        return

    if debug:
        print '        missing alignments for %d %s genes' % (len(genes_without_alignments), region)
        if len(aligned_seqs) > 0:
            print '      existing alignments:'
            for g, seq in aligned_seqs.items():
                print '    %s   %s' % (seq, utils.color_gene(g))

    # find the longest aligned sequence, so we can pad everybody else with dots on the right out to that length
    biggest_length = None
    for gene in genes_with_alignments:
        if biggest_length is None or len(aligned_seqs[gene]) > biggest_length:
            biggest_length = len(aligned_seqs[gene])

    tmpdir = tempfile.mkdtemp()
    already_aligned_fname = tmpdir + '/already-aligned.fasta'
    not_aligned_fname = tmpdir + '/not-aligned.fasta'
    msa_table_fname = tmpdir + '/msa-table.txt'
    aligned_and_not_fnamefname = tmpdir + '/aligned-and-not.fasta'
    mafft_outfname = tmpdir + '/everybody-aligned.fasta'
    with open(already_aligned_fname, 'w') as tmpfile, open(msa_table_fname, 'w') as msafile:
        mysterious_index = 1
        msa_str = ''
        for gene in genes_with_alignments:
            dotstr = '.' * (biggest_length - len(aligned_seqs[gene]))
            alistr = aligned_seqs[gene] + dotstr
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
        print '          RUN %s' % cmd
    proc = Popen(cmd, shell=True, stderr=PIPE)
    _, err = proc.communicate()  # debug info goes to err

    if debug and False:  # aw, screw it, I don't even know what any of mafft's output means
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
        aligned_seqs[gene] = seq  # overwrite the old alignment with the new one
    if debug and False:  # too damn verbose with all the v genes
        print '  new alignments:'
        for g, seq in aligned_seqs.items():
            print '    %s   %s  %s' % (seq, utils.color_gene(g), '<--- new' if g in genes_without_alignments else '')

    os.remove(already_aligned_fname)
    os.remove(not_aligned_fname)
    os.remove(msa_table_fname)
    os.remove(aligned_and_not_fnamefname)
    os.remove(mafft_outfname)
    os.rmdir(tmpdir)

    return aligned_seqs

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
        assert utils.codon_ok(codon, seq, pos, debug=debug)  # this only gets called on the gene with the *known* position, so it shouldn't fail
        pos_in_alignment = pos + get_n_gaps_up_to_pos(aligned_seq, pos)
        assert utils.codon_ok(codon, aligned_seq, pos_in_alignment, debug=debug)
        return pos_in_alignment

    for region, codon in utils.conserved_codons[glfo['chain']].items():
        missing_genes = set(glfo['seqs'][region]) - set(glfo[codon + '-positions'])
        if len(missing_genes) == 0:
            if debug:
                print '      no missing %s info' % codon
            continue

        if debug:
            print '      missing %d %s positions' % (len(missing_genes), codon)

        aligned_seqs = get_new_alignments(glfo, region, debug=debug)
        # for g, s in aligned_seqs.items():
        #     print s, utils.color_gene(g)

        # if region == 'j':
        #     raise Exception('missing tryp position for %s, and we can\'t infer it because tryp positions don\'t reliably align to the same position' % ' '.join(missing_genes))

        # existing codon position (this assumes that once aligned, all genes have the same codon position -- which is only really true for the imgt-gapped alignment)
        if len(glfo[codon + '-positions']) > 0:
            known_gene, known_pos = None, None
            for gene, pos in glfo[codon + '-positions'].items():  # take the first one for which we have the sequence (NOTE it would be safer to check that they're all the same)
                if gene in glfo['seqs'][region] and gene in aligned_seqs and utils.codon_ok(codon, glfo['seqs'][region][gene], pos):
                    known_gene, known_pos = gene, pos
                    break
            if known_gene is None:
                raise Exception('couldn\'t find a known %s position' % codon)
            # NOTE for cyst, should be 309 if alignments are imgt [which they used to usually be, but now probably aren't] (imgt says 104th codon --> subtract 1 to get zero-indexing, then multiply by three 3 * (104 - 1) = 309
            known_pos_in_alignment = get_pos_in_alignment(codon, aligned_seqs[known_gene], glfo['seqs'][region][known_gene], known_pos)
            if debug:
                print '        using known position %d (aligned %d) in %s' % (known_pos, known_pos_in_alignment, known_gene)
        elif codon == 'cyst':
            known_pos_in_alignment = 309
            print '      assuming aligned %s position is %d (this will %s work if you\'re using imgt alignments)' % (codon, known_pos_in_alignment, utils.color('red', 'only'))
            raise Exception('not really using imgt alignments much any more, so this isn\'t really going to work')
        else:
            raise Exception('no existing %s info, and couldn\'t guess it, either' % codon)

        n_added = 0
        seqons = []  # (seq, pos) pairs
        for gene in missing_genes:
            unaligned_pos = known_pos_in_alignment - utils.count_gaps(aligned_seqs[gene], istop=known_pos_in_alignment)
            seq_to_check = glfo['seqs'][region][gene]
            seqons.append((seq_to_check, unaligned_pos))
            glfo[codon + '-positions'][gene] = unaligned_pos
            n_added += 1
            # if debug:
            #     tmpseq = aligned_seqs[gene]
            #     tmppos = known_pos_in_alignment
            #     print '    %s%s%s   %s    (new)' % (tmpseq[:tmppos], utils.color('reverse_video', tmpseq[tmppos : tmppos + 3]), tmpseq[tmppos + 3:], utils.color_gene(gene))

        utils.check_a_bunch_of_codons(codon, seqons, extra_str='          ', debug=debug)
        if debug:
            print '      added %d %s positions' % (n_added, codon)

#----------------------------------------------------------------------------------------
def remove_extraneouse_info(glfo, debug=False):
    """ remove codon info corresponding to genes that aren't in 'seqs' """
    for region, codon in utils.conserved_codons[glfo['chain']].items():
        genes_to_remove = set(glfo[codon + '-positions']) - set(glfo['seqs'][region])
        if debug:
            print '    removing %s info for %d genes (leaving %d)' % (codon, len(genes_to_remove), len(glfo[codon + '-positions']) - len(genes_to_remove))
        for gene in genes_to_remove:
                del glfo[codon + '-positions'][gene]

# ----------------------------------------------------------------------------------------
def read_extra_info(glfo, gldir):
    for codon in utils.conserved_codons[glfo['chain']].values():
        glfo[codon + '-positions'] = {}
    with open(gldir + '/' + glfo['chain'] + '/' + extra_fname) as csvfile:
        reader = csv.DictReader(csvfile)
        for line in reader:
            for codon in utils.conserved_codons[glfo['chain']].values():
                if line[codon + '_position'] != '':
                    glfo[codon + '-positions'][line['gene']] = int(line[codon + '_position'])

#----------------------------------------------------------------------------------------
def read_glfo(gldir, chain, only_genes=None, skip_pseudogenes=True, debug=False):
    if not os.path.exists(gldir + '/' + chain):
        raise Exception('germline set directory \'%s\' does not exist (maybe --parameter-dir is corrupted, maybe crashed while writing parameters?)' % (gldir + '/' + chain))

    if debug:
        print '  reading %s chain glfo from %s' % (chain, gldir)
    glfo = {'chain' : chain}
    glfo['seqs'] = read_germline_seqs(gldir, chain, skip_pseudogenes)
    read_extra_info(glfo, gldir)
    get_missing_codon_info(glfo, debug=debug)
    restrict_to_genes(glfo, only_genes, debug=debug)

    for region, codon in utils.conserved_codons[glfo['chain']].items():
        seqons = [(seq, glfo[codon + '-positions'][gene]) for gene, seq in glfo['seqs'][region].items()]  # (seq, pos) pairs
        utils.check_a_bunch_of_codons(codon, seqons, extra_str='      ', debug=debug)

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
def generate_snpd_gene(gene, cpos, seq, positions):
    assert utils.get_region(gene) == 'v'  # others not yet handled
    def choose_position():
        snp_pos = None
        while snp_pos is None or snp_pos in snpd_positions or not utils.codon_ok('cyst', tmpseq, cpos, debug=True):
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

    assert utils.codon_ok('cyst', seq, cpos, debug=True)  # this is probably unnecessary
    snpd_name, mutfo = get_new_allele_name_and_change_mutfo(gene, mutfo)
    return {'template-gene' : gene, 'gene' : snpd_name, 'seq' : seq}

# ----------------------------------------------------------------------------------------
def restrict_to_genes(glfo, only_genes, debug=False):
    """
    Remove from <glfo> any genes which are not in <only_genes>.
    Any regions in <unrestricted_regions> will be left entirely unrestricted (i.e. there shouldn't be anybody from those regions in <only_genes>.
    """
    if only_genes is None:
        return

    unrestricted_regions = copy.deepcopy(utils.regions)
    for g in only_genes:
        if utils.get_region(g) in unrestricted_regions:
            unrestricted_regions.remove(utils.get_region(g))

    restricted_regions = copy.deepcopy(utils.regions)
    if len(unrestricted_regions) > 0:
        restricted_regions = [r for r in utils.regions if r not in unrestricted_regions]

    only_genes_not_in_glfo = set(only_genes) - set([g for r in restricted_regions for g in glfo['seqs'][r]])
    if len(only_genes_not_in_glfo) > 0:
        print '  %s genes %s in <only_genes> aren\'t in glfo to begin with' % (utils.color('red', 'warning'), ' '.join(only_genes_not_in_glfo))

    genes_to_remove = set([g for r in restricted_regions for g in glfo['seqs'][r]]) - set(only_genes)
    if debug:
        print '    removing %d genes from glfo' % len(genes_to_remove)
    remove_genes(glfo, genes_to_remove)

# ----------------------------------------------------------------------------------------
def remove_v_genes_with_bad_cysteines(glfo, debug=False):
    prelength = len(glfo['seqs']['v'])
    for gene in glfo['seqs']['v'].keys():
        # if len(glfo['seqs']['v'][gene]) < glfo['cyst-positions'][gene] + 3:
        if not utils.codon_ok('cyst', glfo['seqs']['v'][gene], glfo['cyst-positions'][gene]):
            remove_gene(glfo, gene, debug=debug)
    if True:  # debug:
        print '  removed %d / %d v genes with bad cysteines' % (prelength - len(glfo['seqs']['v']), len(glfo['seqs']['v']))

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
    region = utils.get_region(gene)
    if gene in glfo['seqs'][region]:
        if debug:
            print '  removing %s from glfo' % utils.color_gene(gene)
        del glfo['seqs'][region][gene]
        if region in utils.conserved_codons[glfo['chain']]:
            del glfo[utils.conserved_codons[glfo['chain']][region] + '-positions'][gene]
    else:
        if debug:
            print '  can\'t remove %s from glfo, it\'s not there' % utils.color_gene(gene)

# ----------------------------------------------------------------------------------------
def add_new_alleles(glfo, newfos, remove_template_genes=False, debug=False):
    for newfo in newfos:
        add_new_allele(glfo, newfo, remove_template_genes=remove_template_genes, debug=debug)

# ----------------------------------------------------------------------------------------
def add_new_allele(glfo, newfo, remove_template_genes=False, debug=False):
    """
    Add a new allele to <glfo>, specified by <newfo> which is of the
    form: {'template-gene' : 'IGHV3-71*01', 'gene' : 'IGHV3-71*01+C35T.T47G', 'seq' : 'ACTG yadda yadda CGGGT'}
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

    added_snp_names = []
    for isnp in range(len(snps_to_add)):
        snpinfo = snps_to_add[isnp]
        gene, positions = snpinfo['gene'], snpinfo['positions']
        print '    adding %d %s to %s' % (len(positions), utils.plural_str('snp', len(positions)), gene)
        seq = glfo['seqs'][utils.get_region(gene)][gene]
        assert utils.get_region(gene) == 'v'
        cpos = glfo['cyst-positions'][gene]

        snpfo = None
        itry = 0
        while snpfo is None or snpfo['gene'] in glfo['seqs'][utils.get_region(gene)]:
            if itry > 0:
                print '      already in glfo, try again'
                if itry > 99:
                    raise Exception('too many tries while trying to generate new snps -- did you specify a lot of snps on the same position?')
            snpfo = generate_snpd_gene(gene, cpos, seq, positions)
            itry += 1

        if remove_template_genes:
            templates_to_remove.add(gene)
        add_new_allele(glfo, snpfo, remove_template_genes=False, debug=debug)  # *don't* remove the templates here, since we don't know if there's another snp later that needs them
        added_snp_names.append(snpfo['gene'])

    remove_the_stupid_godamn_template_genes_all_at_once(glfo, templates_to_remove)  # works fine with zero-length <templates_to_remove>

    return added_snp_names

# ----------------------------------------------------------------------------------------
def write_glfo(output_dir, glfo, only_genes=None, debug=False):
    if debug:
        print '  writing glfo to %s%s' % (output_dir, '' if only_genes is None else ('  (restricting to %d genes)' % len(only_genes)))
    if os.path.exists(output_dir + '/' + glfo['chain']):
        remove_glfo_files(output_dir, glfo['chain'])  # also removes output_dir
    os.makedirs(output_dir + '/' + glfo['chain'])

    for fname in glfo_fasta_fnames(glfo['chain']):
        with open(output_dir + '/' + glfo['chain'] + '/' + fname, 'w') as outfile:
            for gene in glfo['seqs'][utils.get_region(fname)]:
                if only_genes is not None and gene not in only_genes:
                    continue
                outfile.write('>' + gene + '\n')
                outfile.write(glfo['seqs'][utils.get_region(fname)][gene] + '\n')

    with open(output_dir + '/' + glfo['chain'] + '/' + extra_fname, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, csv_headers)
        writer.writeheader()
        for region, codon in utils.conserved_codons[glfo['chain']].items():
            for gene, istart in glfo[codon + '-positions'].items():
                if only_genes is not None and gene not in only_genes:
                    continue
                writer.writerow({'gene' : gene, codon + '_position' : istart})

    # make sure there weren't any files lingering in the output dir when we started
    # NOTE this will ignore the dirs corresponding to any *other* chains (which is what we want now, I think)
    unexpected_files = set(glob.glob(output_dir + '/' + glfo['chain'] + '/*')) - set([output_dir + '/' + glfo['chain'] + '/' + fn for fn in glfo_fnames(glfo['chain'])])
    if len(unexpected_files) > 0:
        raise Exception('unexpected file(s) while writing germline set: %s' % (' '.join(unexpected_files)))

# ----------------------------------------------------------------------------------------
def remove_glfo_files(gldir, chain):
    for fname in glfo_fnames(chain):
        os.remove(gldir + '/' + chain + '/' + fname)
    os.rmdir(gldir + '/' + chain)
    os.rmdir(gldir)  # at the moment, we should only be running on single-chain stuff, so the only dir with info for more than one chain should be data/germlines

# ----------------------------------------------------------------------------------------
def choose_some_alleles(region, genes_to_use, allelic_groups, n_alleles_per_gene, debug=False):
    """ choose a gene (i.e. a primary and sub-version) from <allelic_groups>, and its attendant alleles """
    # NOTE also modifies <allelic_groups>

    if len(allelic_groups[region]) == 0:
        raise Exception('ran out of %s alleles (either --n-genes-per-region or --n-alleles-per-gene are probably too big)' % region)

    available_versions = None
    while available_versions is None or len(available_versions) == 0:
        if available_versions is not None:
            print '  %s couldn\'t find any versions that have %d alleles, so trying again' % (utils.color('red', 'warning'), n_alleles)
        n_alleles = numpy.random.choice(n_alleles_per_gene[region])
        available_versions = [(pv, subv) for pv in allelic_groups[region] for subv in allelic_groups[region][pv] if len(allelic_groups[region][pv][subv]) >= n_alleles]
    ichoice = numpy.random.randint(0, len(available_versions) - 1) if len(available_versions) > 1 else 0  # numpy.random.choice() can't handle list of tuples (and barfs if you give it only one thing to choose from)
    primary_version, sub_version = available_versions[ichoice]
    new_alleles = set(numpy.random.choice(list(allelic_groups[region][primary_version][sub_version]), size=n_alleles, replace=False))
    if debug:
        print '      %8s %5s   %s' % (primary_version, sub_version, ' '.join([utils.color_gene(g, width=15) for g in new_alleles]))

    assert len(new_alleles & genes_to_use) == 0  # make sure none of the new alleles are already in <genes_to_use>
    genes_to_use |= new_alleles  # actually add them to the final set

    # remove stuff we've used from <allelic_groups>
    del allelic_groups[region][primary_version][sub_version]  # remove this sub-version (we don't want any more alleles from it)
    if len(allelic_groups[region][primary_version]) == 0:
        del allelic_groups[region][primary_version]

# ----------------------------------------------------------------------------------------
def write_allele_prevalence_freqs(allele_prevalence_freqs, fname):
    # NOTE kinda weird to mash all the regions into one file here (as compared to parametercounter), but it seems to make more sense
    if not os.path.exists(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    with open(fname, 'w') as pfile:
        writer = csv.DictWriter(pfile, ('gene', 'freq'))
        writer.writeheader()
        for region in utils.regions:
            for gene, freq in allele_prevalence_freqs[region].items():
                writer.writerow({'gene' : gene, 'freq' : freq})

# ----------------------------------------------------------------------------------------
def read_allele_prevalence_freqs(fname, debug=False):
    # NOTE kinda weird to mash all the regions into one file here (as compared to parametercounter), but it seems to make more sense
    allele_prevalence_freqs = {r : {} for r in utils.regions}
    with open(fname) as pfile:
        reader = csv.DictReader(pfile)
        for line in reader:
            allele_prevalence_freqs[utils.get_region(line['gene'])][line['gene']] = float(line['freq'])
    for region in utils.regions:
        if len(allele_prevalence_freqs[region]) == 0:
            continue
        if debug:
            for gene, freq in allele_prevalence_freqs[region].items():
                print '%14.8f   %s' % (freq, utils.color_gene(gene))
        assert utils.is_normed(allele_prevalence_freqs[region])
    return allele_prevalence_freqs

# ----------------------------------------------------------------------------------------
def choose_allele_prevalence_freqs(glfo, allele_prevalence_freqs, region, min_allele_prevalence_freq, debug=False):
    n_alleles = len(glfo['seqs'][region])
    prevalence_counts = numpy.random.randint(1, int(1. / min_allele_prevalence_freq), size=n_alleles)  # ensures that each pair of alleles has a prevalence ratio between <min_allele_prevalence_freq> and 1. NOTE it's inclusive
    prevalence_freqs = [float(c) / sum(prevalence_counts) for c in prevalence_counts]
    allele_prevalence_freqs[region] = {g : f for g, f in zip(glfo['seqs'][region].keys(), prevalence_freqs)}
    assert utils.is_normed(allele_prevalence_freqs[region])
    if debug:
        print '      counts %s' % ' '.join([('%5d' % c) for c in prevalence_counts])
        print '       freqs %s' % ' '.join([('%5.3f' % c) for c in prevalence_freqs])
        print '   min ratio %.3f' % (min(prevalence_freqs) / max(prevalence_freqs))

# ----------------------------------------------------------------------------------------
def generate_germline_set(glfo, n_genes_per_region, n_alleles_per_gene, min_allele_prevalence_freq, allele_prevalence_fname, debug=True):
    """ NOTE removes genes from  <glfo> """
    if debug:
        print '    choosing germline set'
    allelic_groups = utils.separate_into_allelic_groups(glfo)
    allele_prevalence_freqs = {r : {} for r in utils.regions}
    for region in utils.regions:
        if debug:
            print '  %s' % region
        genes_to_use = set()
        for _ in range(n_genes_per_region[region]):
            choose_some_alleles(region, genes_to_use, allelic_groups, n_alleles_per_gene, debug=debug)
        if debug:
            print '      chose %d alleles' % len(genes_to_use)
        remove_genes(glfo, set(glfo['seqs'][region].keys()) - genes_to_use)  # NOTE would use glutils.restrict_to_genes() but it isn't on a regional basis
        choose_allele_prevalence_freqs(glfo, allele_prevalence_freqs, region, min_allele_prevalence_freq, debug=debug)
    write_allele_prevalence_freqs(allele_prevalence_freqs, allele_prevalence_fname)  # NOTE lumps all the regions together, unlike in the parameter dirs

# ----------------------------------------------------------------------------------------
def check_allele_prevalence_freqs(outfname, glfo, allele_prevalence_fname, only_region=None):
    allele_prevalence_freqs = read_allele_prevalence_freqs(allele_prevalence_fname)
    counts = {r : {g : 0 for g in glfo['seqs'][r]} for r in utils.regions}
    with open(outfname) as outfile:
        reader = csv.DictReader(outfile)
        for line in reader:
            for region in utils.regions:
                counts[region][line[region + '_gene']] += 1
    print '   checking allele prevalence freqs'
    for region in utils.regions:
        if only_region is not None and region != only_region:
            continue
        total = sum(counts[region].values())
        print '       %s   obs / tot  =  freq    expected' % region
        for gene in glfo['seqs'][region]:
            print '          %4d / %-4d = %.3f    %.3f   %s' % (counts[region][gene], total, float(counts[region][gene]) / total, allele_prevalence_freqs[region][gene], utils.color_gene(gene, width=15))

