from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import traceback
import tempfile
import operator
import copy
import sys
import os
import random
import numpy
import glob
from collections import OrderedDict
import csv
from subprocess import check_call, Popen, PIPE

from . import utils
from . import indelutils
from io import open

# ----------------------------------------------------------------------------------------
glfo_dir = 'germline-sets'  # always put germline info into a subdir with this name

# setting defaults here so that bin/test-germline-inference.py and bin/partis don't have to both have defaults in them
# these numbers should really be updated at some point with numbers from corey's and gur's papers
default_n_genes_per_region = {'igh' : '42:18:6',
                              'igk' : '11:1:3',  # light chain d has to be 1, but it's just for the dummy d gene
                              'igl' : '11:1:2'}
default_n_alleles_per_gene = {'igh' : '1.33:1.2:1.2',
                              'igk' : '1.1:1:1.1',
                              'igl' : '1.1:1:1.1'}
default_min_allele_prevalence_freq = 0.1

dummy_d_genes = {l : l.upper() + 'Dx-x*x' if not utils.has_d_gene(l) else None for l in utils.loci}  # e.g. IGKDx-x*x for igk, None for igh

# single-locus file names
def get_fname(gldir, locus, region):
    return gldir + '/' + locus + '/' + locus + region + '.fasta'
def get_extra_fname(gldir, locus):
    return gldir + '/' + locus + '/extras.csv'
def glfo_fasta_fnames(gldir, locus):
    return [get_fname(gldir, locus, r) for r in utils.getregions(locus)]
def glfo_fnames(gldir, locus):
    return [get_extra_fname(gldir, locus), ] + glfo_fasta_fnames(gldir, locus)
def functionality_fname(species=None, gldir=None):  # _not_ generally present (but needs to be present at the moment)
    if species is None:
        species = 'human'
        print('  %s using default species %s for functionality file' % (utils.color('yellow', 'warning'), species))
    if gldir is None:
        gldir = utils.get_partis_dir()
    return '%s/data/germlines/%s/functionalities.csv' % (gldir, species)

csv_headers = ['gene', 'cyst_position', 'tryp_position', 'phen_position', 'aligned_seq']

imgt_info_indices = ('accession-number', 'gene', 'species', 'functionality', 'region')  # , '', '', '', '', '', '', '', '', '', '', '', '')  # here we ignore a bunch of ones at the end that we either don't care about, are blank, or sometimes aren't there
separators = ['()', '[]']  # not actually sure what the parentheses and brackets mean
possible_functionalities = ['F', 'ORF', 'P']
def strip_functionality(funcstr):
    return funcstr.strip(''.join(separators))
def get_imgt_info(infostrs, key):
    return infostrs[imgt_info_indices.index(key)]

duplicate_names = {
    'v' : [
        set(['IGHV3-23*04', 'IGHV3-23D*02']),
        set(['IGHV3-30*18', 'IGHV3-30-5*01']),
        set(['IGHV1-69*01', 'IGHV1-69D*01']),
        set(['IGHV3-30*02', 'IGHV3-30-5*02']),
        set(['IGHV3-30*04', 'IGHV3-30-3*03']),
        set(['IGHV3-23*01', 'IGHV3-23D*01']),
    ],
    'd' : [
        set(['IGHD1/OR15-1a*01', 'IGHD1/OR15-1b*01']),
        set(['IGHD2/OR15-2a*01', 'IGHD2/OR15-2b*01']),
        set(['IGHD4/OR15-4a*01', 'IGHD4/OR15-4b*01']),
        set(['IGHD3/OR15-3a*01', 'IGHD3/OR15-3b*01']),
        set(['IGHD5/OR15-5a*01', 'IGHD5/OR15-5b*01']),
        set(['IGHD5-18*01', 'IGHD5-5*01']),
        set(['IGHD4-11*01', 'IGHD4-4*01']),
    ],
    'j' : []
}

#----------------------------------------------------------------------------------------
def is_snpd(gene):
    primary_version, sub_version, allele = utils.split_gene(gene)
    if '+' in allele or '_' in allele:
        return True
    else:
        return False

#----------------------------------------------------------------------------------------
def is_novel(gene):
    primary_version, sub_version, allele = utils.split_gene(gene)
    return is_snpd(gene) or sub_version == 'x' or (sub_version is not None and len(sub_version) > 4)

#----------------------------------------------------------------------------------------
def convert_to_duplicate_name(glfo, gene):
    for equivalence_class in duplicate_names[utils.get_region(gene)]:
        if gene in equivalence_class:
            for alternate_name in equivalence_class:
                if alternate_name != gene and alternate_name in glfo['seqs'][utils.get_region(gene)]:
                    # print 'converting %s --> %s' % (gene, alternate_name)
                    return alternate_name
    raise Exception('couldn\'t find alternate name for %s (and we\'re probably looking for an alternate name because it wasn\'t in glfo to start with) choices: %s' % (utils.color_gene(gene), utils.color_genes(list(glfo['seqs'][utils.get_region(gene)].keys()))))

#----------------------------------------------------------------------------------------
def check_a_bunch_of_codons(codon, seqons, extra_str='', debug=False):  # seqons: list of (seq, pos) pairs
    """ check a list of sequences, and keep track of some statistics """
    n_total, n_ok, n_too_short, n_bad_codons, n_out_of_frame = 0, 0, 0, 0, 0
    for seq, pos in seqons:
        n_total += 1
        if len(seq) < pos + 3:
            n_too_short += 1
        elif not utils.codon_unmutated(codon, seq, pos):
            n_bad_codons += 1
        elif codon == 'cyst' and not utils.in_frame_germline_v(seq, pos):
            n_out_of_frame += 1
        else:
            n_ok += 1

    if debug:
        print('%s%d %s positions:' % (extra_str, n_total, codon), end=' ')
        if n_ok > 0:
            print('  %d ok' % n_ok, end=' ')
        if n_too_short > 0:
            print('  %d too short' % n_too_short, end=' ')
        if n_bad_codons > 0:
            print('  %d mutated' % n_bad_codons, end=' ')
        if n_out_of_frame > 0:
            print('  %d out of frame' % n_out_of_frame, end=' ')
        print('')

# ----------------------------------------------------------------------------------------
def get_is_imgt_file(infostr):
    if infostr[:3].lower() in utils.loci:  # imgt files have another field (like X97051) before the gene name
        return False
    info_str_list = infostr.split('|')
    if len(info_str_list) < len(imgt_info_indices):  # has to at least have the ones we expect all the time
        return False
    if get_imgt_info(info_str_list, 'gene')[:3].lower() not in utils.loci:
        return False
    region_str = get_imgt_info(info_str_list, 'region')
    if region_str[0].lower() not in utils.regions or region_str[1:] != '-REGION':
        return False
    return True

#----------------------------------------------------------------------------------------
def read_fasta_file(glfo, region, fname, skip_pseudogenes, skip_orfs, aligned=False, add_dummy_name_components=False, locus=None, skip_other_region=False, dont_warn_about_duplicates=False, debug=False):
    glseqs = glfo['seqs'][region]  # shorthand
    n_skipped_pseudogenes, n_skipped_orfs = 0, 0
    seq_to_gene_map = {}
    renamed_genes = []
    is_imgt_file = None
    for seqfo in utils.read_fastx(fname, dont_split_infostrs=True, sanitize_seqs=True):
        if is_imgt_file is None:
            is_imgt_file = get_is_imgt_file(seqfo['infostrs'])  # if the fasta lines aren't all formatted the same, who cares it, should crash somewhere
        # first get gene name
        functy = None
        if (not add_dummy_name_components) and is_imgt_file:  # if it's an imgt file, with a bunch of header info, and the accession number first (if we're adding dummy name components, we just take whatever garbage is first and smash on IGH or whatever)
            seqfo['infostrs'] = seqfo['infostrs'].split('|')
            seqfo['name'] = get_imgt_info(seqfo['infostrs'], 'gene')
            if len(seqfo['infostrs']) < len(imgt_info_indices):
                raise Exception('info str is too short (len %d) to correspond to imgt info indices (len %d):\n  %s\n  %s' % (len(seqfo['infostrs']), len(imgt_info_indices), seqfo['infostrs'], imgt_info_indices))
            gene = seqfo['name']
            functy = strip_functionality(get_imgt_info(seqfo['infostrs'], 'functionality'))
            if functy not in possible_functionalities:
                raise Exception('unexpected functionality %s in %s' % (functy, fname))
            if skip_pseudogenes and functy == 'P':
                n_skipped_pseudogenes += 1
                continue
            if skip_orfs and functy == 'ORF':
                n_skipped_orfs += 1
                continue
        else:  # plain fasta with just the gene name after the '>'
            seqfo['infostrs'] = [s3.strip() for s1 in seqfo['infostrs'].split(' ') for s2 in s1.split('\t') for s3 in s2.split('|')]  # just doing what the old/default is in read_fastx(), in case something depends on that
            seqfo['name'] = seqfo['infostrs'][0]
            gene = seqfo['name']

        if add_dummy_name_components:
            old_name = gene
            gene = utils.construct_valid_gene_name(gene, locus=locus, region=region)
            if gene != old_name:
                cpositions = utils.cdn_positions(glfo, region)  # returns None for d
                if cpositions is not None and old_name in cpositions:
                    cpositions[gene] = cpositions[old_name]
                    del cpositions[old_name]
                if debug:
                    renamed_genes.append((old_name, gene))
        try:
            utils.split_gene(gene)  # just to check if it's a valid gene name
        except:
            lines = traceback.format_exception(*sys.exc_info())
            print(utils.pad_lines(''.join(lines)))
            raise Exception('Unhandled gene name \'%s \' in %s (see above). If you don\'t mind us renaming your genes (we just add locus and dummy allele, e.g. IGH<your stuff>*x), you can set --sanitize-input-germlines.' % (gene, fname))
        if skip_other_region and utils.get_region(gene) != region:
            continue
        if utils.get_region(gene) != region:
            raise Exception('region %s from gene name %s doesn\'t match input region %s' % (utils.get_region(gene), gene, region))
        if gene in glseqs:
            utils.color_mutants(glseqs[gene], utils.remove_gaps(seqfo['seq']), align=True, print_result=True)
            raise Exception('gene name %s appears twice in %s (see the two seqs above)' % (gene, fname))

        # then the sequence
        seq = seqfo['seq']
        if not aligned:
            seq = utils.remove_gaps(seq)
        if 'Y' in seq:
            print('      replacing Y --> N (%d of \'em) in %s' % (seq.count('Y'), utils.color_gene(gene)))
            seq = seq.replace('Y', 'N')
        if len(seq.strip(''.join(utils.expected_characters))) > 0:  # return the empty string if it only contains expected characters
            raise Exception('unexpected character %s in %s (expected %s)' % (seq.strip(''.join(utils.expected_characters)), seq, ' '.join(utils.expected_characters)))
        if seq not in seq_to_gene_map:
            seq_to_gene_map[seq] = []
        seq_to_gene_map[seq].append(gene)

        glseqs[gene] = seq
        if functy is not None and glfo['functionalities'] is not None:
            glfo['functionalities'][gene] = functy  # kind of dumb to split seqs by region but not functionalities, but I kind of wish I hadn't done seqs that way

    if debug and len(renamed_genes) > 0:
        n_max_print = 10
        tmpstr = ',  '.join(('%s --> %s' % (og, ng)) for og, ng in renamed_genes[:n_max_print])
        print('    renamed %d genes from %s: %s%s' % (len(renamed_genes), fname, tmpstr, '  [...]' if len(renamed_genes) > n_max_print else ''))

    tmpcounts = [len(gl) for gl in seq_to_gene_map.values()]  # number of names corresponding to each sequence (should all be ones)
    if tmpcounts.count(1) != len(tmpcounts) and not dont_warn_about_duplicates:
        print('  %s multiple names in %s for the following sequences:' % (utils.color('red', 'warning'), fname))
        for seq, genelist in seq_to_gene_map.items():
            if len(genelist) > 1:
                print('    %-50s   %s' % (' '.join(genelist), seq))
        print('  it is highly recommended to remove the duplicates to avoid ambiguity')

    if n_skipped_pseudogenes > 0:
        print('    skipped %d %s pseudogenes (leaving %d genes)' % (n_skipped_pseudogenes, region, len(glseqs)))
    if n_skipped_orfs > 0:
        print('    skipped %d %s orfs (leaving %d genes)' % (n_skipped_orfs, region, len(glseqs)))

# ----------------------------------------------------------------------------------------
def get_empty_glfo(locus):  # adding this very late, so probably some more places I could use it
    glfo = {'locus' : locus, 'functionalities' : {}, 'seqs' : {r : OrderedDict() for r in utils.regions}}
    for codon in utils.conserved_codons[locus].values():
        glfo[codon + '-positions'] = {}
    return glfo

#----------------------------------------------------------------------------------------
def read_seqs_and_metafo(gldir, locus, skip_pseudogenes, skip_orfs, add_dummy_name_components=False, debug=False):
    glfo = get_empty_glfo(locus)
    read_extra_info(glfo, gldir)
    for region in utils.getregions(locus):
        read_fasta_file(glfo, region, get_fname(gldir, locus, region), skip_pseudogenes, skip_orfs, add_dummy_name_components=add_dummy_name_components, locus=locus, debug=debug)
    if not utils.has_d_gene(locus):  # choose a sequence for the dummy d
        glfo['seqs']['d'][dummy_d_genes[locus]] = 'A'  # this (arbitrary) choice is also made in packages/ham/src/bcrutils.cc
    return glfo

# ----------------------------------------------------------------------------------------
def read_aligned_gl_seqs(fname, glfo, locus, dont_warn_about_duplicates=False):  # only used in partitiondriver with --aligned-germline-fname (which is only used for presto output)
    tmpglfo = {'locus' : locus, 'functionalities' : {}, 'seqs' : {r : OrderedDict() for r in utils.regions}}  # HACK HACK HACK
    for region in utils.regions:
        read_fasta_file(tmpglfo, region, fname, skip_pseudogenes=False, skip_orfs=False, aligned=True, skip_other_region=True, dont_warn_about_duplicates=dont_warn_about_duplicates)
    aligned_gl_seqs = {r : tmpglfo['seqs'][r] for r in utils.regions}
    add_missing_alignments(glfo, aligned_gl_seqs, debug=True)

    # pad all the Vs to the same length (imgt fastas just leave them all unequal lengths on the right of the alignment)
    max_aligned_length = max([len(seq) for seq in aligned_gl_seqs['v'].values()])
    for gene in aligned_gl_seqs['v']:
        n_extra_gaps = max_aligned_length - len(aligned_gl_seqs['v'][gene])
        aligned_gl_seqs['v'][gene] += n_extra_gaps * utils.gap_chars[0]

    # check that we got all the genes
    glfo_genes = set([g for r in utils.regions for g in glfo['seqs'][r]]) - set([dummy_d_genes[glfo['locus']], ])
    aligned_genes = set([g for r in utils.regions for g in aligned_gl_seqs[r]])
    if len(glfo_genes - aligned_genes) > 0:
        raise Exception('missing alignments for %s' % ' '.join(glfo_genes - aligned_genes))

    return aligned_gl_seqs

# ----------------------------------------------------------------------------------------
def add_missing_alignments(glfo, aligned_gl_seqs, debug=False):
    for region in utils.regions:
        _ = get_new_alignments(glfo, region, aligned_seqs=aligned_gl_seqs[region], debug=debug)  # don't need the returned one, since the one we pass in is modified (and it's the same)

# ----------------------------------------------------------------------------------------
def get_new_alignments(glfo, region, aligned_seqs=None, use_old_mafft_merge_method=False, debug=False):
    if aligned_seqs is None:
        aligned_seqs = {}

    genes_with_alignments = set(aligned_seqs)  # used to already have some sequences aligned, and may as well keep around the code to handle that case UPDATE and now we've re-added it!
    genes_without_alignments = set(glfo['seqs'][region]) - set(aligned_seqs)
    if len(genes_without_alignments) == 0:
        if debug:
            print('        no missing %s alignments' % region)
        return

    if debug:
        print('        missing alignments for %d %s gene%s: %s' % (len(genes_without_alignments), region, utils.plural(len(genes_without_alignments)), utils.color_genes(genes_without_alignments)))
        if debug > 1 and len(aligned_seqs) > 0:
            print('      existing alignments:')
            for g, seq in aligned_seqs.items():
                print('            %s   %s' % (seq, utils.color_gene(g)))

    if use_old_mafft_merge_method:  # NOTE could really probably remove this, I doubt I'll really use it any more (I think I only used --merge because I didn't know about --add)
        # find the longest aligned sequence, so we can pad everybody else with dots on the right out to that length
        biggest_length = None
        for gene in genes_with_alignments:
            if biggest_length is None or len(aligned_seqs[gene]) > biggest_length:
                biggest_length = len(aligned_seqs[gene])

        tmpdir = tempfile.mkdtemp()
        already_aligned_fname = tmpdir + '/already-aligned.fasta'
        not_aligned_fname = tmpdir + '/not-aligned.fasta'
        msa_table_fname = tmpdir + '/msa-table.txt'
        aligned_and_not_fname = tmpdir + '/aligned-and-not.fasta'
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

        check_call('cat ' + already_aligned_fname + ' ' + not_aligned_fname + ' >' + aligned_and_not_fname, shell=True)

        # actually run mafft
        cmd = 'mafft --merge ' + msa_table_fname + ' ' + aligned_and_not_fname + ' >' + mafft_outfname  # options=  # "--localpair --maxiterate 1000"
        if debug > 1:
            print('          RUN %s' % cmd)
        proc = Popen(cmd, shell=True, stderr=PIPE, universal_newlines=True)
        _, err = proc.communicate()  # debug info goes to err
        aligned_seqfos = utils.read_fastx(mafft_outfname)
    else:
        aligned_seqfos = utils.align_many_seqs([{'name' : g, 'seq' : glfo['seqs'][region][g]} for g in genes_without_alignments],
                                               existing_aligned_seqfos=[{'name' : g, 'seq' : s} for g, s in aligned_seqs.items()],
                                               ignore_extra_ids=True)

    # overwrite the old alignment with the new one
    for seqfo in aligned_seqfos:
        aligned_seqs[seqfo['name']] = seqfo['seq']

    genes_still_missing = set(genes_without_alignments) - set(aligned_seqs)
    if len(genes_still_missing) > 0:
        raise Exception('missing alignment for %d genes: %s' % (len(genes_still_missing), utils.color_genes(genes_still_missing)))
    if debug:
        print('          added %d new alignments' % len(set(aligned_seqs) & set(genes_without_alignments)))
    if debug > 1:
        print('  new alignments:')
        for g in genes_without_alignments:
            print('            %s   %s  %s' % (aligned_seqs[g], utils.color_gene(g, width=12 if region == 'v' else 8), '<--- new' if g in genes_without_alignments else ''))
        print('')

    if use_old_mafft_merge_method:
        os.remove(already_aligned_fname)
        os.remove(not_aligned_fname)
        os.remove(msa_table_fname)
        os.remove(aligned_and_not_fname)
        os.remove(mafft_outfname)
        os.rmdir(tmpdir)

    return aligned_seqs


#----------------------------------------------------------------------------------------
def get_missing_codon_info(glfo, template_glfo=None, remove_bad_genes=False, debug=False):
    # debug = 2

    for region, codon in utils.conserved_codons[glfo['locus']].items():
        missing_genes = set(glfo['seqs'][region]) - set(glfo[codon + '-positions'])
        if len(missing_genes) == 0:
            if debug:
                print('      no missing %s info' % codon)
            continue

        if debug:
            print('      missing %d %s positions' % (len(missing_genes), codon))

        if template_glfo is not None:  # add one of the genes from the template glfo
            renamed_template_gene = None
            template_gene = list(template_glfo[codon + '-positions'].keys())[0]
            if template_gene in glfo['seqs'][region]:
                if template_gene not in glfo[codon + '-positions']:  # if it's already in the new glfo, but for some reason not in extras.csv (presumably because we don't have an extras.csv)
                    glfo[codon + '-positions'][template_gene] = template_glfo[codon + '-positions'][template_gene]
            else:
                renamed_template_gene = template_gene + '_TEMPLATE'
                if debug:
                    print('    adding gene with codon position %s from template glfo, to align against new genes' % utils.color_gene(template_gene))
                if renamed_template_gene in glfo['seqs'][region]:  # ok there's not really any way that could happen
                    raise Exception('%s already in glfo' % renamed_template_gene)
                newfo = {'gene' : renamed_template_gene, 'seq' : template_glfo['seqs'][region][template_gene], 'cpos' : template_glfo[codon + '-positions'][template_gene]}
                if newfo['seq'] in list(glfo['seqs'][region].values()):  # if it's in there under a different name, just tweak the seq a bit
                    if debug:
                        print('    had to add an N to the right side of template glfo\'s gene %s, since it\'s in <glfo> under a different name %s' % (utils.color_gene(template_gene), ' '.join([utils.color_gene(g) for g, s in glfo['seqs'][region].items() if s == newfo['seq']])))
                    newfo['seq'] += 'N'  # this is kind of hackey, but we don't have Ns in current germline seqs, so it at least shouldn't clash with anybody
                add_new_allele(glfo, newfo, use_template_for_codon_info=False, debug=debug)

        aligned_seqs = get_new_alignments(glfo, region, debug=debug)

        # if region == 'j':
        #     raise Exception('missing tryp position for %s, and we can\'t infer it because tryp positions don\'t reliably align to the same position' % ' '.join(missing_genes))

        # existing codon position (this assumes that once aligned, all genes have the same codon position -- which is only really true for the imgt-gapped alignment)
        if len(glfo[codon + '-positions']) > 0:
            known_gene, known_pos = None, None
            known_but_not_in_glfo, known_but_unaligned, known_but_mutated = [], [], []
            for gene, pos in glfo[codon + '-positions'].items():  # take the first one for which we have the sequence (NOTE it would be safer to check that they're all the same)
                if gene not in glfo['seqs'][region]:
                    known_but_not_in_glfo.append(gene)
                    continue
                if gene not in aligned_seqs:
                    known_but_unaligned.append(gene)
                    continue
                if not utils.codon_unmutated(codon, glfo['seqs'][region][gene], pos):
                    known_but_mutated.append(gene)
                    continue
                known_gene, known_pos = gene, pos
                break
            if known_gene is None:
                print('genes with missing codon info: %s' % ' '.join(sorted(missing_genes)))
                print('        known but not in glfo: %s' % ' '.join(sorted(known_but_not_in_glfo)))
                print('          known but unaligned: %s' % ' '.join(sorted(known_but_unaligned)))
                print('            known but mutated: %s' % ' '.join(sorted(known_but_mutated)))
                raise Exception('couldn\'t find a known %s position (see above)' % codon)
            # NOTE for cyst, should be 309 if alignments are imgt [which they used to usually be, but now probably aren't] (imgt says 104th codon --> subtract 1 to get zero-indexing, then multiply by three 3 * (104 - 1) = 309
            known_pos_in_alignment = utils.get_codon_pos_in_alignment(codon, aligned_seqs[known_gene], glfo['seqs'][region][known_gene], known_pos, known_gene)
            if debug:
                print('  using known position %d (aligned %d) from %s' % (known_pos, known_pos_in_alignment, known_gene))
        elif codon == 'cyst':
            known_pos_in_alignment = 309
            print('      assuming aligned %s position is %d (this will %s work if you\'re using imgt alignments)' % (codon, known_pos_in_alignment, utils.color('red', 'only')))
            raise Exception('not really using imgt alignments much any more, so this isn\'t really going to work')
        else:
            raise Exception('No existing conserved codon position info for \'%s\', and don\'t have any genes with info from which to guess it. If you set --sanitize-input-germlines, though, we\'ll use the default germlines in data/germlines/ as a template.' % codon)

        n_added = 0
        seqons = []  # (seq, pos) pairs
        bad_codons = []
        for gene in [known_gene] + list(missing_genes):
            unaligned_pos = known_pos_in_alignment - utils.count_gap_chars(aligned_seqs[gene], aligned_pos=known_pos_in_alignment)
            seq_to_check = glfo['seqs'][region][gene]
            seqons.append((seq_to_check, unaligned_pos))
            glfo[codon + '-positions'][gene] = unaligned_pos
            n_added += 1

            tmpseq = aligned_seqs[gene]  # NOTE this is aligned
            tmppos = known_pos_in_alignment  # NOTE this is aligned
            if not utils.codon_unmutated(codon, tmpseq, tmppos) or (region == 'v' and not utils.in_frame_germline_v(tmpseq, tmppos)):
                bad_codons.append(gene)
            if debug > 1:
                print('       %3s  %s%s%s   %s %3s %5s' % (utils.color('red', 'bad') if gene in bad_codons else '',
                                                           tmpseq[:tmppos],
                                                           utils.color('reverse_video', tmpseq[tmppos : tmppos + 3]),
                                                           tmpseq[tmppos + 3:],
                                                           utils.color_gene(gene, width=12 if region == 'v' else 8),
                                                           utils.color('red', 'bad') if gene in bad_codons else '',
                                                           'new' if gene != known_gene else ''))

        check_a_bunch_of_codons(codon, seqons, extra_str='          ', debug=debug)  # kind of redundant with the check that happens in the loop above, but prints some summary info
        if len(bad_codons) > 0:
            print('    %s %d bad %s positions (%s)' % (utils.wrnstr(), len(bad_codons), codon, 'removing those genes, since <remove_bad_genes>' if remove_bad_genes else 'leaving, since <remove_bad_genes> isn\'t set'))
            if remove_bad_genes:
                remove_genes(glfo, bad_codons, debug=True)
        if debug:
            print('      added %d %s positions' % (n_added, codon))

        if template_glfo is not None and renamed_template_gene is not None:
            remove_gene(glfo, renamed_template_gene)

# ----------------------------------------------------------------------------------------
def get_merged_glfo(glfo_a, glfo_b, debug=False):  # doesn't modify either of the arguments
    if debug:
        print('  merging glfos')
    assert set(glfo_a['seqs']) == set(glfo_b['seqs'])
    merged_glfo = copy.deepcopy(glfo_a)
    name_mapping = {region: {} for region in utils.regions}  # maps dropped names (from glfo_b) to retained names (from glfo_a)
    for region in utils.regions:
        duplicate_genes, duplicate_seqs = [], []
        inconsistent_names = []  # names corresponding to <duplicate_seqs>
        merged_seqs = set(merged_glfo['seqs'][region].values())
        for gene, seq in glfo_b['seqs'][region].items():
            if gene in merged_glfo['seqs'][region]:
                duplicate_genes.append(gene)
                continue
            if seq in merged_seqs:
                duplicate_seqs.append(seq)
                a_names = [g for g, s in glfo_a['seqs'][region].items() if s == seq]
                assert len(a_names) == 1  # I really don't want to deal with duplicate sequences here. They should've been cleaned up elsewhere
                inconsistent_names.append((a_names[0], gene))
                name_mapping[region][gene] = a_names[0]  # track that gene name from glfo_b maps to a_names[0] from glfo_a
                continue
            add_new_allele(merged_glfo, {'gene' : gene, 'seq' : seq, 'cpos' : utils.cdn_pos(glfo_b, region, gene)}, use_template_for_codon_info=False)
            merged_seqs.add(seq)
        if debug:
            print('   %s: %d + %d --> %d' % (region, len(glfo_a['seqs'][region]), len(glfo_b['seqs'][region]), len(merged_glfo['seqs'][region])))
            genes_from_a = set(merged_glfo['seqs'][region]) - set(glfo_b['seqs'][region])
            genes_from_b = set(merged_glfo['seqs'][region]) - set(glfo_a['seqs'][region])
            print('     %d only from a: %s' % (len(genes_from_a), utils.color_genes(genes_from_a)))
            print('     %d only from b: %s' % (len(genes_from_b), utils.color_genes(genes_from_b)))
            if len(duplicate_genes) > 0:
                print('     %d gene names in both: %s' % (len(duplicate_genes), utils.color_genes(duplicate_genes)))
            # assert set(duplicate_genes) | genes_from_a | genes_from_b == set(merged_glfo['seqs'][region])
        for dgene in duplicate_genes:  # check for inconsistent sequences for the same name
            if glfo_a['seqs'][region][dgene] != glfo_b['seqs'][region][dgene]:
                print('      %s different seqs for name %s' % (utils.color('red', 'warning'), utils.color_gene(dgene)))
                utils.color_mutants(glfo_a['seqs'][region][dgene], glfo_b['seqs'][region][dgene], align=True, print_result=True, extra_str='        ')
        if len(duplicate_seqs) > 0:
            print('     %d seqs in both, but with different names (we use the name from glfo_a, the first arg): %s' % (len(duplicate_seqs), '   '.join([('%s %s' % (utils.color_gene(ga), utils.color_gene(gb))) for ga, gb in inconsistent_names])))

    return merged_glfo, name_mapping

#----------------------------------------------------------------------------------------
def remove_extraneouse_info(glfo, debug=False):
    """ remove codon info corresponding to genes that aren't in 'seqs' """
    for region, codon in utils.conserved_codons[glfo['locus']].items():
        genes_to_remove = set(glfo[codon + '-positions']) - set(glfo['seqs'][region])
        if debug:
            print('    removing %s info for %d genes (leaving %d)' % (codon, len(genes_to_remove), len(glfo[codon + '-positions']) - len(genes_to_remove)))
        for gene in genes_to_remove:
                del glfo[codon + '-positions'][gene]

# ----------------------------------------------------------------------------------------
def read_extra_info(glfo, gldir):
    if not os.path.exists(get_extra_fname(gldir, glfo['locus'])):
        print('  %s extra file %s not found (i.e. no codon positions for glfo)' % (utils.color('yellow', 'warning'), get_extra_fname(gldir, glfo['locus'])))
        return
    with open(get_extra_fname(gldir, glfo['locus'])) as csvfile:
        reader = csv.DictReader(csvfile)
        for line in reader:
            for codon in utils.conserved_codons[glfo['locus']].values():
                if line[codon + '_position'] != '':
                    glfo[codon + '-positions'][line['gene']] = int(line[codon + '_position'])

#----------------------------------------------------------------------------------------
# groups by gene family unless <use_primary_version> is set
def print_glfo(glfo, use_primary_version=False, gene_groups=None, print_separate_cons_seqs=False, only_region=None, input_groupfcn=None):  # NOTE kind of similar to bin/cf-alleles.py
    if gene_groups is None:
        gene_groups = {}
        for region in (utils.regions if only_region is None else [only_region]):
            if input_groupfcn is None:
                groupfcn = utils.primary_version if use_primary_version else utils.gene_family
            else:
                groupfcn = input_groupfcn
            group_labels = sorted(set([groupfcn(g) for g in glfo['seqs'][region]]))
            gene_groups[region] = [(glabel, {g : glfo['seqs'][region][g] for g in glfo['seqs'][region] if groupfcn(g) == glabel}) for glabel in group_labels]

    for region in [r for r in utils.regions if r in gene_groups]:
        print('%s' % utils.color('reverse_video', utils.color('green', region)))
        for group_label, group_seqs in gene_groups[region]:
            print('  %s' % utils.color('blue', group_label))
            workdir = tempfile.mkdtemp()
            with tempfile.NamedTemporaryFile() as tmpfile:  # kind of hilarious that i use vsearch here, but mafft up there... oh well it shouldn't matter
                _ = utils.run_vsearch('cluster', group_seqs, workdir, threshold=0.3, minseqlength=5, msa_fname=tmpfile.name)  # <threshold> is kind of random, i just set it to something that seems to group all the V genes with the same ggroup together
                # NOTE I"m not sure why the stupid thing sometimes splits the group apart no matter what I set <threshold> to
                msa_seqs = utils.read_fastx(tmpfile.name, ftype='fa')
            msa_info = []
            for seqfo in msa_seqs:
                # print '    %s    %s' % (seqfo['seq'], seqfo['name'])
                if seqfo['name'][0] == '*':  # start of new cluster (centroid is first, and is marked with a '*')
                    centroid = seqfo['name'].lstrip('*')
                    msa_info.append({'centroid' : centroid, 'seqfos' : [{'name' : centroid, 'seq' : seqfo['seq']}]})
                elif seqfo['name'] == 'consensus':
                    msa_info[-1]['cons_seq'] = seqfo['seq'].replace('+', '-')  # gaaaaah not sure what the +s mean
                else:
                    msa_info[-1]['seqfos'].append(seqfo)

            first_cons_seq = None
            for clusterfo in msa_info:
                align = region == 'd' or input_groupfcn is not None
                if first_cons_seq is None:  # shenanigans to account for vsearch splitting up my groups
                    first_cons_seq = clusterfo['cons_seq']
                    print('    %s    consensus (first cluster)' % clusterfo['cons_seq'])
                else:  # these are mostly necessary for d
                    if print_separate_cons_seqs:
                        # print '    %s    extra consensus' % utils.color_mutants(first_cons_seq, clusterfo['cons_seq'], align=align)  # not very informative to print this, since the first-cluster-consensus was already printed without the proper gaps to make this make sense
                        print('    %s    extra consensus' % clusterfo['cons_seq'])
                for seqfo in clusterfo['seqfos']:
                    emphasis_positions = None
                    extra_str = ''

                    if region in utils.conserved_codons[glfo['locus']]:
                        codon = utils.conserved_codons[glfo['locus']][region]
                        aligned_cpos = utils.get_codon_pos_in_alignment(codon, seqfo['seq'], group_seqs[seqfo['name']], utils.cdn_pos(glfo, region, seqfo['name']), seqfo['name'])
                        emphasis_positions = [aligned_cpos + i for i in range(3)]
                        if region == 'v':
                            if utils.cdn_pos(glfo, region, seqfo['name']) % 3 != 0:  # flag out of frame cysteines
                                extra_str += '%s %s frame (%d)' % (utils.color('red', 'bad'), codon, utils.cdn_pos(glfo, region, seqfo['name']))
                            if not utils.codon_unmutated(codon, glfo['seqs'][region][seqfo['name']], utils.cdn_pos(glfo, region, seqfo['name'])) or not utils.in_frame_germline_v(glfo['seqs'][region][seqfo['name']], utils.cdn_pos(glfo, region, seqfo['name'])):
                                extra_str += '   %s codon' % utils.color('red', 'bad')

                    # NOTE if you <align> below here the codon info is wrong (and I _think_ I've fixed it so I no longer ever need to align) UPDATE nope, just happened, so adding the continue below

                    cons_seq = clusterfo['cons_seq'] if print_separate_cons_seqs else first_cons_seq
                    leftpad = len(cons_seq) - len(cons_seq.lstrip('-'))
                    if not align and len(cons_seq) != len(seqfo['seq']):
                        print('%s cons seq and seq different lengths, so skipping %s' % (utils.color('red', 'error'), utils.color_gene(seqfo['name'])))
                        continue
                    print('    %s%s    %s      %s' % (' ' * leftpad, utils.color_mutants(cons_seq, seqfo['seq'], align=align, emphasis_positions=emphasis_positions), utils.color_gene(seqfo['name']), extra_str))

#----------------------------------------------------------------------------------------
def read_glfo(gldir, locus, only_genes=None, skip_pseudogenes=True, skip_orfs=True, remove_orfs=False, template_glfo=None, remove_bad_genes=False, add_dummy_name_components=False, dont_crash=False, debug=False):  # <skip_orfs> is for use when reading just-downloaded imgt files, while <remove_orfs> tells us to look for a separate functionality file
    # NOTE <skip_pseudogenes> and <skip_orfs> only have an effect with just-downloaded imgt files (otherwise we don't in general have the functionality info)
    if not os.path.exists(gldir + '/' + locus):  # NOTE doesn't re-link it if we already made the link before
        if locus[:2] == 'ig' and os.path.exists(gldir + '/' + locus[2]):  # backwards compatibility
            print('    note: linking new germline dir name to old name in %s' % gldir)
            check_call('cd '  + gldir + ' && ln -sfv ' + locus[2] + ' ' + locus, shell=True)
        else:
            if dont_crash:
                return None
            else:
                raise Exception('germline set directory \'%s\' does not exist (maybe --parameter-dir is corrupted or mis-specified?)' % (gldir + '/' + locus))

    if debug:
        print('  reading %s locus glfo from %s' % (locus, gldir))
    glfo = read_seqs_and_metafo(gldir, locus, skip_pseudogenes, skip_orfs, add_dummy_name_components=add_dummy_name_components, debug=debug)
    get_missing_codon_info(glfo, template_glfo=template_glfo, remove_bad_genes=remove_bad_genes, debug=debug)
    restrict_to_genes(glfo, only_genes, debug=debug)

    for region, codon in utils.conserved_codons[glfo['locus']].items():
        seqons = [(seq, glfo[codon + '-positions'][gene]) for gene, seq in glfo['seqs'][region].items()]  # (seq, pos) pairs
        check_a_bunch_of_codons(codon, seqons, extra_str='      ', debug=debug)

    if debug:
        print('  read %s' % '  '.join([('%s: %d' % (r, len(glfo['seqs'][r]))) for r in utils.regions]))

    if remove_orfs:  # shouldn't need this any more (ended up removing 'em from the default glfos), but I don't feel like removing it at the moment
        orfs_removed = {r : [] for r in utils.regions}
        assert len(glfo['functionalities']) == 0
        func_info = {}  # keep them in a spaerate dict so it's easier to loop over the genes in the glfo (to ensure they're all in the functionality file)
        with open(functionality_fname(species=None, gldir=None)) as ffile:  # using basename() to get the species is hackey... but it'll typically work in the (rare) situations where this code'll get called
            reader = csv.DictReader(ffile)
            for line in reader:
                func_info[line['gene']] = line['functionality']
        for gene in [g for r in utils.regions for g in glfo['seqs'][r]]:
            if is_novel(gene):
                func_info[gene] = 'F'
            if gene not in func_info:
                raise Exception('no func info for %s' % utils.color_gene(gene))
            if func_info[gene] == 'ORF':
                remove_gene(glfo, gene)
                orfs_removed[utils.get_region(gene)].append(gene)
                continue
            assert func_info[gene] == 'F'  # should've already removed all the pseudogenes (and there shouldn't be any other functionality)
            glfo['functionalities'][gene] = func_info[gene]
        if len([ors for ors in orfs_removed.values()]) > 0:
            print('     removed %2d ORFs: %s' % (sum(len(ors) for ors in orfs_removed.values()), '  '.join(['%s %s' % (r, str(len(orfs_removed[r])) if len(orfs_removed[r]) > 0 else ' ') for r in utils.regions])), end=' ')
            if len(orfs_removed['v']) > 0:
                print('   (%s)' % ' '.join([utils.color_gene(g) for g in sorted(orfs_removed['v'])]), end=' ')
            print('')

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
def split_inferred_allele_name(gene_name, debug=False):
    if '+' in gene_name:  # partis snp'd alleles (e.g. IGHVx-y*z+G35T.A77C)
        if len(gene_name.split('+')) != 2:
            raise Exception('couldn\'t split \'%s\' into two pieces from \'+\'' % gene_name)
        method = 'partis'
        template_name, mutstrs = gene_name.split('+')
        if mutstrs[0] not in utils.nukes:  # hashed name
            int(mutstrs)  # make sure it's actually a hash
            mutstrs = None
        else:
            mutstrs = mutstrs.split('.')
    elif utils.sub_version(gene_name) == 'x' or len(utils.sub_version(gene_name)) > 4:  # partis indel'd alleles
        method = 'partis'
        template_name = gene_name
        mutstrs = None
    elif '_' in gene_name:
        template_name = gene_name.split('_')[0]
        mutstrs = gene_name.split('_')[1:]  # always a list of strings
        if len(mutstrs) == 1 and mutstrs[0][0] == 'S':  # igdiscover (e.g. IGHVx-y*z_S1098)
            method = 'igdiscover'
            mutstrs = None
        else:  # tigger (e.g. IGHVx-y*z_G35T_A77C)
            method = 'tigger'
    else:
        raise Exception('couldn\'t figure out template gene and snp info for %s' % gene_name)

    return method, template_name, mutstrs

# ----------------------------------------------------------------------------------------
def get_template_gene(gene_name, debug=False):
    _, template_name, _ = split_inferred_allele_name(gene_name, debug=debug)
    return template_name

# ----------------------------------------------------------------------------------------
def try_to_get_name_and_mutfo_from_seq(gene_name, new_seq, glfo, consider_indels=True, debug=False):  # it's way slower with <consider_indels>, since find_nearest_gene_in_glfo() aligns everybody first
    assert is_novel(gene_name)
    method, _, _ = split_inferred_allele_name(gene_name)
    if method != 'igdiscover':
        raise Exception('unhandled method %s' % method)

    if consider_indels:
        nearest_gene, n_snps, n_indels, realigned_new_seq, realigned_nearest_seq = find_nearest_gene_in_glfo(glfo, new_seq, new_name=gene_name, exclusion_3p=3, debug=debug)
        if n_indels > 0:
            return None, None, None
        _, positions = utils.hamming_distance(realigned_new_seq, realigned_nearest_seq, return_mutated_positions=True)
        snpfo = {pos : {'original' : realigned_nearest_seq[pos], 'new' : realigned_new_seq[pos]} for pos in positions}
    else:
        n_snps, nearest_gene, nearest_seq = find_nearest_gene_with_same_cpos(glfo, new_seq, new_name=gene_name, exclusion_3p=3, debug=debug)
        if n_snps is None:
            return None, None, None
        tmp_len = min(len(new_seq), len(nearest_seq))  # this should be ok, since find_nearest_gene_with_same_cpos() will have ensure that their cyst positions are the same
        _, positions = utils.hamming_distance(new_seq[:tmp_len], nearest_seq[:tmp_len], return_mutated_positions=True)
        snpfo = {pos : {'original' : nearest_seq[pos], 'new' : new_seq[pos]} for pos in positions}

    new_name, snpfo = choose_new_allele_name(nearest_gene, new_seq, snpfo=snpfo)  # NOTE <nearest_distance> is aligned, so if there's gaps to the nearest gene it's not really right

    return new_name, snpfo, nearest_gene

# ----------------------------------------------------------------------------------------
def try_to_get_mutfo_from_name(gene_name, aligned_seq=None, unaligned_seq=None, debug=False):
    method, _, mutstrs = split_inferred_allele_name(gene_name)

    if mutstrs is None:
        return None

    mutfo = {}
    for mutstr in mutstrs:
        if len(mutstr) < 3:
            print('couldn\'t extract mutation info from \'%s\' in gene %s' % (mutstr, gene_name))
            return None
        original, new = mutstr[0], mutstr[-1]
        if original not in utils.nukes or new not in utils.nukes:
            print('couldn\'t extract mutation info from \'%s\' in gene %s' % (mutstr, gene_name))
            return None
        try:
            position = int(mutstr[1:-1])
        except ValueError:
            print('couldn\'t convert \'%s\' to a position in gene %s' % (mutstr[1:-1], gene_name))
            return None
        if method == 'tigger':  # convert from 1-indexed aligned to 0-indexed unaligned
            assert aligned_seq is not None
            unaligned_seq = utils.remove_gaps(aligned_seq)
            imgt_aligned_pos = position - 1  # convert to 0-indexed
            if debug:
                print('    imgt aligned %s%d%s: %s%s%s' % (original, imgt_aligned_pos, new, aligned_seq[:imgt_aligned_pos], utils.color('red', aligned_seq[imgt_aligned_pos]), aligned_seq[imgt_aligned_pos + 1:]))
            assert aligned_seq[imgt_aligned_pos] == new
            n_gaps = utils.count_gap_chars(aligned_seq, aligned_pos=imgt_aligned_pos)
            unaligned_pos = imgt_aligned_pos - n_gaps
            if debug:
                print('       unaligned %s%d%s: %s%s%s' % (original, unaligned_pos, new, unaligned_seq[:unaligned_pos], utils.color('red', unaligned_seq[unaligned_pos]), unaligned_seq[unaligned_pos + 1:]))
            assert unaligned_seq[unaligned_pos] == new
            position = unaligned_pos
        if position not in mutfo:
            mutfo[position] = {}
            mutfo[position]['original'] = original  # if it *is* already there, we want to *keep* the old 'original'
        mutfo[position]['new'] = new
        if mutfo[position]['new'] == mutfo[position]['original']:  # reverted back to the original base
            del mutfo[position]

    return mutfo

# ----------------------------------------------------------------------------------------
def simplify_snpfo(template_gene, mutfo):
    """ convert snp info in <mutfo> to the snpd/inferred name (with '+'s and whatnot), accounting for <template_gene> being possibly already snpd/inferred """
    old_mutfo = try_to_get_mutfo_from_name(template_gene)  # start from the template gene's snpd positions
    if old_mutfo is None:
        return None

    for position, posfo in mutfo.items():  # loop over the new gene's snpd positions
        if position not in old_mutfo:  # ...adding the new positions that aren't already there
            old_mutfo[position] = {'original' : posfo['original']}
        old_mutfo[position]['new'] = posfo['new']
        if old_mutfo[position]['new'] == old_mutfo[position]['original']:  # the new mutation reverted the position back to the original base
            del old_mutfo[position]
    final_mutfo = old_mutfo  # <old_mutfo> is modified at this point, so I should really call it something else here
    return final_mutfo

# ----------------------------------------------------------------------------------------
def generate_single_new_allele(template_gene, template_cpos, template_seq, snp_positions, indel_positions):
    assert utils.get_region(template_gene) == 'v'  # others not yet handled
    mean_indel_length = 3  # I really can't think of a reason that the indel lengths are all that important, so it's staying hard-coded here for now (could of course use the mean shm indel length from recombinator, but that's confusig since we're not modeling shm here)
    already_snpd_positions = set()  # only used if a position wasn't specified (i.e. was None) in <snps_to_add> (also not used for indels, which probably makes sense)

    def choose_position(seq, cpos):
        chosen_pos = None
        n_tries = 0
        codon_already_mutated = not utils.codon_unmutated('cyst', seq, cpos)
        while chosen_pos is None or chosen_pos in already_snpd_positions or (not codon_already_mutated and not utils.codon_unmutated('cyst', tmpseq, cpos)):
            chosen_pos = random.randint(0, cpos - 1)  # len(seq) - 1)  # note that randint() is inclusive
            tmpseq = seq[: chosen_pos] + 'X' + seq[chosen_pos + 1 :]  # only used for checking cyst position
            n_tries += 1
            if n_tries > 10000:
                raise Exception('too many tries %d while choosing position in generate_single_new_allele()' % n_tries)
        return chosen_pos

    new_cpos = template_cpos
    new_seq = template_seq

    # first do indels (it's kind of arbitrary to do indels before snps, but it'd be pretty common for a deletion to annihilate a snp, and doing indels first avoids that (at the "cost" of making the snp positions more obtuse))
    codon_positions = {'v' : new_cpos}  # do *not* use <new_cpos> itself from this point until it's re-set after the loop
    indel_location = None
    if None in indel_positions:
        assert indel_positions.count(None) == len(indel_positions)  # doesn't really make sense to specify one, but not all,  of them
        indel_location = 'v'
    new_seq, indelfo = indelutils.add_indels(len(indel_positions), new_seq, new_seq, mean_indel_length, codon_positions, indel_location=indel_location, indel_positions=indel_positions, keep_in_frame=True)
    new_cpos = codon_positions['v']  # ick

    # then do snps
    snpfo = OrderedDict()  # this is only really used for figuring out the final name
    for snp_pos in snp_positions:
        if snp_pos is None:
            snp_pos = choose_position(new_seq, new_cpos)
        already_snpd_positions.add(snp_pos)
        new_base = None
        while new_base is None or new_base == new_seq[snp_pos]:
            new_base = utils.nukes[random.randint(0, len(utils.nukes) - 1)]
        print('        %3d   %s --> %s' % (snp_pos, new_seq[snp_pos], new_base))
        snpfo[snp_pos] = {'original' : new_seq[snp_pos], 'new' : new_base}

        new_seq = new_seq[: snp_pos] + new_base + new_seq[snp_pos + 1 :]

    if not utils.codon_unmutated('cyst', new_seq, new_cpos):
        print('  reverting messed up codon')
        new_seq = new_seq[:new_cpos] + utils.codon_table['cyst'][0] + new_seq[new_cpos + 3 :]
        assert utils.codon_unmutated('cyst', new_seq, new_cpos, debug=True)

    new_name, snpfo = choose_new_allele_name(template_gene, new_seq, snpfo=snpfo, indelfo=indelfo)  # shouldn't actually change <snpfo>

    return {'template-gene' : template_gene, 'gene' : new_name, 'seq' : new_seq, 'cpos' : new_cpos}

# ----------------------------------------------------------------------------------------
def restrict_to_genes(glfo, only_genes, debug=False):
    """
    Remove from <glfo> any genes which are not in <only_genes>.
    Any regions which are not represented in a non-None <only_genes> will be *un*restricted (i.e. any gene from that region is fair game).
    """
    if only_genes is None:
        return

    restricted_regions = set([utils.get_region(g) for g in only_genes])
    unrestricted_regions = set(utils.regions) - restricted_regions

    only_genes_not_in_glfo = set(only_genes) - set([g for r in restricted_regions for g in glfo['seqs'][r]])
    if len(only_genes_not_in_glfo) > 0:
        print('  %s genes %s in <only_genes> aren\'t in glfo to begin with' % (utils.wrnstr(), ' '.join(only_genes_not_in_glfo)))

    genes_to_remove = set([g for r in restricted_regions for g in glfo['seqs'][r]]) - set(only_genes)
    remove_genes(glfo, genes_to_remove)
    if debug:
        print('    removed %-2d genes from glfo (leaving %s)' % (len(genes_to_remove), '  '.join(['%s %-d' % (r, len(glfo['seqs'][r])) for r in utils.regions])))

# ----------------------------------------------------------------------------------------
def restrict_to_observed_genes(glfo, parameter_dir, debug=False):  # remove from <glfo> any genes that were not observed in <parameter_dir>
    only_genes = set([g for genes in utils.read_overall_gene_probs(parameter_dir).values() for g in genes])
    restrict_to_genes(glfo, only_genes, debug=debug)

# ----------------------------------------------------------------------------------------
def remove_v_genes_with_bad_cysteines(glfo, debug=False):
    prelength = len(glfo['seqs']['v'])
    n_mutated, n_out_of_frame = 0, 0
    for gene in glfo['seqs']['v'].keys():  # have to use a copy of the keys, since we modify the dict in the loop
        mutated = not utils.codon_unmutated('cyst', glfo['seqs']['v'][gene], glfo['cyst-positions'][gene])
        in_frame = utils.in_frame_germline_v(glfo['seqs']['v'][gene], glfo['cyst-positions'][gene])
        if mutated:
            n_mutated += 1
        if not in_frame:
            n_out_of_frame += 1
        if mutated or not in_frame:
            remove_gene(glfo, gene, debug=debug)
    if True:  # debug:
        print('  removed %d / %d v genes with bad cysteines (%d mutated, %d out of frame)' % (prelength - len(glfo['seqs']['v']), len(glfo['seqs']['v']), n_mutated, n_out_of_frame))

# ----------------------------------------------------------------------------------------
def remove_genes(glfo, genes, debug=False):
    """ remove <genes> from <glfo> """
    if debug:
        print('  removing %d gene%s from glfo: %s' % (len(genes), utils.plural(len(genes)), ' '.join([utils.color_gene(g) for g in genes])))

    for gene in genes:
        remove_gene(glfo, gene)

    for region in utils.regions:
        if len(glfo['seqs'][region]) == 0:
            print('%s removed all %s genes from glfo' % (utils.color('yellow', 'warning'), region))

# ----------------------------------------------------------------------------------------
def remove_gene(glfo, gene, debug=False):
    """ remove <gene> from <glfo> """
    region = utils.get_region(gene)
    if gene in glfo['seqs'][region]:
        if debug:
            print('  removing %s from glfo' % utils.color_gene(gene))
        del glfo['seqs'][region][gene]
        if region in utils.conserved_codons[glfo['locus']]:
            del glfo[utils.conserved_codons[glfo['locus']][region] + '-positions'][gene]
    else:
        print('  %s tried to remove %s from glfo, but it isn\'t there' % (utils.color('yellow', 'warning'), utils.color_gene(gene)))

# ----------------------------------------------------------------------------------------
def add_new_alleles(glfo, newfos, remove_template_genes=False, use_template_for_codon_info=True, simglfo=None, debug=False):
    for newfo in newfos:
        add_new_allele(glfo, newfo, remove_template_genes=remove_template_genes, use_template_for_codon_info=use_template_for_codon_info, simglfo=simglfo, debug=debug)

# ----------------------------------------------------------------------------------------
def rename_gene(glfo, old_name, new_name, debug=False):
    if debug:
        print('    renaming %s --> %s' % (utils.color_gene(old_name), utils.color_gene(new_name)))
    region = utils.get_region(old_name)
    codon = utils.cdn(glfo, region)
    glfo[codon + '-positions'][new_name] = glfo[codon + '-positions'][old_name]
    glfo['seqs'][region][new_name] = glfo['seqs'][region][old_name]
    remove_gene(glfo, old_name)

# ----------------------------------------------------------------------------------------
def add_new_allele(glfo, newfo, remove_template_genes=False, use_template_for_codon_info=True, simglfo=None, debug=False):
    # NOTE <use_template_for_codon_info> has default True
    """
    Add an allele to <glfo>, specified by <newfo> which has mandatory keys: ['gene', 'seq'] and optional keys: ['cpos', 'template-gene', 'remove-template-gene'] (note c in 'cpos' is for 'codon', not 'cysteinie', i.e. it's not just for v)
    If <remove_template_genes>, or if 'remove-template-gene' in <newfo>, we also remove 'template-gene' from <glfo>.
    Fails if name or seq is already in glfo.
    """
    region = utils.get_region(newfo['gene'])

    if len(set(newfo['seq']) - utils.alphabet) > 0:
        raise Exception('unexpected characters %s in new gl seq %s' % (set(newfo['seq']) - utils.alphabet, newfo['seq']))
    if newfo['gene'] in glfo['seqs'][region]:
        print('  %s attempted to add name %s that already exists in glfo (returning without doing anything)' % (utils.color('yellow', 'warning'), utils.color_gene(newfo['gene'])))
        return
    if newfo['seq'] in list(glfo['seqs'][region].values()):  # it's nice to allow duplicate seqs e.g. if you want to rename a gene
        names_with_this_seq = [utils.color_gene(g) for g, seq in glfo['seqs'][region].items() if seq == newfo['seq']]
        print('  %s attempted to add sequence (with name %s) that\'s already in glfo with name%s %s (returning without doing anything)\n    %s' % (utils.color('yellow', 'warning'), utils.color_gene(newfo['gene']), utils.plural(len(names_with_this_seq)), ' '.join(names_with_this_seq), newfo['seq']))
        return

    glfo['seqs'][region][newfo['gene']] = newfo['seq']

    if region in utils.conserved_codons[glfo['locus']]:
        codon = utils.cdn(glfo, region)
        if use_template_for_codon_info:
            if 'cpos' in newfo and newfo['cpos'] is not None:
                raise Exception('cpos set in newfo, but <use_template_for_codon_info> also set to True')
            if 'template-gene' not in newfo:
                raise Exception('no template gene specified in newfo')
            if newfo['template-gene'] not in glfo[codon + '-positions']:
                raise Exception('template gene %s not found in codon info' % newfo['template-gene'])
            glfo[codon + '-positions'][newfo['gene']] = glfo[codon + '-positions'][newfo['template-gene']]
        elif 'cpos' in newfo and newfo['cpos'] is not None:  # codon position explicitly set in newfo
            glfo[codon + '-positions'][newfo['gene']] = newfo['cpos']
        else:
            get_missing_codon_info(glfo)

    if debug:
        simstr = ''
        if simglfo is not None:
            if newfo['gene'] in simglfo['seqs'][region]:  # exact name is in simglfo
                simstr = 'same name'
                if newfo['seq'] == simglfo['seqs'][region][newfo['gene']]:  # I mean if they have the same name they should be identical, but maybe not?
                    simstr += ' and seq'
                    simstr = utils.color('green', simstr)
                else:
                    simstr += utils.color('red', ' different seq')
                simstr += ' as in simulation'
            else:  # see if an equivalent gene is in simglfo
                sim_name, sim_seq = find_equivalent_gene_in_glfo(simglfo, newfo['seq'], glfo[codon + '-positions'][newfo['gene']], new_name=newfo['gene'], glfo_str='sim glfo', debug=True)
                if sim_name is not None:
                    if sim_seq == newfo['seq']:
                        simstr = 'same sequence found in simulation, but with name %s' % utils.color_gene(sim_name)
                    else:
                        simstr = '%s (%s) found in simulation' % (utils.color('green', 'equivalent sequence'), utils.color_gene(sim_name))
                else:
                    simstr = 'doesn\'t seem to correspond to any simulation genes'
            simstr = '(' + simstr + ')'
        print('    adding new allele to glfo: %s' % simstr)
        if 'template-gene' in newfo:
            if newfo['template-gene'] in glfo['seqs'][region]:  # it should be in there, unless it was a removed template gene corresponding to a previous new gene
                aligned_new_seq, aligned_template_seq = utils.color_mutants(glfo['seqs'][region][newfo['template-gene']], newfo['seq'], align=True, return_ref=True)
                print('      template %s   %s' % (aligned_template_seq, utils.color_gene(newfo['template-gene'])))
                print('           new %s   %s' % (aligned_new_seq, utils.color_gene(newfo['gene'])))
            else:
                print('      template %s not in glfo -- this _should_ be because it was a removed template gene a few lines up ^' % utils.color_gene(newfo['template-gene']))
                print('           new %s   %s' % (newfo['seq'], utils.color_gene(newfo['gene'])))
        else:
            print('               %s   %s (no template)' % (newfo['seq'], utils.color_gene(newfo['gene'])))

    if remove_template_genes or ('remove-template-gene' in newfo and newfo['remove-template-gene']):
        # assert newfo['remove-template-gene']  # damnit, i have two ways to specify this and it pisses me off. But too hard to change it just now
        if 'template-gene' not in newfo:
            raise Exception('no template gene specified in newfo')
        print('    %s template gene %s' % (utils.color('red', 'removing'), utils.color_gene(newfo['template-gene'])))
        remove_gene(glfo, newfo['template-gene'])

# ----------------------------------------------------------------------------------------
def remove_the_stupid_godamn_template_genes_all_at_once(glfo, templates_to_remove):
    for gene in templates_to_remove:
        remove_gene(glfo, gene, debug=True)

# ----------------------------------------------------------------------------------------
def generate_new_alleles(glfo, new_allele_info, remove_template_genes=False, debug=False):
    """
    Generate some snp'd genes and add them to glfo, specified with <new_allele_info>.
    e.g. [{'gene' : 'IGHV3-71*01', 'positions' : (35, None)}, ] will add a snp at position 35 and at a random location.
    The resulting snp'd gene will have a name like IGHV3-71*01+C35T.T47G
    """
    templates_to_remove = set()  # NOTE they're not removed if you don't actually add a new gene

    added_names = []
    for genefo in new_allele_info:
        if len(genefo['snp-positions']) == 0 and len(genefo['indel-positions']) == 0:
            continue
        template_gene = genefo['gene']

        print('    generating new allele from %s: ' % utils.color_gene(template_gene), end=' ')
        for mtype in ['snp', 'indel']:
            if len(genefo[mtype + '-positions']) > 0:
                print('%d %s%s' % (len(genefo[mtype + '-positions']), mtype, utils.plural(len(genefo[mtype + '-positions']))), end=' ')
        print('')

        assert utils.get_region(template_gene) == 'v'

        newfo = None
        itry = 0
        while newfo is None or newfo['gene'] in glfo['seqs'][utils.get_region(template_gene)]:
            if itry > 0:
                print('      already in glfo, try again')
                if itry > 99:
                    raise Exception('too many tries while trying to generate new snps -- did you specify a lot of snps on the same position?')
            newfo = generate_single_new_allele(template_gene, glfo['cyst-positions'][template_gene], glfo['seqs'][utils.get_region(template_gene)][template_gene], genefo['snp-positions'], genefo['indel-positions'])
            itry += 1

        if remove_template_genes:
            templates_to_remove.add(template_gene)
        add_new_allele(glfo, newfo, remove_template_genes=False, use_template_for_codon_info=False, debug=debug)  # *don't* remove the templates here, since we don't know if there's another snp later that needs them
        added_names.append(newfo['gene'])

    remove_the_stupid_godamn_template_genes_all_at_once(glfo, templates_to_remove)  # works fine with zero-length <templates_to_remove>

    return added_names  # need the order of the names so we can get allele prevalence freqs from the command line right

# ----------------------------------------------------------------------------------------
def write_glfo(output_dir, glfo, only_genes=None, debug=False):
    if debug:
        print('  writing glfo to %s%s' % (output_dir, '' if only_genes is None else ('  (restricting to %d genes)' % len(only_genes))))
    if os.path.exists(output_dir + '/' + glfo['locus']):
        remove_glfo_files(output_dir, glfo['locus'])  # also removes output_dir
    os.makedirs(output_dir + '/' + glfo['locus'])

    for region in utils.getregions(glfo['locus']):
        with open(get_fname(output_dir, glfo['locus'], region), 'w') as outfile:
            for gene in glfo['seqs'][region]:
                if only_genes is not None and gene not in only_genes:
                    continue
                outfile.write('>' + gene + '\n')
                outfile.write(glfo['seqs'][region][gene] + '\n')

    with open(get_extra_fname(output_dir, glfo['locus']), utils.csv_wmode()) as csvfile:
        writer = csv.DictWriter(csvfile, csv_headers)
        writer.writeheader()
        for region, codon in utils.conserved_codons[glfo['locus']].items():
            for gene, istart in glfo[codon + '-positions'].items():
                if only_genes is not None and gene not in only_genes:
                    continue
                writer.writerow({'gene' : gene, codon + '_position' : istart})

    # make sure there weren't any files lingering in the output dir when we started
    # NOTE this will ignore the dirs corresponding to any *other* loci (which is what we want now, I think)
    unexpected_files = set(glob.glob(output_dir + '/' + glfo['locus'] + '/*')) - set(glfo_fnames(output_dir, glfo['locus']))
    if len(unexpected_files) > 0:
        raise Exception('unexpected file(s) while writing germline set: %s' % (' '.join(unexpected_files)))

# ----------------------------------------------------------------------------------------
def remove_glfo_files(gldir, locus, print_warning=True):
    locusdir = gldir + '/' + locus
    if not os.path.exists(locusdir):
        # print '    %s tried to remove nonexistent glfo dir %s' % (utils.color('yellow', 'warning'), locusdir)  # this seems to only happen when something crashed in the middle, i.e. this code is old/mature enough I'm commenting the check
        return
    if os.path.islink(locusdir):  # linked to original, for backwards compatibility (see read_glfo())
        print('    note: removing link new germline dir name %s' % locusdir)
        os.remove(locusdir)
        if os.path.exists(gldir + '/' + locus[2]):  # presumably the link's target also exists and needs to be removed
            locusdir = gldir + '/' + locus[2]
            print('    note: also removing old germline dir name (i.e. link target) %s' % locusdir)
    for fname in glfo_fnames(gldir, locus):
        if os.path.exists(fname):
            os.remove(fname)
        else:
            if print_warning:
                print('    %s tried to remove non-existent glfo file %s' % (utils.color('yellow', 'warning'), fname))
    os.rmdir(locusdir)
    if len(os.listdir(gldir)) == 0:  # if there aren't any other locus dirs in here, remove the parent dir as well
        os.rmdir(gldir)

# ----------------------------------------------------------------------------------------
def get_alleles_per_gene_weights(n_alleles_per_gene):  # given desired mean alleles per gene, figure out the required probability for 1 and 2 alleles
    if n_alleles_per_gene == 1.:
        return 1., 0.
    elif n_alleles_per_gene == 2.:
        return 0., 1.
    assert n_alleles_per_gene > 1. and n_alleles_per_gene < 2.
    tmpfrac = (n_alleles_per_gene - 2.) / (1. - n_alleles_per_gene)
    a = tmpfrac / (1. + tmpfrac)
    b = 1. - a
    return a, b

# ----------------------------------------------------------------------------------------
def choose_some_alleles(region, genes_to_use, allelic_groups, n_alleles_per_gene, n_max_alleles_per_gene=2, debug=False):
    """ choose a gene (i.e. a primary and sub-version) from <allelic_groups>, and its attendant alleles """
    # NOTE also modifies <allelic_groups>

    if n_alleles_per_gene[region] >= n_max_alleles_per_gene:
        raise Exception('requested mean number of alleles per gene has to be less than the max alleles per gene, but %f >= %f' % (n_alleles_per_gene[region], n_max_alleles_per_gene))

    if len(allelic_groups[region]) == 0:
        raise Exception('ran out of %s alleles (either --n-genes-per-region or --n-sim-alleles-per-gene are probably too big)' % region)  # note that we don't reuse pv/sv pairs (the idea being such a pair represents an actual gene), and we don't directly control how many alleles are chosen from each such pair, so there isn't really a way to make sure you get every single allele in the germline set.
    available_versions = None
    while available_versions is None or len(available_versions) == 0:
        # if available_versions is not None:
        #     print '  %s couldn\'t find any versions that have %d alleles, so trying again' % (utils.color('red', 'warning'), n_alleles)
        n_alleles = numpy.random.choice(list(range(1, n_max_alleles_per_gene + 1)), p=get_alleles_per_gene_weights(n_alleles_per_gene[region]))
        available_versions = [(pv, subv) for pv in allelic_groups[region] for subv in allelic_groups[region][pv] if len(allelic_groups[region][pv][subv]) >= n_alleles]
    ichoice = numpy.random.randint(0, len(available_versions) - 1) if len(available_versions) > 1 else 0  # numpy.random.choice() can't handle list of tuples (and barfs if you give it only one thing to choose from)
    primary_version, sub_version = available_versions[ichoice]
    new_alleles = set(numpy.random.choice(list(allelic_groups[region][primary_version][sub_version]), size=n_alleles, replace=False))
    if debug:
        print('      %8s %5s   %s' % (primary_version, sub_version, ' '.join([utils.color_gene(g, width=15) for g in new_alleles])))

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
    with open(fname, utils.csv_wmode()) as pfile:
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
                print('%14.8f   %s' % (freq, utils.color_gene(gene)))
        assert utils.is_normed(allele_prevalence_freqs[region])
    return allele_prevalence_freqs

# ----------------------------------------------------------------------------------------
def choose_allele_prevalence_freqs(glfo, allele_prevalence_freqs, region, min_allele_prevalence_freq, debug=False):
    n_alleles = len(glfo['seqs'][region])
    prevalence_counts = numpy.random.randint(1, int(1. / min_allele_prevalence_freq), size=n_alleles)  # ensures that any two alleles have a prevalence ratio between <min_allele_prevalence_freq> and 1. NOTE it's inclusive
    prevalence_freqs = [float(c) / sum(prevalence_counts) for c in prevalence_counts]
    allele_prevalence_freqs[region] = {g : f for g, f in zip(list(glfo['seqs'][region].keys()), prevalence_freqs)}
    assert utils.is_normed(allele_prevalence_freqs[region])
    if debug:
        print('    choosing %s allele prevalence freqs:' % region)
        print('      counts %s' % ' '.join([('%5d' % c) for c in prevalence_counts]))
        print('      freqs  %s' % ' '.join([('%5.3f' % c) for c in prevalence_freqs]))
        print('      min ratio: %.3f' % (min(prevalence_freqs) / max(prevalence_freqs)))

# ----------------------------------------------------------------------------------------
def generate_germline_set(glfo, args, new_allele_info=None, dont_remove_template_genes=False, debug=False):
    """ NOTE removes genes from  <glfo> """
    allelic_groups = utils.separate_into_allelic_groups(glfo)  # NOTE by design, these are not the same as the groups created by alleleremover
    allele_prevalence_freqs = {r : {} for r in utils.regions}
    for region in utils.regions:
        genes_to_use = set()
        for _ in range(args.n_genes_per_region[region]):
            choose_some_alleles(region, genes_to_use, allelic_groups, args.n_sim_alleles_per_gene)
        remove_genes(glfo, set(glfo['seqs'][region].keys()) - genes_to_use)  # NOTE would use glutils.restrict_to_genes() but it isn't on a regional basis

        if region == 'v' and new_allele_info is not None:
            # assert len(new_allele_info) <= len(glfo['seqs'][region])  # hmm, not actually sure why I had this here
            template_genes = numpy.random.choice(list(glfo['seqs'][region].keys()), size=len(new_allele_info), replace=True)  # (template) genes in <new_allele_info> should all be None
            for ig in range(len(new_allele_info)):
                new_allele_info[ig]['gene'] = template_genes[ig]
            _ = generate_new_alleles(glfo, new_allele_info, remove_template_genes=not dont_remove_template_genes, debug=debug)

        choose_allele_prevalence_freqs(glfo, allele_prevalence_freqs, region, args.min_sim_allele_prevalence_freq)

    print('  generated germline set with n alleles (n genes):', end=' ')
    allelic_groups = utils.separate_into_allelic_groups(glfo)
    for region in utils.regions:
        n_genes = sum([len(list(allelic_groups[region][pv].keys())) for pv in allelic_groups[region]])
        print(' %s %2d (%d) ' % (region, len(glfo['seqs'][region]), n_genes), end=' ')
    print('')

    utils.separate_into_allelic_groups(glfo, allele_prevalence_freqs=allele_prevalence_freqs, debug=debug)

    write_allele_prevalence_freqs(allele_prevalence_freqs, args.allele_prevalence_fname)  # NOTE lumps all the regions together, unlike in the parameter dirs

# ----------------------------------------------------------------------------------------
def check_allele_prevalence_freqs(outfname, glfo, allele_prevalence_fname, only_region=None, debug=False):
    allele_prevalence_freqs = read_allele_prevalence_freqs(allele_prevalence_fname)
    counts = {r : {g : 0 for g in glfo['seqs'][r]} for r in utils.regions}
    _, annotation_list, _ = utils.read_output(outfname, dont_add_implicit_info=True)
    for line in annotation_list:
        for region in utils.regions:
            counts[region][line[region + '_gene']] += 1
    if debug:
        print('   checking allele prevalence freqs')
    for region in utils.regions:
        if only_region is not None and region != only_region:
            continue
        total = sum(counts[region].values())
        width = str(max(3, len(str(max(counts[region].values())))))
        if debug:
            print('       %s  (%d total)' % (region, total))
            print(('           %' + width + 's    freq    expected') % 'obs')
            for gene in sorted(glfo['seqs'][region], key=lambda g: counts[region][g], reverse=True):
                print(('          %' + width + 'd     %.3f    %.3f   %s') % (counts[region][gene], float(counts[region][gene]) / total, allele_prevalence_freqs[region][gene], utils.color_gene(gene, width=15)))

# ----------------------------------------------------------------------------------------
def choose_new_allele_name(template_gene, new_seq, snpfo=None, indelfo=None):  # may modify <snpfo>, if it's passed
    new_name = template_gene
    hashstr = utils.uidhashstr(new_seq, max_len=5)

    if indelfo is not None and len(indelfo['indels']) > 0:  # call it a new sub-family (er, sub-version, depending on nomenclature)
        new_name = utils.rejoin_gene(utils.get_locus(template_gene), utils.get_region(template_gene), utils.gene_family(template_gene), hashstr, '01')
    elif snpfo is not None:
        if '+' in utils.allele(template_gene) and len(template_gene.split('+')) == 2:  # if template was snpd we need to fix up <snpfo>
            simplified_snpfo = simplify_snpfo(template_gene, snpfo)  # returns None if it failed to parse mutation info in template name
            if simplified_snpfo is not None:
                new_name = template_gene.split('+')[0]  # original template name, before we did any snp'ing
                snpfo = simplified_snpfo
        if '+' not in new_name:  # if there's a plus in <new_name> at this juncture, it means we failed to simplify snpfo, i.e. probably the template gene has a hash name
            if len(snpfo) > 0 and len(snpfo) <= 5:
                new_name += '+' + stringify_mutfo(snpfo)
            else:
                new_name += '+' + hashstr
        else:  # just give up and use a new hash str
            new_name = new_name.split('+')[0] + '+' + hashstr
    else:
        new_name = template_gene + '+' + hashstr

    if new_name.count('+') > 1:
        raise Exception('somehow ended up with too many \'+\'s in %s' % new_name)

    return new_name, snpfo

# ----------------------------------------------------------------------------------------
def find_nearest_gene_in_glfo(glfo, new_seq, new_name=None, exclusion_3p=None, region='v', debug=False):  # NOTE should really be merged with find_nearest_gene_with_same_cpos()
    if new_seq in list(glfo['seqs'][region].values()):
        raise Exception('exact sequence already in glfo')
    seqfos = [{'name' : g, 'seq' : s} for g, s in glfo['seqs'][region].items()]
    seqfos.append({'name' : 'new', 'seq' : new_seq})
    aligned_seqs = {sfo['name'] : sfo['seq'] for sfo in utils.align_many_seqs(seqfos)}
    hdists = [(name, utils.hamming_distance(seq, aligned_seqs['new'])) for name, seq in aligned_seqs.items() if name != 'new']
    hdists = sorted(hdists, key=operator.itemgetter(1))
    if len(hdists) == 0:
        raise Exception('also no nearby genes (should only happen if the gl set is pretty trivial)')
    nearest_gene, nearest_distance = hdists[0]
    n_snps = nearest_distance

    realigned_new_seq, realigned_nearest_seq = utils.align_seqs(aligned_seqs['new'], aligned_seqs[nearest_gene])  # have to re-align 'em in order to get rid of extraneous gaps from other seqs in the previous alignment
    n_indels = utils.count_n_separate_gaps(realigned_new_seq, exclusion_3p=exclusion_3p) + utils.count_n_separate_gaps(realigned_nearest_seq, exclusion_3p=exclusion_3p)

    if debug:
        colored_realigned_new_seq, colored_realigned_nearest_seq = utils.color_mutants(realigned_new_seq, realigned_nearest_seq, return_ref=True)  # have to re-align 'em in order to get rid of extraneous gaps from other seqs in the previous alignment
        print('      %s: nearest is %s (%d snp%s, %d indel%s):' % (utils.color_gene(new_name), utils.color_gene(nearest_gene), n_snps, utils.plural(n_snps), n_indels, utils.plural(n_indels)))
        print('            %s %s' % (colored_realigned_new_seq, utils.color_gene(new_name) if new_name is not None else ''))
        print('            %s %s' % (colored_realigned_nearest_seq, utils.color_gene(nearest_gene)))

    return nearest_gene, n_snps, n_indels, realigned_new_seq, realigned_nearest_seq

# ----------------------------------------------------------------------------------------
def find_nearest_gene_using_names(glfo, gene, debug=False):
    def handle_return(genelist, label):
        if len(genelist) > 0:
            if debug:
                print('  found %d %s%s, returning %s' % (len(genelist), label, utils.plural(len(genelist)), utils.color_gene(genelist[0])))
            return genelist[0]
        else:
            if debug:
                print('  no %s' % label)
            return None

    region = utils.get_region(gene)
    if debug:
        print('looking for nearest gene to %s (%d choices)' % (utils.color_gene(gene), len(glfo['seqs'][region])))
    family = utils.gene_family(gene)
    primary_version, sub_version, allele = utils.split_gene(gene)
    alleles = [g for g in glfo['seqs'][region] if utils.are_alleles(g, gene)]
    ngene = handle_return(alleles, 'allele')
    if ngene is not None:
        return ngene
    family_genes = [g for g in glfo['seqs'][region] if utils.gene_family(g) == family]
    sv_genes = [g for g in family_genes if utils.sub_version(g) == sub_version]
    ngene = handle_return(sv_genes, 'subversion')
    if ngene is not None:
        return ngene
    ngene = handle_return(family_genes, 'family gene')
    if ngene is not None:
        return ngene
    raise Exception('couldn\'t find any similar genes for %s among %s' % (gene, ' '.join(glfo['seqs'][region])))

# ----------------------------------------------------------------------------------------
def find_nearest_gene_with_same_cpos(glfo, new_seq, new_cpos=None, new_name=None, exclusion_5p=0, exclusion_3p=3, glfo_str='glfo', debug=False, debug_only_zero_distance=False):  # NOTE should really be merged with find_nearest_gene_in_glfo()
    region = 'v'
    if new_cpos is None:  # try to guess it...
        if new_name is None:
            raise Exception('have to specify either <new_cpos> or <new_name>')
        if new_name in glfo['seqs'][region]:
            new_cpos = utils.cdn_pos(glfo, region, new_name)
        elif get_template_gene(new_name) in glfo['seqs'][region]:  # that'll fail if it isn't an inferred allele... but then again if it's not an inferred allele hopefully it was in there as itself
            new_cpos = utils.cdn_pos(glfo, region, get_template_gene(new_name))
        else:  # it's not very safe to use a different allele... but it's probably ok (could do the darn alignment thing, but then that's slow)
            other_alleles = [g for g in glfo['seqs'][region] if utils.are_alleles(g, get_template_gene(new_name))]
            if len(other_alleles) > 0:
                new_cpos = utils.cdn_pos(glfo, region, other_alleles[0])
            else:
                other_alleles = [g for g in glfo['seqs'][region] if utils.gene_family(g) == utils.gene_family(new_name) and len(glfo['seqs'][region][g]) == len(new_seq)]
                if len(other_alleles) > 0:
                    print('    arg, giving up and using cpos from %s for %s' % (utils.color_gene(g), utils.color_gene(new_name)))
                else:  # this is a terrible hack, but it hardly ever comes up (I think only when tigger/igdiscover infer a new allele that's really dissimilar to anything in the glfo)
                    glfo['seqs'][region][new_name] = new_seq  # I should really add an option to do the alignment and stuff from get_missing_codon_info() without havving to add it to glfo... but this is ok for now
                    get_missing_codon_info(glfo)
                    new_cpos = utils.cdn_pos(glfo, region, new_name)
                    remove_gene(glfo, new_name)
                    # raise Exception('couldn\'t guess a codon position for %s (glfo has: %s)' % (utils.color_gene(new_name), utils.color_genes(glfo['seqs'][region].keys())))

    min_distance, nearest_gene, nearest_seq = None, None, None
    for oldname_gene, oldname_seq in glfo['seqs'][region].items():  # NOTE <oldname_{gene,seq}> is the old *name* corresponding to the new (snp'd) allele, whereas <old_seq> is the allele from which we inferred the new (snp'd) allele
        oldpos = utils.cdn_pos(glfo, region, oldname_gene)
        if oldpos != new_cpos:
            continue

        # not sure why this happens, but whatever
        if len(oldname_seq[exclusion_5p : oldpos + 3]) != len(new_seq[exclusion_5p : new_cpos + 3]):
            # print '%s' % utils.color('red', 'wtf')
            # print 'oldname_seq[%2d : %2d + 3] = %s' % (exclusion_5p, oldpos, oldname_seq[exclusion_5p : oldpos + 3])
            # print '    new_seq[%2d : %2d + 3] = %s' % (exclusion_5p, new_cpos, new_seq[exclusion_5p : new_cpos + 3])
            continue

        distance = 0

        # snps up through cysteine
        distance += utils.hamming_distance(oldname_seq[exclusion_5p : oldpos + 3], new_seq[exclusion_5p : new_cpos + 3])

        if abs(len(oldname_seq) - len(new_seq)) > exclusion_3p:  # allow differences in length, but only if they're <= the number of 3' excluded bases
            continue

        # distance to right of cysteine
        bases_to_right_of_cysteine = min(len(oldname_seq) - (oldpos + 3), len(new_seq) - exclusion_3p - (new_cpos + 3))  # NOTE this is kind of dumb, it excludes <exclusion_3p> *more* bases, even if we've already excluded <exclusion_3p> bases due to length differences. But, oh, well, it's just equivalent to a somewhat larger exclusion, anyway
        if bases_to_right_of_cysteine > 0:
            distance += utils.hamming_distance(oldname_seq[oldpos + 3 : oldpos + 3 + bases_to_right_of_cysteine], new_seq[new_cpos + 3 : new_cpos + 3 + bases_to_right_of_cysteine])

        if min_distance is None or distance < min_distance:
            min_distance = distance
            nearest_gene = oldname_gene
            nearest_seq = oldname_seq

    if min_distance is None:  # nobody had the same cpos
        return None, None, None

    if debug and (not debug_only_zero_distance or min_distance == 0):
        new_seq_str = new_seq + ' ' * max(0, len(nearest_seq) - len(new_seq))
        nearest_seq_str = nearest_seq + ' ' * max(0, len(new_seq) - len(nearest_seq))
        new_print_strs, nearest_print_strs = [], []
        for istart, istop, color in ((0, exclusion_5p, 'blue'), (exclusion_5p, new_cpos, None), (new_cpos, new_cpos + 3, 'reverse_video'), (new_cpos + 3, new_cpos + 3 + bases_to_right_of_cysteine, None), (new_cpos + 3 + bases_to_right_of_cysteine, 999999, 'blue')):  # arg, I don't like the 999999, but can't figure out a better way
            new_print_strs += [utils.color(color, new_seq_str[istart : istop])]
            if color == 'blue':  # don't color mutated bases in the blue (excluded) parts (tried to fix this in commented line below but it doesn't quite work)
                nearest_print_strs += [nearest_seq_str[istart : istop]]
            else:
                nearest_print_strs += [utils.color_mutants(new_seq_str[istart : istop], nearest_seq_str[istart : istop])]
            nearest_print_strs[-1] = utils.color(color, nearest_print_strs[-1])
            # nearest_print_strs += [utils.color(color, utils.color_mutants(new_seq_str[istart : istop], nearest_seq_str[istart : istop], red_bkg=color=='blue'))]
        print('        %s gene %s with same cpos in %s for %s (blue bases are not considered):' % (utils.color('blue', 'equivalent') if min_distance == 0 else 'nearest', utils.color_gene(nearest_gene), glfo_str, ' ' if new_name is None else utils.color_gene(new_name)))
        print('            %s   %s' % (''.join(new_print_strs), 'new' if new_name is None else utils.color_gene(new_name)))
        print('            %s   %s' % (''.join(nearest_print_strs), utils.color_gene(nearest_gene)))

    return min_distance, nearest_gene, nearest_seq

# ----------------------------------------------------------------------------------------
def find_equivalent_gene_in_glfo(glfo, new_seq, new_cpos=None, new_name=None, exclusion_5p=0, exclusion_3p=3, glfo_str='glfo', debug=False):
    # if <new_seq> likely corresponds to an allele that's already in <glfo>, return that name and its sequence, otherwise return (None, None).
    # NOTE that the calling code, in general, is determining whether we want this sequence in the calling code's glfo.
    # Here, we're trying to find any existing names in <glfo>, which is typically a different glfo (it's usually either default_initial or simulation)
    region = 'v'  # conserved codon stuff below will have to be changed for j

    if new_name is not None and new_name in glfo['seqs'][region]:
        raise Exception('you have to check for new name in glfo before calling this (%s)' % utils.color_gene(new_name))

    # first see if the exact sequence is in there
    if new_seq in list(glfo['seqs'][region].values()):
        names_for_this_seq = [g for g in glfo['seqs'][region] if glfo['seqs'][region][g] == new_seq]
        assert len(names_for_this_seq) == 1  # this should have already been verified in glutils
        new_name = names_for_this_seq[0]
        if debug:
            print('        exact sequence in %s under name %s' % (glfo_str, utils.color_gene(new_name)))
        return new_name, new_seq

    min_distance, nearest_gene, nearest_seq = find_nearest_gene_with_same_cpos(glfo, new_seq, new_cpos=new_cpos, new_name=new_name, exclusion_5p=exclusion_5p, exclusion_3p=exclusion_3p, glfo_str=glfo_str, debug=debug, debug_only_zero_distance=True)
    if min_distance is not None and min_distance == 0:
        return nearest_gene, nearest_seq

    if debug:
        print('        no equivalent gene for %s' % utils.color_gene(new_name) if new_name is not None else '', end=' ')
        _, _, _, _, _ = find_nearest_gene_in_glfo(glfo, new_seq, new_name=new_name, debug=True)

    return None, None

# ----------------------------------------------------------------------------------------
def synchronize_glfos(ref_glfo, new_glfo, region, ref_label='ref glfo', debug=False):
    # note that this modifies not only the names in <new_glfo>, but also the sequences (since we want equivalent genes to have identical sequences after synchronization)
    assert region == 'v'  # cysteine stuff would need to be generalized
    genes_in_common, equivalent_gene_pairs = [], []
    for new_name, new_seq in new_glfo['seqs'][region].items():
        if new_name in ref_glfo['seqs'][region]:
            if new_seq != ref_glfo['seqs'][region][new_name]:
                print('%s different sequences for %s:' % (utils.color('red', 'error'), utils.color_gene(new_name)))  # this should really probably be an exception
                utils.color_mutants(ref_glfo['seqs'][region][new_name], new_seq, align=True, print_result=True, extra_str='  ', ref_label=ref_label + '  ')
            genes_in_common.append(new_name)
            continue
        if debug:
            print('     %s:' % utils.color_gene(new_name))
        equiv_name, equiv_seq = find_equivalent_gene_in_glfo(ref_glfo, new_seq, utils.cdn_pos(new_glfo, region, new_name), new_name=new_name, glfo_str=ref_label, debug=debug)
        if equiv_name is not None:
            if equiv_name in new_glfo['seqs'][region]:
                if debug:
                    print('        (already in new glfo [probably not a snpd allele, i.e. you\'ve got two alleles in your gl set that are equivalent])')
                continue

            if debug:
                print('        %s --> %s' % (utils.color_gene(new_name), utils.color_gene(equiv_name)))
            equivalent_gene_pairs.append((new_name, equiv_name))
            remove_gene(new_glfo, new_name)
            add_new_allele(new_glfo, {'gene' : equiv_name, 'seq' : equiv_seq, 'cpos' : utils.cdn_pos(ref_glfo, region, equiv_name)}, use_template_for_codon_info=False)

    if debug:
        print('  %d genes in common' % len(genes_in_common))
        print('  %d (%d) different genes in %s (new)' % (len(ref_glfo['seqs'][region]) - len(genes_in_common) - len(equivalent_gene_pairs), len(new_glfo['seqs'][region]) - len(genes_in_common) - len(equivalent_gene_pairs), ref_label))
        print('  %d equivalent genes: ' % len(equivalent_gene_pairs))
        gene_str_width = max([utils.len_excluding_colors(utils.color_gene(g)) for g, _ in equivalent_gene_pairs])  # only care about the first one
        def novelstr(g):
            return utils.color('blue', 'x') if is_novel(g) else ' '
        for name, ename in equivalent_gene_pairs:
            print('            %s %s  %s %s' % (novelstr(name), utils.color_gene(name, width=gene_str_width), novelstr(ename), utils.color_gene(ename, width=gene_str_width)))

    # try to convert igdiscover novel allele names (randomish string at end) to partis notation (e.g. +A78C)
    assert 'igdiscover' not in ref_label # partis has to be the <ref_label> for this to work
    if debug:
        print('  attempting to convert igdiscover new-allele names')
    for new_name, new_seq in new_glfo['seqs'][region].items():
        if not is_novel(new_name):
            continue
        method, _, _ = split_inferred_allele_name(new_name)
        if method != 'igdiscover':
            continue
        if debug:
            print('     %s:' % utils.color_gene(new_name))
        better_name, better_snpfo, nearest_gene = try_to_get_name_and_mutfo_from_seq(new_name, new_seq, ref_glfo, consider_indels=False, debug=debug)
        if better_name is not None:
            remove_gene(new_glfo, new_name)
            newfo = {'gene' : better_name, 'seq' : new_seq}
            if nearest_gene in new_glfo['seqs'][region]:
                newfo['template-gene'] = nearest_gene
                use_template_for_codon_info = True
            elif is_novel(nearest_gene) and get_template_gene(nearest_gene) in new_glfo['seqs'][region]:
                newfo['template-gene'] = get_template_gene(nearest_gene)
                use_template_for_codon_info = True
            else:
                # newfo['cpos'] = nearest_gene
                use_template_for_codon_info = False
                pass  # dammit, I guess just let it do the alignment
            add_new_allele(new_glfo, newfo, use_template_for_codon_info=use_template_for_codon_info)
            if debug:
                print('        %s --> %s' % (utils.color_gene(new_name), utils.color_gene(better_name)))

    return new_glfo

# ----------------------------------------------------------------------------------------
def create_glfo_from_fasta(fastafname, locus, region, template_germline_dir, simulation_germline_dir=None, debug=False):
    glfo = read_glfo(template_germline_dir, locus)
    simglfo = None
    if simulation_germline_dir is not None:
        simglfo = read_glfo(simulation_germline_dir, locus)
    fasta_alleles = set()
    for seqfo in utils.read_fastx(fastafname):
        new_name, aligned_seq = seqfo['name'], seqfo['seq']  # well, tigger or igdiscover
        unaligned_seq = utils.remove_gaps(aligned_seq)

        if new_name in glfo['seqs'][region]:
            if glfo['seqs'][region][new_name] != unaligned_seq:
                print('%s different sequences in template glfo and fasta output for %s:\n    %s\n    %s' % (utils.color('red', 'error'), new_name, glfo['seqs'][region][new_name], aligned_seq))
            fasta_alleles.add(new_name)
            continue

        equiv_name, equiv_seq = find_equivalent_gene_in_glfo(glfo, unaligned_seq, new_cpos=None, new_name=new_name, debug=debug)
        if equiv_name is not None:
            new_name = equiv_name
            unaligned_seq = equiv_seq  # this has no effect, but it's a nice reminder that <unaligned_seq> isn't correct any more
            fasta_alleles.add(new_name)
        else:
            print('  %s allele %s' % (utils.color('red', 'new'), utils.color_gene(new_name)))
            method, template_gene, mutstrs = split_inferred_allele_name(new_name, debug=debug)  # will fail (and we want it to) if it's not an inferred allele
            snpfo = None
            if mutstrs is not None:
                snpfo = try_to_get_mutfo_from_name(new_name, aligned_seq=aligned_seq, debug=debug)
            if snpfo is not None:
                new_name, snpfo = choose_new_allele_name(template_gene, unaligned_seq, snpfo=snpfo)  # reminder: can modify <snpfo>
            newfo = {'gene' : new_name, 'seq' : unaligned_seq, 'template-gene' : template_gene}
            add_new_allele(glfo, newfo, use_template_for_codon_info=True, simglfo=simglfo, debug=True)
            fasta_alleles.add(new_name)

    # remove alleles that *aren't* in fasta's gl set
    for gene in glfo['seqs'][region]:  # can't do it before, since we want to use existing ones to get codon info
        if gene not in fasta_alleles:
            remove_gene(glfo, gene)

    return glfo

# ----------------------------------------------------------------------------------------
# NOTE this is pretty similar to AlleleRemover, but it's different enough that it doesn't make sense to merge them (yes, I tried)
def remove_extra_alleles_per_gene(glfo, n_max_alleles_per_gene, gene_obs_counts, only_regions=None, debug=False):
    debug = True
    if only_regions is None:
        only_regions = utils.regions
    allelic_groups = utils.separate_into_allelic_groups(glfo)
    for region in only_regions:
        if debug:
            print('  %s' % utils.color('green', region))
        for pv in allelic_groups[region]:
            for sv, allele_list in allelic_groups[region][pv].items():
                # if len(allele_list) <= n_max_alleles_per_gene:  # nothing to do
                #     continue
                sorted_alleles = sorted(allele_list, key=lambda q: gene_obs_counts[region][q], reverse=True)
                rm_dbg_str = ''
                if len(allele_list) > n_max_alleles_per_gene:
                    remove_genes(glfo, sorted_alleles[n_max_alleles_per_gene:])  #, extrastr='      ', debug=debug)
                    rm_dbg_str = '   removed: %s' % ' '.join([utils.color_gene(g) for g in sorted_alleles[n_max_alleles_per_gene:]])
                if debug:
                    print('    %s   %8s     %s%s' % (utils.color('purple', pv + ('-' + sv if sv is not None else ''), width=6), ' '.join([('%.0f' % gene_obs_counts[region][g]) for g in  sorted_alleles]), '   '.join([utils.color_gene(g) for g in sorted_alleles]), rm_dbg_str))
