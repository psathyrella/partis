""" A few utility functions. At the moment simply functions used in recombinator which do not
require member variables. """

import sys
import os
import ast
import random
import ast
import re
import math
import glob
from collections import OrderedDict
import itertools
import csv
from subprocess import check_output, CalledProcessError
# from sklearn.metrics.cluster import adjusted_mutual_info_score
import sklearn.metrics.cluster
import numpy
import multiprocessing
import shutil
import copy

from opener import opener

from Bio import SeqIO

#----------------------------------------------------------------------------------------
# NOTE I also have an eps defined in hmmwriter. Simplicity is the hobgoblin of... no, wait, that's just plain ol' stupid to have two <eps>s defined
eps = 1.0e-10  # if things that should be 1.0 are this close to 1.0, blithely keep on keepin on. kinda arbitrary, but works for the moment
def is_normed(probs, this_eps=eps):
    if hasattr(probs, 'keys'):  # if it's a dict, call yourself with a list of the dict's values
        return is_normed([val for key, val in probs.items()])
    elif hasattr(probs, '__getitem__'):  # if it's a list call yourself with their sum
        return is_normed(sum(probs))
    else:  # and if it's a float actually do what you're supposed to do
        return math.fabs(probs - 1.0) < this_eps

# ----------------------------------------------------------------------------------------
def pass_fcn(val):  # dummy function for conversions (see beloww)
    return val

# ----------------------------------------------------------------------------------------
def get_arg_list(arg, intify=False, floatify=False, translation=None, list_of_pairs=False):  # make lists from args that are passed as strings of colon-separated values
    if arg is None:
        return None

    convert_fcn = pass_fcn
    if intify:
        convert_fcn = int
    elif floatify:
        convert_fcn = float

    arglist = arg.strip().split(':')  # to allow ids with minus signs, need to add a space, which you then have to strip() off
    if list_of_pairs:
        arglist = [pairstr.split(',') for pairstr in arglist]
        arglist = [[convert_fcn(p) for p in pair] for pair in arglist]
    else:
        arglist = [convert_fcn(x) for x in arglist]

    if translation is not None:
        for ia in range(len(arglist)):
            if arglist[ia] in translation:
                arglist[ia] = translation[arglist[ia]]

    return arglist

# # ----------------------------------------------------------------------------------------
# hackey_default_gene_versions = {'v':'IGHV3-23*04', 'd':'IGHD3-10*01', 'j':'IGHJ4*02_F'}
# ----------------------------------------------------------------------------------------
regions = ['v', 'd', 'j']
real_erosions = ['v_3p', 'd_5p', 'd_3p', 'j_5p']
# NOTE since we now handle v_5p and j_3p deletions by padding with Ns, the hmm does *not* allow actual v_5p and j_3p deletions.
# This means that while we write parameters for v_5p and j_3p deletions to the parameter dir, these are *not* used in making the
# hmm yamels -- which is what we want, because we want to be able to read in short data reads but make full-length simulation.
effective_erosions = ['v_5p', 'j_3p']
boundaries = ['vd', 'dj']
effective_boundaries = ['fv', 'jf']
humans = ['A', 'B', 'C']
nukes = ['A', 'C', 'G', 'T']
ambiguous_bases = ['N', ]
gap_chars = ['.', '-']
naivities = ['M', 'N']
conserved_codon_names = {'v':'cyst', 'd':'', 'j':'tryp'}
# Infrastrucure to allow hashing all the columns together into a dict key.
# Uses a tuple with the variables that are used to index selection frequencies
# NOTE fv and jf insertions are *effective* (not real) insertions between v or j and the framework. They allow query sequences that extend beyond the v or j regions
index_columns = ('v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'v_5p_del', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del', 'j_3p_del', 'fv_insertion', 'vd_insertion', 'dj_insertion', 'jf_insertion')
# not_used_for_simulation = ('fv_insertion', 'jf_insertion', 'v_5p_del')
index_keys = {}
for i in range(len(index_columns)):  # dict so we can access them by name instead of by index number
    index_keys[index_columns[i]] = i

# ----------------------------------------------------------------------------------------
# Info specifying which parameters are assumed to correlate with which others. Taken from mutual
# information plot in bcellap repo

# key is parameter of interest, and associated list gives the parameters (other than itself) which are necessary to predict it
# TODO am I even using these any more? I think I might just be using the all-probs file
column_dependencies = {}
column_dependencies['v_gene'] = [] # NOTE v choice actually depends on everything... but not super strongly, so a.t.m. I ignore it
column_dependencies['v_5p_del'] = ['v_gene']
column_dependencies['v_3p_del'] = ['v_gene']
column_dependencies['d_gene'] = []
column_dependencies['d_5p_del'] = ['d_gene']  # NOTE according to the mebcell mutual information plot d_5p_del is also correlated to d_3p_del (but we have no way to model that a.t.m. in the hmm)
column_dependencies['d_3p_del'] = ['d_gene']  # NOTE according to the mebcell mutual information plot d_3p_del is also correlated to d_5p_del (but we have no way to model that a.t.m. in the hmm)
column_dependencies['j_gene'] = []
column_dependencies['j_5p_del'] = ['j_gene']  # NOTE mebcell plot showed this correlation as being small, but I'm adding it here for (a possibly foolish) consistency
column_dependencies['j_3p_del'] = ['j_gene']
column_dependencies['fv_insertion'] = []
column_dependencies['vd_insertion'] = ['d_gene']
column_dependencies['dj_insertion'] = ['j_gene']
column_dependencies['jf_insertion'] = []

# column_dependencies['vd_insertion_content'] = []
# column_dependencies['dj_insertion_content'] = []

# definitions here: http://clip.med.yale.edu/changeo/manuals/Change-O_Data_Format.pdf
presto_headers = {
    'unique_id' : 'SEQUENCE_ID',
    'v_gene' : 'V_CALL',
    'd_gene' : 'D_CALL',
    'j_gene' : 'J_CALL',
    'aligned_v_seq' : 'SEQUENCE_IMGT',  # TODO is this supposed to be the whole sequence, or just the v bit?
    'cdr3_length' : 'JUNCTION_LENGTH'
}

# test partition
# reco_info = {'a' : {'reco_id' : '1'}, 'b' : {'reco_id' : '0'}, 'c' : {'reco_id' : '1'}, 'd' : {'reco_id' : '2'}}
# partition = [['a'], ['b', 'c', 'd']]
# true_partition = [['b'], ['a', 'c'], ['d']]

# ----------------------------------------------------------------------------------------
column_configs = {
    'ints' : ('nth_best', 'v_5p_del', 'd_5p_del', 'cdr3_length', 'j_5p_del', 'j_3p_del', 'd_3p_del', 'v_3p_del'),
    'floats' : ('logprob'),
    'literals' : ('indels'),
    'lists' : ('unique_ids', 'seqs', 'aligned_seqs', 'aligned_v_seqs', 'aligned_d_seqs', 'aligned_j_seqs')
}

# ----------------------------------------------------------------------------------------
def convert_to_presto(glfo, line):
    """ convert <line> to presto csv format """
    if len(line['unique_ids']) > 1:
        print line['unique_ids']
        raise Exception('multiple seqs not handled in convert_to_presto')

    if len(line['indelfo']['indels']) > 0:
        raise Exception('passing indel info to presto requires some more thought')
    else:
        del line['indelfo']

    single_info = synthesize_single_seq_line(glfo, line, iseq=0)

    presto_line = {}
    for head, phead in presto_headers.items():
        presto_line[phead] = single_info[head]

    return presto_line

# these are the top 10 v and d genes, and top six js, from mishmash.csv. Restricting to these should make testing much more stable and much faster.
test_only_genes = 'IGHV4-61*08:IGHV3-48*01:IGHV5-51*02:IGHV3-69-1*02:IGHV1/OR15-1*04:IGHV3-66*03:IGHV3-23D*01:IGHV3-71*03:IGHV1-2*04:IGHV1-2*02:IGHD3-16*02:IGHD2-2*03:IGHD2-8*01:IGHD3-22*01:IGHD6-13*01:IGHD4-17*01:IGHD6-19*01:IGHD3-10*01:IGHD2-15*01:IGHD2-21*02:IGHJ5*02:IGHJ3*02:IGHJ2*01:IGHJ1*01:IGHJ6*03:IGHJ4*02'

# tuples with the column and its dependencies mashed together
# (first entry is the column of interest, and it depends upon the following entries)
column_dependency_tuples = []
for column, deps in column_dependencies.iteritems():
    tmp_list = [column]
    tmp_list.extend(deps)
    column_dependency_tuples.append(tuple(tmp_list))

def get_parameter_fname(column=None, deps=None, column_and_deps=None):
    """ return the file name in which we store the information for <column>. Either pass in <column> and <deps> *or* <column_and_deps> """
    if column == 'all':
        return 'all-probs.csv'
    if column_and_deps is None:
        assert column is not None and deps is not None
    if column_and_deps == None:
        column_and_deps = [column]
        column_and_deps.extend(deps)
    outfname = 'probs.csv'
    for ic in column_and_deps:
        outfname = ic + '-' + outfname
    return outfname

# ----------------------------------------------------------------------------------------
def rewrite_germline_fasta(input_dir, output_dir, only_genes):
    """ rewrite the germline set files in <input_dir> to <output_dir>, only keeping the genes in <only_genes> """
    print 'rewriting germlines from %s to %s (using %d genes)' % (input_dir, output_dir, len(only_genes))
    glfo = read_germline_set(input_dir)
    input_germlines = glfo['seqs']
    input_aligned_v_genes = glfo['aligned-v-genes']
    expected_files = []  # list of files that we write here -- if anything else is in the output dir, we barf

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    def write_gl_file(fname, region, igls):
        expected_files.append(fname)
        with open(fname, 'w') as outfile:
            for gene in igls[region]:
                if gene not in only_genes:
                    continue
                outfile.write('>' + gene + '\n')
                outfile.write(igls[region][gene] + '\n')

    for region in regions:
        write_gl_file(output_dir + '/igh' + region + '.fasta', region, input_germlines)
    write_gl_file(output_dir + '/ighv-aligned.fasta', 'v', input_aligned_v_genes)

    for fname in ['cyst-positions.csv', 'tryp-positions.csv']:
        expected_files.append(output_dir + '/' + fname)
        shutil.copyfile(input_dir + '/' + fname, output_dir + '/' + fname)

    # make sure there weren't any files lingering in the output dir when we started
    final_file_list = glob.glob(output_dir + '/*')
    for fname in final_file_list:
        if fname not in expected_files:
            raise Exception('unexpected file %s (expected %s)' % (fname, ' '.join(expected_files)))

    return expected_files  # return the rewritten files so they can be deleted if desired

# ----------------------------------------------------------------------------------------
def from_same_event(reco_info, query_names):
    if len(query_names) > 1:
        # reco_id = reco_info[query_names[0]]['reco_id']  # the first one's reco id
        # for iq in range(1, len(query_names)):  # then loop through the rest of 'em to see if they're all the same
        #     if reco_id != reco_info[query_names[iq]]['reco_id']:
        #         return False
        # return True

        # darn it, this isn't any faster
        reco_ids = [reco_info[query]['reco_id'] for query in query_names]
        return reco_ids.count(reco_ids[0]) == len(reco_ids)  # they're clonal if every reco id is the same as the first one

    else:  # one or zero sequences
        return True

# ----------------------------------------------------------------------------------------
# bash color codes
Colors = {}
Colors['head'] = '\033[95m'
Colors['bold'] = '\033[1m'
Colors['purple'] = '\033[95m'
Colors['blue'] = '\033[94m'
Colors['light_blue'] = '\033[1;34m'
Colors['green'] = '\033[92m'
Colors['yellow'] = '\033[93m'
Colors['red'] = '\033[91m'
Colors['end'] = '\033[0m'

def color(col, seq):
    assert col in Colors
    return Colors[col] + seq + Colors['end']

# ----------------------------------------------------------------------------------------
def color_chars(chars, col, seq):
    return_str = ''
    for ch in seq:
        if ch in chars:
            return_str += color(col, ch)
        else:
            return_str += ch
    return return_str

# ----------------------------------------------------------------------------------------
def color_mutants(ref_seq, seq, print_result=False, extra_str='', ref_label='', post_str=''):
    assert len(ref_seq) == len(seq)
    return_str = ''
    for inuke in range(len(seq)):
        if inuke >= len(ref_seq) or seq[inuke] == ref_seq[inuke]:
            return_str += seq[inuke]
        else:
            return_str += color('red', seq[inuke])
    if print_result:
        print '%s%s%s' % (extra_str, ref_label, ref_seq)
        print '%s%s%s%s' % (extra_str, ' '*len(ref_label), return_str, post_str)
    return return_str

# ----------------------------------------------------------------------------------------
def summarize_gene_name(gene):
    region = get_region(gene)
    primary_version, sub_version, allele = split_gene(gene)
    return ' '.join([region, primary_version, sub_version, allele])

# ----------------------------------------------------------------------------------------
def color_gene(gene):
    """ color gene name (and remove extra characters), eg IGHV3-h*01 --> v3-h1 """
    region = get_region(gene)
    primary_version, sub_version, allele = split_gene(gene)

    return_str =  color('red', region) + color('purple', primary_version)
    if region != 'j':
        return_str += '-' + color('purple', sub_version)
    return_str += color('yellow', allele)
    return return_str

#----------------------------------------------------------------------------------------
def int_to_nucleotide(number):
    """ Convert between (0,1,2,3) and (A,C,G,T) """
    if number == 0:
        return 'A'
    elif number == 1:
        return 'C'
    elif number == 2:
        return 'G'
    elif number == 3:
        return 'T'
    else:
        print 'ERROR nucleotide number not in [0,3]'
        sys.exit()

# ----------------------------------------------------------------------------------------                    
def check_conserved_cysteine(seq, cyst_position, debug=False, extra_str='', assert_on_fail=True):
    """ Ensure there's a cysteine at <cyst_position> in <seq>. """
    if len(seq) < cyst_position+3:
        if debug:
            print '%sseq not long enough in cysteine checker %d %s' % (extra_str, cyst_position, seq)
        if assert_on_fail:
            assert False
        else:
            return False
    cyst_word = str(seq[cyst_position:cyst_position+3])
    if cyst_word != 'TGT' and cyst_word != 'TGC':
        if debug:
            print '%scysteine in v is messed up: %s (%s %d)' % (extra_str, cyst_word, seq, cyst_position)
        if assert_on_fail:
            assert False
        else:
            return False

    return True

# ----------------------------------------------------------------------------------------
def check_conserved_tryptophan(seq, tryp_position, debug=False, extra_str='', assert_on_fail=True):
    """ Ensure there's a tryptophan at <tryp_position> in <seq>. """
    if len(seq) < tryp_position+3:
        if debug:
            print '%sseq not long enough in tryp checker %d %s' % (extra_str, tryp_position, seq)
        if assert_on_fail:
            assert False
        else:
            return False
    tryp_word = str(seq[tryp_position:tryp_position+3])
    if tryp_word != 'TGG':
        if debug:
            print '%stryptophan in j is messed up: %s (%s %d)' % (extra_str, tryp_word, seq, tryp_position)
        if assert_on_fail:
            assert False
        else:
            return False

    return True

# ----------------------------------------------------------------------------------------
def check_both_conserved_codons(seq, cyst_position, tryp_position, debug=False, extra_str='', assert_on_fail=True):
    """ Double check that we conserved the cysteine and the tryptophan. """
    cyst_ok = check_conserved_cysteine(seq, cyst_position, debug, extra_str=extra_str, assert_on_fail=assert_on_fail)
    tryp_ok = check_conserved_tryptophan(seq, tryp_position, debug, extra_str=extra_str, assert_on_fail=assert_on_fail)
    return cyst_ok and tryp_ok

# ----------------------------------------------------------------------------------------
def are_conserved_codons_screwed_up(reco_event):
    """ Version that checks all the final seqs in reco_event.

    Returns True if codons are screwed up, or if no sequences have been added.
    """
    if len(reco_event.final_seqs) == 0:
        return True
    for seq in reco_event.final_seqs:
        codons_ok = check_both_conserved_codons(seq, reco_event.final_cyst_position, reco_event.final_tryp_position, assert_on_fail=False)
        if not codons_ok:
            return True

    return False

#----------------------------------------------------------------------------------------
def stop_codon_check(seq, cyst_position, debug=False):
    """ 
    Make sure there is no in-frame stop codon, where frame is inferred from <cyst_position>.
    Returns True if no stop codon is found
    """
    if cyst_position >= len(seq):
        if debug:
            print '      not sure if there\'s a stop codon (invalid cysteine position)'
        return False  # not sure if there is one, since we have to way to establish the frame
    # jump leftward in steps of three until we reach the start of the sequence
    ipos = cyst_position
    while ipos > 2:
        ipos -= 3
    # ipos should now bet the index of the start of the first complete codon
    while ipos + 2 < len(seq):  # then jump forward in steps of three bases making sure none of them are stop codons
        codon = seq[ipos:ipos+3]
        if codon == 'TAG' or codon == 'TAA' or codon == 'TGA':  # found a stop codon
            if debug:
                print '      stop codon %s at %d in %s' % (codon, ipos, seq)
            return False
        ipos += 3

    return True  # no stop codon

#----------------------------------------------------------------------------------------
def is_position_protected(protected_positions, prospective_position):
    """ Would a mutation at <prospective_position> screw up a protected codon? """
    for position in protected_positions:
        if (prospective_position == position or
            prospective_position == (position + 1) or
            prospective_position == (position + 2)):
            return True
    return False

#----------------------------------------------------------------------------------------
def would_erode_conserved_codon(reco_event):
    """ Would any of the erosion <lengths> delete a conserved codon? """
    lengths = reco_event.erosions
    # check conserved cysteine
    if len(reco_event.seqs['v']) - lengths['v_3p'] <= reco_event.cyst_position + 2:
        print '      about to erode cysteine (%d), try again' % lengths['v_3p']
        return True  # i.e. it *would* screw it up
    # check conserved tryptophan
    if lengths['j_5p'] - 1 >= reco_event.tryp_position:
        print '      about to erode tryptophan (%d), try again' % lengths['j_5p']
        return True

    return False  # *whew*, it won't erode either of 'em

#----------------------------------------------------------------------------------------
def is_erosion_longer_than_seq(reco_event):
    """ Are any of the proposed erosion <lengths> longer than the seq to be eroded? """
    lengths = reco_event.erosions
    if lengths['v_3p'] > len(reco_event.seqs['v']):  # NOTE not actually possible since we already know we didn't erode the cysteine
        print '      v_3p erosion too long (%d)' % lengths['v_3p']
        return True
    if lengths['d_5p'] + lengths['d_3p'] > len(reco_event.seqs['d']):
        print '      d erosions too long (%d)' % (lengths['d_5p'] + lengths['d_3p'])
        return True
    if lengths['j_5p'] > len(reco_event.seqs['j']):  # NOTE also not possible for the same reason
        print '      j_5p erosion too long (%d)' % lengths['j_5p']
        return True
    return False

#----------------------------------------------------------------------------------------
def find_tryp_in_joined_seq(gl_tryp_position_in_j, v_seq, vd_insertion, d_seq, dj_insertion, j_seq, j_erosion, debug=False):
    """ Find the <end> tryptophan in a joined sequence.

    Given local tryptophan position in the j region, figure
    out what position it's at in the final sequence.
    NOTE gl_tryp_position_in_j is the position *before* the j was eroded,
    but this fcn assumes that the j *has* been eroded.
    also NOTE <[vdj]_seq> are assumed to already be eroded
    """
    if debug:
        print 'checking tryp with: %s, %d - %d = %d' % (j_seq, gl_tryp_position_in_j, j_erosion, gl_tryp_position_in_j - j_erosion)
    check_conserved_tryptophan(j_seq, gl_tryp_position_in_j - j_erosion)  # make sure tryp is where it's supposed to be
    length_to_left_of_j = len(v_seq + vd_insertion + d_seq + dj_insertion)
    if debug:
        print '  finding tryp position as'
        print '    length_to_left_of_j = len(v_seq + vd_insertion + d_seq + dj_insertion) = %d + %d + %d + %d' % (len(v_seq), len(vd_insertion), len(d_seq), len(dj_insertion))
        print '    result = gl_tryp_position_in_j - j_erosion + length_to_left_of_j = %d - %d + %d = %d' % (gl_tryp_position_in_j, j_erosion, length_to_left_of_j, gl_tryp_position_in_j - j_erosion + length_to_left_of_j)
    return gl_tryp_position_in_j - j_erosion + length_to_left_of_j

# ----------------------------------------------------------------------------------------
def is_mutated(original, final, n_muted=-1, n_total=-1, also_allow=['*', ]):
    alphabet = nukes + ambiguous_bases + also_allow
    if original not in alphabet or final not in alphabet:
        raise Exception('bad base (%s or %s) in utils.is_mutated()' % (original, final))

    # print '%s %s' % (original, final)
    return_str = final
    if original in ambiguous_bases or final in ambiguous_bases:  # don't count Ns in the total
        return return_str, n_muted, n_total
    # if original in also_allow:
    #     return color('light_blue', return_str), n_muted, n_total
        
    n_total += 1
    if original != final:
        return_str = color('red', final)
        n_muted += 1
    return return_str, n_muted, n_total

# ----------------------------------------------------------------------------------------
def get_v_5p_del(original_seqs, line):
    # deprecated
    assert False  # this method will no longer work when I need to get v left *and* j right 'deletions'
    original_length = len(original_seqs['v']) + len(original_seqs['d']) + len(original_seqs['j'])
    total_deletion_length = int(line['v_3p_del']) + int(line['d_5p_del']) + int(line['d_3p_del']) + int(line['j_5p_del'])
    total_insertion_length = len(line['vd_insertion']) + len(line['dj_insertion'])
    return original_length - total_deletion_length + total_insertion_length - len(line['seq'])

# ----------------------------------------------------------------------------------------
def get_reco_event_seqs(germlines, line, original_seqs, lengths, eroded_seqs):
    """
    get original and eroded germline seqs
    NOTE does not modify line
    """
    for region in regions:
        del_5p = int(line[region + '_5p_del'])
        del_3p = int(line[region + '_3p_del'])
        original_seqs[region] = germlines[region][line[region + '_gene']]
        lengths[region] = len(original_seqs[region]) - del_5p - del_3p
        # assert del_5p + lengths[region] != 0
        eroded_seqs[region] = original_seqs[region][del_5p : del_5p + lengths[region]]

# ----------------------------------------------------------------------------------------
def add_cdr3_info(glfo, line, debug=False):
    """
    Add the cyst_position, tryp_position, and cdr3_length to <line> based on the information already in <line>.
    If info is already there, make sure it's the same as what we calculate here
    """

    original_seqs = {}  # original (non-eroded) germline seqs
    lengths = {}  # length of each match (including erosion)
    eroded_seqs = {}  # eroded germline seqs
    get_reco_event_seqs(glfo['seqs'], line, original_seqs, lengths, eroded_seqs)

    # if len(line['fv_insertion']) > 0:
    #     print line
    # if len(line['jf_insertion']) > 0:
    #     print line

    # NOTE see get_conserved_codon_position -- they do similar things, but start from different information
    eroded_gl_cpos = glfo['cyst-positions'][line['v_gene']]  - int(line['v_5p_del']) + len(line['fv_insertion'])  # cysteine position in eroded germline sequence. EDIT darn, actually you *don't* want to subtract off the v left deletion, because that (deleted) base is presumably still present in the query sequence
    # if debug:
    #     print '  cysteine: cpos - v_5p_del + fv_insertion = %d - %d + %d = %d' % (glfo['cyst-positions'][line['v_gene']], int(line['v_5p_del']), len(line['fv_insertion']), eroded_gl_cpos)
    eroded_gl_tpos = glfo['tryp-positions'][line['j_gene']] - int(line['j_5p_del'])
    values = {}
    values['cyst_position'] = eroded_gl_cpos
    tpos_in_joined_seq = eroded_gl_tpos + len(line['fv_insertion']) + len(eroded_seqs['v']) + len(line['vd_insertion']) + len(eroded_seqs['d']) + len(line['dj_insertion'])
    values['tryp_position'] = tpos_in_joined_seq
    values['cdr3_length'] = tpos_in_joined_seq - eroded_gl_cpos + 3

    for key, val in values.items():
        if key in line:  # if <key> was previously added to <line> make sure we got the same value this time
            if int(line[key]) != int(val):
                # raise Exception('previously calculated value of ' + key + ' (' + str(line[key]) + ') not equal to new value (' + str(val) + ') for ' + str(line['unique_id']))
                print 'ERROR previously calculated value of ' + key + ' (' + str(line[key]) + ') not equal to new value (' + str(val) + ') for ' + str(line['unique_id'])
        else:
            line[key] = val

    cyst_ok = check_conserved_cysteine(line['fv_insertion'] + eroded_seqs['v'], eroded_gl_cpos, debug=debug, extra_str='      ', assert_on_fail=False)
    tryp_ok = check_conserved_tryptophan(eroded_seqs['j'], eroded_gl_tpos, debug=debug, extra_str='      ', assert_on_fail=False)
    if not cyst_ok or not tryp_ok:
        if debug:
            print '    bad codon[s] (%s %s) in %s' % ('cyst' if not cyst_ok else '', 'tryp' if not tryp_ok else '', ':'.join(line['unique_ids']) if 'unique_ids' in line else line)

    line['naive_seq'] = get_full_naive_seq(glfo['seqs'], line)

# ----------------------------------------------------------------------------------------
def disambiguate_effective_insertions(bound, line, seq, unique_id, debug=False):
    # These are kinda weird names, but the distinction is important
    # If an insert state with "germline" N emits one of [ACGT], then the hmm will report this as an inserted N. Which is what we want -- we view this as a germline N which "mutated" to [ACGT].
    # This concept of insertion germline state is mostly relevant for simultaneous inference on several sequences, i.e. in ham we don't want to just say the inserted base was the base in the query sequence.
    # But here, we're trimming off the effective insertions and we have to treat the inserted germline N to [ACGT] differently than we would an insertion which was "germline" [ACGT] which emitted an N,
    # and also differently to a real germline [VDJ] state that emitted an N.
    naive_insertion = line[bound + '_insertion']  # reminder: ham gets this from the last character in the insertion state name, e.g. 'insert_left_A' or 'insert_right_N'
    if bound == 'fv':  # ...but to accomodate multiple sequences, the insert states can emit non-germline states, so the mature bases might be different.
        mature_insertion = seq[ : len(line[bound + '_insertion'])]
    elif bound == 'jf':
        if len(line[bound + '_insertion']) > 0:
            mature_insertion = seq[-len(line[bound + '_insertion']) : ]
        else:
            mature_insertion = ''
    else:
        assert False
    if naive_insertion == mature_insertion:  # all is simple and hunky-dory: no insertion 'mutations'
        final_insertion = ''  # leave this bit as an insertion in the final <line>
        insertion_to_remove = naive_insertion  # this bit we'll remove -- it's just Ns (note that this is only equal to the N padding if we correctly inferred the right edge of the J [for jf bound])
    else:
        if len(naive_insertion) != len(mature_insertion):
            raise Exception('naive and mature insertions not the same length\n   %s\n   %s\n' % (naive_insertion, mature_insertion))
        assert naive_insertion.count('N') == len(naive_insertion)  # should be e.g. naive: NNN   mature: ANN
        if bound == 'fv':  # ...but to accomodate multiple sequences, the insert states can emit non-germline states, so the mature bases might be different.
            i_first_non_N = find_first_non_ambiguous_base(mature_insertion)
            final_insertion = mature_insertion[i_first_non_N : ]
            insertion_to_remove = mature_insertion[ : i_first_non_N]
        elif bound == 'jf':
            i_first_N = find_last_non_ambiguous_base_plus_one(mature_insertion)
            final_insertion = mature_insertion[ : i_first_N]
            insertion_to_remove = mature_insertion[i_first_N : ]
        else:
            assert False
        if debug:
            print 'naive and mature %s insertions differ for %s' % (bound, unique_id)
            color_mutants(naive_insertion, mature_insertion, print_result=True, extra_str='          ')
            print '   removing %s and leaving %s' % (insertion_to_remove, final_insertion)

    # remove the insertion that we want to remove
    if bound == 'fv':  # ...but to accomodate multiple sequences, the insert states can emit non-germline states, so the mature bases might be different.
        trimmed_seq = seq[len(insertion_to_remove) : ]
    elif bound == 'jf':
        if len(insertion_to_remove) > 0:
            trimmed_seq = seq[ : -len(insertion_to_remove)]
        else:
            trimmed_seq = seq
    if debug:
        print '    %s insertion   final %s   to_remove %s    trimmed_seq %s' % (bound, final_insertion, insertion_to_remove, trimmed_seq)

    return trimmed_seq, final_insertion, insertion_to_remove

# ----------------------------------------------------------------------------------------
def reset_effective_erosions_and_effective_insertions(line, debug=False):
    """ 
    Ham does not allow (well, no longer allows) v_5p and j_3p deletions -- we instead pad sequences with Ns.
    This means that the info we get from ham always has these effective erosions set to zero, but for downstream
    things we sometimes want to know where the reads stopped (e.g. if we want to mimic them in simulation).
    Note that these effective erosion values will be present in the parameter dir, but are *not* incorporated into
    the hmm yaml files.
    """

    assert line['v_5p_del'] == 0  # just to be safe
    assert line['j_3p_del'] == 0

    if debug:
        print 'resetting effective erosions'
        print '     %s' % line['seqs'][0]

    # first remove effective (fv and jf) insertions
    trimmed_seqs = []
    final_insertions, insertions_to_remove = [], []
    for iseq in range(len(line['seqs'])):
        trimmed_seq = line['seqs'][iseq]
        final_insertions.append({})
        insertions_to_remove.append({})
        for bound in effective_boundaries:
            trimmed_seq, final_insertion, insertion_to_remove = disambiguate_effective_insertions(bound, line, trimmed_seq, line['unique_ids'][iseq], debug)
            final_insertions[-1][bound] = final_insertion
            insertions_to_remove[-1][bound] = insertion_to_remove
        trimmed_seqs.append(trimmed_seq)

    # arbitrarily use the zeroth sequence
    if len(trimmed_seqs) > 1:
        print 'TODO don\'t just use the zeroth sequence'
    trimmed_seq = trimmed_seqs[0]  # TODO right now I'm setting these to the same values for the entire clonal family, but at some point we should allow different sequences to have different read lengths/start positions
    final_fv_insertion = final_insertions[0]['fv']
    final_jf_insertion = final_insertions[0]['jf']
    fv_insertion_to_remove = insertions_to_remove[0]['fv']
    jf_insertion_to_remove = insertions_to_remove[0]['jf']
    line['v_5p_del'] = find_first_non_ambiguous_base(trimmed_seq)
    line['j_3p_del'] = len(trimmed_seq) - find_last_non_ambiguous_base_plus_one(trimmed_seq)
    line['cyst_position'] -= line['v_5p_del'] + len(fv_insertion_to_remove)
    line['tryp_position'] -= line['v_5p_del'] + len(fv_insertion_to_remove)

    for iseq in range(len(line['seqs'])): # TODO note, though, that this trims *all* the seqs according to the read truncation from the zeroth sequence
        line['seqs'][iseq] = trimmed_seqs[iseq][line['v_5p_del'] : ]
        if line['j_3p_del'] > 0:
            line['seqs'][iseq] = line['seqs'][iseq][ : -line['j_3p_del']]

    if 'naive_seq' in line:
        line['naive_seq'] = line['naive_seq'][len(fv_insertion_to_remove) : len(line['naive_seq']) - len(jf_insertion_to_remove)]
        line['naive_seq'] = line['naive_seq'][line['v_5p_del'] : len(line['naive_seq']) - line['j_3p_del']]
        if len(line['naive_seq']) != len(line['seqs'][0]):
            raise Exception('didn\'t trim naive seq to proper length:\n  %s\n  %s' % (line['naive_seq'], line['seqs'][0]))
        # color_mutants(line['naive_seq'], line['seqs'][0], print_result=True)

    if debug:
        print '     fv %d   v_5p %d   j_3p %d   jf %d    %s' % (len(fv_insertion_to_remove), line['v_5p_del'], line['j_3p_del'], len(jf_insertion_to_remove), line['seqs'][0])

    line['fv_insertion'] = final_fv_insertion
    line['jf_insertion'] = final_jf_insertion

# ----------------------------------------------------------------------------------------
def get_full_naive_seq(germlines, line):
    for erosion in real_erosions + effective_erosions:
        if line[erosion + '_del'] < 0:  # wow, I have no idea why I thought this would happen
            print 'ERROR %s less than zero %d' % (erosion, line[erosion + '_del'])
        assert line[erosion + '_del'] >= 0
    original_seqs = {}  # original (non-eroded) germline seqs
    lengths = {}  # length of each match (including erosion)
    eroded_seqs = {}  # eroded germline seqs
    get_reco_event_seqs(germlines, line, original_seqs, lengths, eroded_seqs)
    return line['fv_insertion'] + eroded_seqs['v'] + line['vd_insertion'] + eroded_seqs['d'] + line['dj_insertion'] + eroded_seqs['j'] + line['jf_insertion']

# ----------------------------------------------------------------------------------------
def is_this_a_validate_event(germlines, line):
    # this is just a copy of stuff in get_regional_naive_seq_bounds, but so we can run this first and see if it fails and then skip the event

    original_seqs = {}  # original (non-eroded) germline seqs
    lengths = {}  # length of each match (including erosion)
    eroded_seqs = {}  # eroded germline seqs
    get_reco_event_seqs(germlines, line, original_seqs, lengths, eroded_seqs)

    start, end = {}, {}
    start['v'] = int(line['v_5p_del'])
    end['v'] = start['v'] + len(line['fv_insertion'] + eroded_seqs['v'])  # base just after the end of v
    start['d'] = end['v'] + len(line['vd_insertion'])
    end['d'] = start['d'] + len(eroded_seqs['d'])
    start['j'] = end['d'] + len(line['dj_insertion'])
    end['j'] = start['j'] + len(eroded_seqs['j'] + line['jf_insertion'])

    for tmpreg in regions:  # subtract_unphysical_erosions
        start[tmpreg] -= int(line['v_5p_del'])
        end[tmpreg] -= int(line['v_5p_del'])

    try:
        for chkreg in regions:
            assert start[chkreg] >= 0
            assert end[chkreg] >= 0
            assert end[chkreg] >= start[chkreg]
            assert end[chkreg] <= len(line['seq'])
        assert end['j'] == len(line['seq'])
        return False
    except:
        return True

# ----------------------------------------------------------------------------------------
def get_regional_naive_seq_bounds(return_reg, germlines, line):
    # NOTE it's kind of a matter of taste whether unphysical deletions (v left and j right) should be included in the 'naive sequence'.
    # Unless <subtract_unphysical_erosions>, here we assume the naive sequence has *no* unphysical deletions

    subtract_unphysical_erosions = True

    original_seqs = {}  # original (non-eroded) germline seqs
    lengths = {}  # length of each match (including erosion)
    eroded_seqs = {}  # eroded germline seqs
    get_reco_event_seqs(germlines, line, original_seqs, lengths, eroded_seqs)

    start, end = {}, {}
    start['v'] = int(line['v_5p_del'])
    end['v'] = start['v'] + len(line['fv_insertion'] + eroded_seqs['v'])  # base just after the end of v
    start['d'] = end['v'] + len(line['vd_insertion'])
    end['d'] = start['d'] + len(eroded_seqs['d'])
    start['j'] = end['d'] + len(line['dj_insertion'])
    end['j'] = start['j'] + len(eroded_seqs['j'] + line['jf_insertion'])

    if subtract_unphysical_erosions:
        for tmpreg in regions:
            start[tmpreg] -= int(line['v_5p_del'])
            end[tmpreg] -= int(line['v_5p_del'])
        # end['j'] -= line['j_3p_del']  # ARG.ARG.ARG

    def elegantishfail():
        for k, v in line.items():
            print '%30s %s' % (k, v)
        raise Exception('end of j %d not equal to sequence length %d in %s (or maybe something else is wrong...)' % (end['j'], len(line['seq']), line['unique_id']))

    try:
        for chkreg in regions:
            assert start[chkreg] >= 0
            assert end[chkreg] >= 0
            assert end[chkreg] >= start[chkreg]
            assert end[chkreg] <= len(line['seq'])
        assert end['j'] == len(line['seq'])
    except:
        elegantishfail()

    return (start[return_reg], end[return_reg])

# ----------------------------------------------------------------------------------------
def add_match_info(glfo, line, debug=False):
    """
    add to <line> the query match seqs (sections of the query sequence that are matched to germline) and their corresponding germline matches.

    """

    original_seqs = {}  # original (non-eroded) germline seqs
    lengths = {}  # length of each match (including erosion)
    eroded_seqs = {}  # eroded germline seqs
    get_reco_event_seqs(glfo['seqs'], line, original_seqs, lengths, eroded_seqs)
    add_cdr3_info(glfo, line, debug=debug)  # add cyst and tryp positions, and cdr3 length

    # add the <eroded_seqs> to <line> so we can find them later
    for region in regions:
        line[region + '_gl_seq'] = eroded_seqs[region]

    # the sections of the query sequence which are assigned to each region
    v_start = len(line['fv_insertion'])
    d_start = v_start + len(eroded_seqs['v']) + len(line['vd_insertion'])
    j_start = d_start + len(eroded_seqs['d']) + len(line['dj_insertion'])
    line['v_qr_seq'] = line['seq'][v_start : v_start + len(eroded_seqs['v'])]
    line['d_qr_seq'] = line['seq'][d_start : d_start + len(eroded_seqs['d'])]
    line['j_qr_seq'] = line['seq'][j_start : j_start + len(eroded_seqs['j'])]

    # if line['cdr3_length'] == -1:
    #     print '      ERROR %s failed to add match info' % ' '.join([str(i) for i in ids])

# ----------------------------------------------------------------------------------------
def print_reco_event(germlines, line, one_line=False, extra_str='', return_string=False, label='', indelfo=None, indelfos=None):
    """ if <line> describes a whole event, with multiple sequences, it'll have <unique_ids> and <seqs> fields; otherwise it'll have <unique_id> and <seq> """
    event_str_list = []
    if 'unique_ids' in line:
        for iseq in range(len(line['unique_ids'])):
            tmpline = copy.deepcopy(line)
            del tmpline['unique_ids']  # not that they'd cause any harm, but may as well be tidy
            del tmpline['seqs']
            tmpline['seq'] = line['seqs'][iseq]
            this_extra_str = ''
            if indelfos[iseq] is not None:  # for now, just print the reversed seq, i.e. the seq with the indels undone
                # if indelfos[iseq]['indels'] is not None:
                #     for ii in range(len(indelfos[iseq]['indels'])):
                #         idl = indelfos[iseq]['indels'][ii]
                #         if ii > 0:
                #             this_extra_str += '\nxxx'
                #         this_extra_str += ' %10s: %2d bases at %3d'  % (idl['type'], idl['len'], idl['pos'])
                # # if 'indels' not in extra_str:
                # #     extra_str += color('yellow', 'indels')
                if indelfos[iseq] is not None and len(indelfos[iseq]['indels']) > 0:
                    tmpline['seq'] = indelfos[iseq]['reversed_seq']
                    if find_first_non_ambiguous_base(tmpline['seq']) > 0:
                        tmpline['seq'] = tmpline['seq'][find_first_non_ambiguous_base(tmpline['seq']) : ]
                    if find_last_non_ambiguous_base_plus_one(tmpline['seq']) < len(tmpline['seq']):
                        tmpline['seq'] = tmpline['seq'][ : find_last_non_ambiguous_base_plus_one(tmpline['seq'])]
                    # if len(tmpline['fv_insertion']) + tmpline['v_5p_del'] == find_first_non_ambiguous_base(tmpline['seq']):  # if we hack in effective erosions (for the hmm) based on Ns at either end, we have to hack that info into the indel info as well and it doesn't happen until here
                    #     tmpline['seq'] = tmpline['seq'][len(tmpline['fv_insertion']) + tmpline['v_5p_del'] : ]
                    # if len(tmpline['jf_insertion']) + tmpline['j_3p_del'] == len(tmpline['seq']) - find_last_non_ambiguous_base_plus_one(tmpline['seq']):  # if we hack in effective erosions (for the hmm) based on Ns at either end, we have to hack that info into the indel info as well and it doesn't happen until here
                    #     if tmpline['j_3p_del'] > 0:
                    #         tmpline['seq'] = tmpline['seq'][ : -tmpline['j_3p_del']]
            else:
                pass
                # this_extra_str = ' %10s  %2s          %3s'  % ('', '', '')
            event_str = print_seq_in_reco_event(germlines, tmpline, extra_str=extra_str + this_extra_str, return_string=return_string,
                                                      label=(label if iseq==0 else ''),
                                                      one_line=(iseq>0),
                                                      indelfo=None)
            event_str_list.append(event_str)
    else:
        tmpline = dict(line)
        if indelfo is not None and indelfo['reversed_seq'] != '':  # and len(indelfo['indels']) > 1:  # TODO allow printing with more than one indel
            tmpline['seq'] = indelfo['reversed_seq']
            indelfo = None
            if 'indels' not in extra_str:
                extra_str += color('yellow', 'indels')
        event_str = print_seq_in_reco_event(germlines, tmpline, extra_str=extra_str, return_string=return_string, label=label, one_line=one_line, indelfo=indelfo)
        event_str_list.append(event_str)

    if return_string:
        return ''.join(event_str_list)

# ----------------------------------------------------------------------------------------
def print_seq_in_reco_event(germlines, line, extra_str='', return_string=False, label='', indelfo=None, one_line=False):
    """ Print ascii summary of recombination event and mutation.

    If <one_line>, then only print out the final_seq line.
    """
    reverse_indels = True  # for inferred sequences, we want to un-reverse the indels that we previously reversed in smith-waterman
    if 'indels' in line:
        # raise Exception('')
        pass  # need to implement printing for multiple indels and multiple sequences
        # if len(line['indels']) == 0:
        #     indelfo = None
        # else:
        #     assert indelfo is None  # don't want indel info from two places
        #     indelfo = line['indels']  #{'reversed_seq' : None, 'indels' : ast.literal_eval(line['indels'])}
        #     reverse_indels = False  # ...whereas for simulation, we indels were not reverse, so we just want to color insertions

    v_5p_del = int(line['v_5p_del'])
    v_3p_del = int(line['v_3p_del'])
    d_5p_del = int(line['d_5p_del'])
    d_3p_del = int(line['d_3p_del'])
    j_5p_del = int(line['j_5p_del'])
    j_3p_del = int(line['j_3p_del'])

    original_seqs = {}  # original (non-eroded) germline seqs
    lengths = {}  # length of each match (including erosion)
    eroded_seqs = {}  # eroded germline seqs
    get_reco_event_seqs(germlines, line, original_seqs, lengths, eroded_seqs)

    if indelfo is not None:
        if len(indelfo['indels']) == 0:  # TODO make this less hackey
            indelfo = None
        else:
            add_indels_to_germline_strings(line, indelfo, original_seqs, lengths, eroded_seqs, reverse_indels)

    # build up the query sequence line, including colors for mutations and conserved codons
    final_seq = ''
    n_muted, n_total = 0, 0
    j_right_extra = ''  # portion of query sequence to right of end of the j match
    n_inserted = 0
    for inuke in range(len(line['seq'])):
        if indelfo is not None:
            lastfo = indelfo['indels'][-1]  # if the "last" (arbitrary but necessary ordering) indel starts here
              # if we're at the position that the insertion started at (before we removed it)
        if indelfo is not None and lastfo['type'] == 'insertion':
            if reverse_indels and inuke == lastfo['pos']:
                final_seq += lastfo['seqstr']  # put the insertion back into the query sequence
                n_inserted += len(lastfo['seqstr'])
            elif not reverse_indels and inuke - lastfo['pos'] >= 0 and inuke - lastfo['pos'] < lastfo['len']:
                final_seq += line['seq'][inuke]
                n_inserted += 1
                continue
            # if inuke > lastfo['pos'] + lastfo['len']:
            #     inuke -= lastfo['len']
        if indelfo is not None and lastfo['type'] == 'deletion':
            if reverse_indels and inuke - lastfo['pos'] >= 0 and inuke - lastfo['pos'] < lastfo['len']:  # if we're within the bases that we added to make up for the deletionlen
                final_seq += color('light_blue', '*')
                # j_right_extra += ' '
                continue
            elif not reverse_indels and inuke == lastfo['pos']:
                final_seq += lastfo['len'] * color('light_blue', '*')
                n_inserted = - lastfo['len']
                # j_right_extra += ' ' * lastfo['len']

        new_nuke = ''  # now not necessarily a single base
        ilocal = inuke
        if indelfo is not None and reverse_indels:
            ilocal += n_inserted
        if indelfo is not None and not reverse_indels and lastfo['type'] == 'deletion':
            ilocal -= n_inserted
        if ilocal < len(line['fv_insertion']):  # haven't got to start of v match yet, so just add on the query seq nuke
            new_nuke, n_muted, n_total = line['seq'][inuke], n_muted, n_total + 1
        else:
            ilocal -= len(line['fv_insertion'])
            if ilocal < lengths['v']:
                new_nuke, n_muted, n_total = is_mutated(eroded_seqs['v'][ilocal], line['seq'][inuke], n_muted, n_total)
            else:
                ilocal -= lengths['v']
                if ilocal < len(line['vd_insertion']):
                    new_nuke, n_muted, n_total = is_mutated(line['vd_insertion'][ilocal], line['seq'][inuke], n_muted, n_total)
                else:
                    ilocal -= len(line['vd_insertion'])
                    if ilocal < lengths['d']:
                        new_nuke, n_muted, n_total = is_mutated(eroded_seqs['d'][ilocal], line['seq'][inuke], n_muted, n_total)
                    else:
                        ilocal -= lengths['d']
                        if ilocal < len(line['dj_insertion']):
                            new_nuke, n_muted, n_total = is_mutated(line['dj_insertion'][ilocal], line['seq'][inuke], n_muted, n_total)
                        else:
                            ilocal -= len(line['dj_insertion'])
                            if ilocal < lengths['j']:
                                new_nuke, n_muted, n_total = is_mutated(eroded_seqs['j'][ilocal], line['seq'][inuke], n_muted, n_total)
                            else:
                                new_nuke, n_muted, n_total = line['seq'][inuke], n_muted, n_total + 1
                                j_right_extra += ' '

        if 'cyst_position' in line and 'tryp_position' in line:
            for pos in (line['cyst_position'], line['tryp_position']):  # reverse video for the conserved codon positions
                if indelfo is not None and not reverse_indels:
                    pos += n_inserted
                if inuke >= pos and inuke < pos + 3:
                    new_nuke = '\033[7m' + new_nuke + '\033[m'

        final_seq += new_nuke

    # check if there isn't enough space for dots in the vj line
    no_space = False
    interior_length = len(line['vd_insertion']) + len(eroded_seqs['d']) + len(line['dj_insertion'])  # length of the portion of the vj line that is normally taken up by dots (and spaces)
    if v_3p_del + j_5p_del > interior_length:
        no_space = True

    if no_space:
        v_3p_del_str = '.' + str(v_3p_del) + '.'
        j_5p_del_str = '.' + str(j_5p_del) + '.'
        extra_space_because_of_fixed_nospace = max(0, interior_length - len(v_3p_del_str + j_5p_del_str))
        if len(v_3p_del_str + j_5p_del_str) <= interior_length:  # ok, we've got space now
            no_space = False
    else:
        v_3p_del_str = '.'*v_3p_del
        j_5p_del_str = '.'*j_5p_del
        extra_space_because_of_fixed_nospace = 0

    eroded_seqs_dots = {}
    eroded_seqs_dots['v'] = eroded_seqs['v'] + v_3p_del_str
    eroded_seqs_dots['d'] = '.'*d_5p_del + eroded_seqs['d'] + '.'*d_3p_del
    eroded_seqs_dots['j'] = j_5p_del_str + eroded_seqs['j'] + '.'*j_3p_del

    v_5p_del_str = '.'*v_5p_del
    if v_5p_del > 50:
        v_5p_del_str = '...' + str(v_5p_del) + '...'

    insert_line = ' '*len(line['fv_insertion']) + ' '*lengths['v']
    insert_line += ' '*len(v_5p_del_str)
    insert_line += line['vd_insertion']
    insert_line += ' ' * lengths['d']
    insert_line += line['dj_insertion']
    insert_line += ' ' * lengths['j']
    insert_line += j_right_extra
    insert_line += ' ' * j_3p_del  # no damn idea why these need to be commented out for some cases in the igblast parser...
    # insert_line += ' '*len(line['jf_insertion'])

    germline_d_start = len(line['fv_insertion']) + lengths['v'] + len(line['vd_insertion']) - d_5p_del
    germline_d_end = germline_d_start + len(original_seqs['d'])
    d_line = ' ' * germline_d_start
    d_line += ' '*len(v_5p_del_str)
    d_line += eroded_seqs_dots['d']
    d_line += ' ' * (len(eroded_seqs['j']) + len(line['dj_insertion']) - d_3p_del)
    d_line += j_right_extra
    d_line += ' ' * j_3p_del
    # d_line += ' '*len(line['jf_insertion'])

    vj_line = ' ' * len(line['fv_insertion'])
    vj_line += v_5p_del_str
    vj_line += eroded_seqs_dots['v'] + '.'*extra_space_because_of_fixed_nospace
    germline_v_end = len(line['fv_insertion']) + len(eroded_seqs['v']) + v_3p_del - 1  # position in the query sequence at which we find the last base of the v match. NOTE we subtract off the v_5p_del because we're *not* adding dots for that deletion (it's just too long)
    germline_j_start = germline_d_end + 1 - d_3p_del + len(line['dj_insertion']) - j_5p_del
    vj_line += ' ' * (germline_j_start - germline_v_end - 2)
    vj_line += eroded_seqs_dots['j']
    vj_line += j_right_extra
    # vj_line += ' '*len(line['jf_insertion'])

    # if len(insert_line) != len(d_line) or len(insert_line) != len(vj_line):
    #     # print '\nERROR lines unequal lengths in event printer -- insertions %d d %d vj %d' % (len(insert_line), len(d_line), len(vj_line)),
    #     # assert no_space
    #     if not no_space:
    #         print 'ERROR no space'
    #     # print ' ...but we\'re out of space so it\'s expected'

    insert_line = color_chars(ambiguous_bases + ['*', ], 'light_blue', insert_line)
    d_line = color_chars(ambiguous_bases + ['*', ], 'light_blue', d_line)
    vj_line = color_chars(ambiguous_bases + ['*', ], 'light_blue', vj_line)

    out_str_list = []
    # insert, d, and vj lines
    if not one_line:
        out_str_list.append('%s    %s   inserts\n' % (extra_str, insert_line))
        if label != '':
            out_str_list[-1] = extra_str + label + out_str_list[-1][len(extra_str + label) :]
        out_str_list.append('%s    %s   %s\n' % (extra_str, d_line, color_gene(line['d_gene'])))
        out_str_list.append('%s    %s   %s %s\n' % (extra_str, vj_line, color_gene(line['v_gene']), color_gene(line['j_gene'])))

    # then query sequence
    v_5p_del_space_str = ' '*len(v_5p_del_str)
    j_3p_del_space_str = ' ' * j_3p_del
    if 'chops' in line:
        assert False  # deprecated
        if line['chops']['left'] > 0:
            v_5p_del_space_str = v_5p_del_space_str[ : -line['chops']['left']]
            final_seq = color('green', '.'*line['chops']['left']) + final_seq
        if line['chops']['right'] > 0:
            j_3p_del_space_str = j_3p_del_space_str[line['chops']['right'] : ]
            final_seq = final_seq + color('green', '.'*line['chops']['right'])
    final_seq = v_5p_del_space_str + final_seq + j_3p_del_space_str
    final_seq = color_chars(ambiguous_bases + ['*', ], 'light_blue', final_seq)

    out_str_list.append('%s    %s' % (extra_str, final_seq))
    # and finally some extra info
    out_str_list.append('   %4.2f mut' % (0. if n_total == 0. else float(n_muted) / n_total))
    # out_str_list.append(' %4.0f%%' % (100 * float(n_muted) / n_total))
    # if 'logprob' in line:
    #     out_str_list.append('  score: %s' % line['logprob'])
    # if 'cdr3_length' in line:
    #     out_str_list.append('   cdr3: %d' % int(line['cdr3_length']))
    out_str_list.append('\n')

    # if print_width is not None:
    #     a_str_list, b_str_list = [], []
    #     for line in out_str_list:
    #         a_str_list.append(line[ : print_width] + '\n')
    #         b_str_list.append(line[print_width : ])
    #     out_str_list = a_str_list + b_str_list

    if return_string:
        return ''.join(out_str_list)
    else:
        print ''.join(out_str_list),

    assert '.' not in line['seq']  # make sure I'm no longer altering line['seq']
    assert ' ' not in line['seq']
    assert '[' not in line['seq']
    assert ']' not in line['seq']

#----------------------------------------------------------------------------------------
def sanitize_name(name):
    """ Replace characters in gene names that make crappy filenames. """
    saniname = name.replace('*', '_star_')
    saniname = saniname.replace('/', '_slash_')
    return saniname

#----------------------------------------------------------------------------------------
def unsanitize_name(name):
    """ Re-replace characters in gene names that make crappy filenames. """
    unsaniname = name.replace('_star_', '*')
    unsaniname = unsaniname.replace('_slash_', '/')
    return unsaniname

#----------------------------------------------------------------------------------------
def read_germline_set(datadir):
    glfo = {}
    glfo['seqs'] = read_germline_seqs(datadir)
    glfo['aligned-v-genes'] = read_germline_seqs(datadir, only_region='v', aligned=True)
    for codon in ['cyst', 'tryp']:
        glfo[codon + '-positions'] = read_codon_positions(datadir + '/' + codon + '-positions.csv')
    return glfo

#----------------------------------------------------------------------------------------
def read_germline_seqs(datadir, only_region=None, aligned=False):
    glseqs = {}
    for region in regions:
        if only_region is not None and region != only_region:
            continue
        fname = datadir + '/igh' + region + '.fasta'
        if aligned:
            fname = fname.replace('.fasta', '-aligned.fasta')
        glseqs[region] = OrderedDict()
        for seq_record in SeqIO.parse(fname, 'fasta'):
            gene_name = seq_record.name
            seq_str = str(seq_record.seq).upper()
            glseqs[region][gene_name] = seq_str
    return glseqs

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

# ----------------------------------------------------------------------------------------
def get_region(gene):
    """ return v, d, or j of gene"""
    region = gene[3:4].lower()
    if 'IGH' not in gene or region not in regions:
        raise Exception('faulty gene name %s' % gene)
    return region

# ----------------------------------------------------------------------------------------
def maturity_to_naivety(maturity):
    if maturity == 'memory':
        return 'M'
    elif maturity == 'naive':
        return 'N'
    else:
        assert False

# ----------------------------------------------------------------------------------------
def are_alleles(gene1, gene2):
    return primary_version(gene1) == primary_version(gene2) and sub_version(gene1) == sub_version(gene2)

# ----------------------------------------------------------------------------------------
def split_gene(gene):
    """ returns (primary version, sub version, allele) """
    # make sure IGH[VDJ] is at the start, and - and * are as expected
    if gene[:4] != 'IGH' + get_region(gene).upper():
        raise Exception('unexpected string in gene name %s' % gene)
    if '*' not in gene:
        raise Exception('no \'*\' found in %s' % gene)

    if get_region(gene) == 'j':  # js don't have sub versions
        primary_version = gene[4 : gene.find('*')]  # the bit between the IGH[VDJ] and the star
        sub_version = None
        allele = gene[gene.find('*') + 1 : ]  # the bit after the star
        if gene != 'IGH' + get_region(gene).upper() + primary_version + '*' + allele:
            raise Exception('could build gene name %s from %s %s' % (gene, primary_version, allele))
    else:
        if '-' not in gene:
            raise Exception('no \'-\' found in %s' % gene)
        if gene.find('-') >= gene.find('*'):
            raise Exception('found \'*\' before \'-\' in %s' % gene)
        primary_version = gene[4 : gene.find('-')]  # the bit between the IGH[VDJ] and the first dash (sometimes there's a second dash as well)
        sub_version = gene[gene.find('-') + 1 : gene.find('*')]  # the bit between the first dash and the star
        allele = gene[gene.find('*') + 1 : ]  # the bit after the star
        if gene != 'IGH' + get_region(gene).upper() + primary_version + '-' + sub_version + '*' + allele:
            raise Exception('could build gene name %s from %s %s %s' % (gene, primary_version, sub_version, allele))

    return primary_version, sub_version, allele

# ----------------------------------------------------------------------------------------
def primary_version(gene):
    return split_gene(gene)[0]

# ----------------------------------------------------------------------------------------
def sub_version(gene):
    return split_gene(gene)[1]

# ----------------------------------------------------------------------------------------
def allele(gene):
    return split_gene(gene)[2]

# ----------------------------------------------------------------------------------------
def are_same_primary_version(gene1, gene2):
    """
    Return true if the bit up to the dash is the same.
    """
    if get_region(gene1) != get_region(gene2):
        return False
    if primary_version(gene1) != primary_version(gene2):
        return False
    return True

# ----------------------------------------------------------------------------------------
def read_overall_gene_probs(indir, only_gene='', normalize=True):
    """
    Return the observed counts/probabilities of choosing each gene version.
    If <normalize> then return probabilities
    If <only_gene> is specified, just return the prob/count for that gene
    """
    counts = { region:{} for region in regions }
    probs = { region:{} for region in regions }
    for region in regions:
        total = 0
        with opener('r')(indir + '/' + region + '_gene-probs.csv') as infile:  # NOTE note this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
            reader = csv.DictReader(infile)
            for line in reader:
                line_count = int(line['count'])
                gene = line[region + '_gene']
                total += line_count
                if gene not in counts[region]:
                    counts[region][gene] = 0
                counts[region][gene] += line_count
        if total < 1:
            assert total == 0
            print 'ERROR zero counts in %s' % indir + '/' + region + '_gene-probs.csv'
            assert False
        for gene in counts[region]:
            probs[region][gene] = float(counts[region][gene]) / total

    if only_gene not in counts[get_region(only_gene)]:
        print '      WARNING %s not found in overall gene probs, returning zero' % only_gene
        if normalize:
            return 0.0
        else:
            return 0

    if only_gene == '':
        if normalize:
            return probs
        else:
            return counts
    else:
        if normalize:
            return probs[get_region(only_gene)][only_gene]
        else:
            return counts[get_region(only_gene)][only_gene]

# ----------------------------------------------------------------------------------------
def find_replacement_genes(indir, min_counts, gene_name=None, single_gene=False, debug=False, all_from_region=''):
    if gene_name != None:
        assert all_from_region == ''
        region = get_region(gene_name)
    else:
        assert all_from_region in regions
        assert single_gene == False
        assert min_counts == -1
        region = all_from_region
    lists = OrderedDict()  # we want to try alleles first, then primary versions, then everything and it's mother
    lists['allele'] = []  # list of genes that are alleles of <gene_name>
    lists['primary_version'] = []  # same primary version as <gene_name>
    lists['all'] = []  # give up and return everything
    with opener('r')(indir + '/' + region + '_gene-probs.csv') as infile:  # NOTE note this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
        reader = csv.DictReader(infile)
        for line in reader:
            gene = line[region + '_gene']
            count = int(line['count'])
            vals = {'gene':gene, 'count':count}
            if all_from_region == '':
                if are_alleles(gene, gene_name):
                    lists['allele'].append(vals)
                if are_same_primary_version(gene, gene_name):
                    lists['primary_version'].append(vals)
            lists['all'].append(vals)

    if single_gene:
        for list_type in lists:
            # return the first which has at least <min_counts> counts
            lists[list_type].sort(reverse=True, key=lambda vals: vals['count'])  # sort by score
            for vals in lists[list_type]:
                if vals['count'] >= min_counts:
                    if debug:
                        print '    return replacement %s %s' % (list_type, color_gene(vals['gene']))
                    return vals['gene']

        print 'ERROR didn\'t find any genes with at least %d for %s in %s' % (min_counts, gene_name, indir)
        assert False
    else:
        # return the whole list NOTE we're including here <gene_name>
        if all_from_region != '':
            return [vals['gene'] for vals in lists['all']]
        for list_type in lists:
            total_counts = sum([vals['count'] for vals in lists[list_type]])
            if total_counts >= min_counts:
                return_list = [vals['gene'] for vals in lists[list_type]]
                if debug:
                    print '      returning all %s for %s (%d genes, %d total counts)' % (list_type + 's', gene_name, len(return_list), total_counts)
                return return_list
            else:
                if debug:
                    print '      not enough counts in %s' % (list_type + 's')

        print 'ERROR couldn\'t find genes for %s in %s' % (gene_name, indir)
        assert False
    
    # print '    \nWARNING return default gene %s \'cause I couldn\'t find anything remotely resembling %s' % (color_gene(hackey_default_gene_versions[region]), color_gene(gene_name))
    # return hackey_default_gene_versions[region]

# ----------------------------------------------------------------------------------------
def hamming_fraction(seq1, seq2, return_len_excluding_ambig=False, extra_bases=None):
    assert len(seq1) == len(seq2)
    if len(seq1) == 0:
        if return_len_excluding_ambig:
            return 0., 0
        else:
            return 0.

    alphabet = nukes + ambiguous_bases

    if extra_bases is not None:
        alphabet += extra_bases

    for ch in seq1 + seq2:  # check that all characters are in the expected alphabet
        if ch not in alphabet:
            raise Exception('unexpected character \'%s\' not among %s in hamming_fraction() with input:\n  %s\n  %s' % (ch, alphabet, seq1, seq2))

    distance, len_excluding_ambig = 0, 0
    for ch1, ch2 in zip(seq1, seq2):
        if ch1 in ambiguous_bases or ch2 in ambiguous_bases:
            continue

        len_excluding_ambig += 1
        if ch1 != ch2:
            distance += 1

    fraction = 0.
    if len_excluding_ambig > 0:
        fraction = distance / float(len_excluding_ambig)
    if return_len_excluding_ambig:
        return fraction, len_excluding_ambig
    else:
        return fraction

# ----------------------------------------------------------------------------------------
def get_key(names):
    """
    Return a hashable combination of the two query names that's the same if we reverse their order.
    """
    return '.'.join(sorted([str(name) for name in names]))

# ----------------------------------------------------------------------------------------
def split_key(key):
    """ 
    Reverse the action of get_key(). 
    NOTE does not necessarily give a_ and b_ in the same order, though
    NOTE also that b_name may not be the same (if 0), and this just returns strings, even if original names were ints
    """
    # assert len(re.findall('.', key)) == 1  # make sure none of the keys had a dot in it
    return key.split('.')

# ----------------------------------------------------------------------------------------
def prep_dir(dirname, wildling=None, multilings=None):
    """ make <dirname> if it d.n.e., and if shell glob <wildling> is specified, remove existing files which are thereby matched """
    assert dirname is not None  # would typically happen if <dirname> is from a command line argument
    if os.path.exists(dirname):
        if wildling is not None:
            for fname in glob.glob(dirname + '/' + wildling):
                os.remove(fname)
        if multilings is not None:  # allow multiple file name suffixes
            for wild in multilings:
                for fname in glob.glob(dirname + '/' + wild):
                    os.remove(fname)
    else:
        os.makedirs(dirname)
    if wildling is not None or multilings is not None:
        if len([fname for fname in os.listdir(dirname) if os.path.isfile(dirname + '/' + fname)]) != 0:  # make sure there's no other files in the dir
            print 'ERROR files remain in',dirname,'despite wildling',
            if wildling != None:
                print wildling
            if multilings != None:
                print multilings
            assert False

# ----------------------------------------------------------------------------------------
def process_input_line(info):
    """ 
    Attempt to convert all the keys and values in <info> from str to int.
    The keys listed in <list_columns> will be split as colon-separated lists before intification.
    """

    ccfg = column_configs  # shorten the name a bit

    for key, val in info.items():
        if key is None:
            continue

        convert_fcn = pass_fcn  # dummy fcn, just returns the argument
        if key in ccfg['ints']:
            convert_fcn = int
        elif key in ccfg['floats']:
            convert_fcn = float
        elif key in ccfg['literals']:
            convert_fcn = ast.literal_eval

        if key in ccfg['lists']:
            info[key] = [convert_fcn(val) for val in info[key].split(':')]
        else:
            info[key] = convert_fcn(info[key])

# ----------------------------------------------------------------------------------------
def get_line_for_output(info):
    """ Reverse the action of process_input_line() """ 
    outfo = {}
    for key, val in info.items():
        if key in column_configs['lists']:
            outfo[key] = ':'.join([str(v) for v in info[key]])
        else:
            outfo[key] = str(info[key])
    return outfo

# ----------------------------------------------------------------------------------------
def merge_csvs(outfname, csv_list, cleanup=True):
    """ NOTE copy of merge_hmm_outputs in partitiondriver, I should really combine the two functions """
    header = None
    outfo = []
    # print 'merging'
    for infname in csv_list:
        # print '  ', infname
        workdir = os.path.dirname(infname)
        with opener('r')(infname) as sub_outfile:
            reader = csv.DictReader(sub_outfile)
            header = reader.fieldnames
            for line in reader:
                outfo.append(line)
        if cleanup:
            os.remove(infname)
            os.rmdir(workdir)

    outdir = '.' if os.path.dirname(outfname) == '' else os.path.dirname(outfname)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with opener('w')(outfname) as outfile:
        writer = csv.DictWriter(outfile, header)
        writer.writeheader()
        for line in outfo:
            writer.writerow(line)

# ----------------------------------------------------------------------------------------
def get_mutation_rate(germlines, line, restrict_to_region=''):
    naive_seq = get_full_naive_seq(germlines, line)  # NOTE this includes the fv and jf insertions
    muted_seq = line['seq']
    if restrict_to_region == '':  # NOTE this is very similar to code in performanceplotter. I should eventually cut it out of there and combine them, but I'm nervous a.t.m. because of all the complications there of having the true *and* inferred sequences so I'm punting
        mashed_naive_seq = ''
        mashed_muted_seq = ''
        for region in regions:  # can't use the full sequence because we have no idea what the mutations were in the inserts. So have to mash together the three regions
            bounds = get_regional_naive_seq_bounds(region, germlines, line)
            mashed_naive_seq += naive_seq[bounds[0] : bounds[1]]
            mashed_muted_seq += muted_seq[bounds[0] : bounds[1]]
    else:
        bounds = get_regional_naive_seq_bounds(restrict_to_region, germlines, line)
        naive_seq = naive_seq[bounds[0] : bounds[1]]
        muted_seq = muted_seq[bounds[0] : bounds[1]]


    # print 'restrict %s' % restrict_to_region
    # color_mutants(naive_seq, muted_seq, print_result=True, extra_str='  ')
    return hamming_fraction(naive_seq, muted_seq)

# ----------------------------------------------------------------------------------------
def print_linsim_output(outstr):
    import ast
    linsim_out = ast.literal_eval(outstr)
    print '       true clusters %d' % linsim_out['true_cluster_count']
    print '   inferred clusters %d' % linsim_out['inferred_cluster_count']
    print '  mutual information %f' % linsim_out['metrics']['mi']
    print '         adjusted mi %f' % linsim_out['metrics']['adjusted_mi']
    print '       normalized mi %f' % linsim_out['metrics']['normalized_mi']
    print '  completeness score %f' % linsim_out['metrics']['completeness_score']
    print '   homogeneity score %f' % linsim_out['metrics']['homogeneity_score']

# ----------------------------------------------------------------------------------------
def process_out_err(out, err, extra_str='0', info=None, subworkdir=None):
    """ NOTE something in this chain seems to block or truncate or some such nonsense if you make it too big """
    if subworkdir is not None:
        def readfile(fname):
            ftmp = open(fname)
            fstr = ''.join(ftmp.readlines())
            ftmp.close()
            os.remove(fname)
            return fstr
        out = readfile(subworkdir + '/out')
        err = readfile(subworkdir + '/err')

    print_str = ''
    for line in err.split('\n'):
        if 'srun: job' in line and 'queued and waiting for resources' in line:
            continue
        if 'srun: job' in line and 'has been allocated resources' in line:
            continue
        if 'GSL_RNG_TYPE=' in line or 'GSL_RNG_SEED=' in line:
            continue
        if '[ig_align] Read' in line or '[ig_align] Aligned' in line:
            continue
        if len(line.strip()) > 0:
            print_str += line + '\n'

    for line in out.split('\n'):
        if info is not None and 'calculated' in line:  # keep track of how many vtb and fwd calculations the process made
            info['vtb'], info['fwd'] = 0, 0
            words = line.split()
            if words[1] == 'vtb' and words[3] == 'fwd':
                info['vtb'] = int(words[2])
                info['fwd'] = int(words[4])
            else:
                print 'ERROR bad calculated line: %s' % line
    # if info is not None and ('vtb' not in info or 'fwd' not in info):
    #     print 'weird, didnt find anything for info:'
    #     print 'out x', out, 'x'
    #     print 'err x', err, 'x'

    print_str += out

    if print_str != '':
        print '      --> proc %s' % extra_str
        print print_str

# ----------------------------------------------------------------------------------------
def find_first_non_ambiguous_base(seq):
    """ return index of first non-ambiguous base """
    for ib in range(len(seq)):
        if seq[ib] not in ambiguous_bases:
            return ib
    assert False  # whole sequence was ambiguous... probably shouldn't get here

# ----------------------------------------------------------------------------------------
def find_last_non_ambiguous_base_plus_one(seq):
    for ib in range(len(seq) - 1, -1, -1):  # count backwards from the end
        if seq[ib] not in ambiguous_bases:  # find first non-ambiguous base
            return ib + 1  # add one for easy slicing

    assert False  # whole sequence was ambiguous... probably shouldn't get here

# ----------------------------------------------------------------------------------------
def remove_ambiguous_ends(seq):
    """ remove ambiguous bases from the left and right ends of <seq> """
    i_seq_start = find_first_non_ambiguous_base(seq)
    i_seq_end = find_last_non_ambiguous_base_plus_one(seq)
    return seq[i_seq_start : i_seq_end]

# ----------------------------------------------------------------------------------------
def get_true_partition(reco_info, ids=None):
    """ 
    Group ids into their true clonal families.
    If <ids> is specified, only do those, otherwise do all of the ones in <reco_info>.
    """
    if ids is None:
        id_list = reco_info.keys()
    else:
        id_list = ids
    true_partition = {}
    for uid in id_list:
        rid = reco_info[uid]['reco_id']
        if rid not in true_partition:
            true_partition[rid] = []
        true_partition[rid].append(uid)
    return true_partition.values()

# ----------------------------------------------------------------------------------------
def get_partition_from_str(partition_str):
    """ NOTE there's code in some other places that do the same thing """
    clusters = partition_str.split(';')
    partition = [cl.split(':') for cl in clusters]
    return partition

# ----------------------------------------------------------------------------------------
def get_str_from_partition(partition):
    """ NOTE there's code in some other places that do the same thing """
    clusters = [':'.join(cl) for cl in partition]
    partition_str = ';'.join(clusters)
    return partition_str

# ----------------------------------------------------------------------------------------
def new_ccfs_that_need_better_names(partition, true_partition, reco_info, debug=False):
    check_intersection_and_complement(partition, true_partition)

    reco_ids = {uid : reco_info[uid]['reco_id'] for cluster in partition for uid in cluster}  # just a teensy lil' optimization

    def get_clonal_fraction(uid, inferred_cluster):
        """ Return the fraction of seqs in <uid>'s inferred cluster which are really clonal. """
        n_clonal = 0
        for tmpid in inferred_cluster:  # NOTE this includes the case where tmpid equal to uid
            if reco_ids[tmpid] == reco_ids[uid]:
                n_clonal += 1
        return float(n_clonal) / len(inferred_cluster)

    def get_fraction_present(uid, inferred_cluster, true_cluster):
        """ Return the fraction of <uid>'s true clonemates which appear in <uid>'s inferred cluster. """
        n_present = 0
        for tmpid in true_cluster:  # NOTE this includes the case where tmpid equal to uid
            if tmpid in inferred_cluster:
                n_present += 1
        return float(n_present) / len(true_cluster)

    mean_clonal_fraction, mean_fraction_present = 0., 0.
    n_uids = 0
    for true_cluster in true_partition:
        for uid in true_cluster:
            inferred_cluster = partition[find_uid_in_partition(uid, partition)]
            mean_clonal_fraction += get_clonal_fraction(uid, inferred_cluster)
            mean_fraction_present +=  get_fraction_present(uid, inferred_cluster, true_cluster)
            n_uids += 1

    if n_uids > 1e6:
        raise Exception('you should start worrying about numerical precision if you\'re going to run on this many queries')

    return mean_clonal_fraction / n_uids, mean_fraction_present / n_uids

# ----------------------------------------------------------------------------------------
def correct_cluster_fractions(partition, true_partition, debug=False):
    # return new_ccfs_that_need_better_names(partition, true_partition, debug)  # hey, look, I'm a hack! Seriously, though, the new ccfs above are pretty similar, except they're per-sequence rather than per-cluster, so they don't get all scatterbrained and shit when a sample's only got a few clusters. Also, you still get partial credit for how good your cluster is, it's not just all-or-nothing.
    raise Exception('deprecated!')

    def find_clusters_with_ids(ids, partition):
        """ find all clusters in <partition> that contain at least one of <ids> """
        clusters = []
        for cluster in partition:
            for uid in ids:
                if uid in cluster:
                    clusters.append(cluster)
                    break
        return clusters

    check_intersection_and_complement(partition, true_partition)

    n_under_merged, n_over_merged = 0, 0
    for trueclust in true_partition:
        if debug:
            print ''
            print '   true %s' % (len(trueclust) if len(trueclust) > 15 else trueclust)
        infclusters = find_clusters_with_ids(trueclust, partition)  # list of inferred clusters that contain any ids from the true cluster
        if debug and len(infclusters) > 1:
            print '  infclusters %s' % infclusters
        assert len(infclusters) > 0
        under_merged = len(infclusters) > 1  # ids in true cluster are not all in the same inferred cluster
        over_merged = False  # at least one inferred cluster with an id in true cluster also contains an id not in true cluster
        for iclust in infclusters:
            if debug:
                print '   inferred %s' % (len(iclust) if len(iclust) > 15 else iclust)
            for uid in iclust:
                if uid not in trueclust:
                    over_merged = True
                    break
            if over_merged:
                break
        if debug:
            print '  under %s   over %s' % (under_merged, over_merged)
        if under_merged:
            n_under_merged += 1
        if over_merged:
            n_over_merged += 1

    under_frac = float(n_under_merged) / len(true_partition)
    over_frac = float(n_over_merged) / len(true_partition)
    if debug:
        print '  under %.2f   over %.2f' % (under_frac, over_frac)
    return (1. - under_frac, 1. - over_frac)

# ----------------------------------------------------------------------------------------
def partition_similarity_matrix(meth_a, meth_b, partition_a, partition_b, n_biggest_clusters, debug=False):
    """ Return matrix whose ij^th entry is the size of the intersection between <partition_a>'s i^th biggest cluster and <partition_b>'s j^th biggest """
    def intersection_size(cl_1, cl_2):
        isize = 0
        for uid in cl_1:
            if uid in cl_2:
                isize += 1
        return isize

    # n_biggest_clusters = 10
    def sort_within_clusters(part):
        for iclust in range(len(part)):
            part[iclust] = sorted(part[iclust])

    # a_clusters = sorted(partition_a, key=len, reverse=True)[ : n_biggest_clusters]  # i.e. the n biggest clusters
    # b_clusters = sorted(partition_b, key=len, reverse=True)[ : n_biggest_clusters]
    sort_within_clusters(partition_a)
    sort_within_clusters(partition_b)
    a_clusters = sorted(sorted(partition_a), key=len, reverse=True)[ : n_biggest_clusters]  # i.e. the n biggest clusters
    b_clusters = sorted(sorted(partition_b), key=len, reverse=True)[ : n_biggest_clusters]

    smatrix = []
    pair_info = []  # list of full pair info (e.g. [0.8, ick)
    max_pair_info = 5
    for clust_a in a_clusters:
        # if debug:
        #     print clust_a
        smatrix.append([])
        for clust_b in b_clusters:
            # norm_factor = 1.  # don't normalize
            norm_factor = 0.5 * (len(clust_a) + len(clust_b))  # mean size
            # norm_factor = min(len(clust_a), len(clust_b))  # smaller size
            intersection = intersection_size(clust_a, clust_b)
            isize = float(intersection) / norm_factor
            # if debug:
            #     print '    %.2f  %5d   %5d %5d' % (isize, intersection, len(clust_a), len(clust_b))
            if isize == 0.:
                isize = None
            smatrix[-1].append(isize)
            if isize is not None:
                if len(pair_info) < max_pair_info:
                    pair_info.append([intersection, [clust_a, clust_b]])
                    pair_info = sorted(pair_info, reverse=True)
                elif intersection > pair_info[-1][0]:
                    pair_info[-1] = [intersection, [clust_a, clust_b]]
                    pair_info = sorted(pair_info, reverse=True)

    if debug:
        print 'intersection     a_rank  b_rank'
        for it in pair_info:
            print '%-4d  %3d %3d   %s  %s' % (it[0], a_clusters.index(it[1][0]), b_clusters.index(it[1][1]), ':'.join(it[1][0]), ':'.join(it[1][1]))

    # with open('erick/' + meth_a + '_' + meth_b + '.csv', 'w') as csvfile:
    #     writer = csv.DictWriter(csvfile, ['a_meth', 'b_meth', 'intersection', 'a_rank', 'b_rank', 'a_cluster', 'b_cluster'])
    #     writer.writeheader()
    #     for it in pair_info:
    #         writer.writerow({'a_meth' : meth_a, 'b_meth' : meth_b,
    #                          'intersection' : it[0],
    #                          'a_rank' : a_clusters.index(it[1][0]), 'b_rank' : b_clusters.index(it[1][1]),
    #                          'a_cluster' : ':'.join(it[1][0]), 'b_cluster' : ':'.join(it[1][1])})

    a_cluster_lengths, b_cluster_lengths = [len(c) for c in a_clusters], [len(c) for c in b_clusters]
    return a_cluster_lengths, b_cluster_lengths, smatrix

# ----------------------------------------------------------------------------------------
def find_uid_in_partition(uid, partition):
    found = False
    for iclust in range(len(partition)):
        if uid in partition[iclust]:
            found = True
            break
    if found:
        return iclust
    else:
        raise Exception('couldn\'t find %s in %s\n' % (uid, partition))

# ----------------------------------------------------------------------------------------
def check_intersection_and_complement(part_a, part_b):
    """ make sure two partitions have identical uid lists """
    for cluster in part_a:
        for uid in cluster:
            find_uid_in_partition(uid, part_b)
    for cluster in part_b:  # NOTE we could avoid looping over some of these if we were so inclined
        for uid in cluster:
            find_uid_in_partition(uid, part_a)

# ----------------------------------------------------------------------------------------
def get_cluster_list_for_sklearn(part_a, part_b):
    # convert from partition format {cl_1 : [seq_a, seq_b], cl_2 : [seq_c]} to [cl_1, cl_1, cl_2]
    # NOTE this will be really slow for larger partitions

    # first make sure that <part_a> has every uid in <part_b> (the converse is checked below)
    for jclust in range(len(part_b)):
        for uid in part_b[jclust]:
            find_uid_in_partition(uid, part_a)  # raises exception if not found

    # then make the cluster lists
    clusts_a, clusts_b = [], []
    for iclust in range(len(part_a)):
        for uid in part_a[iclust]:
            clusts_a.append(iclust)
            clusts_b.append(find_uid_in_partition(uid, part_b))

    return clusts_a, clusts_b

# ----------------------------------------------------------------------------------------
def adjusted_mutual_information(partition_a, partition_b):
    return -1.  # not using it any more, and it's really slow
    # clusts_a, clusts_b = get_cluster_list_for_sklearn(partition_a, partition_b)
    # return sklearn.metrics.cluster.adjusted_mutual_info_score(clusts_a, clusts_b)

# ----------------------------------------------------------------------------------------
def add_missing_uids_as_singletons_to_inferred_partition(partition_with_missing_uids, true_partition=None, all_ids=None, debug=True):
    """ return a copy of <partition_with_missing_uids> which has had any uids which were missing inserted as singletons (i.e. uids which were in <true_partition>) """

    if true_partition is None:  # it's less confusing than the alternatives, I swear
        true_partition = [[uid, ] for uid in all_ids]

    partition_with_uids_added = copy.deepcopy(partition_with_missing_uids)
    missing_ids = []
    for cluster in true_partition:
        for uid in cluster:
            try:
                find_uid_in_partition(uid, partition_with_missing_uids)
            except:
                partition_with_uids_added.append([uid, ])
                missing_ids.append(uid)
    if debug:
        print '  %d (of %d) ids missing from partition (%s)' % (len(missing_ids), sum([len(c) for c in true_partition]), ' '.join(missing_ids))
    return partition_with_uids_added

# ----------------------------------------------------------------------------------------
def remove_missing_uids_from_true_partition(true_partition, partition_with_missing_uids, debug=True):
    """ return a copy of <true_partition> which has had any uids which do not occur in <partition_with_missing_uids> removed """
    true_partition_with_uids_removed = []
    missing_ids = []
    for cluster in true_partition:
        new_cluster = []
        for uid in cluster:
            try:
                find_uid_in_partition(uid, partition_with_missing_uids)
                new_cluster.append(uid)
            except:
                missing_ids.append(uid)
        if len(new_cluster) > 0:
            true_partition_with_uids_removed.append(new_cluster)
    if debug:
        print '  %d (of %d) ids missing from partition (%s)' % (len(missing_ids), sum([len(c) for c in true_partition]), ' '.join(missing_ids))
    return true_partition_with_uids_removed

# ----------------------------------------------------------------------------------------
def generate_incorrect_partition(true_partition, misassign_fraction, error_type, debug=False):
    """ 
    Generate an incorrect partition from <true_partition>.
    We accomplish this by removing <n_misassigned> seqs at random from their proper cluster, and putting each in either a
    cluster chosen at random from the non-proper clusters (<error_type> 'reassign') or in its own partition (<error_type> 'singleton').
    """
    # true_partition = [['b'], ['a', 'c', 'e', 'f'], ['d', 'g'], ['h', 'j']]
    # debug = True

    new_partition = copy.deepcopy(true_partition)
    nseqs = sum([len(c) for c in true_partition])
    n_misassigned = int(misassign_fraction * nseqs)
    if debug:
        print '  misassigning %d / %d seqs (should be clsoe to %.3f)' % (n_misassigned, nseqs, misassign_fraction)
        print '  before', new_partition

    uids = [uid for cluster in true_partition for uid in cluster]
    for _ in range(n_misassigned):
        uid = uids[random.randint(0, len(uids) - 1)]  # choose a uid to misassign (note that randint() is inclusive)
        uids.remove(uid)
        iclust = find_uid_in_partition(uid, new_partition)
        new_partition[iclust].remove(uid)  # remove it
        if [] in new_partition:
            new_partition.remove([])
        if error_type == 'singletons':  # put the sequence in a cluster by itself
            new_partition.append([uid, ])
            if debug:
                print '    %s: %d --> singleton' % (uid, iclust)
        elif error_type == 'reassign':  # choose a different cluster to add it to
            inewclust = iclust
            while inewclust == iclust:  # hm, this won't work if there's only one cluster in the partition. Oh, well, that probably won't happen
                inewclust = random.randint(0, len(new_partition) - 1)
            new_partition[inewclust].append(uid)
            if debug:
                print '    %s: %d --> %d' % (uid, iclust, inewclust)
        else:
            raise Exception('%s not among %s' % (error_type, 'singletons, reassign'))
    if debug:
        print '  after', new_partition
    return new_partition

# ----------------------------------------------------------------------------------------
def subset_files(uids, fnames, outdir, uid_header='Sequence ID', delimiter='\t', debug=False):
    """ rewrite csv files <fnames> to <outdir>, removing lines with uid not in <uids> """
    for fname in fnames:
        with open(fname) as infile:
            reader = csv.DictReader(infile, delimiter=delimiter)
            with open(outdir + '/' + os.path.basename(fname), 'w') as outfile:
                writer = csv.DictWriter(outfile, reader.fieldnames, delimiter=delimiter)
                writer.writeheader()
                for line in reader:
                    if line[uid_header] in uids:
                        writer.writerow(line)

# ----------------------------------------------------------------------------------------
def add_indels_to_germline_strings(line, indelfo, original_seqs, lengths, eroded_seqs, reverse_indels):
    lastfo = indelfo['indels'][-1]
    if lastfo['type'] == 'insertion':
        chunks = [line['fv_insertion'], eroded_seqs['v'], line['vd_insertion'], eroded_seqs['d'], line['dj_insertion'], eroded_seqs['j'], line['jf_insertion']]
        chunknames =  ['fv_insertion', 'v', 'vd_insertion', 'd', 'dj_insertion', 'j', 'jf_insertion']
        weirdolist = []
        for ichunk in range(len(chunks)):
            for inuke in range(len(chunks[ichunk])):
                weirdolist.append(chunknames[ichunk])
        thischunk = weirdolist[lastfo['pos']]
        # offset = 0
        # if reverse_indels:
        offset = weirdolist.index(thischunk)  # index of first occurence
        if thischunk in eroded_seqs:
            # original_seqs[thischunk] = original_seqs[thischunk][ : lastfo['pos'] - offset] + '*' * lastfo['len'] + eroded_seqs[thischunk][lastfo['pos'] - offset : ]
            lengths[thischunk] += lastfo['len']
            eroded_seqs[thischunk] = eroded_seqs[thischunk][ : lastfo['pos'] - offset] + '*' * lastfo['len'] + eroded_seqs[thischunk][lastfo['pos'] - offset : ]
        else:
            print '     unhandled indel in NTIs'
            pass
    else:
        pass

# ----------------------------------------------------------------------------------------
def undo_indels(indelfo):
    """ not finished """
    rseq = indelfo['reversed_seq']
    oseq = rseq  # original sequence
    for i_indel in range(len(indelfo['indels']) - 1, 0, -1):
    # for i_indel in range(len(indelfo['indels'])):
        idl = indelfo['indels'][i_indel]
        if idl['type'] == 'insertion':
            oseq = oseq[ : idl['pos']] + idl['seqstr'] + oseq[idl['pos'] : ]
        elif idl['type'] == 'deletion':
            if rseq[idl['pos'] : idl['pos'] + idl['len']] != idl['seqstr']:
                raise Exception('found %s instead of expected insertion (%s)' % (rseq[idl['pos'] : idl['pos'] + idl['len']], idl['seqstr']))
            oseq = oseq[ : idl['pos']] + oseq[idl['pos'] + idl['len'] : ]
    # print '              reversed %s' % rseq
    # print '              original %s' % oseq

# ----------------------------------------------------------------------------------------
def csv_to_fasta(infname, outfname=None, name_column='unique_id', seq_column='seq', n_max_lines=None):
    if not os.path.exists(infname):
        raise Exception('input file %s d.n.e.' % infname)
    if outfname is None:
        assert '.csv' in infname
        outfname = infname.replace('.csv', '.fa')
    if os.path.exists(outfname):
        print '  csv --> fasta: overwriting %s' % outfname

    if '.csv' in infname:
        delimiter = ','
    elif '.tsv' in infname:
        delimiter = '\t'
    else:
        assert False
    
    with open(infname) as infile:
        reader = csv.DictReader(infile, delimiter=delimiter)
        with open(outfname, 'w') as outfile:
            n_lines = 0
            for line in reader:
                if name_column not in line:
                    name_column = 'name'
                    seq_column = 'nucleotide'
                n_lines += 1
                if n_max_lines is not None and n_lines > n_max_lines:
                    break
                outfile.write('>%s\n' % line[name_column])
                outfile.write('%s\n' % line[seq_column])

# ----------------------------------------------------------------------------------------
def print_heapy(extrastr, heap):
    'Partition of a set of 1511530 objects. Total size = 188854824 bytes.'
    heapstr = heap.__str__()
    total = None
    for line in heapstr.split('\n'):
        if 'Total size' in line:
            total = int(line.split()[10])
    if total is None:
        print 'oops'
        print heapstr
        sys.exit()
    print 'mem total %.3f MB    %s' % (float(total) / 1e6, extrastr)

# ----------------------------------------------------------------------------------------
def auto_slurm(n_procs):
    """ Return true if we want to force slurm usage, e.g. if there's more processes than cores """
    try:
        check_output(['which', 'srun'])
        slurm_exists = True
    except CalledProcessError:
        slurm_exists = False
    ncpu = multiprocessing.cpu_count()
    if slurm_exists and n_procs > ncpu:
        return True
    return False

# ----------------------------------------------------------------------------------------
def synthesize_single_seq_line(glfo, line, iseq):
    """ without modifying <line>, make a copy of it corresponding to a single-sequence event with the <iseq>th sequence """
    hmminfo = copy.deepcopy(line)  # make a copy of the info, into which we'll insert the sequence-specific stuff
    hmminfo['seq'] = line['seqs'][iseq]
    hmminfo['unique_id'] = line['unique_ids'][iseq]
    del hmminfo['unique_ids']
    del hmminfo['seqs']
    if 'aligned_v_seqs' in hmminfo:
        hmminfo['aligned_v_seq'] = line['aligned_v_seqs'][iseq]
        del hmminfo['aligned_v_seqs']
    add_match_info(glfo, hmminfo)
    return hmminfo

# ----------------------------------------------------------------------------------------                    
def count_gaps(seq, istop=None):
    """ return number of gap characters up to, but not including <istop> """
    if istop is not None:
        seq = seq[ : istop]
    return sum([seq.count(gc) for gc in gap_chars])

# ----------------------------------------------------------------------------------------
def add_v_alignments(glfo, line, debug=False):
    """ add dots according to the imgt gapping scheme """

    def n_gaps(seq):
        return sum([seq.count(gc) for gc in gap_chars])

    aligned_v_seqs = []
    for iseq in range(len(line['seqs'])):
        hmminfo = synthesize_single_seq_line(glfo, line, iseq)

        v_qr_seq = hmminfo['v_qr_seq']
        v_gl_seq = hmminfo['v_gl_seq']
        aligned_v_gl_seq = glfo['aligned-v-genes']['v'][hmminfo['v_gene']]
        if len(v_qr_seq) != len(v_gl_seq):
            raise Exception('v bits not same length for %s\n%s' % (line['unique_ids'], line))

        if debug:
            print 'before alignment'
            print '   qr   ', v_qr_seq
            print '   gl   ', v_gl_seq
            print '   al gl', aligned_v_gl_seq

        if len(aligned_v_gl_seq) != line['v_5p_del'] + len(v_gl_seq) + hmminfo['v_3p_del'] + n_gaps(aligned_v_gl_seq):
            raise Exception('lengths don\'t match up\n%s\n%s + %d' % (aligned_v_gl_seq, v_gl_seq + gap_chars[0] * n_gaps(aligned_v_gl_seq), hmminfo['v_3p_del']))

        v_qr_seq = 'N' * line['v_5p_del'] + v_qr_seq + 'N' * line['v_3p_del']
        v_gl_seq = 'N' * line['v_5p_del'] + v_gl_seq + 'N' * line['v_3p_del']

        for ibase in range(len(aligned_v_gl_seq)):
            if aligned_v_gl_seq[ibase] in gap_chars:
                v_qr_seq = v_qr_seq[ : ibase] + gap_chars[0] + v_qr_seq[ibase : ]
                v_gl_seq = v_gl_seq[ : ibase] + gap_chars[0] + v_gl_seq[ibase : ]
            else:
                if v_gl_seq[ibase] != 'N' and v_gl_seq[ibase] != aligned_v_gl_seq[ibase]:
                    raise Exception('bases don\'t match at position %d in\n%s\n%s' % (ibase, v_gl_seq, aligned_v_gl_seq))

        if debug:
            print 'after alignment'
            print '   qr   ', v_qr_seq
            print '   gl   ', v_gl_seq
            print '   al gl', aligned_v_gl_seq

        if len(v_qr_seq) != len(v_gl_seq) or len(v_qr_seq) != len(aligned_v_gl_seq):
            raise Exception('lengths don\'t match up:\n%s\n%s\n%s' % (v_qr_seq, v_gl_seq, aligned_v_gl_seq))
        aligned_v_seqs.append(v_qr_seq)  # TODO is this supposed to be just the v section of the query sequence, or the whole sequence? (if it's the latter, I don't know what to do about alignments)

    line['aligned_v_seqs'] = aligned_v_seqs

# ----------------------------------------------------------------------------------------
def intexterpolate(x1, y1, x2, y2, x):
    """ interpolate/extrapolate linearly based on two points in 2-space, returning y-value corresponding to <x> """
    m = (y2 - y1) / (x2 - x1);
    b = 0.5 * (y1 + y2 - m*(x1 + x2));
    # if debug:
    #     for x in [x1, x2]:
    #         print '%f x + %f = %f' % (m, b, m*x + b)
    return m * x + b


# ----------------------------------------------------------------------------------------
def get_padding_parameters(queries, info, glfo, debug=False):
    all_v_matches, all_j_matches = set(), set()
    for query in queries:
        swfo = info[query]
        for match in swfo['all'].split(':'):
            if get_region(match) == 'v':
                all_v_matches.add(match)
            elif get_region(match) == 'j':
                all_j_matches.add(match)

    maxima = {'gl_cpos' : None, 'gl_cpos_to_j_end' : None}  #, 'fv_insertion_len' : None, 'jf_insertion_len' : None}
    for query in queries:
        fvstuff = max(0, len(swfo['fv_insertion']) - swfo['v_5p_del'])  # we always want to pad out to the entire germline sequence, so don't let this go negative
        jfstuff = max(0, len(swfo['jf_insertion']) - swfo['j_3p_del'])

        for v_match in all_v_matches:  # NOTE have to loop over all gl matches, even ones for other sequences, because we want bcrham to be able to compare any sequence to any other
            gl_cpos = glfo['cyst-positions'][v_match] + fvstuff
            if maxima['gl_cpos'] is None or gl_cpos > maxima['gl_cpos']:
                maxima['gl_cpos'] = gl_cpos

        swfo = info[query]
        seq = swfo['seq']
        cpos = swfo['cyst_position']  # cyst position in query sequence (as opposed to gl_cpos, which is in germline allele)
        for j_match in all_j_matches:  # NOTE have to loop over all gl matches, even ones for other sequences, because we want bcrham to be able to compare any sequence to any other
            # TODO this is totally wrong -- I'm only storing j_3p_del for the best match... but hopefully it'll give enough padding for the moment
            gl_cpos_to_j_end = len(seq) - cpos + swfo['j_3p_del'] + jfstuff
            if maxima['gl_cpos_to_j_end'] is None or gl_cpos_to_j_end > maxima['gl_cpos_to_j_end']:
                maxima['gl_cpos_to_j_end'] = gl_cpos_to_j_end

        # if maxima['fv_insertion_len'] is None or len(swfo['fv_insertion']) > maxima['fv_insertion_len']:
        #     maxima['fv_insertion_len'] = len(swfo['fv_insertion'])
        # if maxima['jf_insertion_len'] is None or len(swfo['jf_insertion']) > maxima['jf_insertion_len']:
        #     maxima['jf_insertion_len'] = len(swfo['jf_insertion'])

    if debug:
        print '    maxima:',
        for k, v in maxima.items():
            print '%s %d    ' % (k, v),
        print ''
    return maxima

# ----------------------------------------------------------------------------------------
def pad_seqs_to_same_length(queries, info, glfo, indelfo, debug=False):
    # TODO it makes no sense to move this out of waterer.py, I should put it back
    """
    Pad all sequences in <seqinfo> to the same length to the left and right of their conserved cysteine positions.
    Next, pads all sequences further out (if necessary) such as to eliminate all v_5p and j_3p deletions.
    """

    maxima = get_padding_parameters(queries, info, glfo, debug)

    for query in queries:
        swfo = info[query]
        if 'padded' in swfo:  # already added padded information (we're probably partitioning, and this is not the first step)
            return
        seq = swfo['seq']
        cpos = swfo['cyst_position']
        if cpos < 0 or cpos >= len(seq):
            print 'hm now what do I want to do here?'
        k_v = swfo['k_v']

        # padleft = maxima['fv_insertion_len'] + maxima['gl_cpos'] - cpos  # left padding: biggest germline cpos minus cpos in this sequence
        # padright = maxima['gl_cpos_to_j_end'] + maxima['jf_insertion_len'] - (len(seq) - cpos)
        padleft = maxima['gl_cpos'] - cpos  # left padding: biggest germline cpos minus cpos in this sequence
        padright = maxima['gl_cpos_to_j_end'] - (len(seq) - cpos)
        if padleft < 0 or padright < 0:
            raise Exception('bad padding %d %d for %s' % (padleft, padright, query))

        swfo['padded'] = {}
        padfo = swfo['padded']  # shorthand
        assert len(ambiguous_bases) == 1  # could allow more than one, but it's not implemented a.t.m.
        padfo['seq'] = padleft * ambiguous_bases[0] + seq + padright * ambiguous_bases[0]
        if query in indelfo:
            if debug:
                print '    also padding reversed sequence'
            indelfo[query]['reversed_seq'] = padleft * ambiguous_bases[0] + indelfo[query]['reversed_seq'] + padright * ambiguous_bases[0]
        padfo['k_v'] = {'min' : k_v['min'] + padleft, 'max' : k_v['max'] + padleft}
        padfo['cyst_position'] = swfo['cyst_position'] + padleft
        padfo['padleft'] = padleft
        padfo['padright'] = padright
        if debug:
            print '      pad %d %d   %s' % (padleft, padright, query)
            print '     %d --> %d (%d-%d --> %d-%d)' % (len(seq), len(padfo['seq']),
                                                        k_v['min'], k_v['max'],
                                                        padfo['k_v']['min'], padfo['k_v']['max'])

    if debug:
        for query in queries:
            print '%20s %s' % (query, info[query]['padded']['seq'])

# ----------------------------------------------------------------------------------------
def find_genes_that_have_hmms(parameter_dir):
    yamels = glob.glob(parameter_dir + '/hmms/*.yaml')
    if len(yamels) == 0:
        raise Exception('no yamels in %s' % parameter_dir + '/hmms')

    genes = []
    for yamel in yamels:
        gene = unsanitize_name(os.path.basename(yamel).replace('.yaml', ''))
        genes.append(gene)

    return genes
