""" A few utility functions. At the moment simply functions used in recombinator which do not
require member variables. """

import sys
import os
import random
import ast
import math
import glob
from collections import OrderedDict
import csv
from subprocess import check_output, CalledProcessError, Popen
import multiprocessing
import copy

from opener import opener
import seqfileopener
import glutils

#----------------------------------------------------------------------------------------
# NOTE I also have an eps defined in hmmwriter. Simplicity is the hobgoblin of... no, wait, that's just plain ol' stupid to have two <eps>s defined
eps = 1.0e-10  # if things that should be 1.0 are this close to 1.0, blithely keep on keepin on. kinda arbitrary, but works for the moment
def is_normed(probs, this_eps=eps):
    if hasattr(probs, 'keys'):  # if it's a dict, call yourself with a list of the dict's values
        return is_normed([val for val in probs.values()])
    elif hasattr(probs, '__iter__'):  # if it's a list call yourself with their sum
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

# ----------------------------------------------------------------------------------------
# values used when simulating from scratch
# scratch_mean_mute_freqs = {'v' : 0.03, 'd' : 0.8, 'j' : 0.06}
# scratch_mean_mute_freqs['all'] = numpy.mean([v for v in scratch_mean_mute_freqs.values()])
scratch_mean_erosion_lengths = {'v_3p' : 2, 'd_5p' : 3, 'd_3p' : 3, 'j_5p' : 4}
scratch_mean_insertion_lengths = {
    'h': {'vd' : 4, 'dj' : 4},
    'k': {'vd' : 0, 'dj' : 8},
    'l': {'vd' : 0, 'dj' : 8}
}

# ----------------------------------------------------------------------------------------
regions = ['v', 'd', 'j']
chains = ['h', 'k', 'l']

def getregions(chain):
    if chain == 'h':
        return regions
    else:
        return [r for r in regions if r != 'd']

# def getboundaries(chain):
#     regions = getregions(chain)
#     return [regions[i] + regions[i+1] for i in range(len(regions) - 1)]
# def get_region_pairs(chain):
#     return [{'left' : bound[0], 'right' : bound[1]} for bound in getboundaries(chain)]
def region_pairs():
    return [{'left' : bound[0], 'right' : bound[1]} for bound in boundaries]

real_erosions = ['v_3p', 'd_5p', 'd_3p', 'j_5p']
# NOTE since we now handle v_5p and j_3p deletions by padding with Ns, the hmm does *not* allow actual v_5p and j_3p deletions.
# This means that while we write parameters for v_5p and j_3p deletions to the parameter dir, these are *not* used in making the
# hmm yamels -- which is what we want, because we want to be able to read in short data reads but make full-length simulation.
effective_erosions = ['v_5p', 'j_3p']
boundaries = ['vd', 'dj']
effective_boundaries = ['fv', 'jf']
nukes = ['A', 'C', 'G', 'T']
ambiguous_bases = ['N', ]
gap_chars = ['.', '-']
expected_characters = set(nukes + ambiguous_bases + gap_chars)
conserved_codons = {
    'h' : {'v' : 'cyst', 'j' : 'tryp'},
    'k' : {'v' : 'cyst', 'j' : 'phen'},
    'l' : {'v' : 'cyst', 'j' : 'phen'}
}
codon_table = {
    'cyst' : ['TGT', 'TGC'],
    'tryp' : ['TGG', ],
    'phen' : ['TTT', 'TTC'],
    'stop' : ['TAG', 'TAA', 'TGA']
}
# Infrastrucure to allow hashing all the columns together into a dict key.
# Uses a tuple with the variables that are used to index selection frequencies
# NOTE fv and jf insertions are *effective* (not real) insertions between v or j and the framework. They allow query sequences that extend beyond the v or j regions
index_columns = tuple(['v_gene', 'd_gene', 'j_gene', 'v_5p_del', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del', 'j_3p_del', 'fv_insertion', 'vd_insertion', 'dj_insertion', 'jf_insertion'])

index_keys = {}
for i in range(len(index_columns)):  # dict so we can access them by name instead of by index number
    index_keys[index_columns[i]] = i

# ----------------------------------------------------------------------------------------
def get_codon(fname):
    codon = fname.split('-')[0]
    if codon not in [c for chain in chains for c in conserved_codons[chain].values()]:
        raise Exception('couldn\'t get codon from file name %s' % fname)
    return codon

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
    'unique_ids' : 'SEQUENCE_ID',
    'v_gene' : 'V_CALL',
    'd_gene' : 'D_CALL',
    'j_gene' : 'J_CALL',
    'aligned_v_plus_unaligned_dj' : 'SEQUENCE_IMGT',  # TODO is this supposed to be the whole sequence, or just the v bit? UPDATE looks like aligned v, plus the rest of the sequence without any alignment info
    'cdr3_length' : 'JUNCTION_LENGTH'
}

# ----------------------------------------------------------------------------------------
forbidden_characters = set([':', ';', ','])  # strings that are not allowed in sequence ids

functional_columns = ['mutated_invariants', 'in_frames', 'stops']

column_configs = {
    'ints' : ['cdr3_length', 'padlefts', 'padrights'] + [e + '_del' for e in real_erosions + effective_erosions],
    'floats' : ['logprob', 'mut_freqs'],
    'bools' : functional_columns,
    'literals' : ['indelfo', 'indelfos', 'k_v', 'k_d', 'all_matches'],  # simulation has indelfo[s] singular, annotation output has it plural... and I think it actually makes sense to have it that way
    'lists' : ['unique_ids', 'seqs', 'aligned_seqs', 'mut_freqs', 'padlefts', 'padrights'] + ['aligned_' + r + '_seqs' for r in regions] + [r + '_per_gene_support' for r in regions] + functional_columns,
    'lists-of-string-float-pairs' : [r + '_per_gene_support' for r in regions]
}

# keep track of all the *@*@$!ing different keys that happen in the <line>/<hmminfo>/whatever dictionaries
linekeys = {}
linekeys['per_family'] = ['naive_seq', 'cdr3_length', 'codon_positions', 'lengths', 'regional_bounds'] + \
                         [r + '_gene' for r in regions] + \
                         [e + '_del' for e in real_erosions + effective_erosions] + \
                         [b + '_insertion' for b in boundaries + effective_boundaries] + \
                         [r + '_gl_seq' for r in regions] + \
                         [r + '_per_gene_support' for r in regions]
linekeys['per_seq'] = ['seqs', 'unique_ids', 'indelfos', 'mut_freqs'] + \
                      [r + '_qr_seqs' for r in regions] + \
                      ['aligned_' + r + '_seqs' for r in regions] + \
                      functional_columns
linekeys['hmm'] = ['logprob', 'errors']
linekeys['sw'] = ['k_v', 'k_d', 'all_matches', 'padlefts', 'padrights']
linekeys['extra'] = ['invalid', ]
linekeys['simu'] = ['reco_id', ]
all_linekeys = set([k for cols in linekeys.values() for k in cols])

# keys that are added by add_implicit_info()
implicit_linekeys = set(['naive_seq', 'cdr3_length', 'codon_positions', 'lengths', 'regional_bounds', 'invalid'] + \
                       [r + '_gl_seq' for r in regions] + \
                       ['mut_freqs', ] + functional_columns + [r + '_qr_seqs' for r in regions] + ['aligned_' + r + '_seqs' for r in regions])

# ----------------------------------------------------------------------------------------
annotation_headers = ['unique_ids', 'v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'mut_freqs', 'seqs', 'naive_seq', 'indelfos'] \
                     + ['aligned_' + r + '_seqs' for r in regions] \
                     + [r + '_per_gene_support' for r in regions] \
                     + [e + '_del' for e in real_erosions + effective_erosions] + [b + '_insertion' for b in boundaries + effective_boundaries] \
                     + functional_columns
sw_cache_headers = ['k_v', 'k_d', 'padlefts', 'padrights', 'all_matches', 'mut_freqs']
partition_cachefile_headers = ('unique_ids', 'logprob', 'naive_seq', 'naive_hfrac', 'errors')  # these have to match whatever bcrham is expecting

# ----------------------------------------------------------------------------------------
def convert_to_presto_headers(line):
    """ convert <line> to presto csv format """
    if len(line['unique_ids']) > 1:  # has to happen *before* utils.get_line_for_output()
        print line['unique_ids']
        raise Exception('multiple seqs not handled for presto output')

    presto_line = {}
    for head, phead in presto_headers.items():
        if head == 'aligned_v_plus_unaligned_dj':
            presto_line[phead] = line['aligned_v_seqs'][0] + line['vd_insertion'] + line['d_qr_seqs'][0] + line['dj_insertion'] + line['j_qr_seqs'][0]
        elif head == 'unique_ids':
            presto_line[phead] = line[head][0]
        else:
            presto_line[phead] = line[head]

    return presto_line

# ----------------------------------------------------------------------------------------
# these are the top 10 v and d genes, and top six js, from mishmash.csv. Restricting to these makes testing much more stable and much faster.
test_only_genes = 'IGHV3-53*02:IGHV3-7*01:IGHV3-7*03:IGHV3-23D*01:IGHV4-61*08:IGHV3-66*03:IGHV5-51*02:IGHV3-66*01:IGHV3-23D*02:IGHV1-2*02:' + \
                  'IGHD2-2*03:IGHD6-6*01:IGHD2-8*01:IGHD2-21*02:IGHD3-10*01:IGHD3-16*02:IGHD3-22*01:IGHD6-19*01:IGHD2-15*01:' + \
                  'IGHJ6*03:IGHJ6*02:IGHJ2*01:IGHJ1*01:IGHJ5*02:IGHJ4*02'

# tuples with the column and its dependencies mashed together
# (first entry is the column of interest, and it depends upon the following entries)
column_dependency_tuples = []
for column, deps in column_dependencies.iteritems():
    tmp_list = [column]
    tmp_list.extend(deps)
    column_dependency_tuples.append(tuple(tmp_list))

# ----------------------------------------------------------------------------------------
def get_parameter_fname(column=None, deps=None, column_and_deps=None):
    """ return the file name in which we store the information for <column>. Either pass in <column> and <deps> *or* <column_and_deps> """
    if column == 'all':
        return 'all-probs.csv'
    if column_and_deps is None:
        if column is None or deps is None:
            raise Exception('you have to either pass <column_and_deps>, or else pass both <column> and <deps>')
        column_and_deps = [column]
        column_and_deps.extend(deps)
    outfname = 'probs.csv'
    for ic in column_and_deps:
        outfname = ic + '-' + outfname
    return outfname

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
Colors['reverse_video'] = '\033[7m'
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
def color_mutants(ref_seq, seq, print_result=False, extra_str='', ref_label='', post_str='', print_hfrac=False, print_isnps=False):
    assert len(ref_seq) == len(seq)
    return_str = ''
    isnps = []
    for inuke in range(len(seq)):
        if inuke >= len(ref_seq) or seq[inuke] == ref_seq[inuke]:
            return_str += seq[inuke]
        else:
            return_str += color('red', seq[inuke])
            isnps.append(inuke)
    if print_result:
        print '%s%s%s' % (extra_str, ref_label, ref_seq)
        print '%s%s%s%s' % (extra_str, ' '*len(ref_label), return_str, post_str),
        if print_hfrac:
            print '   hfrac %.3f' % hamming_fraction(ref_seq, seq),
        if print_isnps:
            print '   snps at: %s' % ' '.join([str(i) for i in isnps]),
        print ''
    return return_str

# ----------------------------------------------------------------------------------------
def plural_str(pstr, count):
    if count == 1:
        return pstr
    else:
        return pstr + 's'

# ----------------------------------------------------------------------------------------
def summarize_gene_name(gene):
    region = get_region(gene)
    primary_version, sub_version, allele = split_gene(gene)
    return ' '.join([region, primary_version, sub_version, allele])

# ----------------------------------------------------------------------------------------
def color_gene(gene, width=None):
    """ color gene name (and remove extra characters), eg IGHV3-h*01 --> hv3-h1 """
    chain = get_chain(gene)
    region = get_region(gene)
    primary_version, sub_version, allele = split_gene(gene)

    n_chars = len(chain + region + primary_version)  # number of non-special characters
    return_str = color('purple', chain) + color('red', region) + color('purple', primary_version)
    if sub_version is not None:
        n_chars += 1 + len(sub_version)
        return_str += '-' + color('purple', sub_version)
    n_chars += len(allele)
    return_str += color('yellow', allele)
    if width is not None:
        return_str += (width - n_chars) * ' '
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

#----------------------------------------------------------------------------------------
def remove_gaps(seq):
    return seq.translate(None, ''.join(gap_chars))

# ----------------------------------------------------------------------------------------
def both_codons_ok(chain, seq, positions, debug=False, extra_str=''):
    both_ok = True
    for region, codon in conserved_codons[chain].items():
        both_ok &= codon_ok(codon, seq, positions[region], debug=debug, extra_str=extra_str)
    return both_ok

# ----------------------------------------------------------------------------------------
def codon_ok(codon, seq, position, debug=False, extra_str=''):
    if len(seq) < position + 3:
        if debug:
            print '%ssequence length %d less than %s position %d + 3' % (extra_str, len(seq), codon, position)
        return False

    if seq[position : position + 3] not in codon_table[codon]:
        if debug:
            print '%s%s codon %s not among expected codons (%s)' % (extra_str, codon, seq[position : position + 3], ' '.join(codon_table[codon]))
        return False

    return True

#----------------------------------------------------------------------------------------
def check_a_bunch_of_codons(codon, seqons, extra_str='', debug=False):  # seqons: list of (seq, pos) pairs
    """ check a list of sequences, and keep track of some statistics """
    n_total, n_ok, n_too_short, n_bad_codons = 0, 0, 0, 0
    for seq, pos in seqons:
        n_total += 1
        if len(seq) < pos + 3:
            n_too_short += 1
        elif codon_ok(codon, seq, pos):
            n_ok += 1
        else:
            n_bad_codons += 1

    if debug:
        print '%s%d %s positions:' % (extra_str, n_total, codon),
        if n_ok > 0:
            print '  %d ok' % n_ok,
        if n_too_short > 0:
            print '  %d too short' % n_too_short,
        if n_bad_codons > 0:
            print '  %d mutated' % n_bad_codons,
        print ''

#----------------------------------------------------------------------------------------
def is_there_a_stop_codon(seq, cyst_position, debug=False):
    """
    Make sure there is no in-frame stop codon, where frame is inferred from <cyst_position>.
    Returns True if no stop codon is found
    """
    if cyst_position >= len(seq):
        if debug:
            print '      not sure if there\'s a stop codon (invalid cysteine position)'
        return True  # not sure if there is one, since we have to way to establish the frame
    # jump leftward in steps of three until we reach the start of the sequence
    ipos = cyst_position
    while ipos > 2:
        ipos -= 3
    # ipos should now bet the index of the start of the first complete codon
    while ipos + 2 < len(seq):  # then jump forward in steps of three bases making sure none of them are stop codons
        codon = seq[ipos : ipos + 3]
        if codon in codon_table['stop']:
            if debug:
                print '      stop codon %s at %d in %s' % (codon, ipos, seq)
            return True
        ipos += 3

    return False  # no stop codon

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
def reset_effective_erosions_and_effective_insertions(glfo, padded_line, aligned_gl_seqs=None, debug=False):  # , padfo=None
    """
    Ham does not allow (well, no longer allows) v_5p and j_3p deletions -- we instead pad sequences with Ns.
    This means that the info we get from ham always has these effective erosions set to zero, but for downstream
    things we sometimes want to know where the reads stopped (e.g. if we want to mimic them in simulation).
    Note that these effective erosion values will be present in the parameter dir, but are *not* incorporated into
    the hmm yaml files.
    """

    line = copy.deepcopy(padded_line)
    remove_all_implicit_info(line)

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

    # arbitrarily use the zeroth sequence (in principle v_5p and j_3p should be per-sequence, not per-rearrangement... but that'd be a mess to implement, since the other deletions are per-rearrangement)
    tmpiseq = 0  # NOTE this is pretty hackey: we just use the values from the first sequence. But it's actually not that bad -- we can either have some extra pad Ns showing, or chop of some bases.
    trimmed_seq = trimmed_seqs[tmpiseq]
    final_fv_insertion = final_insertions[tmpiseq]['fv']
    final_jf_insertion = final_insertions[tmpiseq]['jf']
    fv_insertion_to_remove = insertions_to_remove[tmpiseq]['fv']
    jf_insertion_to_remove = insertions_to_remove[tmpiseq]['jf']
    line['v_5p_del'] = find_first_non_ambiguous_base(trimmed_seq)
    line['j_3p_del'] = len(trimmed_seq) - find_last_non_ambiguous_base_plus_one(trimmed_seq)

    for iseq in range(len(line['seqs'])):
        line['seqs'][iseq] = trimmed_seqs[iseq][line['v_5p_del'] : ]
        if line['j_3p_del'] > 0:
            line['seqs'][iseq] = line['seqs'][iseq][ : -line['j_3p_del']]

        if line['indelfos'][iseq]['reversed_seq'] != '':
            rseq = line['indelfos'][iseq]['reversed_seq']
            rseq = rseq[len(fv_insertion_to_remove) + line['v_5p_del'] : ]
            if len(jf_insertion_to_remove) + line['j_3p_del'] > 0:
                rseq = rseq[ : -(len(jf_insertion_to_remove) + line['j_3p_del'])]
            line['indelfos'][iseq]['reversed_seq'] = rseq

    if debug:
        print '     fv %d   v_5p %d   j_3p %d   jf %d    %s' % (len(fv_insertion_to_remove), line['v_5p_del'], line['j_3p_del'], len(jf_insertion_to_remove), line['seqs'][0])

    line['fv_insertion'] = final_fv_insertion
    line['jf_insertion'] = final_jf_insertion

    # if padfo is None:
    #     line['padlefts'], line['padrights'] = [0 for _ in range(len(line['seqs']))], [0 for _ in range(len(line['seqs']))]
    # else:
    #     line['padlefts'], line['padrights'] = [padfo[uid]['padded']['padleft'] for uid in line['unique_ids']], [padfo[uid]['padded']['padright'] for uid in line['unique_ids']]

    add_implicit_info(glfo, line, aligned_gl_seqs=aligned_gl_seqs)

    return line

# ----------------------------------------------------------------------------------------
def add_qr_seqs(line):
    """ Add [vdj]_qr_seq, i.e. the sections of the query sequence which are assigned to each region. """

    starts = {}
    starts['v'] = len(line['fv_insertion'])
    starts['d'] = starts['v'] + len(line['v_gl_seq']) + len(line['vd_insertion'])
    starts['j'] = starts['d'] + len(line['d_gl_seq']) + len(line['dj_insertion'])

    def get_single_qr_seq(region, seq):
        return seq[starts[region] : starts[region] + len(line[region + '_gl_seq'])]

    for region in regions:
        line[region + '_qr_seqs'] = [get_single_qr_seq(region, seq) for seq in line['seqs']]

# ----------------------------------------------------------------------------------------
def add_functional_info(chain, line):
    def get_val(ftype, iseq):
        if ftype == 'mutated_invariants':
            return not both_codons_ok(chain, line['seqs'][iseq], line['codon_positions'])
        elif ftype == 'in_frames':
            return line['cdr3_length'] % 3 == 0
        elif ftype == 'stops':
            return is_there_a_stop_codon(line['seqs'][iseq], line['codon_positions']['v'])
        else:
            assert False

    for fc in functional_columns:
        line[fc] = [get_val(fc, iseq) for iseq in range(len(line['seqs']))]

# ----------------------------------------------------------------------------------------
def remove_all_implicit_info(line):
    for col in implicit_linekeys:
        if col in line:
            del line[col]

# ----------------------------------------------------------------------------------------
def process_per_gene_support(line, debug=False):
    for region in regions:
        if debug:
            print region
        support = OrderedDict()
        logtotal = float('-inf')
        for gene, logprob in line[region + '_per_gene_support'].items():
            support[gene] = logprob
            logtotal = add_in_log_space(logtotal, logprob)

        for gene in support:
            if debug:
                print '   %5.2f     %5.2f   %s' % (support[gene], math.exp(support[gene] - logtotal), color_gene(gene))
            support[gene] = math.exp(support[gene] - logtotal)

        if len(support.keys()) > 0 and support.keys()[0] != line[region + '_gene']:
            print '   WARNING best-supported gene %s not same as viterbi gene %s' % (color_gene(support.keys()[0]), color_gene(line[region + '_gene']))

        line[region + '_per_gene_support'] = support

# ----------------------------------------------------------------------------------------
def add_implicit_info(glfo, line, existing_implicit_keys=None, aligned_gl_seqs=None, debug=False):
    """ Add to <line> a bunch of things that are initially only implicit. """

    # check for existing and unexpected keys
    if existing_implicit_keys is not None:  # remove any existing implicit keys, keeping track of their values to make sure they're the same afterwards
        pre_existing_info = {}
        for ekey in existing_implicit_keys:
            pre_existing_info[ekey] = copy.deepcopy(line[ekey])
            del line[ekey]
    initial_keys = set(line.keys())  # keep track of the keys that are in <line> to start with (so we know which ones we added)
    if len(initial_keys - all_linekeys) > 0:  # make sure there aren't any extra keys to start with
        raise Exception('unexpected keys %s' % ' '.join(initial_keys - all_linekeys))

    # add the regional germline seqs and their lengths
    line['lengths'] = {}  # length of each match (including erosion)
    for region in regions:
        uneroded_gl_seq = glfo['seqs'][region][line[region + '_gene']]
        del_5p = line[region + '_5p_del']
        del_3p = line[region + '_3p_del']
        length = len(uneroded_gl_seq) - del_5p - del_3p  # eroded length
        line[region + '_gl_seq'] = uneroded_gl_seq[del_5p : del_5p + length]
        line['lengths'][region] = length

    # add codon-related stuff
    line['codon_positions'] = {}
    for region, codon in conserved_codons[glfo['chain']].items():
        eroded_gl_pos = glfo[codon + '-positions'][line[region + '_gene']] - line[region + '_5p_del']
        if region == 'v':
            line['codon_positions'][region] = eroded_gl_pos + len(line['f' + region + '_insertion'])
        elif region == 'j':
            line['codon_positions'][region] = eroded_gl_pos + len(line['fv_insertion']) + line['lengths']['v'] + len(line['vd_insertion']) + line['lengths']['d'] + len(line['dj_insertion'])
        else:
            assert False
    line['cdr3_length'] = line['codon_positions']['j'] - line['codon_positions']['v'] + 3  # i.e. first base of cysteine to last base of tryptophan inclusive

    # add naive seq stuff
    line['naive_seq'] = line['fv_insertion'] + line['v_gl_seq'] + line['vd_insertion'] + line['d_gl_seq'] + line['dj_insertion'] + line['j_gl_seq'] + line['jf_insertion']
    start, end = {}, {}  # add naive seq bounds for each region (could stand to make this more concise)
    start['v'] = 0  # holy fuck, look at that, I start at zero here, but at the end of the fv insertion in add_qr_seqs(). Scary!
    end['v'] = start['v'] + len(line['fv_insertion'] + line['v_gl_seq'])  # base just after the end of v
    start['d'] = end['v'] + len(line['vd_insertion'])
    end['d'] = start['d'] + len(line['d_gl_seq'])
    start['j'] = end['d'] + len(line['dj_insertion'])
    end['j'] = start['j'] + len(line['j_gl_seq'] + line['jf_insertion'])
    line['regional_bounds'] = {r : (start[r], end[r]) for r in regions}

    # add regional query seqs
    add_qr_seqs(line)

    add_functional_info(glfo['chain'], line)

    line['mut_freqs'] = [hamming_fraction(line['naive_seq'], mature_seq) for mature_seq in line['seqs']]

    # set validity (alignment addition can also set invalid)  # TODO clean up this checking stuff
    line['invalid'] = False
    seq_length = len(line['seqs'][0])  # they shouldn't be able to be different lengths
    for chkreg in regions:
        if start[chkreg] < 0 or end[chkreg] < 0 or end[chkreg] < start[chkreg] or end[chkreg] > seq_length:
            line['invalid'] = True
    if end['j'] != seq_length:
        line['invalid'] = True
    if line['cdr3_length'] < 6:  # i.e. if cyst and tryp overlap  NOTE six is also hardcoded in waterer
        line['invalid'] = True

    # add alignment info
    if aligned_gl_seqs is None:
        add_dummy_alignments(line)
    else:
        add_alignments(glfo, aligned_gl_seqs, line, debug)

    # make sure we added exactly what we expected to
    new_keys = set(line.keys()) - initial_keys
    if len(new_keys - implicit_linekeys) > 0 or len(implicit_linekeys - new_keys) > 0:
        print ''
        print '           new   %s' % ' '.join(sorted(new_keys))
        print '      implicit   %s' % ' '.join(sorted(implicit_linekeys))
        print 'new - implicit:  %s' % (' '.join(sorted(new_keys - implicit_linekeys)))
        print 'implicit - new:  %s' % (' '.join(sorted(implicit_linekeys - new_keys)))
        print ''
        raise Exception('column/key problems (see above)')

    # make sure that any pre-existing implicit info matches what we just added
    if existing_implicit_keys is not None:
        for ekey in existing_implicit_keys:
            if pre_existing_info[ekey] != line[ekey]:
                print '  WARNING pre-existing info\n    %s\n    doesn\'t match new info\n    %s\n    for %s in %s' % (pre_existing_info[ekey], line[ekey], ekey, line['unique_ids'])

# ----------------------------------------------------------------------------------------
def print_true_events(glfo, reco_info, line, print_uid=False):
    """ print the true events which contain the seqs in <line> """
    # true_naive_seqs = []
    for uids in get_true_partition(reco_info, ids=line['unique_ids']):  # make a multi-seq line that has all the seqs from this clonal family
        multiline = synthesize_multi_seq_line(uids, reco_info)
        print_reco_event(glfo['seqs'], multiline, extra_str='    ', label='true:', print_uid=print_uid)
        # true_naive_seqs.append(multiline['naive_seq'])

    # print '\ntrue vs inferred naive sequences:'
    # for tseq in true_naive_seqs:
    #     color_mutants(tseq, line['naive_seq'], print_result=True, print_hfrac=True, ref_label='true ')
    # print ''

# ----------------------------------------------------------------------------------------
def add_gaps_ignoring_color_characters(colored_seq, ipos, gapstr):
    """ add <gapstr> to <colored_seq> at position <ipos>, where <ipos> is the index ignoring characters for bash colors """
    # NOTE don't really need this any more, I realized I can just use lists of strings to keep track of the indices
    iwithout = 0  # i.e. without bash color characters
    for iwith in range(len(colored_seq)):
        ch = colored_seq[iwith]
        if iwithout == ipos:  # first position after ipos
            break
        if ch in expected_characters:  # i.e. if it isn't a character having to do with bash colors
            # if iwithout == ipos:  # first position after ipos
            #     break
            iwithout += 1

    # print repr(colored_seq[:iwith])
    # print repr(colored_seq[iwith:])
    return colored_seq[:iwith] + gapstr + colored_seq[iwith:]

# ----------------------------------------------------------------------------------------
def print_reco_event(germlines, line, one_line=False, extra_str='', label='', print_uid=False):
    for iseq in range(len(line['unique_ids'])):
        print_seq_in_reco_event(germlines, line, iseq, extra_str=extra_str, label=(label if iseq==0 else ''), one_line=(iseq>0), print_uid=print_uid)

# ----------------------------------------------------------------------------------------
def print_seq_in_reco_event(germlines, original_line, iseq, extra_str='', label='', one_line=False, print_uid=False):
    """
    Print ascii summary of recombination event and mutation.
    If <one_line>, then skip the germline lines, and only print the final_seq line.
    """
    line = copy.deepcopy(original_line)  # copy that we can modify without changing <line>

    lseq = line['seqs'][iseq]
    indelfo = None if line['indelfos'][iseq]['reversed_seq'] == '' else line['indelfos'][iseq]
    reverse_indels = False  # for inferred sequences, we want to un-reverse the indels that we previously reversed in smith-waterman
    if indelfo is not None:
        if indelfo['reversed_seq'] == lseq:  # if <line> has the reversed sequence in it, then this is an inferred <line>, i.e. we removed the info, then passed the reversed sequence to the sw/hmm, so we need to reverse the indels now in order to get a sequence with indels in it
            reverse_indels = True
        if len(indelfo['indels']) > 1:
            print 'WARNING multiple indels not really handled'
        add_indels_to_germline_strings(line, indelfo)

    # ----------------------------------------------------------------------------------------
    def process_position(original, final):
        if original not in expected_characters or final not in expected_characters:
            raise Exception('one of %s %s not among expected characters' % (original, final))

        if original in ambiguous_bases or final in ambiguous_bases:  # don't count Ns in the total
            return final

        if original != final:
            return color('red', final)

        return final

    # build up the query sequence line, including colors for mutations and conserved codons
    j_right_extra = ''  # portion of query sequence to right of end of the j match
    n_inserted = 0
    final_seq_list = []
    for inuke in range(len(lseq)):
        if indelfo is not None:
            lastfo = indelfo['indels'][-1]  # if the "last" (arbitrary but necessary ordering) indel starts here
              # if we're at the position that the insertion started at (before we removed it)
        if indelfo is not None and lastfo['type'] == 'insertion':
            if reverse_indels and inuke == lastfo['pos']:
                final_seq_list.append(lastfo['seqstr'])  # put the insertion back into the query sequence
                n_inserted += len(lastfo['seqstr'])
            elif not reverse_indels and inuke - lastfo['pos'] >= 0 and inuke - lastfo['pos'] < lastfo['len']:
                final_seq_list.append(lseq[inuke])
                n_inserted += 1
                continue
        if indelfo is not None and lastfo['type'] == 'deletion':
            if reverse_indels and inuke - lastfo['pos'] >= 0 and inuke - lastfo['pos'] < lastfo['len']:  # if we're within the bases that we added to make up for the deletionlen
                final_seq_list.append(color('light_blue', '*'))
                continue
            elif not reverse_indels and inuke == lastfo['pos']:
                final_seq_list.append(lastfo['len'] * color('light_blue', '*'))
                n_inserted = - lastfo['len']

        new_nuke = ''
        key = None
        ilocal = inuke
        if indelfo is not None and reverse_indels:
            ilocal += n_inserted
        if indelfo is not None and not reverse_indels and lastfo['type'] == 'deletion':
            ilocal -= n_inserted
        if ilocal < len(line['fv_insertion']):  # haven't got to start of v match yet, so just add on the query seq nuke
            pass
        else:
            ilocal -= len(line['fv_insertion'])
            if ilocal < line['lengths']['v']:
                key = 'v_gl_seq'
            else:
                ilocal -= line['lengths']['v']
                if ilocal < len(line['vd_insertion']):
                    key = 'vd_insertion'
                else:
                    ilocal -= len(line['vd_insertion'])
                    if ilocal < line['lengths']['d']:
                        key = 'd_gl_seq'
                    else:
                        ilocal -= line['lengths']['d']
                        if ilocal < len(line['dj_insertion']):
                            key = 'dj_insertion'
                        else:
                            ilocal -= len(line['dj_insertion'])
                            if ilocal < line['lengths']['j']:
                                key = 'j_gl_seq'
                            else:
                                j_right_extra += ' '

        if key is None:
            original = lseq[inuke]  # dummy value
        else:
            original = line[key][ilocal]
        new_nuke = process_position(original, lseq[inuke])

        for region, pos in line['codon_positions'].items():  # reverse video for the conserved codon positions
            if indelfo is not None and not reverse_indels:
                pos += n_inserted
            if inuke >= pos and inuke < pos + 3:
                new_nuke = '\033[7m' + new_nuke + '\033[m'

        final_seq_list.append(new_nuke)

    # check if there isn't enough space for dots in the vj line
    no_space = False
    interior_length = len(line['vd_insertion']) + len(line['d_gl_seq']) + len(line['dj_insertion'])  # length of the portion of the vj line that is normally taken up by dots (and spaces)
    if line['v_3p_del'] + line['j_5p_del'] > interior_length:
        no_space = True

    gaps_to_add = 0
    if no_space:
        v_3p_del_str = '.' + str(line['v_3p_del']) + '.'  # if line['v_3p_del'] > 0 else ' '
        j_5p_del_str = '.' + str(line['j_5p_del']) + '.'  # if line['j_5p_del'] > 0 else ' '
        extra_space_because_of_fixed_nospace = max(0, interior_length - len(v_3p_del_str + j_5p_del_str))

        gap_insertion_point = len(line['fv_insertion'] + line['v_gl_seq'])
        gaps_to_add = len(v_3p_del_str + j_5p_del_str) - interior_length
        # gapstr = gaps_to_add * color('blue', '-')
        # final_seq = add_gaps_ignoring_color_characters(final_seq, gap_insertion_point, gapstr)
        final_seq_list = final_seq_list[:gap_insertion_point] + gaps_to_add * [color('blue', '-'), ] + final_seq_list[gap_insertion_point:]
    else:
        v_3p_del_str = '.' * line['v_3p_del']
        j_5p_del_str = '.' * line['j_5p_del']
        extra_space_because_of_fixed_nospace = 0

    eroded_seqs_dots = {}
    eroded_seqs_dots['v'] = line['v_gl_seq'] + v_3p_del_str
    eroded_seqs_dots['d'] = '.'*line['d_5p_del'] + line['d_gl_seq'] + '.'*line['d_3p_del']
    eroded_seqs_dots['j'] = j_5p_del_str + line['j_gl_seq'] + '.'*line['j_3p_del']

    v_5p_del_str = '.'*line['v_5p_del']
    if line['v_5p_del'] > 50:
        v_5p_del_str = '...' + str(line['v_5p_del']) + '...'

    insert_line = ' '*len(line['fv_insertion']) + ' '*line['lengths']['v'] + color('blue', '-')*gaps_to_add
    insert_line += ' '*len(v_5p_del_str)
    insert_line += line['vd_insertion']
    insert_line += ' ' * line['lengths']['d']
    insert_line += line['dj_insertion']
    insert_line += ' ' * line['lengths']['j']
    insert_line += j_right_extra
    insert_line += ' ' * line['j_3p_del']  # no damn idea why these need to be commented out for some cases in the igblast parser...
    # insert_line += ' '*len(line['jf_insertion'])

    germline_d_start = len(line['fv_insertion']) + line['lengths']['v'] + len(line['vd_insertion']) - line['d_5p_del']
    germline_d_end = germline_d_start + len(germlines['d'][line['d_gene']])
    d_line = ' ' * germline_d_start + color('blue', '-')*gaps_to_add
    d_line += ' '*len(v_5p_del_str)
    d_line += eroded_seqs_dots['d']
    d_line += ' ' * (len(line['j_gl_seq']) + len(line['dj_insertion']) - line['d_3p_del'])
    d_line += j_right_extra
    d_line += ' ' * line['j_3p_del']
    # d_line += ' '*len(line['jf_insertion'])

    vj_line = ' ' * len(line['fv_insertion'])
    vj_line += v_5p_del_str
    vj_line += eroded_seqs_dots['v'] + '.'*extra_space_because_of_fixed_nospace
    germline_v_end = len(line['fv_insertion']) + len(line['v_gl_seq']) + line['v_3p_del'] - 1  # position in the query sequence at which we find the last base of the v match. NOTE we subtract off the v_5p_del because we're *not* adding dots for that deletion (it's just too long)
    germline_j_start = germline_d_end + 1 - line['d_3p_del'] + len(line['dj_insertion']) - line['j_5p_del']
    vj_line += ' ' * (germline_j_start - germline_v_end - 2)
    vj_line += eroded_seqs_dots['j']
    vj_line += j_right_extra
    # vj_line += ' '*len(line['jf_insertion'])

    insert_line = color_chars(ambiguous_bases + ['*', ], 'light_blue', insert_line)
    d_line = color_chars(ambiguous_bases + ['*', ], 'light_blue', d_line)
    vj_line = color_chars(ambiguous_bases + ['*', ], 'light_blue', vj_line)

    chain = get_chain(line['v_gene'])  # kind of hackey
    dont_show_d_stuff = chain != 'h' and line['lengths']['d'] == 0 and len(line['vd_insertion']) == 0

    if print_uid:
        extra_str += '%20s ' % line['unique_ids'][iseq]
    out_str_list = []
    # insert, d, and vj lines
    if not one_line:
        out_str_list.append('%s    %s   insert%s\n' % (extra_str, insert_line, '' if dont_show_d_stuff else 's'))
        if label != '':
            out_str_list[-1] = extra_str + label + out_str_list[-1][len(extra_str + label) :]
        if not dont_show_d_stuff:
            out_str_list.append('%s    %s   %s\n' % (extra_str, d_line, color_gene(line['d_gene'])))
        out_str_list.append('%s    %s   %s %s\n' % (extra_str, vj_line, color_gene(line['v_gene']), color_gene(line['j_gene'])))

    # if indelfo is not None:
    #     for ii in range(len(indelfo['indels'])):
    #         idl = indelfo['indels'][ii]
    #         if ii > 0:
    #             extra_str += '\nxxx'
    #     extra_str += ' %10s: %2d bases at %3d'  % (idl['type'], idl['len'], idl['pos'])

    # then query sequence
    v_5p_del_space_str = ' '*len(v_5p_del_str)
    j_3p_del_space_str = ' ' * line['j_3p_del']
    final_seq = ''.join(final_seq_list)
    final_seq = v_5p_del_space_str + final_seq + j_3p_del_space_str
    final_seq = color_chars(ambiguous_bases + ['*', ], 'light_blue', final_seq)
    out_str_list.append(extra_str)
    # if print_uid:
    #     extra_str += '%20s' % line['unique_id']
    out_str_list.append('    %s' % final_seq)
    out_str_list.append('   %4.2f mut' % line['mut_freqs'][iseq])
    if 'logprob' in line:
        out_str_list.append('     %8.2f  logprob' % line['logprob'])
    out_str_list.append('\n')

    print ''.join(out_str_list),

    # if copy_of_original_line != original_line:
    #     for k in copy_of_original_line:
    #         if copy_of_original_line[k] != original_line[k]:
    #             print k, copy_of_original_line[k], original_line[k]
    #     raise Exception('')

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

# ----------------------------------------------------------------------------------------
def get_chain(inputstr):
    """ return chain weight/locus given gene of file name """
    if inputstr[:2] == 'IG':  # it's a gene name
        chain = inputstr[2:3].lower()
    elif inputstr in glutils.all_glfo_fasta_fnames():  # it's a file name
        chain = inputstr[2:3]
    else:
        raise Exception('couldn\'t figure out if %s was a gene or file name' % inputstr)
    if chain not in chains:
        raise Exception('couldn\'t get chain from input string %s' % inputstr)
    return chain

# ----------------------------------------------------------------------------------------
def get_region(inputstr):
    """ return v, d, or j of gene or gl fname """
    if inputstr[:2] == 'IG':  # it's a gene name
        region = inputstr[3:4].lower()
    elif inputstr in glutils.all_glfo_fasta_fnames():  # it's a file name
        region = inputstr[3:4]
    else:
        raise Exception('couldn\'t figure out if %s was a gene or file name' % inputstr)
    if region not in regions:
        raise Exception('couldn\'t get region from input string %s' % inputstr)
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
    # make sure IG[HKL][VDJ] is at the start, and there's a *
    if gene[:4] != 'IG' + get_chain(gene).upper() + get_region(gene).upper():
        raise Exception('unexpected string in gene name %s' % gene)
    if gene.count('*') != 1:
        raise Exception('expected exactly 1 \'*\' in %s but found %d' % (gene, gene.count('*')))

    if '-' in gene and gene.find('-') < gene.find('*'):  # Js (and a few Vs) don't have sub versions
        primary_version = gene[4 : gene.find('-')]  # the bit between the IG[HKL][VDJ] and the first dash (sometimes there's a second dash as well)
        sub_version = gene[gene.find('-') + 1 : gene.find('*')]  # the bit between the first dash and the star
        allele = gene[gene.find('*') + 1 : ]  # the bit after the star
        if gene != 'IG' + get_chain(gene).upper() + get_region(gene).upper() + primary_version + '-' + sub_version + '*' + allele:
            raise Exception('couldn\'t build gene name %s from %s %s %s' % (gene, primary_version, sub_version, allele))
    else:
        primary_version = gene[4 : gene.find('*')]  # the bit between the IG[HKL][VDJ] and the star
        sub_version = None
        allele = gene[gene.find('*') + 1 : ]  # the bit after the star
        if gene != 'IG' + get_chain(gene).upper() + get_region(gene).upper() + primary_version + '*' + allele:
            raise Exception('couldn\'t build gene name %s from %s %s' % (gene, primary_version, allele))

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
def separate_into_allelic_groups(glfo, debug=False):
    allelic_groups = {r : {} for r in regions}
    for region in regions:
        for gene in glfo['seqs'][region]:
            primary_version, sub_version, allele = split_gene(gene)
            if primary_version not in allelic_groups[region]:
                allelic_groups[region][primary_version] = {}
            if sub_version not in allelic_groups[region][primary_version]:
                allelic_groups[region][primary_version][sub_version] = set()
            allelic_groups[region][primary_version][sub_version].add(gene)
    if debug:
        for r in allelic_groups:
            print r
            for p in allelic_groups[r]:
                print '    %15s' % p
                for s in allelic_groups[r][p]:
                    print '        %15s      %s' % (s, ' '.join([color_gene(g, width=12) for g in allelic_groups[r][p][s]]))
    return allelic_groups

# ----------------------------------------------------------------------------------------
def read_single_gene_count(indir, gene, expect_zero_counts=False, debug=False):
    region = get_region(gene)
    count = 0
    with opener('r')(indir + '/' + region + '_gene-probs.csv') as infile:  # NOTE this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
        reader = csv.DictReader(infile)
        for line in reader:
            if line[region + '_gene'] == gene:
                count = int(line['count'])
                break

    if count == 0 and not expect_zero_counts:
        print '          %s %s not found in %s_gene-probs.csv, returning zero' % (color('red', 'warning'), gene, region)

    if debug:
        print '      %d observations of %s' % (count, color_gene(gene))

    return count

# ----------------------------------------------------------------------------------------
def read_overall_gene_probs(indir, only_gene=None, normalize=True, expect_zero_counts=False, debug=False):
    """
    Return the observed counts/probabilities of choosing each gene version.
    If <normalize> then return probabilities
    If <only_gene> is specified, just return the prob/count for that gene
    """
    counts, probs = {r : {} for r in regions}, {r : {} for r in regions}
    for region in regions:
        total = 0
        with opener('r')(indir + '/' + region + '_gene-probs.csv') as infile:  # NOTE this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
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

    if only_gene is not None and only_gene not in counts[get_region(only_gene)]:
        if not expect_zero_counts:
            print '      WARNING %s not found in overall gene probs, returning zero' % only_gene
        if normalize:
            return 0.0
        else:
            return 0

    if only_gene is None:
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
def find_replacement_genes(param_dir, min_counts, gene_name=None, debug=False, all_from_region=''):
    if gene_name is not None:  # if you specify <gene_name> you shouldn't specify <all_from_region>
        assert all_from_region == ''
        region = get_region(gene_name)
    else:  # and vice versa
        assert all_from_region in regions
        assert min_counts == -1
        region = all_from_region
    lists = OrderedDict()  # we want to try alleles first, then primary versions, then everything and it's mother
    lists['allele'] = []  # list of genes that are alleles of <gene_name>
    lists['primary_version'] = []  # same primary version as <gene_name>
    lists['all'] = []  # give up and return everything
    with opener('r')(param_dir + '/' + region + '_gene-probs.csv') as infile:  # NOTE note this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
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

    raise Exception('couldn\'t find genes for %s in %s' % (gene_name, param_dir))

    # print '    \nWARNING return default gene %s \'cause I couldn\'t find anything remotely resembling %s' % (color_gene(hackey_default_gene_versions[region]), color_gene(gene_name))
    # return hackey_default_gene_versions[region]

# ----------------------------------------------------------------------------------------
def hamming_distance(seq1, seq2, extra_bases=None, return_len_excluding_ambig=False):
    if len(seq1) != len(seq2):
        raise Exception('unequal length sequences %d %d' % (len(seq1), len(seq2)))
    if len(seq1) == 0:
        if return_len_excluding_ambig:
            return 0, 0
        else:
            return 0

    alphabet = nukes + ambiguous_bases  # I think I don't want to allow gap characters here

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

    if return_len_excluding_ambig:
        return distance, len_excluding_ambig
    else:
        return distance

# ----------------------------------------------------------------------------------------
def hamming_fraction(seq1, seq2, extra_bases=None):
    distance, len_excluding_ambig = hamming_distance(seq1, seq2, extra_bases=extra_bases, return_len_excluding_ambig=True)
    fraction = 0.
    if len_excluding_ambig > 0:
        fraction = distance / float(len_excluding_ambig)
    return fraction

# ----------------------------------------------------------------------------------------
def subset_sequences(line, iseq, restrict_to_region):
    naive_seq = line['naive_seq']  # NOTE this includes the fv and jf insertions
    muted_seq = line['seqs'][iseq]
    if restrict_to_region != '':  # NOTE this is very similar to code in performanceplotter. I should eventually cut it out of there and combine them, but I'm nervous a.t.m. because of all the complications there of having the true *and* inferred sequences so I'm punting
        bounds = line['regional_bounds'][restrict_to_region]
        naive_seq = naive_seq[bounds[0] : bounds[1]]
        muted_seq = muted_seq[bounds[0] : bounds[1]]
    return naive_seq, muted_seq

# ----------------------------------------------------------------------------------------
def get_n_muted(line, iseq, restrict_to_region=''):
    naive_seq, muted_seq = subset_sequences(line, iseq, restrict_to_region)
    return hamming_distance(naive_seq, muted_seq)

# ----------------------------------------------------------------------------------------
def get_mutation_rate(line, iseq, restrict_to_region=''):
    naive_seq, muted_seq = subset_sequences(line, iseq, restrict_to_region)
    return hamming_fraction(naive_seq, muted_seq)

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
def prep_dir(dirname, wildlings=None, subdirs=None):
    """
    Make <dirname> if it d.n.e.
    Also, if shell glob <wildling> is specified, remove existing files which are thereby matched.
    """

    if wildlings is None:
        wildlings = []
    elif isinstance(wildlings, basestring):  # allow to pass in just a string, instead of a list of strings
        wildlings = [wildlings, ]

    if subdirs is not None:  # clean out the subdirs first
        for subdir in subdirs:
            prep_dir(dirname + '/' + subdir, wildlings=wildlings)

    if os.path.exists(dirname):
        for wild in wildlings:
            for fname in glob.glob(dirname + '/' + wild):
                os.remove(fname)
        remaining_files = [fn for fn in os.listdir(dirname) if subdirs is not None and fn not in subdirs]
        if len(remaining_files) > 0:  # make sure there's no other files in the dir
            raise Exception('files (%s) remain in %s despite wildlings %s' % (' '.join(['\'' + fn + '\'' for fn in remaining_files]), dirname, wildlings))
    else:
        os.makedirs(dirname)

# ----------------------------------------------------------------------------------------
def useful_bool(bool_str):
    if bool_str == 'True':
        return True
    elif bool_str == 'False':
        return False
    elif bool_str == '1':
        return True
    elif bool_str == '0':
        return False
    else:
        raise Exception('couldn\'t convert \'%s\' to bool' % bool_str)

# ----------------------------------------------------------------------------------------
def process_input_line(info):
    """
    Attempt to convert all the keys and values in <info> according to the specifications in <column_configs> (e.g. splitting lists, casting to int/float, etc).
    """

    ccfg = column_configs  # shorten the name a bit

    for key in info:
        if key is None:
            raise Exception('none type key')
        if info[key] == '':
            if key in ccfg['lists']:
                info[key] = ['', ]
            continue

        convert_fcn = pass_fcn  # dummy fcn, just returns the argument
        if key in ccfg['ints']:
            convert_fcn = int
        elif key in ccfg['floats']:
            convert_fcn = float
        elif key in ccfg['bools']:
            convert_fcn = useful_bool
        elif key in ccfg['literals']:
            convert_fcn = ast.literal_eval

        if key in ccfg['lists']:
            if key in ccfg['lists-of-string-float-pairs']:  # ok, that's getting a little hackey

                def splitstrpair(pairstr):
                    pairlist = pairstr.split(':')
                    if len(pairlist) != 2:
                        raise Exception('couldn\'t split %s into two pieces with \':\'' % (pairstr))
                    return (pairlist[0], float(pairlist[1]))

                info[key] = [convert_fcn(val) for val in info[key].split(';')]
                info[key] = OrderedDict(splitstrpair(pairstr) for pairstr in info[key])
            else:
                info[key] = [convert_fcn(val) for val in info[key].split(':')]
        else:
            info[key] = convert_fcn(info[key])

# ----------------------------------------------------------------------------------------
def get_line_for_output(info):
    """ Reverse the action of process_input_line() """
    outfo = {}
    for key in info:
        str_fcn = str
        if key in column_configs['floats']:
            str_fcn = repr  # keeps it from losing precision (we only care because we want it to match expectation if we read it back in)
        if key in column_configs['lists']:
            if key in column_configs['lists-of-string-float-pairs']:  # ok, that's getting a little hackey
                outfo[key] = ';'.join([k + ':' + str_fcn(v) for k, v in info[key].items()])
            else:
                outfo[key] = ':'.join([str_fcn(v) for v in info[key]])
        else:
            outfo[key] = str_fcn(info[key])
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
def run_cmd(cmd_str, workdir):
    # print cmd_str
    # sys.exit()
    proc = Popen(cmd_str + ' 1>' + workdir + '/out' + ' 2>' + workdir + '/err', shell=True)
    return proc

# ----------------------------------------------------------------------------------------
# deal with a process once it's finished (i.e. check if it failed, and restart if so)
def finish_process(iproc, procs, n_tries, workdir, outfname, cmd_str, info=None, debug=True):
    procs[iproc].communicate()
    process_out_err('', '', extra_str='' if len(procs) == 1 else str(iproc), info=info, subworkdir=workdir, debug=debug)
    if procs[iproc].returncode == 0 and os.path.exists(outfname):  # TODO also check cachefile, if necessary
        procs[iproc] = None  # job succeeded
    elif n_tries[iproc] > 5:
        raise Exception('exceeded max number of tries for command\n    %s\nlook for output in %s' % (cmd_str, workdir))
    else:
        print '    rerunning proc %d (exited with %d' % (iproc, procs[iproc].returncode),
        if not os.path.exists(outfname):
            print ', output %s d.n.e.' % outfname,
        print ')'
        procs[iproc] = run_cmd(cmd_str, workdir)
        n_tries[iproc] += 1

# ----------------------------------------------------------------------------------------
def process_out_err(out, err, extra_str='', info=None, subworkdir=None, debug=True):
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

    if info is not None:  # keep track of how many vtb and fwd calculations the process made
        for header, variables in {'calcd' : ['vtb', 'fwd'], 'time' : ['bcrham', ]}.items():
            info[header] = {}
            theselines = [ln for ln in out.split('\n') if header + ':' in ln]
            if len(theselines) != 1:
                raise Exception('couldn\'t find \'%s\' line in:\nstdout:\n%s\nstderr:\n%s' % (header, out, err))
            words = theselines[0].split()
            try:
                for var in variables:  # convention: value corresponding to the string <var> is the word immediately vollowing <var>
                    info[header][var] = float(words[words.index(var) + 1])
            except:
                raise Exception('couldn\'t find \'%s\' line in:\nstdout:\n%s\nstderr:\n%s' % (header, out, err))

    print_str += out

    if print_str != '' and debug:
        if extra_str != '':
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
def get_cluster_ids(uids, partition):
    clids = {uid : [] for uid in uids}
    for iclust in range(len(partition)):
        for uid in partition[iclust]:
            if iclust not in clids[uid]:  # in case there's duplicates (from seed unique id)
                clids[uid].append(iclust)
    return clids

# ----------------------------------------------------------------------------------------
def new_ccfs_that_need_better_names(partition, true_partition, reco_info, seed_unique_id=None):
    if seed_unique_id is None:
        check_intersection_and_complement(partition, true_partition)
    reco_ids = {uid : reco_info[uid]['reco_id'] for cluster in partition for uid in cluster}  # just a teensy lil' optimization
    uids = set([uid for cluster in partition for uid in cluster])
    clids = get_cluster_ids(uids, partition)  # inferred cluster ids

    def get_clonal_fraction(uid, inferred_cluster):
        """ Return the fraction of seqs in <uid>'s inferred cluster which are really clonal. """
        n_clonal = 0
        for tmpid in inferred_cluster:  # NOTE this includes the case where tmpid equal to uid
            if reco_ids[tmpid] == reco_ids[uid]:
                n_clonal += 1
        return float(n_clonal) / len(inferred_cluster)

    def get_fraction_present(inferred_cluster, true_cluster):
        """ Return the fraction of the true clonemates in <true_cluster> which appear in <inferred_cluster>. """
        n_present = 0
        for tmpid in true_cluster:  # NOTE this includes the case where tmpid equal to uid
            if tmpid in inferred_cluster:
                n_present += 1
        return float(n_present) / len(true_cluster)

    mean_clonal_fraction, mean_fraction_present = 0., 0.
    n_uids = 0
    for true_cluster in true_partition:
        if seed_unique_id is not None and seed_unique_id not in true_cluster:
            continue
        for uid in true_cluster:
            if seed_unique_id is not None and uid != seed_unique_id:
                continue
            if len(clids[uid]) != 1:
                print 'WARNING %s in multiple clusters' % uid
            inferred_cluster = partition[clids[uid][0]]
            mean_clonal_fraction += get_clonal_fraction(uid, inferred_cluster)
            mean_fraction_present += get_fraction_present(inferred_cluster, true_cluster)
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
    iclust, found = None, False
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

    uids_a = set([uid for cluster in part_a for uid in cluster])
    uids_b = set([uid for cluster in part_b for uid in cluster])
    a_not_b = uids_a - uids_b
    b_not_a = uids_b - uids_a
    if len(a_not_b) > 0 or len(b_not_a) > 0:
        raise Exception('partitions don\'t have the same uids')

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
def add_indels_to_germline_strings(line, indelfo):
    """ Add stars to the germline sequences (for ascii printing) if there were SHM insertions. """

    if len(indelfo['indels']) > 1:
        print '    WARNING found %d indels, but we can only handle 1' % len(indelfo['indels'])

    lastfo = indelfo['indels'][-1]

    if lastfo['type'] == 'deletion':
        return

    # divide it up into its constituent chunks
    chunks = [line['fv_insertion'], line['v_gl_seq'], line['vd_insertion'], line['d_gl_seq'], line['dj_insertion'], line['j_gl_seq'], line['jf_insertion']]
    chunknames = ['fv_insertion', 'v', 'vd_insertion', 'd', 'dj_insertion', 'j', 'jf_insertion']
    weirdolist = []  # list of lists, where each entry is the name of the chunk in which we find ourselves
    for ichunk in range(len(chunks)):
        for inuke in range(len(chunks[ichunk])):
            weirdolist.append(chunknames[ichunk])
    thischunk = weirdolist[lastfo['pos']]
    offset = weirdolist.index(thischunk)  # index of first occurence
    if thischunk in regions:
        line['lengths'][thischunk] += lastfo['len']
        line[thischunk + '_gl_seq'] = line[thischunk + '_gl_seq'][ : lastfo['pos'] - offset] + '*' * lastfo['len'] + line[thischunk + '_gl_seq'][lastfo['pos'] - offset : ]
    else:
        print '     unhandled indel within a non-templated insertion'

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

    def slurm_exists():
        try:
            fnull = open(os.devnull, 'w')
            check_output(['which', 'srun'], stderr=fnull, close_fds=True)
            return True
        except CalledProcessError:
            return False

    ncpu = multiprocessing.cpu_count()
    if n_procs > ncpu and slurm_exists():
        return True
    return False

# ----------------------------------------------------------------------------------------
def synthesize_single_seq_line(line, iseq):
    """ without modifying <line>, make a copy of it corresponding to a single-sequence event with the <iseq>th sequence """
    hmminfo = copy.deepcopy(line)  # make a copy of the info, into which we'll insert the sequence-specific stuff
    for col in linekeys['per_seq']:
        hmminfo[col] = [hmminfo[col][iseq], ]
    return hmminfo

# ----------------------------------------------------------------------------------------
def synthesize_multi_seq_line(uids, reco_info):
    """ assumes you already added all the implicit info """
    assert len(uids) > 0
    outline = copy.deepcopy(reco_info[uids[0]])
    for col in linekeys['per_seq']:
        assert [len(reco_info[uid][col]) for uid in uids].count(1) == len(uids)  # make sure they're all length one
        outline[col] = [copy.deepcopy(reco_info[uid][col][0]) for uid in uids]
    return outline

# ----------------------------------------------------------------------------------------
def count_gaps(seq, istop=None):
    """ return number of gap characters up to, but not including <istop> """
    if istop is not None:
        seq = seq[ : istop]
    return sum([seq.count(gc) for gc in gap_chars])

# ----------------------------------------------------------------------------------------
def add_dummy_alignments(line):
    for region in regions:
        line['aligned_' + region + '_seqs'] = ['' for _ in range(len(line['seqs']))]

# ----------------------------------------------------------------------------------------
def add_regional_alignments(glfo, aligned_gl_seqs, line, region, debug=False):
    aligned_seqs = []
    for iseq in range(len(line['seqs'])):
        qr_seq = line[region + '_qr_seqs'][iseq]
        gl_seq = line[region + '_gl_seq']
        aligned_gl_seq = aligned_gl_seqs[region][line[region + '_gene']]
        if len(qr_seq) != len(gl_seq):
            line['invalid'] = True
            continue

        if debug:
            print 'before alignment'
            print '   qr   ', qr_seq
            print '   gl   ', gl_seq
            print '   al gl', aligned_gl_seq

        n_gaps = sum([aligned_gl_seq.count(gc) for gc in gap_chars])
        if len(aligned_gl_seq) != line[region + '_5p_del'] + len(gl_seq) + line[region + '_3p_del'] + n_gaps:
            print len(aligned_gl_seq), line[region + '_5p_del'], len(gl_seq), line[region + '_3p_del'], n_gaps
            line['invalid'] = True
            continue

        qr_seq = 'N' * line[region + '_5p_del'] + qr_seq + 'N' * line[region + '_3p_del']
        gl_seq = 'N' * line[region + '_5p_del'] + gl_seq + 'N' * line[region + '_3p_del']

        for ibase in range(len(aligned_gl_seq)):
            if aligned_gl_seq[ibase] in gap_chars:
                qr_seq = qr_seq[ : ibase] + gap_chars[0] + qr_seq[ibase : ]
                gl_seq = gl_seq[ : ibase] + gap_chars[0] + gl_seq[ibase : ]
            else:
                if gl_seq[ibase] != 'N' and gl_seq[ibase] != aligned_gl_seq[ibase]:
                    line['invalid'] = True
                    break
        if line['invalid']:
            continue

        if debug:
            print 'after alignment'
            print '   qr   ', qr_seq
            print '   gl   ', gl_seq
            print '   al gl', aligned_gl_seq

        if len(qr_seq) != len(gl_seq) or len(qr_seq) != len(aligned_gl_seq):
            line['invalid'] = True
            continue
        aligned_seqs.append(qr_seq)  # TODO is this supposed to be just the v section of the query sequence, or the whole sequence? (if it's the latter, I don't know what to do about alignments)

    if line['invalid']:  # this seems to only happening when we simulate with debug 1, and in that case we don't actually do anything with the aligned seqs, so screw it
        # aligned_seqs = [None for _ in range(len(line['seqs']))]
        raise Exception('failed adding alignment info for %s' % ' '.join(line['unique_ids']))
    line['aligned_' + region + '_seqs'] = aligned_seqs

# ----------------------------------------------------------------------------------------
def add_alignments(glfo, aligned_gl_seqs, line, debug=False):
    for region in regions:
        add_regional_alignments(glfo, aligned_gl_seqs, line, region, debug)

# ----------------------------------------------------------------------------------------
def intexterpolate(x1, y1, x2, y2, x):
    """ interpolate/extrapolate linearly based on two points in 2-space, returning y-value corresponding to <x> """
    m = (y2 - y1) / (x2 - x1)
    b = 0.5 * (y1 + y2 - m*(x1 + x2))
    # if debug:
    #     for x in [x1, x2]:
    #         print '%f x + %f = %f' % (m, b, m*x + b)
    return m * x + b

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

# ----------------------------------------------------------------------------------------
def get_empty_indel():
    return {'reversed_seq' : '', 'indels' : []}

# ----------------------------------------------------------------------------------------
def choose_seed_unique_id(gldir, chain, simfname, seed_cluster_size_low, seed_cluster_size_high, iseed=None, n_max_queries=-1, debug=True):
    glfo = glutils.read_glfo(gldir, chain)
    _, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False, glfo=glfo, n_max_queries=n_max_queries)
    true_partition = get_true_partition(reco_info)

    nth_seed = 0  # don't always take the first one we find
    for cluster in true_partition:
        if len(cluster) < seed_cluster_size_low or len(cluster) > seed_cluster_size_high:
            continue
        if iseed is not None and int(iseed) > nth_seed:
            nth_seed += 1
            continue
        if debug:
            print '    chose seed %s in cluster %s with size %d' % (cluster[0], reco_info[cluster[0]]['reco_id'], len(cluster))
        return cluster[0], len(cluster)  # arbitrarily use the first member of the cluster as the seed

    raise Exception('couldn\'t find seed in cluster between size %d and %d' % (seed_cluster_size_low, seed_cluster_size_high))

# ----------------------------------------------------------------------------------------
# Takes two logd values and adds them together, i.e. takes (log a, log b) --> log a+b
# i.e. a *or* b
def add_in_log_space(first, second):
    if first == -float('inf'):
        return second
    elif second == -float('inf'):
        return first
    elif first > second:
        return first + math.log(1 + math.exp(second - first))
    else:
        return second + math.log(1 + math.exp(first - second))
