import platform
import resource
import psutil
import numpy
import tempfile
import string
import time
import sys
import os
import random
import itertools
import ast
import math
import glob
from collections import OrderedDict
import csv
import subprocess
import multiprocessing
import copy

import indelutils

# ----------------------------------------------------------------------------------------
def fsdir():
    fsdir = '/fh/fast/matsen_e'
    if not os.path.exists(fsdir):
        fsdir = '/tmp'
    if os.getenv('USER') is not None:
        fsdir += '/' + os.getenv('USER')
    return fsdir

# ----------------------------------------------------------------------------------------
# putting these up here so glutils import doesn't fail... I think I should be able to do it another way, though
regions = ['v', 'd', 'j']
loci = {'igh' : 'vdj',
        'igk' : 'vj',
        'igl' : 'vj',
        'tra' : 'vj',
        'trb' : 'vdj',
        'trg' : 'vj',
        'trd' : 'vdj',
}

def getregions(locus):  # for clarity, don't use the <loci> dictionary directly to access its .values()
    return list(loci[locus])  # doesn't really need to be a list, but it's more clearly analagous to regions & co that way

def has_d_gene(locus):  # for clarity, don't use the <loci> dictionary directly to access its .values()
    return 'd' in loci[locus]

def region_pairs():
    return [{'left' : bound[0], 'right' : bound[1]} for bound in boundaries]

import seqfileopener
import glutils
import prutils

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
def get_arg_list(arg, intify=False, floatify=False, translation=None, list_of_pairs=False, choices=None):  # make lists from args that are passed as strings of colon-separated values
    if arg is None:
        return None

    convert_fcn = pass_fcn
    if intify:
        convert_fcn = int
    elif floatify:
        convert_fcn = float

    arglist = arg.strip().split(':')  # to allow ids with minus signs, you can add a space (if you don't use --name=val), which you then have to strip() off
    if list_of_pairs:
        arglist = [pairstr.split(',') for pairstr in arglist]
        arglist = [[convert_fcn(p) for p in pair] for pair in arglist]
    else:
        arglist = [convert_fcn(x) for x in arglist]

    if translation is not None:
        for ia in range(len(arglist)):
            if arglist[ia] in translation:
                arglist[ia] = translation[arglist[ia]]

    if choices is not None:
        for arg in arglist:
            if arg not in choices:
                raise Exception('unexpected argument \'%s\' (choices: %s)' % (str(arg), [str(c) for c in choices]))

    return arglist

# ----------------------------------------------------------------------------------------
# values used when simulating from scratch
# scratch_mean_mute_freqs = {'v' : 0.03, 'd' : 0.8, 'j' : 0.06}
# scratch_mean_mute_freqs['all'] = numpy.mean([v for v in scratch_mean_mute_freqs.values()])
scratch_mean_erosion_lengths = {'v_3p' : 2, 'd_5p' : 3, 'd_3p' : 3, 'j_5p' : 4}
scratch_mean_insertion_lengths = {l : {'vd' : 4 if has_d_gene(l) else 0,  # e.g. light chain gets no vd insertion
                                       'dj' : 4 if has_d_gene(l) else 8}  # ...but a longer dj insertion
                                  for l in loci}

real_erosions = ['v_3p', 'd_5p', 'd_3p', 'j_5p']
# NOTE since we now handle v_5p and j_3p deletions by padding with Ns, the hmm does *not* allow actual v_5p and j_3p deletions.
# This means that while we write parameters for v_5p and j_3p deletions to the parameter dir, these are *not* used in making the
# hmm yamels -- which is what we want, because we want to be able to read in short data reads but make full-length simulation.
effective_erosions = ['v_5p', 'j_3p']
all_erosions = real_erosions + effective_erosions
boundaries = ['vd', 'dj']  # NOTE needs to be updated to be consistent with erosion names
effective_boundaries = ['fv', 'jf']
all_boundaries = boundaries + effective_boundaries
nukes = ['A', 'C', 'G', 'T']
ambiguous_bases = ['N', ]
alphabet = set(nukes + ambiguous_bases)  # NOTE not the greatest naming distinction, but note difference to <expected_characters>
gap_chars = ['.', '-']
expected_characters = set(nukes + ambiguous_bases + gap_chars)  # NOTE not the greatest naming distinction, but note difference to <alphabet>
conserved_codons = {l : {'v' : 'cyst',
                          'j' : 'tryp' if l == 'igh' else 'phen'}  # e.g. heavy chain has tryp, light chain has phen
                     for l in loci}

def cdn_positions(glfo, region):
    return glfo[conserved_codons[glfo['locus']][region] + '-positions']
def cdn_pos(glfo, region, gene):
    return cdn_positions(glfo, region)[gene]
def gap_len(seq):
    return len(filter(gap_chars.__contains__, seq))
def non_gap_len(seq):
    return len(seq) - gap_len(seq)

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
    if codon not in [c for locus in loci for c in conserved_codons[locus].values()]:
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

# tuples with the column and its dependencies mashed together
# (first entry is the column of interest, and it depends upon the following entries)
column_dependency_tuples = []
for column, deps in column_dependencies.iteritems():
    tmp_list = [column]
    tmp_list.extend(deps)
    column_dependency_tuples.append(tuple(tmp_list))

adaptive_headers = {
    'seqs' : 'nucleotide',
    'v_gene' : 'vMaxResolved',
    'd_gene' : 'dMaxResolved',
    'j_gene' : 'jMaxResolved',
    'v_3p_del' : 'vDeletion',
    'd_5p_del' : 'd5Deletion',
    'd_3p_del' : 'd3Deletion',
    'j_5p_del' : 'jDeletion'
}

# ----------------------------------------------------------------------------------------
forbidden_characters = set([':', ';', ','])  # strings that are not allowed in sequence ids
forbidden_character_translations = string.maketrans(':;,', 'csm')

functional_columns = ['mutated_invariants', 'in_frames', 'stops']

column_configs = {
    'ints' : ['n_mutations', 'cdr3_length', 'padlefts', 'padrights'] + [e + '_del' for e in all_erosions],
    'floats' : ['logprob', 'mut_freqs'],
    'bools' : functional_columns,
    'literals' : ['indelfo', 'indelfos', 'k_v', 'k_d', 'all_matches'],  # simulation has indelfo[s] singular, annotation output has it plural... and I think it actually makes sense to have it that way
    'lists' : ['unique_ids', 'seqs', 'input_seqs', 'indel_reversed_seqs', 'n_mutations', 'mut_freqs', 'padlefts', 'padrights'] + ['aligned_' + r + '_seqs' for r in regions] + functional_columns,  # indelfos is a list, but we can't just split it by colons since it has colons within the dict string
    'lists-of-lists' : ['duplicates'] + [r + '_per_gene_support' for r in regions]
}

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
def get_str_float_pair_dict(strlist):
    def getpair(pairstr):
        pairlist = pairstr.split(':')
        assert len(pairlist) == 2
        return (pairlist[0], float(pairlist[1]))
    return OrderedDict(getpair(pairstr) for pairstr in strlist)

# ----------------------------------------------------------------------------------------
def get_list_of_str_list(strlist):
    if strlist == '':
        return []
    return [substr.split(':') for substr in strlist]

# ----------------------------------------------------------------------------------------
def reconstruct_full_indelfo(indel_list, reversed_seq):
    if 'reversed_seq' in indel_list:  # handle old files
        return indel_list
    indelfo = indelutils.get_empty_indel()
    indelfo['indels'] = indel_list
    if len(indel_list) > 0:
        indelfo['reversed_seq'] = reversed_seq
    return indelfo

conversion_fcns = {}
for key in column_configs['ints']:
    conversion_fcns[key] = int
for key in column_configs['floats']:
    conversion_fcns[key] = float
for key in column_configs['bools']:
    conversion_fcns[key] = useful_bool
for key in column_configs['literals']:
    conversion_fcns[key] = ast.literal_eval
for region in regions:
    conversion_fcns[region + '_per_gene_support'] = get_str_float_pair_dict
conversion_fcns['duplicates'] = get_list_of_str_list

# keep track of all the *@*@$!ing different keys that happen in the <line>/<hmminfo>/whatever dictionaries
linekeys = {}
linekeys['per_family'] = ['naive_seq', 'cdr3_length', 'codon_positions', 'lengths', 'regional_bounds'] + \
                         [r + '_gene' for r in regions] + \
                         [e + '_del' for e in all_erosions] + \
                         [b + '_insertion' for b in all_boundaries] + \
                         [r + '_gl_seq' for r in regions] + \
                         [r + '_per_gene_support' for r in regions]
# used by the synthesize_[] fcns below
linekeys['per_seq'] = ['seqs', 'unique_ids', 'indelfos', 'mut_freqs', 'n_mutations', 'input_seqs', 'indel_reversed_seqs'] + \
                      [r + '_qr_seqs' for r in regions] + \
                      ['aligned_' + r + '_seqs' for r in regions] + \
                      functional_columns
linekeys['hmm'] = ['logprob', 'errors']
linekeys['sw'] = ['k_v', 'k_d', 'all_matches', 'padlefts', 'padrights', 'duplicates', 'flexbounds', 'relpos']  # TODO move 'duplicates' to 'per_seq' (see note in synthesize_multi_seq_line())
linekeys['extra'] = ['invalid', ]
linekeys['simu'] = ['reco_id', ]
all_linekeys = set([k for cols in linekeys.values() for k in cols])

# keys that are added by add_implicit_info()
implicit_linekeys = set(['naive_seq', 'cdr3_length', 'codon_positions', 'lengths', 'regional_bounds', 'invalid', 'indel_reversed_seqs'] + \
                        [r + '_gl_seq' for r in regions] + \
                        ['mut_freqs', 'n_mutations'] + functional_columns + [r + '_qr_seqs' for r in regions] + ['aligned_' + r + '_seqs' for r in regions])

# ----------------------------------------------------------------------------------------
annotation_headers = ['unique_ids', 'v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'mut_freqs', 'n_mutations', 'input_seqs', 'indel_reversed_seqs', 'naive_seq', 'indelfos', 'duplicates'] \
                     + [r + '_per_gene_support' for r in regions] \
                     + [e + '_del' for e in all_erosions] + [b + '_insertion' for b in all_boundaries] \
                     + functional_columns
sw_cache_headers = ['k_v', 'k_d', 'padlefts', 'padrights', 'all_matches', 'mut_freqs']
partition_cachefile_headers = ('unique_ids', 'logprob', 'naive_seq', 'naive_hfrac', 'errors')  # these have to match whatever bcrham is expecting
bcrham_dbgstrs = {  # corresponds to stdout from glomerator.cc
    'read-cache' : ['logprobs', 'naive-seqs'],
    'calcd' : ['vtb', 'fwd', 'hfrac'],
    'merged' : ['hfrac', 'lratio'],
    'time' : ['bcrham', ]
}
bcrham_dbgstr_types = {
    'sum' : ['calcd', 'merged'],  # for these ones, sum over all procs
    'same' : ['read-cache', ],  # check that these are the same for all procs
    'min-max' : ['time', ]
}

# ----------------------------------------------------------------------------------------
def synthesize_single_seq_line(line, iseq):
    """ without modifying <line>, make a copy of it corresponding to a single-sequence event with the <iseq>th sequence """
    singlefo = {}
    for key in line:
        if key in linekeys['per_seq']:
            singlefo[key] = [copy.deepcopy(line[key][iseq]), ]
        else:
            singlefo[key] = copy.deepcopy(line[key])
    return singlefo

# ----------------------------------------------------------------------------------------
def synthesize_multi_seq_line(uids, multifo):  # assumes you already added all the implicit info
    reco_info = {multifo['unique_ids'][iseq] : synthesize_single_seq_line(multifo, iseq) for iseq in range(len(multifo['unique_ids']))}
    return synthesize_multi_seq_line_from_reco_info(uids, reco_info)

# ----------------------------------------------------------------------------------------
def synthesize_multi_seq_line_from_reco_info(uids, reco_info):  # assumes you already added all the implicit info
    assert len(uids) > 0
    multifo = copy.deepcopy(reco_info[uids[0]])
    for col in [c for c in linekeys['per_seq'] if c in multifo]:
        assert [len(reco_info[uid][col]) for uid in uids].count(1) == len(uids)  # make sure every uid's info for this column is of length 1
        multifo[col] = [copy.deepcopy(reco_info[uid][col][0]) for uid in uids]
    return multifo

# ----------------------------------------------------------------------------------------
def generate_dummy_v(d_gene):
    pv, sv, al = split_gene(d_gene)
    return get_locus(d_gene).upper() + 'VxDx' + pv + '-' + sv + '*' + al

# ----------------------------------------------------------------------------------------
def convert_from_adaptive_headers(glfo, line, uid=None, only_dj_rearrangements=False):
    newline = {}
    print_it = False

    for head, ahead in adaptive_headers.items():
        newline[head] = line[ahead]
        if head in column_configs['lists']:
            newline[head] = [newline[head], ]

    if uid is not None:
        newline['unique_ids'] = [uid, ]

    for erosion in real_erosions:
        newline[erosion + '_del'] = int(newline[erosion + '_del'])
    newline['v_5p_del'] = 0
    newline['j_3p_del'] = 0

    for region in regions:
        if newline[region + '_gene'] == 'unresolved':
            newline[region + '_gene'] = None
            continue
        if region == 'j' and 'P' in newline[region + '_gene']:
            newline[region + '_gene'] = newline[region + '_gene'].replace('P', '')

        if '*' not in newline[region + '_gene']:
            # print uid
            # tmpheads = ['dMaxResolved', 'dFamilyName', 'dGeneName', 'dGeneAllele', 'dFamilyTies', 'dGeneNameTies', 'dGeneAlleleTies']
            # for h in tmpheads:
            #     print '   %s: %s' % (h, line[h]),
            # print ''
            if line['dGeneAlleleTies'] == '':
                newline['failed'] = True
                return newline
            d_alleles = line['dGeneAlleleTies'].split(',')
            newline[region + '_gene'] += '*' + d_alleles[0]
        primary_version, sub_version, allele = split_gene(newline[region + '_gene'])
        primary_version, sub_version = primary_version.lstrip('0'), sub_version.lstrip('0')  # alleles get to keep their leading zero (thank you imgt for being consistent)
        if region == 'j':  # adaptive calls every j sub_version 1
            sub_version = None
        gene = rejoin_gene(glfo['locus'], region, primary_version, sub_version, allele)
        if gene not in glfo['seqs'][region]:
            gene = glutils.convert_to_duplicate_name(glfo, gene)
        if gene not in glfo['seqs'][region]:
            raise Exception('couldn\'t rebuild gene name from adaptive data: %s' % gene)
        newline[region + '_gene'] = gene

    seq = newline['seqs'][0]
    boundlist = ['vIndex', 'n1Index', 'dIndex', 'n2Index', 'jIndex']
    qrbounds, glbounds = {}, {}
    for region in regions:
        if only_dj_rearrangements and region == 'v':
            if newline['d_gene'] is None:
                newline['failed'] = True
                return newline
            newline['v_gene'] = generate_dummy_v(newline['d_gene'])
            line['vIndex'] = 0
            line['n1Index'] = int(line['dIndex']) - 1  # or is it without the -1?
            glfo['seqs']['v'][newline['v_gene']] = seq[line['vIndex'] : line['n1Index']]
            glfo['cyst-positions'][newline['v_gene']] = len(glfo['seqs']['v'][newline['v_gene']]) - 3
        if newline[region + '_gene'] is None:
            newline['failed'] = True
            return newline
        qrb = [int(line[region + 'Index']),
                    int(line[boundlist[boundlist.index(region + 'Index') + 1]]) if region != 'j' else len(seq)]
        glseq = glfo['seqs'][region][newline[region + '_gene']]
        glb = [newline[region + '_5p_del'],
                    len(glseq) - newline[region + '_3p_del']]
        if region == 'j' and glb[1] - glb[0] > qrb[1] - qrb[0]:  # extra adaptive stuff on right side of j
            old = glb[1]
            glb[1] = glb[0] + qrb[1] - qrb[0]
            newline['j_3p_del'] = old - glb[1]

        if qrb[0] == -1 or qrb[1] == -1 or qrb[1] < qrb[0]:  # should this also be equals?
            newline['failed'] = True
            return newline
        if qrb[1] - qrb[0] != glb[1] - glb[0]:
            newline['failed'] = True
            return newline
        qrbounds[region] = qrb
        glbounds[region] = glb

    for bound in boundaries:
        newline[bound + '_insertion'] = seq[qrbounds[bound[0]][1] : qrbounds[bound[1]][0]]  # end of lefthand region to start of righthand region

    newline['fv_insertion'] = ''
    newline['jf_insertion'] = seq[qrbounds['j'][1]:]

    # print seq
    # print seq[:qrbounds['d'][0]],
    # print seq[qrbounds['d'][0] : qrbounds['d'][1]],
    # print seq[qrbounds['d'][1] : qrbounds['j'][0]],
    # print seq[qrbounds['j'][0] : qrbounds['j'][1]],
    # print seq[qrbounds['j'][1] :]

    newline['indelfos'] = [indelutils.get_empty_indel(), ]

    if print_it:
        add_implicit_info(glfo, newline)
        print_reco_event(newline, label=uid)

    # still need to convert to integers/lists/whatnot (?)

    newline['failed'] = False

    return newline

# ----------------------------------------------------------------------------------------
# definitions here: http://clip.med.yale.edu/changeo/manuals/Change-O_Data_Format.pdf
presto_headers = OrderedDict([  # enforce this ordering so the output files are easier to read
    ('SEQUENCE_ID', 'unique_ids'),
    ('V_CALL', 'v_gene'),
    ('D_CALL', 'd_gene'),
    ('J_CALL', 'j_gene'),
    ('JUNCTION_LENGTH', None),
    ('SEQUENCE_INPUT', 'input_seqs'),
    ('SEQUENCE_IMGT', 'aligned_v_plus_unaligned_dj'),
])

# ----------------------------------------------------------------------------------------
def get_line_with_presto_headers(line):  # NOTE doesn't deep copy
    """ convert <line> to presto csv format """
    if len(line['unique_ids']) > 1:  # has to happen *before* utils.get_line_for_output()
        raise Exception('multiple seqs not handled for presto output')

    presto_line = {}
    for phead, head in presto_headers.items():
        if head == 'aligned_v_plus_unaligned_dj':
            presto_line[phead] = line['aligned_v_seqs'][0] + line['vd_insertion'] + line['d_qr_seqs'][0] + line['dj_insertion'] + line['j_qr_seqs'][0]
        elif phead == 'JUNCTION_LENGTH':
            presto_line[phead] = line['cdr3_length'] + 6
        elif head == 'unique_ids' or head == 'input_seqs':
            presto_line[phead] = line[head][0]
        else:
            presto_line[phead] = line[head]

    return presto_line

# ----------------------------------------------------------------------------------------
def write_presto_annotations(outfname, glfo, annotations, failed_queries):
    # outstr = subprocess.check_output(['mv', '-v', outfname, outfname + '.partis'])
    # print '    backing up partis output before converting to presto: %s' % outstr.strip()

    with open(outfname, 'w') as outfile:
        writer = csv.DictWriter(outfile, presto_headers.keys(), delimiter='\t')
        writer.writeheader()

        for full_line in annotations.values():
            writer.writerow(get_line_with_presto_headers(full_line))

        # and write empty lines for seqs that failed either in sw or the hmm
        for uid, seq in failed_queries.items():
            writer.writerow({'SEQUENCE_ID' : uid, 'SEQUENCE_INPUT' : seq})

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

def color(col, seq, width=None, padside='left'):
    return_str = [Colors[col], seq, Colors['end']]
    if width is not None:  # make sure final string prints to correct width
        n_spaces = max(0, width - len(seq))  # if specified <width> is greater than uncolored length of <seq>, pad with spaces so that when the colors show up properly the colored sequences prints with width <width>
        if padside == 'left':
            return_str.insert(0, n_spaces * ' ')
        elif padside == 'right':
            return_str.insert(len(return_str), n_spaces * ' ')
        else:
            assert False
    return ''.join(return_str)

def len_excluding_colors(seq):  # NOTE this won't work if you inserted a color code into the middle of another color code
    for color_code in Colors.values():
        seq = seq.replace(color_code, '')
    return len(seq)

def len_only_letters(seq):  # usually the same as len_excluding_colors(), except it doesn't count gap chars or spaces
    return len(filter((alphabet).__contains__, seq))

# ----------------------------------------------------------------------------------------
def color_chars(chars, col, seq):
    if sum([seq.count(c) for c in chars]) == 0:  # if <chars> aren't present, immediately return
        return seq
    return_str = [color(col, c) if c in chars else c for c in seq]
    return ''.join(return_str)

# ----------------------------------------------------------------------------------------
def align_many_seqs(seqfos, outfname=None):  # if <outfname> is specified, we just tell mafft to write to <outfname> and then return None
    def outfile_fcn():
        if outfname is None:
            return tempfile.NamedTemporaryFile()
        else:
            return open(outfname, 'w')

    with tempfile.NamedTemporaryFile() as fin, outfile_fcn() as fout:
        for seqfo in seqfos:
            fin.write('>%s\n%s\n' % (seqfo['name'], seqfo['seq']))
        fin.flush()
        subprocess.check_call('mafft --quiet %s >%s' % (fin.name, fout.name), shell=True)
        if outfname is not None:
            return None
        msa_info = read_fastx(fout.name, ftype='fa')
        x = set([sfo['name'] for sfo in msa_info])
        y = set([sfo['name'] for sfo in seqfos])
        if set([sfo['name'] for sfo in msa_info]) != set([sfo['name'] for sfo in seqfos]):
            subprocess.check_call(['cat', fin.name])
            raise Exception('incoherent mafft output from %s (cat\'d on previous line)' % fin.name)
    return msa_info

# ----------------------------------------------------------------------------------------
def align_seqs(ref_seq, seq):  # should eventually change name to align_two_seqs() or something
    with tempfile.NamedTemporaryFile() as fin, tempfile.NamedTemporaryFile() as fout:
        fin.write('>%s\n%s\n' % ('ref', ref_seq))
        fin.write('>%s\n%s\n' % ('new', seq))
        fin.flush()
        subprocess.check_call('mafft --quiet %s >%s' % (fin.name, fout.name), shell=True)
        msa_info = {sfo['name'] : sfo['seq'] for sfo in read_fastx(fout.name, ftype='fa')}
        if 'ref' not in msa_info or 'new' not in msa_info:
            subprocess.check_call(['cat', fin.name])
            raise Exception('incoherent mafft output from %s (cat\'d on previous line)' % fin.name)
    return msa_info['ref'], msa_info['new']

# ----------------------------------------------------------------------------------------
def cons_seq(threshold, aligned_seqfos=None, unaligned_seqfos=None, debug=False):
    print '%s decide how you want to do these imports' % color('red', 'hey')
    from cStringIO import StringIO
    from Bio.Align import AlignInfo
    # if 'Bio.Align' not in sys.modules:
    #     from Bio.Align import AlignInfo
    # Bio.Align.AlignInfo = sys.modules['Bio.Align']
    import Bio.AlignIO
    # if 'Bio.AlignIO' not in sys.modules:
    #     import Bio.AlignIO
    # Bio.AlignIO = sys.modules['Bio.AlignIO']

    if aligned_seqfos is not None:
        assert unaligned_seqfos is None
        seqfos = aligned_seqfos
    elif unaligned_seqfos is not None:
        assert aligned_seqfos is None
        seqfos = align_many_seqs(unaligned_seqfos)
    else:
        assert False

    fastalist = ['>%s\n%s' % (sfo['name'], sfo['seq']) for sfo in seqfos]
    alignment = Bio.AlignIO.read(StringIO('\n'.join(fastalist) + '\n'), 'fasta')
    cons_seq = str(AlignInfo.SummaryInfo(alignment).gap_consensus(threshold, ambiguous='N'))

    # debug = True
    # if debug:
    #     print 'consensus %s' % cons_seq
    #     for sfo in seqfos:
    #         print '' % (color_mutants(cons_seq, sfo['seq'], align=True))
    #     sys.exit()

    return cons_seq

# ----------------------------------------------------------------------------------------
def color_mutants(ref_seq, seq, print_result=False, extra_str='', ref_label='', seq_label='', post_str='', print_hfrac=False, print_isnps=False, return_isnps=False, emphasis_positions=None, use_min_len=False, only_print_seq=False, align=False, return_ref=False):
    """ default: return <seq> string with colored mutations with respect to <ref_seq> """

    # NOTE now I've got <return_ref>, I can probably remove a bunch of the label/whatever arguments and do all the damn formatting in the caller

    if use_min_len:
        min_len = min(len(ref_seq), len(seq))
        ref_seq = ref_seq[:min_len]
        seq = seq[:min_len]

    if align:  #  and len(ref_seq) != len(seq):  # it would be nice to avoid aligning when we don't need to... but i'm not sure how to identify cases where multiple indels result in the same length
        ref_seq, seq = align_seqs(ref_seq, seq)

    if len(ref_seq) != len(seq):
        raise Exception('unequal lengths in color_mutants()\n    %s\n    %s' % (ref_seq, seq))

    return_str, isnps = [], []
    for inuke in range(len(seq)):  # would be nice to integrate this with hamming_distance()
        rchar = ref_seq[inuke]
        char = seq[inuke]
        if char in ambiguous_bases or rchar in ambiguous_bases:
            char = color('blue', char)
        elif char in gap_chars or rchar in gap_chars:
            char = color('blue', char)
        elif char != rchar:
            char = color('red', char)
            isnps.append(inuke)
        if emphasis_positions is not None and inuke in emphasis_positions:
            char = color('reverse_video', char)
        return_str.append(char)

    isnp_str = ''
    if print_isnps and len(isnps) > 0:
        isnp_str = '   %d snp%s' % (len(isnps), plural(len(isnps)))
        if len(isnps) < 10:
            isnp_str +=  ' at: %s' % ' '.join([str(i) for i in isnps])
    hfrac_str = ''
    if print_hfrac:
        hfrac_str = '   hfrac %.3f' % hamming_fraction(ref_seq, seq)
    if print_result:
        lwidth = max(len_excluding_colors(ref_label), len_excluding_colors(seq_label))
        if not only_print_seq:
            print '%s%s%s' % (extra_str, ('%' + str(lwidth) + 's') % ref_label, ref_seq)
        print '%s%s%s' % (extra_str, ('%' + str(lwidth) + 's') % seq_label, ''.join(return_str) + post_str + isnp_str + hfrac_str)

    return_list = [extra_str + seq_label + ''.join(return_str) + post_str + isnp_str + hfrac_str]
    if return_ref:
        return_list.append(''.join([ch if ch not in gap_chars else color('blue', ch) for ch in ref_seq]))
    if return_isnps:
        return_list.append(isnps)
    return return_list[0] if len(return_list) == 1 else return_list

# ----------------------------------------------------------------------------------------
def plural_str(pstr, count):
    if count == 1:
        return pstr
    else:
        return pstr + 's'

# ----------------------------------------------------------------------------------------
def plural(count):  # TODO should combine these
    if count == 1:
        return ''
    else:
        return 's'

# ----------------------------------------------------------------------------------------
def summarize_gene_name(gene):
    region = get_region(gene)
    primary_version, sub_version, allele = split_gene(gene)
    return ' '.join([region, primary_version, sub_version, allele])

# ----------------------------------------------------------------------------------------
def color_gene(gene, width=None, leftpad=False):
    """ color gene name (and remove extra characters), eg IGHV3-h*01 --> hv3-h1 """
    locus = get_locus(gene)
    locus = locus[2]  # hmm... maybe?
    region = get_region(gene)
    primary_version, sub_version, allele = split_gene(gene)

    n_chars = len(locus + region + primary_version)  # number of non-special characters
    return_str = color('purple', locus) + color('red', region) + color('purple', primary_version)
    if sub_version is not None:
        n_chars += 1 + len(sub_version)
        return_str += color('purple', '-' + sub_version)
    n_chars += len(allele)
    return_str += color('yellow', allele)
    if width is not None:
        if leftpad:
            return_str = (width - n_chars) * ' ' + return_str
        else:
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
def both_codons_unmutated(locus, seq, positions, debug=False, extra_str=''):
    both_ok = True
    for region, codon in conserved_codons[locus].items():
        both_ok &= codon_unmutated(codon, seq, positions[region], debug=debug, extra_str=extra_str)
    return both_ok

# ----------------------------------------------------------------------------------------
def codon_unmutated(codon, seq, position, debug=False, extra_str=''):
    if len(seq) < position + 3:
        if debug:
            print '%ssequence length %d less than %s position %d + 3' % (extra_str, len(seq), codon, position)
        return False
    if seq[position : position + 3] not in codon_table[codon]:  # NOTE this allows it to be mutated to one of the other codons that codes for the same amino acid
        if debug:
            print '%s%s codon %s not among expected codons (%s)' % (extra_str, codon, seq[position : position + 3], ' '.join(codon_table[codon]))
        return False
    return True

#----------------------------------------------------------------------------------------
def in_frame_germline_v(seq, cyst_position):  # NOTE duplication with in_frame() (this is for when all we have is the germline v gene, whereas in_frame() is for when we have the whole rearrangement line)
    return cyst_position <= len(seq) - 3 and cyst_position % 3 == 0

#----------------------------------------------------------------------------------------
def in_frame(seq, codon_positions, fv_insertion, v_5p_del, debug=False):  # NOTE I'm not passing the whole <line> in order to make it more explicit that <seq> and <codon_positions> need to correspond to each other, i.e. are either both for input seqs, or both for indel-reversed seqs
    # NOTE duplication with in_frame_germline_v()
    """ return true if the start and end of the cdr3 are both in frame with respect to the start of the V """
    germline_v_start = len(fv_insertion) - v_5p_del  # position in <seq> (the query sequence) to which the first base of the germline sequence aligns
    v_cpos = codon_positions['v'] - germline_v_start
    j_cpos = codon_positions['j'] - germline_v_start  # NOTE I'm actually not sure how necessary it is that the right side of the J is in frame. I mean, I think it's pretty framework-ey, but I'm not sure.
    if debug:
        print '    in frame:   v codon %d   j codon %d  -->  %s' % (v_cpos % 3 == 0, j_cpos % 3 == 0, v_cpos % 3 == 0 and j_cpos % 3 == 0)
    return v_cpos % 3 == 0 and j_cpos % 3 == 0

#----------------------------------------------------------------------------------------
def is_there_a_stop_codon(seq, fv_insertion, jf_insertion, v_5p_del, debug=False):  # NOTE getting the indexing correct here is extremely non-trivial
    """ true if there's a stop codon in frame with respect to the start of the V """
    germline_v_start = len(fv_insertion) - v_5p_del  # position in <seq> (the query sequence) to which the first base of the germline sequence aligns
    istart = germline_v_start  # start with the first complete codon after <germline_v_start>
    while istart < len(fv_insertion):  # while staying in frame with the start of the v, skootch up to the first base in the query sequence that's actually aligned to the germline (i.e. up to 0 if no fv_insertion, and further if there is one)
        istart += 3
    germline_j_end = len(seq) - len(jf_insertion)  # position immediately after the end of the germline j (if there's a j_3p_del it's accounted for with len(seq))
    istop = germline_j_end - ((germline_j_end - istart) % 3)
    codons = [seq[i : i + 3] for i in range(istart, istop, 3)]
    if debug:
        print '%25s  %3d %3d  %6s   %s' % (seq, istart, istop, len(set(codons) & set(codon_table['stop'])) > 0, ' '.join(codons))
    return len(set(codons) & set(codon_table['stop'])) > 0  # true if any of the stop codons from <codon_table> are present in <codons>

# ----------------------------------------------------------------------------------------
def disambiguate_effective_insertions(bound, line, iseq, debug=False):
    # These are kinda weird names, but the distinction is important
    # If an insert state with "germline" N emits one of [ACGT], then the hmm will report this as an inserted N. Which is what we want -- we view this as a germline N which "mutated" to [ACGT].
    # This concept of insertion germline state is mostly relevant for simultaneous inference on several sequences, i.e. in ham we don't want to just say the inserted base was the base in the query sequence.
    # But here, we're trimming off the effective insertions and we have to treat the inserted germline N to [ACGT] differently than we would an insertion which was "germline" [ACGT] which emitted an N,
    # and also differently to a real germline [VDJ] state that emitted an N.
    naive_insertion = line[bound + '_insertion']  # reminder: ham gets this from the last character in the insertion state name, e.g. 'insert_left_A' or 'insert_right_N'
    insert_len = len(line[bound + '_insertion'])
    if bound == 'fv':  # ...but to accomodate multiple sequences, the insert states can emit non-germline states, so the mature bases might be different.
        mature_insertion = line['seqs'][iseq][ : insert_len]
    elif bound == 'jf':
        mature_insertion = line['seqs'][iseq][len(line['seqs'][iseq]) - insert_len : ]
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
        print '     %s      final: %s' % (color('red', bound), color('purple', final_insertion))
        print '         to_remove: %s' % color('blue', insertion_to_remove)

    return final_insertion, insertion_to_remove

# ----------------------------------------------------------------------------------------
def reset_effective_erosions_and_effective_insertions(glfo, padded_line, aligned_gl_seqs=None, debug=False):  # , padfo=None
    """
    Ham does not allow (well, no longer allows) v_5p and j_3p deletions -- we instead pad sequences with Ns.
    This means that the info we get from ham always has these effective erosions set to zero, but for downstream
    things we sometimes want to know where the reads stopped (e.g. if we want to mimic them in simulation).
    Note that these effective erosion values will be present in the parameter dir, but are *not* incorporated into
    the hmm yaml files.
    """

    if debug:
        print 'resetting effective erosions/insertions for %s' % ' '.join(padded_line['unique_ids'])

    line = {k : copy.deepcopy(padded_line[k]) for k in padded_line if k not in implicit_linekeys}

    assert line['v_5p_del'] == 0  # just to be safe
    assert line['j_3p_del'] == 0
    nseqs = len(line['unique_ids'])  # convenience

    # first disambiguate/remove effective (fv and jf) insertions
    if debug:
        print '   disambiguating effective insertions'
    trimmed_seqs = [line['seqs'][iseq] for iseq in range(nseqs)]
    trimmed_input_seqs = [line['input_seqs'][iseq] for iseq in range(nseqs)]
    final_insertions = [{} for _ in range(nseqs)]  # the effective insertions that will remain in the final info
    insertions_to_remove = [{} for _ in range(nseqs)]  # the effective insertions that we'll be removing, so which won't be in the final info
    for iseq in range(nseqs):
        fin = final_insertions[iseq]
        rem = insertions_to_remove[iseq]
        for bound in effective_boundaries:
            final_insertion, insertion_to_remove = disambiguate_effective_insertions(bound, line, iseq, debug=debug)
            fin[bound] = final_insertion
            rem[bound] = insertion_to_remove
        trimmed_seqs[iseq] = trimmed_seqs[iseq][len(rem['fv']) : len(trimmed_seqs[iseq]) - len(rem['jf'])]
        trimmed_input_seqs[iseq] = trimmed_input_seqs[iseq][len(rem['fv']) : len(trimmed_input_seqs[iseq]) - len(rem['jf'])]
        if debug:
            print '       %s  %s%s%s%s%s' % (' '.join(line['unique_ids']),
                                             color('blue', rem['fv']), color('purple', fin['fv']),
                                             trimmed_seqs[iseq][len(fin['fv']) : len(trimmed_seqs[iseq]) - len(fin['jf'])],
                                             color('purple', fin['jf']), color('blue', rem['jf']))

    # arbitrarily use the zeroth sequence (in principle v_5p and j_3p should be per-sequence, not per-rearrangement... but that'd be a mess to implement, since the other deletions are per-rearrangement)
    TMPiseq = 0  # NOTE this is pretty hackey: we just use the values from the first sequence. But it's actually not that bad -- we can either have some extra pad Ns showing, or chop off some bases.
    trimmed_seq = trimmed_seqs[TMPiseq]
    final_fv_insertion = final_insertions[TMPiseq]['fv']
    final_jf_insertion = final_insertions[TMPiseq]['jf']
    fv_insertion_to_remove = insertions_to_remove[TMPiseq]['fv']
    jf_insertion_to_remove = insertions_to_remove[TMPiseq]['jf']

    def max_effective_erosion(erosion):  # don't "erode" more than there is left to erode
        region = erosion[0]
        gl_len = len(glfo['seqs'][region][line[region + '_gene']])
        if '5p' in erosion:
            other_del = line[region + '_3p_del']
        elif '3p' in erosion:
            other_del = line[region + '_5p_del']
        return gl_len - other_del - 1

    line['v_5p_del'] = min(max_effective_erosion('v_5p'), find_first_non_ambiguous_base(trimmed_seq))
    line['j_3p_del'] = min(max_effective_erosion('j_3p'), len(trimmed_seq) - find_last_non_ambiguous_base_plus_one(trimmed_seq))

    if debug:
        v_5p = line['v_5p_del']
        j_3p = line['j_3p_del']
        print '     %s:  %d' % (color('red', 'v_5p'), v_5p)
        print '     %s:  %d' % (color('red', 'j_3p'), j_3p)
        for iseq in range(nseqs):
            print '       %s  %s%s%s' % (' '.join(line['unique_ids']), color('red', v_5p * '.'), trimmed_seqs[iseq][v_5p : len(trimmed_seqs[iseq]) - j_3p], color('red', j_3p * '.'))

    for iseq in range(nseqs):
        line['seqs'][iseq] = trimmed_seqs[iseq][line['v_5p_del'] : len(trimmed_seqs[iseq]) - line['j_3p_del']]
        line['input_seqs'][iseq] = trimmed_input_seqs[iseq][line['v_5p_del'] : len(trimmed_input_seqs[iseq]) - line['j_3p_del']]

        # also de-pad the indel info
        if indelutils.has_indels(line['indelfos'][iseq]):
            rseq = line['indelfos'][iseq]['reversed_seq']
            rseq = rseq[len(fv_insertion_to_remove) + line['v_5p_del'] : ]
            if len(jf_insertion_to_remove) + line['j_3p_del'] > 0:
                rseq = rseq[ : -(len(jf_insertion_to_remove) + line['j_3p_del'])]
            line['indelfos'][iseq]['reversed_seq'] = rseq
            for indel in line['indelfos'][iseq]['indels']:
                indel['pos'] -= len(fv_insertion_to_remove) + line['v_5p_del']

    line['fv_insertion'] = final_fv_insertion
    line['jf_insertion'] = final_jf_insertion

    # if padfo is None:
    #     line['padlefts'], line['padrights'] = [0 for _ in range(len(line['seqs']))], [0 for _ in range(len(line['seqs']))]
    # else:
    #     line['padlefts'], line['padrights'] = [padfo[uid]['padded']['padleft'] for uid in line['unique_ids']], [padfo[uid]['padded']['padright'] for uid in line['unique_ids']]

    # NOTE fixed the problem we were actually seeing, so this shouldn't fail any more, but I'll leave it in for a bit just in case UPDATE totally saved my ass from an unrelated problem (well, maybe not "saved" -- definitely don't remove the add_implicit_info() call though)
    try:
        add_implicit_info(glfo, line, aligned_gl_seqs=aligned_gl_seqs)
    except:
        print '%s failed adding implicit info to \'%s\'' % (color('red', 'error'), ':'.join(line['unique_ids']))
        print color('red', 'padded:')
        for k, v in padded_line.items():
            print '%20s  %s' % (k, v)
        print color('red', 'eroded:')
        for k, v in line.items():
            print '%20s  %s' % (k, v)
        line['invalid'] = True

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
def is_functional_dbg_str(line):  # NOTE code duplication with is_functional(
    dbg_str_list = []
    if True in line['mutated_invariants']:
        dbg_str_list.append('mutated invariant codon')
    if False in line['in_frames']:
        dbg_str_list.append('out of frame cdr3')
    if True in line['stops']:
        dbg_str_list.append('stop codon')
    return ', '.join(dbg_str_list)

# ----------------------------------------------------------------------------------------
def is_functional(line):  # NOTE code duplication with is_functional_dbg_str(
    if True in line['mutated_invariants']:
        return False
    if False in line['in_frames']:
        return False
    if True in line['stops']:
        return False
    return True

# ----------------------------------------------------------------------------------------
def add_functional_info(locus, line, input_codon_positions):
    nseqs = len(line['seqs'])  # would normally use 'unique_ids', but this gets called during simulation before the point at which we choose the uids
    line['mutated_invariants'] = [not both_codons_unmutated(locus, line['input_seqs'][iseq], input_codon_positions[iseq])
                                  for iseq in range(nseqs)]
    line['in_frames'] = [in_frame(line['input_seqs'][iseq], input_codon_positions[iseq], line['fv_insertion'], line['v_5p_del'])
                         for iseq in range(nseqs)]
    line['stops'] = [is_there_a_stop_codon(line['input_seqs'][iseq], line['fv_insertion'], line['jf_insertion'], line['v_5p_del'])
                     for iseq in range(nseqs)]

# ----------------------------------------------------------------------------------------
def remove_all_implicit_info(line):
    for col in implicit_linekeys:
        if col in line:
            del line[col]

# ----------------------------------------------------------------------------------------
def get_non_implicit_copy(line):  # return a deep copy of <line> with only non-implicit info
    return {col : copy.deepcopy(line[col]) for col in line if col not in implicit_linekeys}

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
            print '   %s best-supported gene %s not same as viterbi gene %s' % (color('yellow', 'warning'), color_gene(support.keys()[0]), color_gene(line[region + '_gene']))

        line[region + '_per_gene_support'] = support

# ----------------------------------------------------------------------------------------
def add_implicit_info(glfo, line, aligned_gl_seqs=None, check_line_keys=False):  # should turn on <check_line_keys> for a bit if you change anything
    """ Add to <line> a bunch of things that are initially only implicit. """
    if line['v_gene'] == '':
        raise Exception('can\'t add implicit info to line with failed annotation:\n%s' % (''.join(['  %+20s  %s\n' % (k, v) for k, v in line.items()])))

    if check_line_keys:
        initial_keys = set(line)
        # first make sure there aren't any unauthorized keys
        if len(initial_keys - all_linekeys) > 0:
            raise Exception('unexpected keys: \'%s\'' % '\' \''.join(initial_keys - all_linekeys))
        # then keep track of the keys we got to start with
        pre_existing_implicit_info = {ek : copy.deepcopy(line[ek]) for ek in implicit_linekeys if ek in line}

    for region in regions:  # backwards compatibility with old simulation files should be removed when you're no longer running on them
        if line[region + '_gene'] not in glfo['seqs'][region]:
            alternate_name = glutils.convert_to_duplicate_name(glfo, line[region + '_gene'])
            # print ' using alternate name %s instead of %s' % (alternate_name, line[region + '_gene'])
            line[region + '_gene'] = alternate_name

    # add the regional germline seqs and their lengths
    line['lengths'] = {}  # length of each match (including erosion)
    for region in regions:
        uneroded_gl_seq = glfo['seqs'][region][line[region + '_gene']]
        del_5p = line[region + '_5p_del']
        del_3p = line[region + '_3p_del']
        length = len(uneroded_gl_seq) - del_5p - del_3p  # eroded length
        if length < 0:
            raise Exception('invalid %s lengths passed to add_implicit_info()\n    gl seq: %d  5p: %d  3p: %d' % (region, len(uneroded_gl_seq), del_5p, del_3p))
        line[region + '_gl_seq'] = uneroded_gl_seq[del_5p : del_5p + length]
        line['lengths'][region] = length

    # add codon-related stuff
    line['codon_positions'] = {}
    for region, codon in conserved_codons[glfo['locus']].items():
        eroded_gl_pos = glfo[codon + '-positions'][line[region + '_gene']] - line[region + '_5p_del']
        if region == 'v':
            line['codon_positions'][region] = eroded_gl_pos + len(line['f' + region + '_insertion'])
        elif region == 'j':
            line['codon_positions'][region] = eroded_gl_pos + len(line['fv_insertion']) + line['lengths']['v'] + len(line['vd_insertion']) + line['lengths']['d'] + len(line['dj_insertion'])
        else:
            assert False
    line['cdr3_length'] = line['codon_positions']['j'] - line['codon_positions']['v'] + 3  # i.e. first base of cysteine to last base of tryptophan inclusive

    # add naive seq stuff
    line['naive_seq'] = len(line['fv_insertion']) * ambiguous_bases[0] + line['v_gl_seq'] + line['vd_insertion'] + line['d_gl_seq'] + line['dj_insertion'] + line['j_gl_seq'] + len(line['jf_insertion']) * ambiguous_bases[0]
    for iseq in range(len(line['seqs'])):
        if len(line['naive_seq']) != len(line['seqs'][iseq]):
            raise Exception('naive and mature sequences different lengths %d %d for %s:\n    %s\n    %s' % (len(line['naive_seq']), len(line['seqs'][iseq]), ' '.join(line['unique_ids']), line['naive_seq'], line['seqs'][iseq]))

    start, end = {}, {}  # add naive seq bounds for each region (could stand to make this more concise)
    start['v'] = len(line['fv_insertion'])  # NOTE this duplicates code in add_qr_seqs()
    end['v'] = start['v'] + len(line['v_gl_seq'])  # base just after the end of v
    start['d'] = end['v'] + len(line['vd_insertion'])
    end['d'] = start['d'] + len(line['d_gl_seq'])
    start['j'] = end['d'] + len(line['dj_insertion'])
    end['j'] = start['j'] + len(line['j_gl_seq'])
    line['regional_bounds'] = {r : (start[r], end[r]) for r in regions}

    input_codon_positions = [indelutils.get_codon_positions_with_indels_reinstated(line, iseq, line['codon_positions']) for iseq in range(len(line['seqs']))]
    if 'indel_reversed_seqs' not in line:  # everywhere internally, we refer to 'indel_reversed_seqs' as simply 'seqs'. For interaction with outside entities, however (i.e. writing files) we use the more explicit 'indel_reversed_seqs'
        line['indel_reversed_seqs'] = line['seqs']

    # add regional query seqs
    add_qr_seqs(line)

    add_functional_info(glfo['locus'], line, input_codon_positions)

    hfracfo = [hamming_fraction(line['naive_seq'], mature_seq, also_return_distance=True) for mature_seq in line['seqs']]
    line['mut_freqs'] = [hfrac for hfrac, _ in hfracfo]
    line['n_mutations'] = [n_mutations for _, n_mutations in hfracfo]

    # set validity (alignment addition [below] can also set invalid)  # TODO clean up this checking stuff
    line['invalid'] = False
    seq_length = len(line['seqs'][0])  # they shouldn't be able to be different lengths
    for chkreg in regions:
        if start[chkreg] < 0 or end[chkreg] < 0 or end[chkreg] < start[chkreg] or end[chkreg] > seq_length:
            line['invalid'] = True
    if end['j'] + len(line['jf_insertion']) != seq_length:
        line['invalid'] = True
    if line['cdr3_length'] < 6:  # i.e. if cyst and tryp overlap  NOTE six is also hardcoded in waterer
        line['invalid'] = True

    # add alignment info (this is only used if presto output has been specified on the command line, which requires specification of your own alignment file)
    if aligned_gl_seqs is None:  # empty/dummy info
        for region in regions:
            line['aligned_' + region + '_seqs'] = ['' for _ in range(len(line['seqs']))]
    else:
        add_alignments(glfo, aligned_gl_seqs, line)

    if check_line_keys:
        new_keys = set(line) - initial_keys
        if len(new_keys - implicit_linekeys) > 0:
            raise Exception('added new keys that aren\'t in implicit_linekeys: %s' % ' '.join(new_keys - implicit_linekeys))
        for ikey in implicit_linekeys:  # make sure every key/value we added is either a) new or b) the same as it was before
            if ikey in initial_keys:
                if pre_existing_implicit_info[ikey] != line[ikey]:
                    print '%s pre-existing info for \'%s\' in %s\n    %s\n    doesn\'t match new info\n    %s' % (color('yellow', 'warning'), ikey, line['unique_ids'], pre_existing_implicit_info[ikey], line[ikey])
            else:
                assert ikey in new_keys  # only really checks the logic of the previous few lines

# ----------------------------------------------------------------------------------------
def print_true_events(glfo, reco_info, line, print_naive_seqs=False, extra_str='    '):
    """ print the true events which contain the seqs in <line> """
    true_naive_seqs = []
    for uids in get_true_partition(reco_info, ids=line['unique_ids']):  # make a multi-seq line that has all the seqs from this clonal family
        multiline = synthesize_multi_seq_line_from_reco_info(uids, reco_info)
        print_reco_event(multiline, extra_str=extra_str, label=color('green', 'true:'))
        true_naive_seqs.append(multiline['naive_seq'])

    if print_naive_seqs:
        print '      naive sequences:'
        for tseq in true_naive_seqs:
            color_mutants(tseq, line['naive_seq'], print_result=True, print_hfrac=True, ref_label='true ', extra_str='          ')

# ----------------------------------------------------------------------------------------
def print_reco_event(line, one_line=False, extra_str='', label='', seed_uid=None):
    if len(line['unique_ids']) > 1:
        label += '%s%d sequences with %.1f mean mutations' % ('' if label == '' else '    ', len(line['unique_ids']), numpy.mean(line['n_mutations']))
    for iseq in range(len(line['unique_ids'])):
        prutils.print_seq_in_reco_event(line, iseq, extra_str=extra_str, label=(label if iseq==0 else ''), one_line=(iseq>0), seed_uid=seed_uid)

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
def get_locus(inputstr):
    """ return locus given gene or gl fname """
    locus = inputstr[:3].lower()  # only need the .lower() if it's a gene name
    if locus not in loci:
        raise Exception('couldn\'t get locus from input string \'%s\'' % inputstr)
    return locus

# ----------------------------------------------------------------------------------------
def get_region(inputstr):
    """ return v, d, or j of gene or gl fname """
    region = inputstr[3].lower()  # only need the .lower() if it's a gene name
    if region not in regions:
        raise Exception('couldn\'t get region from input string \'%s\'' % inputstr)
    return region

# ----------------------------------------------------------------------------------------
def are_alleles(gene1, gene2):
    return primary_version(gene1) == primary_version(gene2) and sub_version(gene1) == sub_version(gene2)

# ----------------------------------------------------------------------------------------
def split_gene(gene):
    """ returns (primary version, sub version, allele) """
    # make sure {IG,TR}{[HKL],[abgd]}[VDJ] is at the start, and there's a *
    if '_star_' in gene or '_slash_' in gene:
        raise Exception('gene name \'%s\' isn\'t entirely unsanitized' % gene)
    if gene[:4] != get_locus(gene).upper() + get_region(gene).upper():
        raise Exception('unexpected string in gene name %s' % gene)
    if gene.count('*') != 1:
        raise Exception('expected exactly 1 \'*\' in %s but found %d' % (gene, gene.count('*')))

    if '-' in gene and gene.find('-') < gene.find('*'):  # Js (and a few Vs) don't have sub versions
        primary_version = gene[4 : gene.find('-')]  # the bit between the IG[HKL][VDJ] and the first dash (sometimes there's a second dash as well)
        sub_version = gene[gene.find('-') + 1 : gene.find('*')]  # the bit between the first dash and the star
        allele = gene[gene.find('*') + 1 : ]  # the bit after the star
        if gene != get_locus(gene).upper() + get_region(gene).upper() + primary_version + '-' + sub_version + '*' + allele:
            raise Exception('couldn\'t build gene name %s from %s %s %s' % (gene, primary_version, sub_version, allele))
    else:
        primary_version = gene[4 : gene.find('*')]  # the bit between the IG[HKL][VDJ] and the star
        sub_version = None
        allele = gene[gene.find('*') + 1 : ]  # the bit after the star
        if gene != get_locus(gene).upper() + get_region(gene).upper() + primary_version + '*' + allele:
            raise Exception('couldn\'t build gene name %s from %s %s' % (gene, primary_version, allele))

    return primary_version, sub_version, allele

# ----------------------------------------------------------------------------------------
def rejoin_gene(locus, region, primary_version, sub_version, allele):
    """ reverse the action of split_gene() """
    return_str = locus.upper() + region.upper() + primary_version
    if sub_version is not None:  # e.g. J genes typically don't have sub-versions
        return_str += '-' + sub_version
    return return_str + '*' + allele

# ----------------------------------------------------------------------------------------
def primary_version(gene):
    return split_gene(gene)[0]

# ----------------------------------------------------------------------------------------
def gene_family(gene):  # same as primary_version(), except ignore stuff after the slash, e.g. 1/OR15 --> 1
    return primary_version(gene).split('/')[0].replace('D', '')

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
    return allelic_groups  # NOTE doesn't return the same thing as separate_into_snp_groups()


# ----------------------------------------------------------------------------------------
def separate_into_snp_groups(glfo, region, n_max_snps, genelist=None):
    """ where each class contains all alleles with the same distance from start to cyst, and within a hamming distance of <n_max_snps> """
    if genelist is None:
        genelist = glfo['seqs'][region].keys()
    assert region == 'v'  # would need to change the up-to-cpos requirement if it isn't v
    snp_groups = []
    for gene in genelist:
        cpos = cdn_pos(glfo, region, gene)
        seq = glfo['seqs'][region][gene][:cpos + 3]  # only go up through the end of the cysteine
        add_new_class = True  # to begin with, assume we'll add a new class for this gene
        for gclass in snp_groups:  # then check if, instead, this gene belongs in any of the existing classes
            for gfo in gclass:
                if len(gfo['seq']) != len(seq):  # everybody in the class has to have the same distance from start of V (rss, I think) to end of cysteine
                    continue
                hdist = hamming_distance(gfo['seq'], seq)
                if hdist < n_max_snps - 2:  # if this gene is close to any gene in the class, add it to this class
                    add_new_class = False
                    snp_groups[snp_groups.index(gclass)].append({'gene' : gene, 'seq' : seq})
                    break
            if not add_new_class:
                break

        if add_new_class:
            snp_groups.append([{'gene' : gene, 'seq' : seq}, ])

    return snp_groups  # NOTE this is a list of lists of dicts, whereas separate_into_allelic_groups() returns a dict of region-keyed dicts

# ----------------------------------------------------------------------------------------
def read_single_gene_count(indir, gene, expect_zero_counts=False, debug=False):
    region = get_region(gene)
    count = 0
    with open(indir + '/' + region + '_gene-probs.csv', 'r') as infile:  # NOTE this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
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
        with open(indir + '/' + region + '_gene-probs.csv', 'r') as infile:  # NOTE this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
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
    with open(param_dir + '/' + region + '_gene-probs.csv', 'r') as infile:  # NOTE note this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
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
def hamming_distance(seq1, seq2, extra_bases=None, return_len_excluding_ambig=False, return_mutated_positions=False, align=False):
    if align:  # way the hell slower, of course
        seq1, seq2 = align_seqs(seq1, seq2)
    if len(seq1) != len(seq2):
        raise Exception('unequal length sequences %d %d:\n  %s\n  %s' % (len(seq1), len(seq2), seq1, seq2))
    if len(seq1) == 0:
        if return_len_excluding_ambig:
            return 0, 0
        else:
            return 0

    skip_chars = set(ambiguous_bases + gap_chars)

    distance, len_excluding_ambig = 0, 0
    mutated_positions = []
    for ich in range(len(seq1)):  # already made sure they're the same length
        if seq1[ich] in skip_chars or seq2[ich] in skip_chars:
            continue
        len_excluding_ambig += 1
        if seq1[ich] != seq2[ich]:
            distance += 1
            if return_mutated_positions:
                mutated_positions.append(ich)

    if return_len_excluding_ambig and return_mutated_positions:
        assert False  # simplifies things to forbid it for the moment

    if return_len_excluding_ambig:
        return distance, len_excluding_ambig
    elif return_mutated_positions:
        return distance, mutated_positions
    else:
        return distance

    # seems like it'd maybe be faster, but it's twice as slow:
    # assert len(ambiguous_bases) == 1  # would just have to update the below if it's longer
    # ambig_base = ambiguous_bases[0]
    # if return_len_excluding_ambig:
    #     spairs = [(ch1, ch2) for ch1, ch2 in zip(seq1, seq2) if ambig_base not in ch1 + ch2]
    #     return sum(ch1 != ch2 for ch1, ch2 in spairs), len(spairs)
    # else:
    #     return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2) if ambig_base not in ch1 + ch2)

# ----------------------------------------------------------------------------------------
def hamming_fraction(seq1, seq2, extra_bases=None, also_return_distance=False):
    distance, len_excluding_ambig = hamming_distance(seq1, seq2, extra_bases=extra_bases, return_len_excluding_ambig=True)

    fraction = 0.
    if len_excluding_ambig > 0:
        fraction = distance / float(len_excluding_ambig)

    if also_return_distance:
        return fraction, distance
    else:
        return fraction

# ----------------------------------------------------------------------------------------
def subset_sequences(line, iseq, restrict_to_region):
    naive_seq = line['naive_seq']  # NOTE this includes the fv and jf insertions
    muted_seq = line['seqs'][iseq]
    if restrict_to_region != '':  # NOTE this is very similar to code in performanceplotter. I should eventually cut it out of there and combine them, but I'm nervous a.t.m. because of all the complications there of having the true *and* inferred sequences so I'm punting
        if restrict_to_region in regions:
            bounds = line['regional_bounds'][restrict_to_region]
        elif restrict_to_region == 'cdr3':
            bounds = (line['codon_positions']['v'], line['codon_positions']['j'] + 3)
        else:
            assert False
        naive_seq = naive_seq[bounds[0] : bounds[1]]
        muted_seq = muted_seq[bounds[0] : bounds[1]]
    return naive_seq, muted_seq

# ----------------------------------------------------------------------------------------
def get_n_muted(line, iseq, restrict_to_region='', return_mutated_positions=False):
    naive_seq, muted_seq = subset_sequences(line, iseq, restrict_to_region)
    return hamming_distance(naive_seq, muted_seq, return_mutated_positions=return_mutated_positions)

# ----------------------------------------------------------------------------------------
def get_mutation_rate(line, iseq, restrict_to_region=''):
    naive_seq, muted_seq = subset_sequences(line, iseq, restrict_to_region)
    return hamming_fraction(naive_seq, muted_seq)

# ----------------------------------------------------------------------------------------
def get_mutation_rate_and_n_muted(line, iseq, restrict_to_region=''):
    naive_seq, muted_seq = subset_sequences(line, iseq, restrict_to_region)
    fraction, distance = hamming_fraction(naive_seq, muted_seq, also_return_distance=True)
    return fraction, distance

# ----------------------------------------------------------------------------------------
def dot_product(naive_seq, seq1, seq2):
    _, imutes1 = hamming_distance(naive_seq, seq1, return_mutated_positions=True)
    _, imutes2 = hamming_distance(naive_seq, seq2, return_mutated_positions=True)
    both_muted = set(imutes1) & set(imutes2)
    both_muted_to_same_thing = [imut for imut in both_muted if seq1[imut] == seq2[imut]]
    dot_product = len(both_muted_to_same_thing)
    # print '    naive  %s' % naive_seq
    # print '           %s' % utils.color_mutants(naive_seq, seq1)
    # print '           %s' % utils.color_mutants(naive_seq, seq2)
    # print '    dot %d' % dot_product
    return dot_product

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
def prep_dir(dirname, wildlings=None, subdirs=None, fname=None, allow_other_files=False):
    """
    Make <dirname> if it d.n.e.
    Also, if shell glob <wildling> is specified, remove existing files which are thereby matched.
    """
    if fname is not None:  # passed in a file name, and we want to prep the file's dir
        assert dirname is None
        dirname = os.path.dirname(fname)
        if dirname == '' or dirname[0] != '/':
            dirname = '/'.join([pn for pn in [os.getcwd(), dirname] if pn != ''])

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
                if os.path.exists(fname):
                    os.remove(fname)
                else:
                    print '%s file %s exists but then it doesn\'t' % (color('red', 'wtf'), fname)
        remaining_files = [fn for fn in os.listdir(dirname) if subdirs is None or fn not in subdirs]  # allow subdirs to still be present
        if len(remaining_files) > 0 and not allow_other_files:  # make sure there's no other files in the dir
            raise Exception('files (%s) remain in %s despite wildlings %s' % (' '.join(['\'' + fn + '\'' for fn in remaining_files]), dirname, wildlings))
    else:
        os.makedirs(dirname)

# ----------------------------------------------------------------------------------------
def process_input_line(info, hmm_cachefile=False):
    """
    Attempt to convert all the keys and values in <info> according to the specifications in <column_configs> (e.g. splitting lists, casting to int/float, etc).
    """

    if 'v_gene' in info and info['v_gene'] == '':
        return

    if 'seq' in info:  # old simulation files
        for key in ['unique_id', 'seq', 'indelfo']:
            if key not in info:
                continue
            info[key + 's'] = info[key]
            del info[key]
        if 'indelfos' not in info:  # hm, at least some old sim files don't have 'indelfo'
            info['indelfos'] = str(indelutils.get_empty_indel())
        info['indelfos'] = '[' + info['indelfos'] + ']'
        info['input_seqs'] = info['seqs']

    for key in info:
        if info[key] == '':  # handle these below, once we know how many seqs in the line
            continue
        convert_fcn = conversion_fcns.get(key, pass_fcn)
        if key in column_configs['lists']:
            info[key] = [convert_fcn(val) for val in info[key].split(':')]
        elif key in column_configs['lists-of-lists']:
            info[key] = convert_fcn(info[key].split(';'))
        else:
            info[key] = convert_fcn(info[key])

    # this column is called 'seqs' internally (for conciseness and to avoid rewriting a ton of stuff) but is called 'indel_reversed_seqs' in the output file to avoid user confusion
    if 'indel_reversed_seqs' in info and 'input_seqs' in info:  # new-style csv output and simulation files, i.e. it stores 'indel_reversed_seqs' instead of 'seqs'
        if info['indel_reversed_seqs'] == '':
            info['indel_reversed_seqs'] = ['' for _ in range(len(info['unique_ids']))]
        info['seqs'] = [info['indel_reversed_seqs'][iseq] if info['indel_reversed_seqs'][iseq] != '' else info['input_seqs'][iseq] for iseq in range(len(info['unique_ids']))]  # if there's no indels, we just store 'input_seqs' and leave 'indel_reversed_seqs' empty
    elif 'seqs' in info:  # old-style csv output file: just copy 'em into the explicit name
        info['indel_reversed_seqs'] = info['seqs']

    # process things for which we first want to know the number of seqs in the line
    for key in info:
        if key == 'indelfos':
            info[key] = [reconstruct_full_indelfo(info['indelfos'][iseq], info['seqs'][iseq]) for iseq in range(len(info['unique_ids']))]
        elif info[key] != '':
            continue

        if key in column_configs['lists']:
            info[key] = ['' for _ in range(len(info['unique_ids']))]
        elif key in column_configs['lists-of-lists']:
            info[key] = [[] for _ in range(len(info['unique_ids']))]

# ----------------------------------------------------------------------------------------
def get_line_for_output(info):
    """ Reverse the action of process_input_line() """
    outfo = {}
    for key in info:
        str_fcn = str
        if key in column_configs['floats']:
            str_fcn = repr  # keeps it from losing precision (we only care because we want it to match expectation if we read it back in)
        elif 'indelfo' in key:  # just write the list of indels -- don't need the reversed seq and debug str
            str_fcn = lambda x: str([sx['indels'] for sx in x])

        if key == 'seqs':  # don't want it to be in the output dict
            continue
        if key == 'indel_reversed_seqs':  # if no indels, it's the same as 'input_seqs', so set indel_reversed_seqs to empty strings
            outfo['indel_reversed_seqs'] = ':'.join(['' if not indelutils.has_indels(info['indelfos'][iseq]) else info['indel_reversed_seqs'][iseq]
                                                     for iseq in range(len(info['unique_ids']))])
        elif key in column_configs['lists']:
            outfo[key] = ':'.join([str_fcn(v) for v in info[key]])
        elif key in column_configs['lists-of-lists']:
            outfo[key] = copy.deepcopy(info[key])
            if '_per_gene_support' in key:
                outfo[key] = outfo[key].items()
            for isl in range(len(outfo[key])):
                outfo[key][isl] = ':'.join([str_fcn(s) for s in outfo[key][isl]])
            outfo[key] = ';'.join(outfo[key])
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
        with open(infname, 'r') as sub_outfile:
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
    with open(outfname, 'w') as outfile:
        writer = csv.DictWriter(outfile, header)
        writer.writeheader()
        for line in outfo:
            writer.writerow(line)

# ----------------------------------------------------------------------------------------
def get_nodelist_from_slurm_shorthand(nodestr, known_nodes, debug=False):
    if debug:
        print '    getting nodelist from \'%s\'' % nodestr

    if '[' not in nodestr and ']' not in nodestr:  # single node (don't really need this, but maybe it's a little faster)
        return [nodestr]

    # first find the indices at which there're square braces
    nodes = []
    bracketfo = []
    ilastcomma = -1  # if the first one has brackets, the "effective" comma is at -1
    thisnodestr = ''
    for ich in range(len(nodestr)):
        ch = nodestr[ich]
        if ch == ',':
            ilastcomma = ich
        if ch == '[':
            bracketfo.append({'comma' : ilastcomma, 'ibrackets' : [ich, None]})
            thisnodestr = ''
            if debug:
                print '      start bracket    %s' % nodestr[ilastcomma + 1 : ich + 1]
        elif ch == ']':
            assert bracketfo[-1]['ibrackets'][1] is None
            bracketfo[-1]['ibrackets'][1] = ich
            thisnodestr = ''
            if debug:
                print '      end bracket      %s' % nodestr[bracketfo[-1]['ibrackets'][0] : bracketfo[-1]['ibrackets'][1] + 1]

        # if we're not within a bracket info
        if len(bracketfo) == 0 or bracketfo[-1]['ibrackets'][1] is not None:
            thisnodestr += ch
            thisnodestr = thisnodestr.strip(',[]')
            if len(thisnodestr) > 1 and ch == ',':  # if we just got to a comma, and there's something worth looking at in <thisnodestr>
                nodes.append(thisnodestr)
                thisnodestr = ''
                if debug:
                    print '      add no-bracket   %s' % nodes[-1]

    if debug:
        if len(nodes) > 0:
            print '    %d bracketless nodes: %s' % (len(nodes), ' '.join(nodes))
        if len(bracketfo) > 0:
            print '      brackets:'

    # the expand the ranges in the brackets
    for bfo in bracketfo:
        ibp = bfo['ibrackets']
        original_str = nodestr[ibp[0] + 1 : ibp[1]]  # NOTE excludes the actual bracket characters
        bracketnodes = []
        for subnodestr in original_str.split(','):
            if '-' in subnodestr:  # range of node numbers
                startstoplist = [int(i) for i in subnodestr.split('-')]
                if len(startstoplist) != 2:
                    raise Exception('wtf %s' % subnodestr)
                istart, istop = startstoplist
                bracketnodes += range(istart, istop + 1)
            else: # single node
                bracketnodes.append(int(subnodestr))
        namestr = nodestr[bfo['comma'] + 1 : ibp[0]]  # the texty bit of the name (i.e. without the numbers)
        if debug:
            print '        %s: \'%s\' --> %s' % (namestr, original_str, ' '.join([str(n) for n in bracketnodes]))
        bracketnodes = [namestr + str(i) for i in bracketnodes]
        nodes += bracketnodes

    unknown_nodes = set(nodes) - set(known_nodes)
    if len(unknown_nodes) > 0:
        print '    %s unknown nodes parsed from \'%s\': %s' % (color('yellow', 'warning'), nodestr, ' '.join(unknown_nodes))

    if debug:
        print '    %d final nodes: %s' % (len(nodes), ' '.join(nodes))

    return nodes

# ----------------------------------------------------------------------------------------
def get_available_node_core_list(batch_config_fname, debug=False):
    if debug:
        print ''
        print '  figuring out slurm config'

    our_nodes = []

    if os.getenv('SLURM_NODELIST') is None:  # not within a slurm allocation
        if debug:
            print '  not inside a slurm allocation'
        return None

    # first get info on all nodes from config file
    nodefo = {}  # node : (that node's specifications in config file)
    with open(batch_config_fname) as bfile:
        for line in bfile:
            linefo = line.strip().split()
            node = None
            if len(linefo) > 0 and linefo[0].find('NodeName=') == 0:  # node config line
                for strfo in linefo:
                    tokenval = strfo.split('=')
                    if len(tokenval) != 2:
                        raise Exception('couldn\'t parse %s into \'=\'-separated key-val pairs' % strfo)
                    key, val = tokenval
                    if ',' in val:
                        val = val.split(',')
                    if key == 'NodeName':
                        node = val
                        # if node not in our_nodes:  # damn, doesn't work
                        #     continue
                        nodefo[node] = {}
                    if node is None or node not in nodefo:
                        raise Exception('first key wasn\'t NodeName')
                    nodefo[node][key] = val
    # multiply sockets times cores/socket
    for node, info in nodefo.items():
        if 'Sockets' not in info or 'CoresPerSocket' not in info:
            raise Exception('missing keys in: %s' % ' '.join(info.keys()))
        info['nproc'] = int(info['Sockets']) * int(info['CoresPerSocket'])
    if debug:
        print '    info for %d nodes in %s' % (len(nodefo), batch_config_fname)

    our_nodes = get_nodelist_from_slurm_shorthand(os.getenv('SLURM_NODELIST'), known_nodes=nodefo.keys())
    if len(our_nodes) == 0:
        return []

    if debug:
        print '    current allocation includes %d nodes' % len(our_nodes)

    # then info on all current allocations
    quefo = {}  # node : (number of tasks allocated to that node, including ours)
    squeue_str = subprocess.check_output(['squeue', '--format', '%.18i %.2t %.6D %R'])
    headers = ['JOBID', 'ST',  'NODES', 'NODELIST(REASON)']
    for line in squeue_str.split('\n'):
        linefo = line.strip().split()
        if len(linefo) == 0:
            continue
        if linefo[0] == 'JOBID':
            assert linefo == headers
        if linefo[headers.index('ST')] != 'R':  # skip jobs that aren't running
            continue
        nodes = get_nodelist_from_slurm_shorthand(linefo[headers.index('NODELIST(REASON)')], known_nodes=nodefo.keys())
        for node in nodes:
            if node not in our_nodes:
                continue
            if node not in nodefo:
                print '  %s node %s in squeue output but not in config file %s' % (color('yellow', 'warning'), node, batch_config_fname)
                continue
            if node not in quefo:
                quefo[node] = 0
            quefo[node] += 1  # NOTE ideally this would be the number of cores slurm gave this task, rather than 1, but I can't figure out how to that info (and the docs make it sound like it might not be possible to)

    if debug:
        print '    %d "total tasks" allocated among nodes in our current allocation' % sum(quefo.values())

    # and finally, decide how many procs we can send to each node
    corelist = []
    for node in our_nodes:
        if node not in nodefo:
            raise Exception('node %s in our allocation not in config file %s' % (node, batch_config_fname))
        if node not in quefo:
            raise Exception('node %s in our allocation not in squeue output' % node)
        n_cores_we_can_use = nodefo[node]['nproc'] - quefo[node] + 1  # add one to account for the fact that quefo[node] includes our allocation
        if n_cores_we_can_use == 0:
            print '  %s huh, n_cores_we_can_use is zero' % color('yellow', 'warning')
            n_cores_we_can_use = 1
        elif n_cores_we_can_use < 0:
            print '  %s more tasks allocated to %s than available cores: %d - %d = %d (setting n_cores_we_can_use to 1 because, uh, not sure what else to do)' % (color('yellow', 'warning'), node, nodefo[node]['nproc'], quefo[node], nodefo[node]['nproc'] - quefo[node])
            n_cores_we_can_use = 1
        corelist += [node for _ in range(n_cores_we_can_use)]

    corelist = sorted(corelist)  # easier to read if it's alphabetical

    if debug:
        print '    %d available cores:' % len(corelist)
        for node in set(corelist):
            print '        %d  %s' % (corelist.count(node), node)

    if len(corelist) == 0:
        return None

    return corelist

# ----------------------------------------------------------------------------------------
def prepare_cmds(cmdfos, batch_system=None, batch_options=None, batch_config_fname=None, debug=False):
    # set cmdfo defaults
    for iproc in range(len(cmdfos)):
        if 'logdir' not in cmdfos[iproc]:  # if logdirs aren't specified, then log files go in the workdirs
            cmdfos[iproc]['logdir'] = cmdfos[iproc]['workdir']
        if 'dbgfo' not in cmdfos[iproc]:
            cmdfos[iproc]['dbgfo'] = None
        if 'env' not in cmdfos[iproc]:
            cmdfos[iproc]['env'] = None

    # get info about any existing slurm allocation
    corelist = None
    if batch_system is not None and batch_system == 'slurm' and batch_config_fname is not None:
        if not os.path.exists(batch_config_fname):
            print '  %s specified --batch-config-fname %s doesn\'t exist' % (color('yellow', 'warning'), batch_config_fname)
        else:
            corelist = get_available_node_core_list(batch_config_fname)  # list of nodes within our current allocation (empty if there isn't one), with each node present once for each core that we've been allocated on that node
            if corelist is not None and len(corelist) < len(cmdfos):
                if 1.5 * len(corelist) < len(cmdfos):
                    print '  %s many fewer cores %d than processes %d' % (color('yellow', 'warning'), len(corelist), len(cmdfos))
                    print '      corelist: %s' % ' '.join(corelist)
                while len(corelist) < len(cmdfos):
                    corelist += sorted(set(corelist))  # add each node once each time through

    if corelist is not None:
        if debug:
            print '    %d final cores for %d procs' % (len(corelist), len(cmdfos))
            print '        iproc     node'
            for iproc in range(len(cmdfos)):
                print '          %-3d    %s' % (iproc, corelist[iproc])
        assert len(corelist) >= len(cmdfos)

    return corelist

# ----------------------------------------------------------------------------------------
def run_r(cmdlines, workdir, dryrun=False, print_time=None, debug=True):
    if not os.path.exists(workdir):
        raise Exception('workdir %s doesn\'t exist' % workdir)
    cmdfname = workdir + '/mds.r'
    with open(cmdfname, 'w') as cmdfile:
        cmdfile.write('\n'.join(cmdlines) + '\n')
    # subprocess.check_call(['cat', cmdfname])
    simplerun('R --slave -f %s' % cmdfname, shell=True, print_time=print_time, swallow_stdout=True, debug=False)
    os.remove(cmdfname)

# ----------------------------------------------------------------------------------------
def simplerun(cmd_str, shell=False, dryrun=False, print_time=None, swallow_stdout=False, cmdfname=None, debug=True):
    if cmdfname is not None:
        with open(cmdfname, 'w') as cmdfile:
            cmdfile.write(cmd_str)
        subprocess.check_call(['chmod', '+x', cmdfname])
        cmd_str = cmdfname

    if debug:
        print '%s %s' % (color('red', 'run'), cmd_str)
    sys.stdout.flush()
    if dryrun:
        return
    if print_time is not None:
        start = time.time()
    runfcn = subprocess.check_output if swallow_stdout else subprocess.check_call
    runfcn(cmd_str if shell else cmd_str.split(), env=os.environ, shell=shell)
    if cmdfname is not None:
        os.remove(cmdfname)
    if print_time is not None:
        print '      %s time: %.1f' % (print_time, time.time() - start)

# ----------------------------------------------------------------------------------------
def memory_usage_fraction(debug=False):  # return fraction of total system memory that this process is using (as always with memory things, this is an approximation)
    if platform.system() != 'Linux':
        print '\n  note: utils.memory_usage_fraction() needs testing on platform \'%s\' to make sure unit conversions don\'t need changing' % platform.system()
    current_usage = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)  # kb
    total = float(psutil.virtual_memory().total) / 1000.  # returns bytes, then convert to kb
    if debug:
        print '  using %.0f / %.0f MB = %.4f' % (current_usage / 1000, total / 1000, current_usage / total)
    return current_usage / total

# ----------------------------------------------------------------------------------------
def auto_n_procs():  # for running on the local machine
    n_procs = multiprocessing.cpu_count()
    if n_procs > 10: # if it's a huge server, we probably shouldn't use all the cores
        n_procs = int(float(n_procs) / 2.)
    return n_procs

# ----------------------------------------------------------------------------------------
def run_proc_functions(procs, n_procs=None, debug=False):  # <procs> is a list of multiprocessing.Process objects
    if n_procs is None:
        n_procs = auto_n_procs()
    if debug:
        print '    running %d proc fcns with %d procs' % (len(procs), n_procs)
        sys.stdout.flush()
    while True:
        while len(procs) > 0 and len(multiprocessing.active_children()) < n_procs:
            procs[0].start()
            procs.pop(0)
        if len(multiprocessing.active_children()) == 0 and len(procs) == 0:
            break

# ----------------------------------------------------------------------------------------
def run_cmd(cmdfo, batch_system=None, batch_options=None, nodelist=None):
    cmd_str = cmdfo['cmd_str']  # don't want to modify the str in <cmdfo>
    # print cmd_str
    # sys.exit()

    fout = cmdfo['logdir'] + '/out'
    ferr = cmdfo['logdir'] + '/err'
    prefix = None
    if batch_system is not None:
        if batch_system == 'slurm':
            prefix = 'srun --nodes 1 --ntasks 1'  # --exclude=data/gizmod.txt'
            if 'threads' in cmdfo:
                prefix += ' --cpus-per-task %d' % cmdfo['threads']
            if nodelist is not None:
                prefix += ' --nodelist ' + ','.join(nodelist)
        elif batch_system == 'sge':
            prefix = 'qsub -sync y -b y -V -o ' + fout + ' -e ' + ferr
            fout = None
            ferr = None
        else:
            assert False
        if batch_options is not None:
            prefix += ' ' + batch_options
    if prefix is not None:
        cmd_str = prefix + ' ' + cmd_str

    if not os.path.exists(cmdfo['logdir']):
        os.makedirs(cmdfo['logdir'])

    # print cmd_str
    proc = subprocess.Popen(cmd_str.split(),
                            stdout=None if fout is None else open(fout, 'w'),
                            stderr=None if ferr is None else open(ferr, 'w'),
                            env=cmdfo['env'])
    return proc

# ----------------------------------------------------------------------------------------
def run_cmds(cmdfos, sleep=True, batch_system=None, batch_options=None, batch_config_fname=None, debug=None, ignore_stderr=False, n_max_tries=3):  # set sleep to False if your commands are going to run really really really quickly
    corelist = prepare_cmds(cmdfos, batch_system=batch_system, batch_options=batch_options, batch_config_fname=batch_config_fname)
    procs, n_tries = [], []
    per_proc_sleep_time = 0.01 / len(cmdfos)
    for iproc in range(len(cmdfos)):
        procs.append(run_cmd(cmdfos[iproc], batch_system=batch_system, batch_options=batch_options, nodelist=[corelist[iproc]] if corelist is not None else None))
        n_tries.append(1)
        if sleep:
            time.sleep(per_proc_sleep_time)
    while procs.count(None) != len(procs):  # we set each proc to None when it finishes
        for iproc in range(len(cmdfos)):
            if procs[iproc] is None:  # already finished
                continue
            if procs[iproc].poll() is not None:  # it just finished
                finish_process(iproc, procs, n_tries, cmdfos[iproc], n_max_tries, dbgfo=cmdfos[iproc]['dbgfo'], batch_system=batch_system, batch_options=batch_options, debug=debug, ignore_stderr=ignore_stderr)
        sys.stdout.flush()
        if sleep:
            time.sleep(per_proc_sleep_time)

# ----------------------------------------------------------------------------------------
def pad_lines(linestr, padwidth=8):
    lines = [padwidth * ' ' + l for l in linestr.split('\n')]
    return '\n'.join(lines)

# ----------------------------------------------------------------------------------------
# deal with a process once it's finished (i.e. check if it failed, and restart if so)
def finish_process(iproc, procs, n_tries, cmdfo, n_max_tries, dbgfo=None, batch_system=None, batch_options=None, debug=None, ignore_stderr=False):
    procs[iproc].communicate()
    if procs[iproc].returncode == 0:
        if not os.path.exists(cmdfo['outfname']):
            print '      proc %d succeded but its output isn\'t there, so sleeping for a bit...' % iproc
            time.sleep(0.5)
        if os.path.exists(cmdfo['outfname']):
            process_out_err(extra_str='' if len(procs) == 1 else str(iproc), dbgfo=dbgfo, logdir=cmdfo['logdir'], debug=debug, ignore_stderr=ignore_stderr)
            procs[iproc] = None  # job succeeded
            return

    # handle failure
    if n_tries[iproc] > n_max_tries:
        failstr = 'exceeded max number of tries for cmd\n    %s\nlook for output in %s and %s' % (cmdfo['cmd_str'], cmdfo['workdir'], cmdfo['logdir'])
        print failstr
        raise Exception(failstr)
    else:
        print '    proc %d try %d' % (iproc, n_tries[iproc]),
        if procs[iproc].returncode == 0 and not os.path.exists(cmdfo['outfname']):  # don't really need both the clauses
            print 'succeded but output is missing'
        else:
            print 'failed with %d (output %s)' % (procs[iproc].returncode, 'exists' if os.path.exists(cmdfo['outfname']) else 'is missing')
        for strtype in ['out', 'err']:
            if os.path.exists(cmdfo['logdir'] + '/' + strtype) and os.stat(cmdfo['logdir'] + '/' + strtype).st_size > 0:
                print '        %s tail:           (%s)' % (strtype, cmdfo['logdir'] + '/' + strtype)
                logstr = subprocess.check_output(['tail', '-n30', cmdfo['logdir'] + '/' + strtype])
                print '\n'.join(['            ' + l for l in logstr.split('\n')])
        print '    restarting proc %d' % iproc
        procs[iproc] = run_cmd(cmdfo, batch_system=batch_system, batch_options=batch_options)
        n_tries[iproc] += 1

# ----------------------------------------------------------------------------------------
def process_out_err(extra_str='', dbgfo=None, logdir=None, debug=None, ignore_stderr=False):
    """ NOTE something in this chain seems to block or truncate or some such nonsense if you make it too big """
    out, err = '', ''
    if logdir is not None:
        def read_and_delete_file(fname):
            fstr = ''
            if os.stat(fname).st_size > 0:
                ftmp = open(fname)
                fstr = ''.join(ftmp.readlines())
                ftmp.close()
            os.remove(fname)
            return fstr
        out = read_and_delete_file(logdir + '/out')
        err = read_and_delete_file(logdir + '/err')

    for line in out.split('\n'):  # temporarily (maybe) print debug info realted to --n-final-clusters/force merging
        if 'force' in line:
            print '    %s %s' % (color('yellow', 'force info:'), line)

    err_str = ''
    for line in err.split('\n'):
        if 'stty: standard input: Inappropriate ioctl for device' in line:
            continue
        if 'srun: job' in line and 'queued and waiting for resources' in line:
            continue
        if 'srun: job' in line and 'has been allocated resources' in line:
            continue
        if 'srun: Required node not available (down, drained or reserved)' in line:
            continue
        if 'GSL_RNG_TYPE=' in line or 'GSL_RNG_SEED=' in line:
            continue
        if '[ig_align] Read' in line or '[ig_align] Aligned' in line:
            continue
        if len(line.strip()) > 0:
            err_str += line + '\n'

    if dbgfo is not None:  # keep track of how many vtb and fwd calculations the process made
        for header, variables in bcrham_dbgstrs.items():
            dbgfo[header] = {var : None for var in variables}
            theselines = [ln for ln in out.split('\n') if header + ':' in ln]
            if len(theselines) == 0:
                continue
            if len(theselines) > 1:
                raise Exception('too many lines with dbgfo for \'%s\' in:\nstdout:\n%s\nstderr:\n%s' % (header, out, err))
            words = theselines[0].split()
            for var in variables:  # convention: value corresponding to the string <var> is the word immediately vollowing <var>
                if var in  words:
                    dbgfo[header][var] = float(words[words.index(var) + 1])

    if debug is None:
        if not ignore_stderr and err_str != '':
            print err_str
    elif err_str + out != '':
        if debug == 'print':
            if extra_str != '':
                print '      --> proc %s' % extra_str
            print err_str + out
        elif debug == 'write':
            logfile = logdir + '/log'
            # print 'writing dbg to %s' % logfile
            with open(logfile, 'w') as dbgfile:
                dbgfile.write(err_str + out)
        else:
            assert False

# ----------------------------------------------------------------------------------------
def summarize_bcrham_dbgstrs(dbgfos):
    def defval(dbgcat):
        if dbgcat in bcrham_dbgstr_types['sum']:
            return 0.
        elif dbgcat in bcrham_dbgstr_types['same']:
            return None
        elif dbgcat in bcrham_dbgstr_types['min-max']:
            return []
        else:
            assert False

    summaryfo = {dbgcat : {vtype : defval(dbgcat) for vtype in tlist} for dbgcat, tlist in bcrham_dbgstrs.items()}  # fill summaryfo with default/initial values
    for procfo in dbgfos:
        for dbgcat in bcrham_dbgstr_types['same']:  # loop over lines in output for which every process should have the same values (e.g. cache-read)
            for vtype in bcrham_dbgstrs[dbgcat]:  # loop over values in that line (e.g. logprobs and naive-seqs)
                if summaryfo[dbgcat][vtype] is None:  # first one
                    summaryfo[dbgcat][vtype] = procfo[dbgcat][vtype]
                if procfo[dbgcat][vtype] != summaryfo[dbgcat][vtype]:  # make sure all subsequent ones are the same
                    print '        %s bcrham procs had different \'%s\' \'%s\' info: %d vs %d' % (color('red', 'warning'), vtype, dbgcat, procfo[dbgcat][vtype], summaryfo[dbgcat][vtype])
        for dbgcat in bcrham_dbgstr_types['sum']:  # lines for which we want to add up the values
            for vtype in bcrham_dbgstrs[dbgcat]:
                summaryfo[dbgcat][vtype] += procfo[dbgcat][vtype]
        for dbgcat in bcrham_dbgstr_types['min-max']:  # lines for which we want to keep track of the smallest and largest values (e.g. time required)
            for vtype in bcrham_dbgstrs[dbgcat]:
                summaryfo[dbgcat][vtype].append(procfo[dbgcat][vtype])

    for dbgcat in bcrham_dbgstr_types['min-max']:
        for vtype in bcrham_dbgstrs[dbgcat]:
            summaryfo[dbgcat][vtype] = min(summaryfo[dbgcat][vtype]), max(summaryfo[dbgcat][vtype])

    return summaryfo

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
    Group queries into their true clonal families.
    If <ids> is specified, only do those, otherwise do all of the ones in <reco_info>.
    """
    if ids is None:
        ids = reco_info.keys()
    def keyfunc(q):
        return reco_info[q]['reco_id']
    return [list(group) for _, group in itertools.groupby(sorted(ids, key=keyfunc), key=keyfunc)]  # sort 'em beforehand so all the ones with the same reco id are consecutive (if there are non-consecutive ones with the same reco id, it means there were independent events with the same rearrangment parameters)

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
    clids = {uid : [] for uid in uids}  # almost always list of length one with index (in <partition>) of the uid's cluster
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
    clids = get_cluster_ids(uids, partition)  # map of {uid : (index of cluster in <partition> in which that uid occurs)} (well, list of indices, in case there's duplicates)

    def get_clonal_fraction(uid, inferred_cluster):
        """ Return the fraction of seqs in <uid>'s inferred cluster which are really clonal. """
        n_clonal = 0
        for tmpid in inferred_cluster:  # NOTE this includes the case where tmpid equal to uid
            if reco_ids[tmpid] == reco_ids[uid]:  # reminder (see event.py) reco ids depend only on rearrangement parameters, i.e. two different rearrangement events with the same rearrangement parameters have the same reco id
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
            inferred_cluster = partition[clids[uid][0]]  # we only look at the first cluster in which it appears
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
def csv_to_fasta(infname, outfname=None, name_column='unique_ids', seq_column='input_seqs', n_max_lines=None, overwrite=True, remove_duplicates=False, debug=True):

    if not os.path.exists(infname):
        raise Exception('input file %s d.n.e.' % infname)
    if outfname is None:
        assert '.csv' in infname
        outfname = infname.replace('.csv', '.fa')
    if os.path.exists(outfname):
        if overwrite:
            if debug:
                print '  csv --> fasta: overwriting %s' % outfname
        else:
            if debug:
                print '  csv --> fasta: leaving existing outfile %s' % outfname
            return

    if '.csv' in infname:
        delimiter = ','
    elif '.tsv' in infname:
        delimiter = '\t'
    else:
        assert False

    uid_set = set()
    n_duplicate_ids = 0
    with open(infname) as infile:
        reader = csv.DictReader(infile, delimiter=delimiter)
        with open(outfname, 'w') as outfile:
            n_lines = 0
            for line in reader:
                if seq_column not in line:
                    raise Exception('specified <seq_column> \'%s\' not in line (keys in line: %s)' % (seq_column, ' '.join(line.keys())))
                if name_column is not None:
                    if name_column not in line:
                        raise Exception('specified <name_column> \'%s\' not in line (keys in line: %s)' % (name_column, ' '.join(line.keys())))
                    uid = line[name_column]
                else:
                    uid = str(abs(hash(line[seq_column])))
                if remove_duplicates:
                    if uid in uid_set:
                        n_duplicate_ids += 1
                        continue
                    uid_set.add(uid)
                n_lines += 1
                if n_max_lines is not None and n_lines > n_max_lines:
                    break
                outfile.write('>%s\n' % uid)
                outfile.write('%s\n' % line[seq_column])
    if debug and n_duplicate_ids > 0:
        print '   skipped %d / %d duplicate uids' % (n_duplicate_ids, len(uid_set))

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
            subprocess.check_output(['which', 'srun'], stderr=fnull, close_fds=True)
            return True
        except subprocess.CalledProcessError:
            return False

    ncpu = multiprocessing.cpu_count()
    if n_procs > ncpu and slurm_exists():
        return True
    return False

# ----------------------------------------------------------------------------------------
def add_regional_alignments(glfo, aligned_gl_seqs, line, region, debug=False):
    if debug:
        print ' %s' % region

    aligned_seqs = [None for _ in range(len(line['unique_ids']))]
    for iseq in range(len(line['seqs'])):
        qr_seq = line[region + '_qr_seqs'][iseq]
        gl_seq = line[region + '_gl_seq']
        aligned_gl_seq = aligned_gl_seqs[region][line[region + '_gene']]
        if len(qr_seq) != len(gl_seq):
            if debug:
                print '    qr %d and gl %d seqs different lengths for %s, setting invalid' % (len(qr_seq), len(gl_seq), ' '.join(line['unique_ids']))
            line['invalid'] = True
            continue

        n_gaps = gap_len(aligned_gl_seq)
        if n_gaps == 0:
            if debug:
                print '   no gaps'
            aligned_seqs[iseq] = qr_seq
            continue

        if debug:
            print '   before alignment'
            print '      qr   ', qr_seq
            print '      gl   ', gl_seq
            print ' aligned gl', aligned_gl_seq

        # add dots for 5p and 3p deletions
        qr_seq = gap_chars[0] * line[region + '_5p_del'] + qr_seq + gap_chars[0] * line[region + '_3p_del']
        gl_seq = gap_chars[0] * line[region + '_5p_del'] + gl_seq + gap_chars[0] * line[region + '_3p_del']

        if len(aligned_gl_seq) - n_gaps != len(gl_seq):
            if debug:
                print '    aligned germline seq without gaps (%d - %d = %d) not the same length as unaligned gl/qr seqs %d' % (len(aligned_gl_seq), n_gaps, len(aligned_gl_seq) - n_gaps, len(gl_seq))
            line['invalid'] = True
            continue

        qr_seq = list(qr_seq)
        gl_seq = list(gl_seq)
        for ibase in range(len(aligned_gl_seq)):
            if aligned_gl_seq[ibase] in gap_chars:  # add gap to the qr and gl seq lists
                qr_seq.insert(ibase, gap_chars[0])
                gl_seq.insert(ibase, gap_chars[0])
            elif gl_seq[ibase] == aligned_gl_seq[ibase] or gl_seq[ibase] in gap_chars:  # latter is 5p or 3p deletion that we filled in above
                pass  # all is well
            else:  # all is not well, don't know why
                line['invalid'] = True
                break
        if line['invalid']:
            if debug:
                print '    unknown error during alignment process'
            continue
        qr_seq = ''.join(qr_seq)
        gl_seq = ''.join(gl_seq)

        if debug:
            print '   after alignment'
            print '      qr   ', qr_seq
            print '      gl   ', gl_seq
            print ' aligned gl', aligned_gl_seq

        if len(qr_seq) != len(gl_seq) or len(qr_seq) != len(aligned_gl_seq):  # I don't think this is really possible as currently written
            if debug:
                print '    lengths qr %d gl %d and aligned gl %d not all the same after alignment' % (len(qr_seq), len(gl_seq), len(aligned_gl_seq))
            line['invalid'] = True
            continue

        aligned_seqs[iseq] = qr_seq

    if line['invalid']:
        print '%s failed adding alignment info for %s' % (color('red', 'error'),' '.join(line['unique_ids']))  # will print more than once if it doesn't fail on the last region
        aligned_seqs = [None for _ in range(len(line['seqs']))]

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
def choose_seed_unique_id(gldir, locus, simfname, seed_cluster_size_low, seed_cluster_size_high, iseed=None, n_max_queries=-1, debug=True):
    glfo = glutils.read_glfo(gldir, locus)
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

# ----------------------------------------------------------------------------------------
def remove_from_arglist(clist, argstr, has_arg=False):
    if argstr not in clist:
        return
    if has_arg:
        clist.pop(clist.index(argstr) + 1)
    clist.remove(argstr)

# ----------------------------------------------------------------------------------------
def replace_in_arglist(clist, argstr, replace_with):
    if argstr not in clist:
        clist.append(argstr)
        clist.append(replace_with)
    else:
        clist[clist.index(argstr) + 1] = replace_with

# ----------------------------------------------------------------------------------------
def kbound_str(kbounds):
    return_str = []
    for region in ['v', 'd']:
        rkb = kbounds[region]
        return_str.append('k_%s %d' % (region, rkb['best']))
        if 'min' in rkb and 'max' in rkb:
            return_str.append(' [%s-%s)' % (str(rkb.get('min', ' ')), str(rkb.get('max', ' '))))
        return_str.append('  ')
    return ''.join(return_str).strip()

# ----------------------------------------------------------------------------------------
def split_partition_with_criterion(partition, criterion_fcn):
    true_cluster_indices = [ic for ic in range(len(partition)) if criterion_fcn(partition[ic])]  # indices of clusters for which <criterion_fcn()> is true
    true_clusters = [partition[ic] for ic in true_cluster_indices]
    false_clusters = [partition[ic] for ic in range(len(partition)) if ic not in true_cluster_indices]
    return true_clusters, false_clusters

# ----------------------------------------------------------------------------------------
def collapse_naive_seqs(swfo, queries=None):
    if queries is None:
        queries = swfo['queries']  # don't modify this
    def keyfunc(q):
        return swfo[q]['naive_seq']  # while this is no longer happening before fwk insertion trimming (which was bad), it is still happening on N-padded sequences, which should be kept in mind
    return [list(group) for _, group in itertools.groupby(sorted(queries, key=keyfunc), key=keyfunc)]

# ----------------------------------------------------------------------------------------
def collapse_naive_seqs_with_hashes(naive_seq_list, sw_info):
    naive_seq_map = {}  # X[cdr3][hash(naive_seq)] : naive_seq
    naive_seq_hashes = {}  # X[hash(naive_seq)] : [uid1, uid2, uid3...]
    for uid, naive_seq in naive_seq_list:
        hashstr = str(hash(naive_seq))
        if hashstr not in naive_seq_hashes:  # first sequence that has this naive
            cdr3_length = sw_info[uid]['cdr3_length']
            if cdr3_length not in naive_seq_map:
                naive_seq_map[cdr3_length] = {}
            naive_seq_map[cdr3_length][hashstr] = naive_seq  # i.e. vsearch gets a hash of the naive seq (which maps to a list of uids with that naive sequence) instead of the uid
            naive_seq_hashes[hashstr] = []
        naive_seq_hashes[hashstr].append(uid)
    print '        collapsed %d sequences into %d unique naive sequences' % (len(naive_seq_list), len(naive_seq_hashes))
    return naive_seq_map, naive_seq_hashes

# ----------------------------------------------------------------------------------------
def read_fastx(fname, name_key='name', seq_key='seq', add_info=True, sanitize=False, queries=None, n_max_queries=-1, istartstop=None, ftype=None, n_random_queries=None):  # Bio.SeqIO takes too goddamn long to import
    if ftype is None:
        suffix = getsuffix(fname)
        if suffix == '.fa' or suffix == '.fasta':
            ftype = 'fa'
        elif suffix == '.fq' or suffix == '.fastq':
            ftype = 'fq'
        else:
            raise Exception('unhandled file type: %s' % suffix)

    finfo = []
    iline = -1  # index of the query/seq that we're currently reading in the fasta
    n_fasta_queries = 0  # number of queries so far added to <finfo> (I guess I could just use len(finfo) at this point)
    missing_queries = set(queries) if queries is not None else None
    already_printed_forbidden_character_warning = False
    with open(fname) as fastafile:
        startpos = None
        while True:
            if startpos is not None:  # rewind since the last time through we had to look to see when the next header line appeared
                fastafile.seek(startpos)
            headline = fastafile.readline()
            if not headline:
                break
            if headline.strip() == '':  # skip a blank line
                headline = fastafile.readline()

            if ftype == 'fa':
                if headline[0] != '>':
                    raise Exception('invalid fasta header line in %s:\n    %s' % (fname, headline))
                headline = headline.lstrip('>')

                seqlines = []
                nextline = fastafile.readline()
                while True:
                    if not nextline:
                        break
                    if nextline[0] == '>':
                        break
                    else:
                        startpos = fastafile.tell()  # i.e. very line that doesn't begin with '>' increments <startpos>
                    seqlines.append(nextline)
                    nextline = fastafile.readline()
                seqline = ''.join([l.strip() for l in seqlines]) if len(seqlines) > 0 else None
            elif ftype == 'fq':
                if headline[0] != '@':
                    raise Exception('invalid fastq header line in %s:\n    %s' % (fname, headline))
                headline = headline.lstrip('@')

                seqline = fastafile.readline()  # NOTE .fq with multi-line entries isn't supported, since delimiter characters are allowed to occur within the quality string
                plusline = fastafile.readline().strip()
                if plusline[0] != '+':
                    raise Exception('invalid fastq quality header in %s:\n    %s' % (fname, plusline))
                qualityline = fastafile.readline()
            else:
                raise Exception('unhandled ftype %s' % ftype)

            if not seqline:
                break

            iline += 1
            if istartstop is not None:
                if iline < istartstop[0]:
                    continue
                elif iline >= istartstop[1]:
                    continue

            infostrs = [ss.strip() for s in headline.split(' ') for ss in s.split('|')]  # NOTE the uid is left untranslated in here
            uid = infostrs[0]
            if sanitize and any(fc in uid for fc in forbidden_characters):
                if not already_printed_forbidden_character_warning:
                    print '  %s: found a forbidden character (one of %s) in sequence id \'%s\'. This means we\'ll be replacing each of these forbidden characters with a single letter from their name (in this case %s). If this will cause problems you should replace the characters with something else beforehand.' % (color('yellow', 'warning'), ' '.join(["'" + fc + "'" for fc in forbidden_characters]), uid, uid.translate(forbidden_character_translations))
                    already_printed_forbidden_character_warning = True
                uid = uid.translate(forbidden_character_translations)

            if queries is not None:
                if uid not in queries:
                    continue
                missing_queries.remove(uid)

            seqfo = {name_key : uid, seq_key : seqline.upper()}
            if add_info:
                seqfo['infostrs'] = infostrs
            finfo.append(seqfo)

            n_fasta_queries += 1
            if n_max_queries > 0 and n_fasta_queries >= n_max_queries:
                break
            if queries is not None and len(missing_queries) == 0:
                break

    if n_random_queries is not None:
        finfo = numpy.random.choice(finfo, n_random_queries, replace=False)

    return finfo

# ----------------------------------------------------------------------------------------
def output_exists(args, outfname, offset=22):
    if os.path.exists(outfname):
        if os.stat(outfname).st_size == 0:
            print '%sdeleting zero length %s' % (offset * ' ', outfname)
            os.remove(outfname)
            return False
        elif args.overwrite:
            print '%soverwriting %s' % (offset * ' ', outfname)
            if os.path.isdir(outfname):
                raise Exception('output %s is a directory, rm it by hand' % outfname)
            else:
                os.remove(outfname)
            return False
        else:
            print '%soutput exists, skipping (%s)' % (offset * ' ', outfname)
            return True
    else:
        return False

# ----------------------------------------------------------------------------------------
def getprefix(fname):  # basename before the dot
    if len(os.path.splitext(fname)) != 2:
        raise Exception('couldn\'t split %s into two pieces using dot' % fname)
    return os.path.splitext(fname)[0]

# ----------------------------------------------------------------------------------------
def getsuffix(fname):  # basename before the dot
    if len(os.path.splitext(fname)) != 2:
        raise Exception('couldn\'t split %s into two pieces using dot' % fname)
    return os.path.splitext(fname)[1]

# ----------------------------------------------------------------------------------------
def read_vsearch_cluster_file(fname):
    id_clusters = {}
    with open(fname) as clusterfile:
        reader = csv.DictReader(clusterfile, fieldnames=['type', 'cluster_id', '3', '4', '5', '6', '7', 'crap', 'query', 'morecrap'], delimiter='\t')
        for line in reader:
            if line['type'] == 'C':  # some lines are a cluster, and some are a query sequence. Skip the cluster ones.
                continue
            cluster_id = int(line['cluster_id'])
            if cluster_id not in id_clusters:
                id_clusters[cluster_id] = []
            uid = line['query']
            id_clusters[cluster_id].append(uid)
    partition = id_clusters.values()
    return partition

# ----------------------------------------------------------------------------------------
def run_vsearch(action, seqs, workdir, threshold, consensus_fname=None, msa_fname=None, glfo=None, print_time=False, vsearch_binary=None):
    # single-pass, greedy, star-clustering algorithm with
    #  - add the target to the cluster if the pairwise identity with the centroid is higher than global threshold <--id>
    #  - pairwise identity definition <--iddef> defaults to: number of (matching columns) / (alignment length - terminal gaps)
    #  - the search process sorts sequences in decreasing order of number of k-mers in common
    #    - the search process stops after --maxaccept matches (default 1), and gives up after --maxreject non-matches (default 32)
    #    - If both are zero, it searches the whole database
    #    - I do not remember why I set both to zero. I just did a quick test, and on a few thousand sequences, it seems to be somewhat faster with the defaults, and a tiny bit less accurate.
    region = 'v'
    start = time.time()

    prep_dir(workdir)
    infname = workdir + '/input.fa'

    # write input
    with open(infname, 'w') as fastafile:
        for name, seq in seqs.items():
            fastafile.write('>' + name + '\n' + seq + '\n')

    # run
    if vsearch_binary is None:
        vsearch_binary = os.path.dirname(os.path.realpath(__file__)).replace('/python', '') + '/bin'
        if platform.system() == 'Linux':
            vsearch_binary += '/vsearch-2.4.3-linux-x86_64'
        elif platform.system() == 'Darwin':
            vsearch_binary += '/vsearch-2.4.3-macos-x86_64'
        else:
            raise Exception('%s no vsearch binary in bin/ for platform \'%s\' (you can specify your own full vsearch path with --vsearch-binary)' % (color('red', 'error'), platform.system()))

    cmd = vsearch_binary
    cmd += ' --id ' + str(1. - threshold)
    # cmd += ' --match '  # default 2
    # cmd += ' --mismatch '  # default -4
    if action == 'cluster':
        outfname = workdir + '/vsearch-clusters.txt'
        cmd += ' --cluster_fast ' + infname
        cmd += ' --uc ' + outfname
        if consensus_fname is not None:  # workdir cleanup below will fail if you put it in this workdir
            cmd += ' --consout ' + consensus_fname  # note: can also output a file with msa and consensus
        if msa_fname is not None:  # workdir cleanup below will fail if you put it in this workdir
            cmd += ' --msaout ' + msa_fname
        # cmd += ' --maxaccept 0 --maxreject 0'  # see note above
    elif action == 'search':
        outfname = workdir + '/aln-info.tsv'
        userfields = ['query', 'target', 'qilo', 'qihi', 'ids']
        dbdir = workdir + '/' + glutils.glfo_dir
        glutils.write_glfo(dbdir, glfo)
        cmd += ' --usearch_global ' + infname
        cmd += ' --maxaccepts 5'  # it's sorted by number of k-mers in common, so this needs to be large enough that we'll almost definitely get the exact best gene match
        cmd += ' --db ' + glutils.get_fname(dbdir, glfo['locus'], region)
        cmd += ' --userfields %s --userout %s' % ('+'.join(userfields), outfname)  # note that --sizeout --dbmatched <fname> adds up all the matches from --maxaccepts, i.e. it's not what we want
    else:
        assert False

    # cmd += ' --threads ' + str(n_procs)  # a.t.m. just let vsearch use all the cores (it'd be nice to be able to control it a little, but it's hard to keep it separate from the number of slurm procs, and we pretty much always want it to be different to that)
    cmd += ' --quiet'
    cmdfos = [{
        'cmd_str' : cmd,
        'outfname' : outfname,
        'workdir' : workdir,
        # 'threads' : n_procs},  # NOTE that this does something very different (adjusts slurm command) to the line above ^ (talks to vsearch)
    }]
    run_cmds(cmdfos)

    # read output
    if action == 'cluster':
        returnfo = read_vsearch_cluster_file(outfname)
    elif action == 'search':
        query_info = {}
        with open(outfname) as alnfile:
            reader = csv.DictReader(alnfile, fieldnames=userfields, delimiter='\t')  # NOTE start/end positions are 1-indexed
            for line in reader:
                query = line['query']
                istart = int(line['qilo']) - 1  # qilo/qihi: first/last nucleotide of query aligned with target (1-indexed, ignores initial gaps)
                length = int(line['qihi']) - int(line['qilo']) + 1
                qr_seq = seqs[query][istart : istart + length]
                id_score = int(line['ids'])
                if query not in query_info:  # note that a surprisingly large number of targets give the same score, and you seem to get a little closer to what sw does if you sort alphabetically, but it don't matter
                    query_info[query] = []
                query_info[query].append({'ids' : id_score,
                                          'gene' : line['target'],
                                          'qr_seq' : qr_seq})

        glutils.remove_glfo_files(dbdir, glfo['locus'])
        if len(query_info) == 0:
            raise Exception('didn\'t read anything from vsearch output file %s from cmd %s' % (outfname, cmd))

        for query in query_info:
            if len(query_info[query]) == 0:
                print '%s zero vsearch matches for query %s' % (color('yellow', 'warning'), query)
                del query_info[query]  # uh... need to handle failures better than this
                continue
            query_info[query] = sorted(query_info[query], key=lambda d: d['ids'], reverse=True)
            best_score = query_info[query][0]['ids']
            query_info[query] = [qinfo for qinfo in query_info[query] if qinfo['ids'] == best_score]  # keep all the matches that have the same score as the best match

        gene_counts = {}
        for query in query_info:
            counts_per_match = 1. / len(query_info[query])  # e.g. if there's four matches with the same score, give 'em each 0.25 counts
            for qinfo in query_info[query]:
                if qinfo['gene'] not in gene_counts:
                    gene_counts[qinfo['gene']] = 0.
                gene_counts[qinfo['gene']] += counts_per_match
        regional_mute_freq = numpy.mean([float(query_info[q][0]['ids']) / len(query_info[q][0]['qr_seq']) for q in query_info])
        # TODO I think I don't need the 'mute-freqs' info any more (?)
        returnfo = {'gene-counts' : gene_counts, 'queries' : query_info, 'mute-freqs' : {region : regional_mute_freq}}  # NOTE <gene_counts> is note integers
    else:
        assert False
    os.remove(infname)
    os.remove(outfname)
    os.rmdir(workdir)
    if print_time:
        print '      vsearch: %d %s annotations in %.1f sec' % (len(seqs), region, time.time() - start)

    return returnfo

# ----------------------------------------------------------------------------------------
def run_swarm(seqs, workdir, differences=1, n_procs=1):
    # groups together all sequence pairs that have <d> or fewer differences (--differences, default 1)
    #  - if d=1, uses algorithm of linear complexity (d=2 or greater uses quadratic algorithm)
    #  - --fastidious (only for d=1) extra pass to reduce the number of small OTUs

    prep_dir(workdir)
    infname = workdir + '/input.fa'
    outfname = workdir + '/clusters.txt'

    dummy_abundance = 1
    with open(infname, 'w') as fastafile:
        for name, seq in seqs.items():
            fastafile.write('>%s_%d\n%s\n' % (name, dummy_abundance, remove_ambiguous_ends(seq).replace('N', 'A')))

    cmd = os.path.dirname(os.path.realpath(__file__)).replace('/python', '') + '/bin/swarm-2.1.13-linux-x86_64 ' + infname
    cmd += ' --differences ' + str(differences)
    if differences == 1:
        cmd += ' --fastidious'
    cmd += ' --threads ' + str(n_procs)
    cmd += ' --output-file ' + outfname
    simplerun(cmd)

    partition = []
    with open(outfname) as outfile:
        for line in outfile:
            partition.append([uidstr.rstrip('_%s' % dummy_abundance) for uidstr in line.strip().split()])

    os.remove(infname)
    os.remove(outfname)
    os.rmdir(workdir)

    from clusterpath import ClusterPath
    cp = ClusterPath()
    cp.add_partition(partition, logprob=0., n_procs=1)
    cp.print_partitions(abbreviate=True)

    return partition


# # ----------------------------------------------------------------------------------------
# def run_mds(seqfos, workdir, plotdir, reco_info=None, title='', debug=False):
#     # ----------------------------------------------------------------------------------------
#     def read_clusterfile(clusterfname, seqfos):
#         all_uids = set([sfo['name'] for sfo in seqfos])

#         partition = []
#         with open(clusterfname) as clusterfile:
#             lines = [l.strip() for l in clusterfile.readlines()]
#             lines = [lines[i : i + 4] for i in range(0, len(lines), 4)]
#             for clusterlines in lines:
#                 clidline = clusterlines[0]
#                 uidline = clusterlines[1]
#                 intline = clusterlines[2]  # some info about the kmean cluster quality i think? don't care a.t.m.
#                 emptyline = clusterlines[3]

#                 if clidline[0] != '$' or int(clidline.lstrip('$').strip('`')) != len(partition) + 1:
#                     raise Exception('couldn\'t convert %s to the expected cluster id %d' % (clidline, len(partition) + 1))

#                 uids = set([u for u in uidline.split()])
#                 if len(uids - all_uids) > 0:
#                     raise Exception('read unexpected uid[s] \'%s\' from %s' % (' '.join(uids - all_uids), clusterfname))
#                 all_uids -= uids
#                 partition.append(list(uids))

#                 integers = [int(istr) for istr in intline.split()]
#                 if len(integers) != len(uids):
#                     raise Exception('uid line %d and integers line %d have different lengths:\n  %s\n  %s' % (len(uids), len(integers), uidline, intline))

#                 if emptyline != '':
#                     raise Exception('expected empty line but got \'%s\'' % emptyline)

#         os.remove(clusterfname)
#         return partition

#     debug = True
#     region = 'v'

#     if not os.path.exists(workdir):
#         os.makedirs(workdir)
#     if not os.path.exists(plotdir):
#         os.makedirs(plotdir)

#     msafname = workdir + '/msa.fa'
#     group_csv_fname = workdir + '/groups.csv'
#     clusterfname = workdir + '/clusters.txt'

#     # R does some horrible truncation or some bullshit when it reads the group csv
#     chmap = ['0123456789', 'abcdefghij']
#     translations = string.maketrans(*chmap)
#     reverse_translations = string.maketrans(*reversed(chmap))
#     def translate(name):
#         return name.translate(translations)
#     def untranslate(trans_name):
#         return trans_name.translate(reverse_translations)

#     print '%s stop modifying <seqfos>' % color('red', 'hey')
#     for seqfo in seqfos:
#         # print seqfo['name'], translate(seqfo['name']), untranslate(translate(seqfo['name']))
#         seqfo['name'] = translate(seqfo['name'])

#     align_many_seqs(seqfos, outfname=msafname)

#     if reco_info is not None:
#         colors = ['red', 'blue', 'forestgreen', 'grey', 'orange', 'green', 'skyblue4', 'maroon', 'salmon', 'chocolate4', 'magenta']
#         all_genes = list(set([reco_info[untranslate(seqfo['name'])][region + '_gene'] for seqfo in seqfos]))
#         if len(all_genes) > len(colors):
#             print '%s more genes %d than colors %d' % (color('yellow', 'warning'), len(all_genes), len(colors))
#         gene_colors = {all_genes[ig] : colors[ig % len(colors)] for ig in range(len(all_genes))}

#         all_genes = set([reco_info[untranslate(seqfo['name'])][region + '_gene'] for seqfo in seqfos])  # R code crashes if there's only one group
#         with open(group_csv_fname, 'w') as groupfile:
#             for iseq in range(len(seqfos)):
#                 seqfo = seqfos[iseq]
#                 gene = reco_info[untranslate(seqfo['name'])][region + '_gene']
#                 if len(all_genes) == 1 and iseq == 0:  # see note above (!#$@!$Q)
#                     gene += '-dummy'
#                 groupfile.write('"%s","%s","%s"\n' % (seqfo['name'], gene, gene_colors.get(gene, 'black')))

#     cmdlines = [
#         # # functional example:
#         # 'require(bios2mds, quietly=TRUE)',
#         # 'data(gpcr)',
#         # 'human <- import.fasta(system.file("msa/human_gpcr.fa", package="bios2mds"))',
#         # 'active <- gpcr$dif$sapiens.sapiens',
#         # 'mmds_active <- mmds(active, group.file=system.file("csv/human_gpcr_group.csv", package = "bios2mds"))',
#         # # 'mmds_active <- mmds(active, group.file="/home/dralph/work/partis/bios2mds/inst/csv/human_gpcr_group.csv")',
#         # 'layout(matrix(1:6, 2, 3))',
#         # 'scree.plot(mmds_active$eigen.perc, lab = TRUE, title = "Scree plot of metric MDS", pdf.file="%s/scree.pdf")' % plotdir,
#         # 'mmds.2D.plot(mmds_active, title = "Sequence space of human GPCRs", outfile.name="%s/mmds-2d", outfile.type="pdf")' % plotdir,

#         # ----------------------------------------------------------------------------------------
#         'require(bios2mds, quietly=TRUE)',
#         # set.seed(1503941627)
#         'human <- import.fasta("%s")' % msafname, #system.file("msa/human_gpcr.fa", package="bios2mds"))',

#         # mat.dif or mat.dis?
#         'active <- mat.dif(human, human)',
#         'kmeans.run1 <- kmeans.run(active, nb.clus = 3, nb.run = 100)',
#         # 'kmeans.run1$clusters',
#         # 'kmeans.run1$elements',
#         'options(width=10000)',
#         'capture.output(kmeans.run1$clusters, file="%s")' % clusterfname,

#         'mmds_active <- mmds(active, group.file=%s)' % ('NULL' if reco_info is None else '"' + group_csv_fname + '"'),

#         # 'layout(matrix(1:6, 2, 3))',

#         'scree.plot(mmds_active$eigen.perc, lab=TRUE, title="%s", pdf.file="%s/scree.pdf")' % (title, plotdir),
#         'mmds.2D.plot(mmds_active, title="%s", outfile.name="%s/mmds-2d", outfile.type="pdf")' % (title, plotdir),

#         # sil.score(mat, nb.clus = c(2:13), nb.run = 100, iter.max = 1000,  # run for every possible number of clusters (?)
#         #               method = "euclidean")
#         # mmds.plot(mmds_active) #, pdf.file="")  # does several of the above steps in one go
#         # random.msa  # builds a random [...]

#         # write.mmds.pdb(mmds_active)
#     ]

#     cmdfname = workdir + '/mds.r'
#     with open(cmdfname, 'w') as cmdfile:
#         cmdfile.write('\n'.join(cmdlines) + '\n')
#     # subprocess.check_call(['cat', cmdfname])
#     subprocess.check_call('R --slave -f %s' % cmdfname, shell=True)
#     partition = read_clusterfile(clusterfname, seqfos)

#     for seqfo in seqfos:
#         seqfo['name'] = untranslate(seqfo['name'])
#     partition = [[untranslate(uid) for uid in cluster] for cluster in partition]

#     os.remove(cmdfname)
#     os.remove(msafname)
#     if reco_info is not None:
#         os.remove(group_csv_fname)
#     os.rmdir(workdir)
