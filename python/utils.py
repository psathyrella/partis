""" A few utility functions. At the moment simply functions used in recombinator which do not
require member variables. """

import sys
import os
import re
import math
import glob
from collections import OrderedDict
import csv

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
def get_arg_list(arg, intify=False, floatify=False):  # make lists from args that are passed as strings of colon-separated values
    if arg == None:
        return arg
    else:
        arglist = arg.strip().split(':')  # to allow ids with minus signs, need to add a space, which you then have to strip() off
        if intify:
            return [int(x) for x in arglist]
        elif floatify:
            return [float(x) for x in arglist]
        else:
            return arglist

# # ----------------------------------------------------------------------------------------
# hackey_default_gene_versions = {'v':'IGHV3-23*04', 'd':'IGHD3-10*01', 'j':'IGHJ4*02_F'}
# ----------------------------------------------------------------------------------------
regions = ['v', 'd', 'j']
real_erosions = ['v_3p', 'd_5p', 'd_3p', 'j_5p']
effective_erosions = ['v_5p', 'j_3p']
boundaries = ['vd', 'dj']
humans = ['A', 'B', 'C']
nukes = ['A', 'C', 'G', 'T']
ambiguous_bases = ['N', ]
naivities = ['M', 'N']
conserved_codon_names = {'v':'cyst', 'd':'', 'j':'tryp'}
# Infrastrucure to allow hashing all the columns together into a dict key.
# Uses a tuple with the variables that are used to index selection frequencies
# NOTE fv and jf insertions are *effective* (not real) insertions between v or j and the framework. The allow query sequences that extend beyond the v or j regions
index_columns = ('v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'v_5p_del', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del', 'j_3p_del', 'fv_insertion', 'vd_insertion', 'dj_insertion', 'jf_insertion')
# not_used_for_simulation = ('fv_insertion', 'jf_insertion', 'v_5p_del')
index_keys = {}
for i in range(len(index_columns)):  # dict so we can access them by name instead of by index number
    index_keys[index_columns[i]] = i

# ----------------------------------------------------------------------------------------
# Info specifying which parameters are assumed to correlate with which others. Taken from mutual
# information plot in bcellap repo

# key is parameter of interest, and associated list gives the parameters (other than itself) which are necessary to predict it
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

def get_parameter_fname(column=None, deps=None, column_and_deps=None):
    """ return the file name in which we store the information for <column>. Either pass in <column> and <deps> *or* <column_and_deps> """
    if column == 'all':
        return 'all-probs.csv'
    if column_and_deps == None:
        column_and_deps = [column]
        column_and_deps.extend(deps)
    outfname = 'probs.csv'
    for ic in column_and_deps:
        outfname = ic + '-' + outfname
    return outfname

# ----------------------------------------------------------------------------------------
def read_cyst_positions(datadir):
    import json
    jsonfname = datadir + '/v-meta.json'
    csvfname = datadir + '/v-meta.csv'
    with opener('r')(jsonfname) as json_file:
        cyst_positions = json.load(json_file)
    if not os.path.exists(csvfname) or os.path.getmtime(csvfname) < os.path.getmtime(jsonfname):  # if csv hasn't been made, or if it's older than the json file
        print 'rewriting v-meta csv'
        # remake csv file
        with opener('w')(csvfname) as csv_file:
            writer = csv.DictWriter(csv_file, ('gene', 'cyst_start'))
            writer.writeheader()
            for gene in cyst_positions:
                writer.writerow({'gene' : gene, 'cyst_start' : cyst_positions[gene]['cysteine-position']})

    return cyst_positions

# ----------------------------------------------------------------------------------------
def from_same_event(is_data, reco_info, query_names):
    if is_data:
        return None
    if len(query_names) > 1:
        reco_id = reco_info[query_names[0]]['reco_id']  # the first one's reco id
        for iq in range(1, len(query_names)):  # then loop through the rest of 'em to see if they're all the same
            if reco_id != reco_info[query_names[iq]]['reco_id']:
                return False
        return True
    else:
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
def color_gene(gene):
    """ color gene name (and remove extra characters), eg IGHV3-h*01 --> v3-h1 """
    return_str = gene[:3] + color('bold', color('red', gene[3]))
    n_version = gene[4 : gene.find('-')]
    n_subversion = gene[gene.find('-')+1 : gene.find('*')]
    if get_region(gene) == 'j':
        n_version = gene[4 : gene.find('*')]
        n_subversion = ''
        return_str += color('purple', n_version)
    else:
        return_str += color('purple', n_version) + '-' + color('purple', n_subversion)

    if gene.find('*') != -1:
        allele_end = gene.find('_')
        if allele_end < 0:
            allele_end = len(gene)
        allele = gene[gene.find('*')+1 : allele_end]
        return_str += '*' + color('yellow', allele)
        if '_' in gene:  # _F or _P in j gene names
            return_str += gene[gene.find('_') :]

    # now remove extra characters
    return_str = return_str.replace('IGH','').lower()
    return_str = return_str.replace('*','')
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
def is_mutated(original, final, n_muted=-1, n_total=-1):
    alphabet = nukes + ambiguous_bases
    if original not in alphabet or final not in alphabet:
        raise Exception('bad base (%s or %s) in utils.is_mutated()' % (original, final))

    return_str = final
    if original in ambiguous_bases or final in ambiguous_bases:  # don't count Ns in the total
        return return_str, n_muted, n_total
        
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
def get_conserved_codon_position(cyst_positions, tryp_positions, region, gene, glbounds, qrbounds, assert_on_fail=True):
    """
    Find location of the conserved cysteine/tryptophan in a query sequence given a germline match which is specified by
    its germline bounds <glbounds> and its bounds in the query sequence <qrbounds>
    """
    # NOTE see add_cdr3_info -- they do similar things, but start from different information
    if region == 'v':
        gl_pos = cyst_positions[gene]['cysteine-position']  # germline cysteine position
    elif region == 'j':
        gl_pos = int(tryp_positions[gene])
    else:  # return -1 for d
        return -1

    if gl_pos == None:
        raise Exception('none type gl_pos for %s ' % gene)

    if assert_on_fail:
        assert glbounds[0] <= gl_pos  # make sure we didn't erode off the conserved codon
    query_pos = gl_pos - glbounds[0] + qrbounds[0]  # position within original germline gene, minus the position in that germline gene at which the match starts, plus the position in the query sequence at which the match starts
    # if region == 'v':
    #     print '%s  %d = %d - %d + %d' % (color_gene(gene), query_pos, gl_pos, glbounds[0], qrbounds[0])
    return query_pos

# ----------------------------------------------------------------------------------------
def add_cdr3_info(germlines, cyst_positions, tryp_positions, line, debug=False):
    """
    Add the cyst_position, tryp_position, and cdr3_length to <line> based on the information already in <line>.
    If info is already there, make sure it's the same as what we calculate here
    """

    original_seqs = {}  # original (non-eroded) germline seqs
    lengths = {}  # length of each match (including erosion)
    eroded_seqs = {}  # eroded germline seqs
    get_reco_event_seqs(germlines, line, original_seqs, lengths, eroded_seqs)

    # NOTE see get_conserved_codon_position -- they do similar things, but start from different information
    eroded_gl_cpos = cyst_positions[line['v_gene']]['cysteine-position']  - int(line['v_5p_del']) + len(line['fv_insertion'])  # cysteine position in eroded germline sequence. EDIT darn, actually you *don't* want to subtract off the v left deletion, because that (deleted) base is presumably still present in the query sequence
    # if debug:
    #     print '  cysteine: cpos - v_5p_del + fv_insertion = %d - %d + %d = %d' % (cyst_positions[line['v_gene']]['cysteine-position'], int(line['v_5p_del']), len(line['fv_insertion']), eroded_gl_cpos)
    eroded_gl_tpos = int(tryp_positions[line['j_gene']]) - int(line['j_5p_del'])
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

# ----------------------------------------------------------------------------------------
def get_full_naive_seq(germlines, line):  #, restrict_to_region=''):
    for erosion in real_erosions + effective_erosions:
        if line[erosion + '_del'] < 0:
            print 'ERROR %s less than zero %d' % (erosion, line[erosion + '_del'])
        assert line[erosion + '_del'] >= 0
    original_seqs = {}  # original (non-eroded) germline seqs
    lengths = {}  # length of each match (including erosion)
    eroded_seqs = {}  # eroded germline seqs
    get_reco_event_seqs(germlines, line, original_seqs, lengths, eroded_seqs)
    # if restrict_to_region == '':
    return line['fv_insertion'] + eroded_seqs['v'] + line['vd_insertion'] + eroded_seqs['d'] + line['dj_insertion'] + eroded_seqs['j'] + line['jf_insertion']
    # else:
    # return eroded_seqs[restrict_to_region]

# ----------------------------------------------------------------------------------------
def get_regional_naive_seq_bounds(return_reg, germlines, line, subtract_unphysical_erosions=True):
    # NOTE it's kind of a matter of taste whether unphysical deletions (v left and j right) should be included in the 'naive sequence'.
    # Unless <subtract_unphysical_erosions>, here we assume the naive sequence has *no* unphysical deletions

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

    # for key, val in line.items():
    #     print key, val
    for chkreg in regions:
        # print chkreg, start[chkreg], end[chkreg]
        assert start[chkreg] >= 0
        assert end[chkreg] >= 0
        assert end[chkreg] >= start[chkreg]
    # print end['j'], len(line['seq']), line['v_5p_del'], line['j_3p_del']
    assert end['j'] == len(line['seq'])

    return (start[return_reg], end[return_reg])

# ----------------------------------------------------------------------------------------
def add_match_info(germlines, line, cyst_positions, tryp_positions, debug=False):
    """
    add to <line> the query match seqs (sections of the query sequence that are matched to germline) and their corresponding germline matches.

    """

    original_seqs = {}  # original (non-eroded) germline seqs
    lengths = {}  # length of each match (including erosion)
    eroded_seqs = {}  # eroded germline seqs
    get_reco_event_seqs(germlines, line, original_seqs, lengths, eroded_seqs)
    add_cdr3_info(germlines, cyst_positions, tryp_positions, line, debug=debug)  # add cyst and tryp positions, and cdr3 length

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
def print_reco_event(germlines, line, one_line=False, extra_str='', return_string=False, label=''):
    """ Print ascii summary of recombination event and mutation.

    If <one_line>, then only print out the final_seq line.
    """
    
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

    germline_v_end = len(line['fv_insertion']) + len(original_seqs['v']) - v_5p_del - 1  # position in the query sequence at which we find the last base of the v match.
                                                                                         # NOTE we subtract off the v_5p_del because we're *not* adding dots for that deletion (it's just too long)
    germline_d_start = len(line['fv_insertion']) + lengths['v'] + len(line['vd_insertion']) - d_5p_del
    germline_d_end = germline_d_start + len(original_seqs['d'])
    germline_j_start = germline_d_end + 1 - d_3p_del + len(line['dj_insertion']) - j_5p_del

    # build up the query sequence line, including colors for mutations and conserved codons
    final_seq = ''
    n_muted, n_total = 0, 0
    j_right_extra = ''  # portion of query sequence to right of end of the j match
    for inuke in range(len(line['seq'])):
        ilocal = inuke
        new_nuke = ''
        if ilocal < len(line['fv_insertion']):  # haven't got to start of v match yet, so just add on the query seq nuke
            new_nuke = line['seq'][inuke]
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
                                new_nuke, n_muted, n_total = line['seq'][inuke], n_muted, n_total+1
                                j_right_extra += ' '

        if 'cyst_position' in line and 'tryp_position' in line:
            for pos in (line['cyst_position'], line['tryp_position']):  # reverse video for the conserved codon positions
                # adjusted_pos = pos - v_5p_del  # adjust positions to allow for reads not extending all the way to left side of v
                adjusted_pos = pos  # don't need to subtract it for smith-waterman stuff... gr, need to make it general
                if inuke >= adjusted_pos and inuke < adjusted_pos+3:
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
    insert_line = color_chars(ambiguous_bases, 'light_blue', insert_line)

    d_line = ' ' * germline_d_start
    d_line += ' '*len(v_5p_del_str)
    d_line += eroded_seqs_dots['d']
    d_line += ' ' * (len(original_seqs['j']) - j_5p_del - j_3p_del + len(line['dj_insertion']) - d_3p_del)
    d_line += j_right_extra
    d_line += ' ' * j_3p_del
    # d_line += ' '*len(line['jf_insertion'])
    d_line = color_chars(ambiguous_bases, 'light_blue', d_line)

    vj_line = ' ' * len(line['fv_insertion'])
    vj_line += v_5p_del_str
    vj_line += eroded_seqs_dots['v'] + '.'*extra_space_because_of_fixed_nospace
    vj_line += ' ' * (germline_j_start - germline_v_end - 2)
    vj_line += eroded_seqs_dots['j']
    vj_line += j_right_extra
    # vj_line += ' '*len(line['jf_insertion'])
    vj_line = color_chars(ambiguous_bases, 'light_blue', vj_line)

    if len(insert_line) != len(d_line) or len(insert_line) != len(vj_line):
        # print '\nERROR lines unequal lengths in event printer -- insertions %d d %d vj %d' % (len(insert_line), len(d_line), len(vj_line)),
        # assert no_space
        if not no_space:
            print 'ERROR no space'
        # print ' ...but we\'re out of space so it\'s expected'

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
        if line['chops']['left'] > 0:
            v_5p_del_space_str = v_5p_del_space_str[ : -line['chops']['left']]
            final_seq = color('green', '.'*line['chops']['left']) + final_seq
        if line['chops']['right'] > 0:
            j_3p_del_space_str = j_3p_del_space_str[line['chops']['right'] : ]
            final_seq = final_seq + color('green', '.'*line['chops']['right'])
    final_seq = v_5p_del_space_str + final_seq + j_3p_del_space_str
    final_seq = color_chars(ambiguous_bases, 'light_blue', final_seq)

    out_str_list.append('%s    %s' % (extra_str, final_seq))
    # and finally some extra info
    out_str_list.append('   muted: %4.2f' % (float(n_muted) / n_total))
    if 'score' in line:
        out_str_list.append('  score: %s' % line['score'])
    if 'cdr3_length' in line:
        out_str_list.append('   cdr3: %d' % int(line['cdr3_length']))
    out_str_list.append('\n')

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
def read_germlines(data_dir):
    """ <remove_fp> sometimes j names have a redundant _F or _P appended to their name. Set to True to remove this """
    germlines = {}
    for region in regions:
        germlines[region] = OrderedDict()
        for seq_record in SeqIO.parse(data_dir + '/igh'+region+'.fasta', "fasta"):
            gene_name = seq_record.name
            seq_str = str(seq_record.seq)
            germlines[region][gene_name] = seq_str
    return germlines

# ----------------------------------------------------------------------------------------
def get_region(gene_name):
    """ return v, d, or j of gene"""
    try:
        assert 'IGH' in gene_name
        region = gene_name[3:4].lower()
        assert region in regions
        return region
    except:
        print 'ERROR faulty gene name %s ' % gene_name
        assert False
# ----------------------------------------------------------------------------------------
def maturity_to_naivety(maturity):
    if maturity == 'memory':
        return 'M'
    elif maturity == 'naive':
        return 'N'
    else:
        assert False

# # ----------------------------------------------------------------------------------------
# def split_gene_name(name):
#     """
#     split name into region, version, allele, etc.
#     e.g. IGHD7-27*01 --> {'region':'d', 'version':7, 'subversion':27, 'allele':1}
#     """
#     return_vals = {}
#     return_vals['region'] = get_region(name)
#     assert name.count('-') == 1
#     return_vals['version'] = name[4 : name.find('-')]
    
#     assert name.count('*') == 1
    

# ----------------------------------------------------------------------------------------
def are_alleles(gene1, gene2):
    """
    Return true if gene1 and gene2 are alleles of the same gene version.
    Assumes they're alleles if everything left of the asterisk is the same, and everything more than two to the right of the asterisk is the same.
    """
    # gene1 = apply_renaming_scheme(gene1)
    # gene2 = apply_renaming_scheme(gene2)

    left_str_1 = gene1[0 : gene1.find('*')]
    left_str_2 = gene2[0 : gene1.find('*')]
    right_str_1 = gene1[gene1.find('*')+3 :]
    right_str_2 = gene2[gene1.find('*')+3 :]
    return left_str_1 == left_str_2 and right_str_1 == right_str_2

# ----------------------------------------------------------------------------------------
def are_same_primary_version(gene1, gene2):
    """
    Return true if the bit up to the dash is the same.
    There's probably a real name for that bit.
    """
    str_1 = gene1[0 : gene1.find('-')]
    str_2 = gene2[0 : gene2.find('-')]
    return str_1 == str_2

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
def get_hamming_distances(pairs):
    return_info = []
    for info in pairs:
        seq_a = info['seq_a']
        seq_b = info['seq_b']
        if True:  #self.args.truncate_pairs:  # chop off the left side of the longer one if they're not the same length
            min_length = min(len(seq_a), len(seq_b))
            seq_a = seq_a[-min_length : ]
            seq_b = seq_b[-min_length : ]
            chopped_off_left_sides = True
        mutation_frac = hamming_fraction(seq_a, seq_b)
        return_info.append({'id_a':info['id_a'], 'id_b':info['id_b'], 'score':mutation_frac})

    return return_info

# ----------------------------------------------------------------------------------------
def hamming_fraction(seq1, seq2, return_len_excluding_ambig=False):
    assert len(seq1) == len(seq2)
    if len(seq1) == 0:
        return 0.

    alphabet = nukes + ambiguous_bases
    distance, len_excluding_ambig = 0, 0
    for ch1, ch2 in zip(seq1, seq2):
        if alphabet is not None:  # check that both characters are in the expected alphabet
            if ch1 not in alphabet or ch2 not in alphabet:
                raise Exception('unexpected character (%s or %s) not among %s in hamming_fraction()' % (ch1, ch2, alphabet))
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
        return  fraction

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
def process_input_line(info, splitargs=(), int_columns=(), float_columns=()):
    """ 
    Attempt to convert all the keys and values in <info> from str to int.
    The keys listed in <splitargs> will be split as colon-separated lists before intification.
    """
    if len(splitargs) == 0 and len(int_columns) == 0 and len(float_columns) == 0:  # nothing to do
        return

    for key, val in info.items():
        if key in splitargs:
            info[key] = info[key].split(':')
            for i in range(len(info[key])):
                if key in int_columns:
                    info[key][i] = int(info[key][i])
                elif key in float_columns:
                    info[key][i] = float(info[key][i])
        else:
            if key in int_columns:
                info[key] = int(info[key])
            elif key in float_columns:
                info[key] = float(info[key])

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
            bounds = get_regional_naive_seq_bounds(region, germlines, line, subtract_unphysical_erosions=True)
            mashed_naive_seq += naive_seq[bounds[0] : bounds[1]]
            mashed_muted_seq += muted_seq[bounds[0] : bounds[1]]
    else:
        bounds = get_regional_naive_seq_bounds(restrict_to_region, germlines, line, subtract_unphysical_erosions=True)
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
def process_out_err(out, err, extra_str):
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

    print_str += out

    if print_str != '':
        print ' --> proc %s' % extra_str
        print print_str

