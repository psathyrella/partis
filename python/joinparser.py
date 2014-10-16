#!/usr/bin/env python
import sys
import csv
import os
import math
import xml.etree.ElementTree as ET

import utils
import plotting
from opener import opener
from performanceplotter import PerformancePlotter


# ----------------------------------------------------------------------------------------
def add_insertions(line):
    for boundary in utils.boundaries:
        left_region = boundary[0]
        right_region = boundary[1]
        left_match_end = line['seq'].find(line[left_region + '_qr_seq']) + len(line[left_region + '_qr_seq'])  # base *after* last base of left region match
        right_match_start = line['seq'].find(line[right_region + '_qr_seq'])  # first base of right region match
        line[boundary + '_insertion'] = line['seq'][left_match_end : right_match_start]


# ----------------------------------------------------------------------------------------
def cut_matches(seq1, seq2):
    """ if <seq1> has white space on either end, cut both seqs down to size so it doesn't """
    i_first_match = len(seq1) - len(seq1.lstrip())  # position of first non-space character in seq1
    seq1 = seq1[i_first_match:]
    seq2 = seq2[i_first_match:]
    i_last_match = len(seq1.rstrip())
    seq1 = seq1[:i_last_match]
    seq2 = seq2[:i_last_match]
    return (seq1, seq2)

# ----------------------------------------------------------------------------------------
def resolve_overlapping_matches(line, debug=False):
    """
    joinsolver allows d and j matches (and v and d matches) to overlap... which makes no sense, so
    arbitrarily split the disputed territory in two.
    """
    # NOTE this is basically a cut and paste job from waterer.py
    for rpairs in ({'left':'v', 'right':'d'}, {'left':'d', 'right':'j'}):
        left_gene = line[rpairs['left'] + '_gene']
        right_gene = line[rpairs['right'] + '_gene']
        left_match_end = line['seq'].find(line[rpairs['left'] + '_qr_seq']) + len(line[rpairs['left'] + '_qr_seq'])  # base after the last left-gene-matched base
        right_match_start = line['seq'].find(line[rpairs['right'] + '_qr_seq'])  # first base of right-hand match
        overlap = left_match_end - right_match_start
        if overlap > 0:
            if debug:
                print '     WARNING %s apportioning %d overlapping bases between %s and %s matches' % (line['unique_id'], overlap, rpairs['left'], rpairs['right'])
            lefthand_portion = int(math.floor(overlap / 2.0))
            righthand_portion = int(math.ceil(overlap / 2.0))
            assert lefthand_portion <= len(line[rpairs['left'] + '_gl_seq'])
            assert righthand_portion <= len(line[rpairs['right'] + '_gl_seq'])
            line[rpairs['left'] + '_gl_seq'] = line[rpairs['left'] + '_gl_seq'][:-lefthand_portion]
            line[rpairs['left'] + '_qr_seq'] = line[rpairs['left'] + '_qr_seq'][:-lefthand_portion]
            line[rpairs['left'] + '_3p_del'] += lefthand_portion

            line[rpairs['right'] + '_gl_seq'] = line[rpairs['right'] + '_gl_seq'][righthand_portion:]
            line[rpairs['right'] + '_qr_seq'] = line[rpairs['right'] + '_qr_seq'][righthand_portion:]
            line[rpairs['right'] + '_5p_del'] += righthand_portion

#--------------------------------------------------------------------------------------
def figure_out_which_damn_gene(germline_seqs, gene_name, seq, debug=False):
    region = utils.get_region(gene_name)
    seq = seq.replace(' ', '')
    if gene_name in germline_seqs[region]:  # already have it
        if len(seq) > len(germline_seqs[region][gene_name]):
            print '      gl match longer than gl!'
            print '       ', seq
            print '       ', germline_seqs[region][gene_name]
            germline_seqs[region][gene_name] = seq
        return gene_name
    candidates = []

    # if it doesn't specify an allele, see if any of the alleles we've got have the same sequence in the match region
    if gene_name.find('*') == -1:
        for candidate_gene in germline_seqs[region]:
            if candidate_gene.find(gene_name) == 0:
                if seq in germline_seqs[region][candidate_gene]:
                     candidates.append(candidate_gene)

    # if it *does* specify an allele, see if any of the other allele have the same sequence in the match region
    if len(candidates) == 0:  # didn't find anything... try other alleles
        for candidate_gene in germline_seqs[region]:
            if utils.are_alleles(candidate_gene, gene_name):
                if seq in germline_seqs[region][candidate_gene]:
                    candidates.append(candidate_gene)

    # sometimes it's 3-9, but sometimes 3-09. *grrrrrr*.
    if len(candidates) == 0:
        for candidate_gene in germline_seqs[region]:
            if gene_name.replace('-0', '-') == candidate_gene:
                if seq in germline_seqs[region][candidate_gene]:
                    candidates.append(candidate_gene)

    # try adding _F and _P to the end of j names
    if len(candidates) == 0:
        for candidate_gene in germline_seqs[region]:
            if gene_name + '_F' == candidate_gene or gene_name + '_P' == candidate_gene:
                if seq in germline_seqs[region][candidate_gene]:
                    candidates.append(candidate_gene)

    # try removing the darn R at the end (and remove the zero). I hope it doesn't mean anything important
    if len(candidates) == 0:
        for candidate_gene in germline_seqs[region]:
            if gene_name.replace('R', '').replace('-0', '-') == candidate_gene:
                if seq in germline_seqs[region][candidate_gene]:
                    candidates.append(candidate_gene)
        
    if len(candidates) == 0:
        print '    ERROR didn\'t find jack for', gene_name, seq
        assert False
    # elif len(candidates) > 1:
    #     print 'NOTE found',len(candidates),'candidates, just using the first one'

    if debug:
        print '     swapping', gene_name, '-->', candidates[0]

    return candidates[0]

# ----------------------------------------------------------------------------------------
class JoinParser(object):
    def __init__(self, seqfname, joinfname, datadir):  # <seqfname>: input to joinsolver, <joinfname> output from joinsolver (I only need both because they don't seem to put the full query seq in the output)
        self.debug = 0
        n_max_queries = -1
        queries = []

        self.germline_seqs = utils.read_germlines(datadir, remove_N_nukes=False)
        perfplotter = PerformancePlotter(self.germline_seqs, os.getenv('www') + '/partis/joinsolver_performance', 'js')

        # get info that was passed to joinsolver
        self.seqinfo = {}
        with opener('r')(seqfname) as seqfile:
            reader = csv.DictReader(seqfile)
            iline = 0
            for line in reader:
                if len(queries) > 0 and line['unique_id'] not in queries:
                    continue
                self.seqinfo[line['unique_id']] = line
                iline += 1
                if n_max_queries > 0 and iline >= n_max_queries:
                    break

        tree = ET.parse(joinfname)
        root = tree.getroot()

        # self.info = {}
        iline = 0
        for query in root:
            iline += 1
            if n_max_queries > 0 and iline > n_max_queries:
                break

            unique_id = query.attrib['id'].replace('>', '').replace(' ', '')
            print iline, unique_id
            if len(queries) > 0 and  unique_id not in queries:
                continue
            # self.info[unique_id] = {}
            line = {}  # self.info[unique_id] TODO make sure that still works with this stuff commented out
            line['unique_id'] = unique_id
            line['seq'] = self.seqinfo[unique_id]['seq']
            for region in utils.regions:
                if self.debug:
                    print ' ', region
                self.get_region_matches(region, query, line)
            if 'v_gene' not in line or 'd_gene' not in line or 'j_gene' not in line:
                print '  ERROR giving up on %s' % unique_id
                continue

            add_insertions(line)
            resolve_overlapping_matches(line, self.debug)

            perfplotter.evaluate(self.seqinfo[unique_id], line, unique_id)

            if self.debug:
                utils.print_reco_event(self.germline_seqs, line)

        perfplotter.plot()

    # ----------------------------------------------------------------------------------------
    def parse_match_seqs(self, match, region_query_seq):
        gl_match_seq = match.text
        if gl_match_seq == None:
            return (None, None)
        if self.debug > 1:
            print '     query', region_query_seq
            print '        gl', gl_match_seq

        # if gl match extends outside of query seq, strip off that part
        region_query_seq, gl_match_seq = cut_matches(region_query_seq, gl_match_seq)
        # and if region_query_seq extends to left of gl match, strip off that stuff as well
        gl_match_seq, region_query_seq = cut_matches(gl_match_seq, region_query_seq)
        # and remove the rest of the spaces
        region_query_seq = region_query_seq.replace(' ', '')
        gl_match_seq = gl_match_seq.replace(' ', '')
    
        # then replace dots in gl_match_seq, and just remove dashes  TODO do something better with dashes
        assert len(gl_match_seq) == len(region_query_seq)
        new_glseq = []
        for inuke in range(len(region_query_seq)):
            if gl_match_seq[inuke] == '.':
                new_glseq.append(region_query_seq[inuke])
            elif gl_match_seq[inuke] == '-':
                pass
            else:
                assert gl_match_seq[inuke] in utils.nukes
                new_glseq.append(gl_match_seq[inuke])
        gl_match_seq = ''.join(new_glseq)
        
        if '-' in region_query_seq:
            print '    WARNING removing gaps in query seq'
            region_query_seq = region_query_seq.replace('-', '')

        if self.debug > 1:
            print '     query', region_query_seq
            print '        gl', gl_match_seq

        return (region_query_seq, gl_match_seq)

    # ----------------------------------------------------------------------------------------
    def get_region_matches(self, region, query, line):
        """ get info for <region> and add it to <line> """
        if query.find(region + 'matches').find('nomatches') != None:  # no matches!
            print '  ERROR no %s matches' % region
            return

        region_query_seq = query.find(region + 'matches').find('userSeq').find('bases').text
        matches = [ match for match in query.find(region + 'matches') if match.tag == 'germline' ]
        if len(matches) == 0:
            print 'ERROR no matches for',line['unique_id']
            return

        imatch = 0
        match = matches[imatch]
        # just take the first (best) match
        region_query_seq, gl_match_seq = self.parse_match_seqs(match, region_query_seq)
        if region_query_seq == None or gl_match_seq == None:
            return

        match_name = match.attrib['id'].replace(' ', '')
        try:
            match_name = figure_out_which_damn_gene(self.germline_seqs, match_name, gl_match_seq, debug=self.debug)
        except AssertionError:  # couldn't find a decent one, so try again with the second match
            # well, ok, I guess I'll just *add* the damn thing to <germline_seqs>. TODO that is so, so, dirty
            self.germline_seqs[region][match.attrib['id'].replace(' ', '')] = gl_match_seq
            print '   WARNING adding %s to <germline_seqs>' % match.attrib['id'].replace(' ', '')
            match_name = figure_out_which_damn_gene(self.germline_seqs, match.attrib['id'].replace(' ', ''), gl_match_seq, debug=self.debug)

        if self.debug > 1:
            print '     ', match_name
        if region_query_seq not in line['seq']:
            print region_query_seq
            print line['seq']
        assert region_query_seq in line['seq']  # they're not coming from the same file, so may as well make sure
            
        del_5p = self.germline_seqs[region][match_name].find(gl_match_seq)
        del_3p = len(self.germline_seqs[region][match_name]) - del_5p - len(gl_match_seq)
        if del_5p < 0 or del_3p < 0:
            print 'ERROR couldn\'t figure out deletions in', region
            print '    germline', self.germline_seqs[region][match_name], match_name
            print '    germline match', gl_match_seq
            print '   ', del_5p, del_3p
            return

        line[region + '_5p_del'] = del_5p
        line[region + '_3p_del'] = del_3p
        line[region + '_gene'] = match_name
        line[region + '_gl_seq'] = gl_match_seq
        line[region + '_qr_seq'] = region_query_seq


plotting.compare_directories('/var/www/sharing/dralph/partis/performance/plots/', 'hmm', '/var/www/sharing/dralph/partis/joinsolver_performance/plots/', 'jsolver', xtitle='inferred - true', stats='')
# jparser = JoinParser('caches/recombinator/simu.csv', '/home/dralph/Dropbox/multijoin.xml', datadir='./data')
