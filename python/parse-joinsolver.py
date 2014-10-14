#!/usr/bin/env python
import sys
import csv
import utils
from opener import opener

import xml.etree.ElementTree as ET

# ----------------------------------------------------------------------------------------
def figure_out_which_damn_gene(germline_seqs, gene_name, seq, debug=False):
    region = utils.get_region(gene_name)
    seq = seq.replace(' ', '')
    if gene_name in germline_seqs[region]:  # already have it
        return
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
        
    if len(candidates) == 0:
        print 'ERROR didn\'t find jack for', gene_name, seq
        sys.exit()
    # elif len(candidates) > 1:
    #     print 'NOTE found',len(candidates),'candidates, just using the first one'

    if debug:
        print '     swapping', gene_name, '-->', candidates[0]

    return candidates[0]

# ----------------------------------------------------------------------------------------
class JoinParser(object):
    def __init__(self, seqfname, joinfname, datadir):  # <seqfname>: input to joinsolver, <joinfname> output from joinsolver (I only need both because they don't seem to put the full query seq in the output)
        self.debug = 2
        n_max_queries = 10
        queries = ['7940688087618100668']

        self.germline_seqs = utils.read_germlines(datadir)

        # get info that was passed to joinsolver
        self.seqinfo = {}
        with opener('r')(seqfname) as seqfile:
            reader = csv.DictReader(seqfile)
            iline = 0
            for line in reader:
                if line['unique_id'] not in queries:
                    continue
                self.seqinfo[line['unique_id']] = line['seq']
                iline += 1
                if iline >= n_max_queries:
                    break

        tree = ET.parse(joinfname)
        root = tree.getroot()

        self.info = {}
        iline = 0
        for query in root:
            unique_id = query.attrib['id'].replace('>', '').replace(' ', '')
            if unique_id not in queries:
                continue
            # number_header = query.find('vmatches').find('userSeq').find('cdntitle').text  # hopefully don't need to use this
            print unique_id
            self.info[unique_id] = {}
            line = self.info[unique_id]
            line['unique_id'] = unique_id
            line['seq'] = self.seqinfo[unique_id]
            for region in utils.regions:
                if self.debug:
                    print ' ', region
                self.get_region_matches(region, query, line)
            self.add_insertions(line)

            total_match_length = len(line['v_qr_seq'] + line['d_qr_seq'] + line['j_qr_seq'])

            utils.print_reco_event(self.germline_seqs, line)

            iline += 1
            if iline >= n_max_queries:
                break

    # ----------------------------------------------------------------------------------------
    def cut_matches(self, seq1, seq2):
        """ if <seq1> has white space on either end, cut both seqs down to size so it doesn't """
        i_first_match = len(seq1) - len(seq1.lstrip())  # position of first non-space character in seq1
        seq1 = seq1[i_first_match:]
        seq2 = seq2[i_first_match:]
        i_last_match = len(seq1.rstrip())
        seq1 = seq1[:i_last_match]
        seq2 = seq2[:i_last_match]
        return (seq1, seq2)

    # ----------------------------------------------------------------------------------------
    def get_region_matches(self, region, query, line):
        """ get info for <region> and add it to <line> """
        region_query_seq = query.find(region + 'matches').find('userSeq').find('bases').text
        matches = [ match for match in query.find(region + 'matches') if match.tag == 'germline' ]
        if len(matches) == 0:
            print 'ERROR no matches for',line['unique_id']
            return

        # just take the first (best) match
        match = matches[0]
        gl_match_seq = match.text

        if self.debug > 1:
            print '     query', region_query_seq
            print '        gl', gl_match_seq
    
        # if gl match extends outside of query seq, strip off that part
        region_query_seq, gl_match_seq = self.cut_matches(region_query_seq, gl_match_seq)
        if self.debug > 1:
            print '     query', region_query_seq
            print '        gl', gl_match_seq
        # and if region_query_seq extends to left of gl match, strip off that stuff as well
        gl_match_seq, region_query_seq = self.cut_matches(gl_match_seq, region_query_seq)
        if self.debug > 1:
            print '     query', region_query_seq
            print '        gl', gl_match_seq
        # and remove the rest of the spaces
        region_query_seq = region_query_seq.replace(' ', '')
        gl_match_seq = gl_match_seq.replace(' ', '')
        if self.debug > 1:
            print '     query', region_query_seq
            print '        gl', gl_match_seq
    
        # then replace dots in gl_match_seq
        assert len(gl_match_seq) == len(region_query_seq)
        assert region_query_seq in line['seq']  # they're not coming from the same file, so may as well make sure
        new_glseq = []
        for inuke in range(len(region_query_seq)):
            if gl_match_seq[inuke] == '.':
                new_glseq.append(region_query_seq[inuke])
            else:
                new_glseq.append(gl_match_seq[inuke])
        gl_match_seq = ''.join(new_glseq)
        if self.debug > 1:
            print '     query', region_query_seq
            print '        gl', gl_match_seq

        match_name = figure_out_which_damn_gene(self.germline_seqs, match.attrib['id'].replace(' ', ''), gl_match_seq, debug=self.debug)

        line[region + '_gene'] = match_name
        line[region + '_gl_seq'] = gl_match_seq
        line[region + '_qr_seq'] = region_query_seq
        del_5p = self.germline_seqs[region][match_name].find(gl_match_seq)
        assert del_5p >= 0
        del_3p = len(self.germline_seqs[region][match_name]) - del_5p - len(gl_match_seq)
        assert del_3p >= 0
        line[region + '_5p_del'] = del_5p
        line[region + '_3p_del'] = del_3p

    # ----------------------------------------------------------------------------------------
    def add_insertions(self, line):
        for boundary in utils.boundaries:
            left_region = boundary[0]
            right_region = boundary[1]
            left_match_end = line['seq'].find(line[left_region + '_qr_seq']) + len(line[left_region + '_qr_seq'])  # base *after* last base of left region match
            right_match_start = line['seq'].find(line[right_region + '_qr_seq'])  # first base of right region match
            line[boundary + '_insertion'] = line['seq'][left_match_end : right_match_start]

jparser = JoinParser('caches/recombinator/simu.csv', '/home/dralph/multijoin.xml', datadir='./data')
