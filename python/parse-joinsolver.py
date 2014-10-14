#!/usr/bin/env python
import sys

import utils

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
                if germline_seqs[region][candidate_gene].find(seq):
                     candidates.append(candidate_gene)
    # if it *does* specify an allele, see if any of the other allele have the same sequence in the match region
    if len(candidates) == 0:  # didn't find anything... try other alleles
        for candidate_gene in germline_seqs[region]:
            if utils.are_alleles(candidate_gene, gene_name):
                if germline_seqs[region][candidate_gene].find(seq):
                    candidates.append(candidate_gene)

    # sometimes it's 3-9, but sometimes 3-09. *grrrrrr*.
    if len(candidates) == 0:
        for candidate_gene in germline_seqs[region]:
            if gene_name.replace('-0', '-') == candidate_gene:
                if germline_seqs[region][candidate_gene].find(seq):
                    candidates.append(candidate_gene)

    # try adding _F and _P to the end of j names
    if len(candidates) == 0:
        for candidate_gene in germline_seqs[region]:
            if gene_name + '_F' == candidate_gene or gene_name + '_P' == candidate_gene:
                if germline_seqs[region][candidate_gene].find(seq):
                    candidates.append(candidate_gene)
        
    if len(candidates) == 0:
        print 'ERROR didn\'t find jack for',gene_name
        sys.exit()
    # elif len(candidates) > 1:
    #     print 'NOTE found',len(candidates),'candidates, just using the first one'

    if debug:
        print '  swapping', gene_name, '-->', candidates[0]

    return candidates[0]

# ----------------------------------------------------------------------------------------
class JoinParser(object):
    def __init__(self, infname, datadir):
        self.germline_seqs = utils.read_germlines(datadir)
        tree = ET.parse(infname)
        root = tree.getroot()

        self.info = {}
        for query in root:
            unique_id = query.attrib['id'].replace('> ', '')
            # number_header = query.find('vmatches').find('userSeq').find('cdntitle').text  # hopefully don't need to use this
            print unique_id
            self.info[unique_id] = {}
            self.info[unique_id]['unique_id'] = unique_id
            for region in utils.regions:
                self.get_region_matches(region, query, self.info[unique_id])
                print self.info[unique_id]

        # line['v_5p_del'] = self.germline_seqs['v'][match_names[0]].find(gl_match_seqs[0])
        # assert line['v_5p_del'] != -1

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
        query_seq = query.find(region + 'matches').find('userSeq').find('bases').text
        matches = [ match for match in query.find(region + 'matches') if match.tag == 'germline' ]
        if len(matches) == 0:
            print 'ERROR no matches for',line['unique_id']
            return

        # just take the first (best) match
        match = matches[0]
        match_name = figure_out_which_damn_gene(self.germline_seqs, match.attrib['id'].replace(' ', ''), match.text, debug=True)
        gl_match_seq = match.text
    
        # if gl match extends outside of query seq, strip off that part
        query_seq, gl_match_seq = self.cut_matches(query_seq, gl_match_seq)
        # and if query_seq extends to left of gl match, strip off that stuff as well
        gl_match_seq, query_seq = self.cut_matches(gl_match_seq, query_seq)
        # and remove the rest of the spaces
        query_seq = query_seq.replace(' ', '')
        gl_match_seq = gl_match_seq.replace(' ', '')
    
        # then replace dots in gl_match_seq
        assert len(gl_match_seq) == len(query_seq)
        new_glseq = []
        for inuke in range(len(query_seq)):
            if gl_match_seq[inuke] == '.':
                new_glseq.append(query_seq[inuke])
            else:
                new_glseq.append(gl_match_seq[inuke])
        gl_match_seq = ''.join(new_glseq)

        line[region + '_gene'] = match_name
        line[region + '_gl_seq'] = gl_match_seq
        line[region + '_qr_seq'] = query_seq

jparser = JoinParser('/home/dralph/multijoin.xml', datadir='./data')
