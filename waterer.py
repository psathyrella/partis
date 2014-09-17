import time
import sys
import csv
import os
import itertools
import operator
import pysam
import contextlib
from subprocess import check_call

from utils import utils
from utils.opener import opener
from parametercounter import ParameterCounter

# ----------------------------------------------------------------------------------------
class Waterer(object):
    """ Run smith-waterman on the query sequences in <infname> """
    def __init__(self, pdriver, args, bootstrap):
        self.pdriver = pdriver
        self.args = args
        self.info = {}
        self.pcounter = ParameterCounter(self.pdriver, 'data/human-beings/' + self.args.human + '/' + self.args.naivety)
        self.bootstrap = bootstrap  # bootsrap=True if we *don't* have any parameters to start with, i.e. we don't yet know anything about this data. The main effect of bootsrap=True is that gene_choice_probs will *not* be applied (since we of course don't know them)
        if not self.bootstrap:
            self.gene_choice_probs = utils.read_overall_gene_prob(self.pdriver.datadir + '/human-beings/' + self.args.human + '/' + self.args.naivety)
        outfname = self.pdriver.workdir + '/query-seqs.bam'
        self.run_smith_waterman(outfname)
        self.read_output(outfname)

    # ----------------------------------------------------------------------------------------
    def run_smith_waterman(self, outfname):
        """
        Run smith-waterman alignment on the seqs in <infname>, and toss all the top matches into <outfname>.
        Then run through <outfname> to get the top hits and their locations to pass to the hmm.
        Then run the hmm on each gene set.
        """
        infname = self.pdriver.workdir + '/query-seqs.fa'
        with opener('w')(infname) as swinfile:  # write *all* the input seqs to file, i.e. run s-w on all of 'em at once
            for query_name,line in self.pdriver.input_info.iteritems():
                swinfile.write('>' + query_name + ' NUKES\n')
                swinfile.write(line['seq'] + '\n')
        start = time.time()
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        check_call('/home/dralph/.local/bin/vdjalign align-fastq --j-subset adaptive --max-drop 50 --match 5 --mismatch 3 --gap-open 100 ' + infname + ' ' + outfname + ' 2>/dev/null', shell=True)
        os.remove(infname)
        print '    s-w time: %.3f' % (time.time()-start)
    
    # ----------------------------------------------------------------------------------------
    def read_output(self, outfname):
        start = time.time()
        with contextlib.closing(pysam.Samfile(outfname)) as bam:
            grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
            for _, reads in grouped:  # loop over query sequences
                self.process_query(bam, list(reads))
        print '    sw read time: %.3f' % (time.time() - start)
        os.remove(outfname)

    # ----------------------------------------------------------------------------------------
    def get_choice_prob(self, region, gene):
        choice_prob = 1.0
        if not self.bootstrap:
            if gene in self.gene_choice_probs[region]:
                choice_prob = self.gene_choice_probs[region][gene]
            else:
                choice_prob = 0.0  # TODO choose something else?
        return choice_prob

    # ----------------------------------------------------------------------------------------
    def process_query(self, bam, reads):
        primary = next((r for r in reads if not r.is_secondary), None)
        query_seq = primary.seq
        query_name = primary.qname
        raw_best = {}
        all_match_names = {}
        for region in utils.regions:
            all_match_names[region] = []
        all_query_bounds, all_germline_bounds = {}, {}
        for read in reads:  # loop over the matches found for each query sequence
            read.seq = query_seq  # only the first one has read.seq set by default, so we need to set the rest by hand
            gene = bam.references[read.tid]
            region = utils.get_region(gene)

            if region not in raw_best:  # best v, d, and j before multiplying by gene choice probs. needed 'cause *these* are the v and j that get excised
                raw_best[region] = gene

            if 'J1P' in gene or 'J3P' in gene:
                continue

            raw_score = read.tags[0][1]  # raw because they don't include the gene choice probs. TODO oh wait shit this isn't right. raw_score isn't a prob. what the hell is it, anyway?
            score = self.get_choice_prob(region, gene) * raw_score  # multiply by the probability to choose this gene
            all_match_names[region].append((score,gene))
            all_query_bounds[gene] = (read.qstart, read.qend)
            # TODO the s-w allows the j right edge to be chopped off -- I should skip the matches where different amounts are chopped off in the query and germline
            all_germline_bounds[gene] = (read.pos, read.aend)
        try:
            self.summarize_query(query_name, query_seq, raw_best, all_match_names, all_query_bounds, all_germline_bounds)
        except AssertionError:
            print 'processing failed on query',query_name

    # ----------------------------------------------------------------------------------------
    def get_conserved_codon_position(self, region, gene, glbounds, qrbounds):
        if region == 'v':
            gl_cpos = self.pdriver.cyst_positions[gene]['cysteine-position']  # germline cystein position
            query_cpos = gl_cpos - glbounds[0] + qrbounds[0]  # cystein position in query sequence match
            # print '%d - %d ( + %d ) = %d (%d)' % (gl_cpos, glbounds[0], qrbounds[0], query_cpos, query_cpos + qrbounds[0])
            return query_cpos
        elif region == 'j':
            gl_tpos = int(self.pdriver.tryp_positions[gene])
            query_tpos = gl_tpos - glbounds[0]  #+ qrbounds[0]  # NOTE this coordinate is w/ respect to the j query *only*. Adjust below for coord in the whole query sequence
            query_tpos += qrbounds[0]  # add the length of the v match, and the length of query seq we'll be matching to d
            return query_tpos
            # print '%3d %s' % (tpos - glbounds[0], query_match_seq[tpos - glbounds[0]: tpos - glbounds[0] + 3])
        else:
            return -1

    # ----------------------------------------------------------------------------------------
    def check_conserved_codons(self, region, gene, query_name, query_seq, glbounds, qrbounds):  # TODO fix name conflict with fcn in utils
        codon_pos = self.get_conserved_codon_position(region, gene, glbounds, qrbounds)  # position in the query sequence, that is
        try:
            if region == 'v':
                utils.check_conserved_cysteine(query_seq[:qrbounds[1]], codon_pos, debug=False)  # don't use qrbounds[0] at start of slice so query_cpos is the same for all v matches, even if they don't start at the same place in the query sequence
            elif region == 'j':
                codon_pos -= qrbounds[0]  # subtract back of the length of the v match and prospective d match
                utils.check_conserved_tryptophan(query_seq[qrbounds[0]:qrbounds[1]], codon_pos, debug=False)
            return True
        except AssertionError:
            return False


    # ----------------------------------------------------------------------------------------
    def print_match(self, region, gene, query_seq, score, glbounds, qrbounds, tmp_codon_pos, skipping=False):
        if self.args.debug:
            buff_str = (17 - len(gene)) * ' '
            print '%8s%s%s%9.1e * %3.0f = %-6.1f' % (' ', utils.color_gene(gene), buff_str, self.get_choice_prob(region, gene), score / self.get_choice_prob(region, gene), score),
            print '%4d%4d   %s' % (glbounds[0], glbounds[1], self.pdriver.germline_seqs[region][gene][glbounds[0]:glbounds[1]])
            print '%48s  %4d%4d' % ('', qrbounds[0], qrbounds[1]),
            print '  %s (%d)' % (utils.color_mutants(self.pdriver.germline_seqs[region][gene][glbounds[0]:glbounds[1]], query_seq[qrbounds[0]:qrbounds[1]]), tmp_codon_pos),
            if skipping:
                print 'skipping!',
            print ''                

    # ----------------------------------------------------------------------------------------
    def summarize_query(self, query_name, query_seq, raw_best, all_match_names, all_query_bounds, all_germline_bounds):
        # best_scores = {}
        best, match_names = {}, {}
        n_matches, n_used, n_skipped = {'v':0, 'd':0, 'j':0}, {'v':0, 'd':0, 'j':0}, {'v':0, 'd':0, 'j':0}
        k_v_min, k_d_min = 999, 999
        k_v_max, k_d_max = 0, 0
        for region in utils.regions:
            all_match_names[region] = sorted(all_match_names[region], reverse=True)
            match_names[region] = []
        codon_positions = {'v':-1, 'd':-1, 'j':-1}  # conserved codon positions (v:cysteine, d:dummy, j:tryptophan)
        if self.args.debug:
            print '  %s' % query_name
        for region in utils.regions:
            for score,gene in all_match_names[region]:
                n_matches[region] += 1
                glbounds = all_germline_bounds[gene]
                qrbounds = all_query_bounds[gene]
                tmp_codon_pos = self.get_conserved_codon_position(region, gene, glbounds, qrbounds)  # position in the query sequence, that is

                if codon_positions[region] == -1:  # set to the position in the best match
                    codon_positions[region] = tmp_codon_pos

                # only use the best few matches
                if n_matches[region] > self.args.n_max_per_region:  # only take the top few from each region
                    # TODO should use *lots* of d matches, but fewer vs and js
                    break

                # add match to the list
                n_used[region] += 1
                match_names[region].append(gene)
                self.print_match(region, gene, query_seq, score, glbounds, qrbounds, tmp_codon_pos, skipping=False)

                if region == 'v':
                    this_k_v = all_query_bounds[gene][1]
                    k_v_min = min(this_k_v, k_v_min)
                    k_v_max = max(this_k_v, k_v_max)
                if region == 'd':
                    this_k_d = all_query_bounds[gene][1] - all_query_bounds[raw_best['v']][1]  # end of d minus end of v
                    k_d_min = min(this_k_d, k_d_min)
                    k_d_max = max(this_k_d, k_d_max)

                # check consistency with best match (since the best match is excised in s-w code, and because stochhmm is run with *one* k_v k_d set)
                if region not in best:
                    best[region] = gene
                    best[region + '_gl_seq'] = self.pdriver.germline_seqs[region][gene][glbounds[0]:glbounds[1]]
                    best[region + '_qr_seq'] = query_seq[qrbounds[0]:qrbounds[1]]
                    # best_scores[region] = score

                glmatchseq = self.pdriver.germline_seqs[region][gene][glbounds[0]:glbounds[1]]
                assert len(glmatchseq) == len(query_seq[qrbounds[0]:qrbounds[1]])  # neurotic double check (um, I think)
                        
        # print how many of the available matches we used
        if self.args.debug:
            print '         used',
            for region in utils.regions:
                if region != 'v':
                    print '            ',
                print ' %d / %d in %s' % (n_used[region], n_matches[region], region),
                if n_skipped[region] != 0:
                    print ' (skipped %d)' % (n_skipped[region]),
                print ''

        for region in utils.regions:
            if region not in best:
                print ' no',region,'match found for',query_name
                print '    true:'
                if not self.args.is_data:
                    utils.print_reco_event(self.pdriver.germline_seqs, self.pdriver.reco_info[query_name], 0, 0, extra_str='    ')
                assert False

        # write sw info to file for hmm
        # best k_v, k_d:
        k_v = all_query_bounds[best['v']][1]  # end of v match
        k_d = all_query_bounds[best['d']][1] - all_query_bounds[best['v']][1]  # end of d minus end of v

        if k_d_max < 5:  # since the s-w step matches to the longest possible j and then excises it, this sometimes gobbles up the d, resulting in a very short d alignment.
            # if self.args.debug:
            print '  expanding k_d'
            k_d_max = max(8, k_d_max)
            
        v_right_length = len(self.pdriver.germline_seqs['v'][best['v']]) - all_germline_bounds[best['v']][0]  # germline v length minus (germline) start of v match
        if 'IGHJ4*' in best['j'] and self.pdriver.germline_seqs['d'][best['d']][-5:] == 'ACTAC':  # the end of some d versions is the same as the start of some j versions, so the s-w frequently kicks out the 'wrong' alignment
            # print '  doubly expanding k_d'
            if k_d_max-k_d_min < 8:
                k_d_min -= 5
                k_d_max += 2

        k_v_min = max(0, k_v_min - self.pdriver.default_v_fuzz)  # ok, so I don't *actually* want it to be zero... oh, well
        k_v_max += self.pdriver.default_v_fuzz
        k_d_min = max(1, k_d_min - self.pdriver.default_d_fuzz)
        k_d_max += self.pdriver.default_d_fuzz
        assert k_v_min > 0 and k_d_min > 0 and k_v_max > 0 and k_d_max > 0 and v_right_length > 0

        # print '  k_v min %d max %d (best %d)' % (k_v_min, k_v_max, k_v)
        # print '  k_d min %d max %d (best %d)' % (k_d_min, k_d_max, k_d)
        self.pdriver.match_names_for_all_queries |= set(match_names['v'] + match_names['d'] + match_names['j'])

        assert query_name not in self.info
        self.info[query_name] = {}
        self.info[query_name]['k_v'] = {'best':k_v, 'min':k_v_min, 'max':k_v_max}
        self.info[query_name]['k_d'] = {'best':k_d, 'min':k_d_min, 'max':k_d_max}
        self.info[query_name]['v_right_length'] = v_right_length
        self.info[query_name]['best'] = best
        self.info[query_name]['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])

        assert codon_positions['v'] != -1
        assert codon_positions['j'] != -1
        self.info[query_name]['cdr3_length'] = codon_positions['j'] - codon_positions['v'] + 3  #tryp_position_in_joined_seq - self.cyst_position + 3
        # erosion, insertion, mutation info for best match
        self.info[query_name]['v_3p_del'] = len(self.pdriver.germline_seqs['v'][best['v']]) - all_germline_bounds[best['v']][1]  # len(germline v) - gl_match_end
        self.info[query_name]['d_5p_del'] = all_germline_bounds[best['d']][0]
        self.info[query_name]['d_3p_del'] = len(self.pdriver.germline_seqs['d'][best['d']]) - all_germline_bounds[best['d']][1]
        self.info[query_name]['j_5p_del'] = all_germline_bounds[best['j']][0]
        self.info[query_name]['vd_insertion'] = query_seq[all_query_bounds[best['v']][1] : all_query_bounds[best['d']][0]]
        self.info[query_name]['dj_insertion'] = query_seq[all_query_bounds[best['d']][1] : all_query_bounds[best['j']][0]]

        for region in utils.regions:
            self.info[query_name][region + '_gene'] = best[region]
        self.pcounter.increment(self.info[query_name], best)
        # tmp_line = {}
        # tmp_line['seq'] = query_seq
        # for region in utils.regions:
        #     tmp_line[region + '_gene'] = best[region]
        # tmp_line['v_3p_del'] = self.info[query_name]['v_3p_del']
        # tmp_line['d_5p_del'] = self.info[query_name]['d_5p_del']
        # tmp_line['d_3p_del'] = self.info[query_name]['d_3p_del']
        # tmp_line['j_5p_del'] = self.info[query_name]['j_5p_del']
        # tmp_line['vd_insertion'] = self.info[query_name]['vd_insertion']
        # tmp_line['dj_insertion'] = self.info[query_name]['dj_insertion']
        
        # utils.print_reco_event(self.pdriver.germline_seqs, self.pdriver.reco_info[query_name], 0, 0, extra_str='    ')
        # utils.print_reco_event(self.pdriver.germline_seqs, tmp_line, 0, 0, extra_str='    ')
 





