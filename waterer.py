import time
import sys
import os
import itertools
import operator
import pysam
import contextlib
from subprocess import check_call

from utils import utils
from utils.opener import opener

# ----------------------------------------------------------------------------------------
class Waterer(object):
    """ Run smith-waterman on the query sequences in <infname> """
    def __init__(self, pdriver):
        self.pdriver = pdriver
        self.sw_info = {}
        outfname = self.pdriver.workdir + '/query-seqs.bam'
        self.run_smith_waterman(outfname)
        self.read_smith_waterman(outfname)

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
    def read_smith_waterman(self, outfname):
        """
        Read bamfile output by s-w step, and write the info (after a bit of collation) to a csv.
        Note that the only reason to bother writing the csv is so you can avoid rerunning the s-w step every time you run the hmm.
        """
        start = time.time()
        gene_choice_probs = utils.read_overall_gene_prob(self.pdriver.datadir + '/human-beings/' + self.pdriver.args.human + '/' + self.pdriver.args.naivety)
        with contextlib.closing(pysam.Samfile(outfname)) as bam:
            grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
            for _, reads in grouped:  # loop over query sequences
                reads = list(reads)
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
                    choice_prob = 0.0  # set to zero if we didn't see it in data. kinda hacky, I suppose
                    if gene in gene_choice_probs[region]:
                        choice_prob = gene_choice_probs[region][gene]
                    score = choice_prob * raw_score  # multiply by the probability to choose this gene
                    all_match_names[region].append((score,gene))
                    all_query_bounds[gene] = (read.qstart, read.qend)
                    # TODO the s-w allows the j right edge to be chopped off -- I should skip the matches where different amounts are chopped off in the query and germline
                    all_germline_bounds[gene] = (read.pos, read.aend)
    
                best = {}
                match_names = {}
                n_matches = {'v':0, 'd':0, 'j':0}
                n_used = {'v':0, 'd':0, 'j':0}
                k_d_min = 999
                k_d_max = 0
                k_v_min = 999
                k_v_max = 0
                for region in utils.regions:
                    all_match_names[region] = sorted(all_match_names[region], reverse=True)
                    match_names[region] = []
                for region in utils.regions:
                    for score,gene in all_match_names[region]:
                        n_matches[region] += 1
                        if n_matches[region] > self.pdriver.args.n_max_per_region:  # only take the top few from each region
                            # TODO should use *lots* of d matches, but fewer vs and js
                            continue
                        n_used[region] += 1
                        match_names[region].append(gene)

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
    
                        if self.pdriver.args.debug:
                            buff_str = (17 - len(gene)) * ' '
                            print '%8s%s%s%6.1e * %2.0f = %-6.1f' % (' ', utils.color_gene(gene), buff_str, gene_choice_probs[region][gene], score / gene_choice_probs[region][gene], score),
                            glbounds = all_germline_bounds[gene]
                            qrbounds = all_query_bounds[gene]
                            print ' %4d%4d   %s' % (glbounds[0], glbounds[1], self.pdriver.germline_seqs[region][gene][glbounds[0]:glbounds[1]])
                            print '%46s  %4d%4d   %s' % ('', qrbounds[0], qrbounds[1], utils.color_mutants(self.pdriver.germline_seqs[region][gene][glbounds[0]:glbounds[1]], query_seq[qrbounds[0]:qrbounds[1]]))
                                
                # print how many of the available matches we used
                if self.pdriver.args.debug:
                    print '  used',
                    for region in utils.regions:
                        if region != 'v':
                            print '      ',
                        print ' %d / %d in %s' % (n_used[region], n_matches[region], region)

                for region in utils.regions:
                    if region not in best:
                        print ' no',region,'match found for',query_name
                        print '    true:'
                        utils.print_reco_event(self.pdriver.germline_seqs, self.pdriver.input_info[query_name], 0, 0, extra_str='    ')
                # write sw info to file for hmm
                # best k_v, k_d:
                k_v = all_query_bounds[best['v']][1]  # end of v match
                k_d = all_query_bounds[best['d']][1] - all_query_bounds[best['v']][1]  # end of d minus end of v

                if k_d_max < 5:  # since the s-w step matches to the longest possible j and then excises it, this sometimes gobbles up the d, resulting in a very short d alignment.
                    # if self.pdriver.args.debug:
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

                assert query_name not in self.sw_info
                self.sw_info[query_name] = {}
                self.sw_info[query_name]['k_v'] = k_v
                self.sw_info[query_name]['k_d'] = k_d
                self.sw_info[query_name]['k_v_min'] = k_v_min
                self.sw_info[query_name]['k_v_max'] = k_v_max
                self.sw_info[query_name]['k_d_min'] = k_d_min
                self.sw_info[query_name]['k_d_max'] = k_d_max
                self.sw_info[query_name]['v_right_length'] = v_right_length
                self.sw_info[query_name]['best'] = best
                self.sw_info[query_name]['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])
        print '    sw read time: %.3f' % (time.time() - start)
        os.remove(outfname)
