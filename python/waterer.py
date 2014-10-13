import time
import sys
import json
import csv
import os
import itertools
import operator
import pysam
import contextlib
from subprocess import check_call

import utils
from opener import opener
from parametercounter import ParameterCounter

# ----------------------------------------------------------------------------------------
class Waterer(object):
    """ Run smith-waterman on the query sequences in <infname> """
    def __init__(self, args, input_info, reco_info, germline_seqs, from_scratch=True, parameter_dir='', write_parameters=False, plotdir=''):
        if write_parameters or not from_scratch:
            assert parameter_dir != ''
        self.args = args
        self.input_info = input_info
        self.reco_info = reco_info
        self.germline_seqs = germline_seqs
        self.pcounter = None
        if write_parameters:
            self.pcounter = ParameterCounter(self.germline_seqs, parameter_dir, plotdir=plotdir)
            if plotdir != '':
                utils.prep_dir(plotdir + '/plots', '*.svg')
        self.info = {}
        self.info['all_best_matches'] = set()  # set of all the matches we found
        self.from_scratch = from_scratch
        if not self.from_scratch:
            if self.args.debug:
                print '  reading gene choice probs from',parameter_dir
            self.gene_choice_probs = utils.read_overall_gene_prob(parameter_dir)

        with opener('r')(self.args.datadir + '/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
            self.cyst_positions = json.load(json_file)
        with opener('r')(self.args.datadir + '/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region (TGG)
            tryp_reader = csv.reader(csv_file)
            self.tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

    # ----------------------------------------------------------------------------------------
    def run(self):
        outfname = self.args.workdir + '/query-seqs.bam'
        self.run_smith_waterman(outfname)
        self.read_output(outfname, plot_performance=self.args.plot_performance)
        if self.pcounter != None:
            self.pcounter.write_counts()

    # ----------------------------------------------------------------------------------------
    def clean(self):
        if self.pcounter != None:
            self.pcounter.clean()

    # ----------------------------------------------------------------------------------------
    def run_smith_waterman(self, outfname):
        """
        Run smith-waterman alignment on the seqs in <infname>, and toss all the top matches into <outfname>.
        Then run through <outfname> to get the top hits and their locations to pass to the hmm.
        Then run the hmm on each gene set.
        """
        infname = self.args.workdir + '/query-seqs.fa'
        with opener('w')(infname) as swinfile:  # write *all* the input seqs to file, i.e. run s-w on all of 'em at once
            for query_name, line in self.input_info.iteritems():
                swinfile.write('>' + query_name + ' NUKES\n')
                swinfile.write(line['seq'] + '\n')
        start = time.time()
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        check_call(self.args.ighutil_dir + '/bin/vdjalign align-fastq --j-subset adaptive --max-drop 50 --match 5 --mismatch 3 --gap-open 1000 ' + infname + ' ' + outfname, shell=True)
        if not self.args.no_clean:
            os.remove(infname)
        print '    s-w time: %.3f' % (time.time()-start)
    
    # ----------------------------------------------------------------------------------------
    def read_output(self, outfname, plot_performance=False):
        start = time.time()

        perfplotter = None
        if plot_performance:
            assert not self.args.is_data
            from performanceplotter import PerformancePlotter
            perfplotter = PerformancePlotter(self.germline_seqs, os.getenv('www') + '/partis/sw_performance', 'sw')

        with contextlib.closing(pysam.Samfile(outfname)) as bam:
            grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
            for _, reads in grouped:  # loop over query sequences
                self.process_query(bam, list(reads), perfplotter)

        if perfplotter != None:
            perfplotter.plot()

        print '    sw read time: %.3f' % (time.time() - start)
        if not self.args.no_clean:
            os.remove(outfname)

    # ----------------------------------------------------------------------------------------
    def get_choice_prob(self, region, gene):
        choice_prob = 1.0
        if not self.from_scratch:
            if gene in self.gene_choice_probs[region]:
                choice_prob = self.gene_choice_probs[region][gene]
            else:
                choice_prob = 0.0  # TODO choose something else?
        return choice_prob

    # ----------------------------------------------------------------------------------------
    def process_query(self, bam, reads, perfplotter=None):
        primary = next((r for r in reads if not r.is_secondary), None)
        query_seq = primary.seq
        query_name = primary.qname
        raw_best = {}
        all_match_names = {}
        warnings = {}  # ick, this is a messy way to pass stuff around
        for region in utils.regions:
            all_match_names[region] = []
        all_query_bounds, all_germline_bounds = {}, {}
        for read in reads:  # loop over the matches found for each query sequence
            read.seq = query_seq  # only the first one has read.seq set by default, so we need to set the rest by hand
            gene = bam.references[read.tid]
            region = utils.get_region(gene)
            warnings[gene] = ''

            if region not in raw_best:  # best v, d, and j before multiplying by gene choice probs. needed 'cause *these* are the v and j that get excised
                raw_best[region] = gene

            raw_score = read.tags[0][1]  # raw because they don't include the gene choice probs. TODO oh wait shit this isn't right. raw_score isn't a prob. what the hell is it, anyway?
            score = self.get_choice_prob(region, gene) * raw_score  # multiply by the probability to choose this gene
            # set bounds  TODO the s-w allows the j right edge to be chopped off -- I should skip the matches where different amounts are chopped off in the query and germline. EDIT well, maybe. I dunno
            qrbounds = (read.qstart, read.qend)
            glbounds = (read.pos, read.aend)

            # check for left-side v erosion  TODO it would make sense to change the score when you expand the boundaries
            if region == 'v' and read.qstart != 0:  # s-w allows the v match to start to the right of the lefthand base in the query sequence, which makes little sense in most cases. So, we artificially expand the left v boundaries, presumably adding a few mutated bases. TODO think about whether this is the proper long-term solution
                qrbounds = (0, read.qend)
                glbounds = (read.pos - read.qstart, read.aend)
                warnings[gene] += 'v left expanded %d -> %d, %d -> %d' % (read.qstart, 0, read.pos, read.pos - read.qstart)
            # check for right-side j erosion  TODO it would make sense to change the score when you expand the boundaries
            gl_length = len(self.germline_seqs[region][gene])
            if region == 'j' and read.aend != gl_length:
                # print 'whoa!', gene, read.pos, read.aend, gl_length
                length_to_end_of_query_seq = len(query_seq) - read.qstart
                new_match_length = min(gl_length, length_to_end_of_query_seq + read.pos)  # don't want to expand beyond the number of query bases we have left
                # print gene, read.qstart, '-', read.qend, '   ', read.pos, '-', read.aend, '    ', gl_length, length_to_end_of_query_seq, new_match_length
                qrbounds = (read.qstart, read.qstart + new_match_length - read.pos)
                glbounds = (read.pos, new_match_length)
                warnings[gene] += 'j right expanded %d -> %d, %d -> %d' % (read.qend, read.qstart + new_match_length - read.pos, read.aend, new_match_length)
            assert qrbounds[1] <= len(query_seq)
            assert glbounds[1] <= len(self.germline_seqs[region][gene])

            # if region == 'j' and qrbounds[1] != len(query_seq):  # TODO it's debatable whether this is a good idea. As far as I can tell, this is certainly going to be a crap match
            #     print 'WARNING couldn\'t expand to right side of query seq, so skipping %s ' % utils.color_gene(gene)
            #     continue
            
            all_match_names[region].append((score,gene))
            all_query_bounds[gene] = qrbounds
            all_germline_bounds[gene] = glbounds

        self.summarize_query(query_name, query_seq, raw_best, all_match_names, all_query_bounds, all_germline_bounds, perfplotter, warnings)

    # ----------------------------------------------------------------------------------------
    def get_conserved_codon_position(self, region, gene, glbounds, qrbounds):
        """
        Find location of the conserved cysteine/tryptophan in a query sequence given a germline match which is specified by
        its germline bounds <glbounds> and its bounds in the query sequence <qrbounds>
        """
        if region == 'v':
            gl_cpos = self.cyst_positions[gene]['cysteine-position']  # germline cysteine position
            query_cpos = gl_cpos - glbounds[0] + qrbounds[0]  # cysteine position in query sequence match
            # print '%d - %d ( + %d ) = %d (%d)' % (gl_cpos, glbounds[0], qrbounds[0], query_cpos, query_cpos + qrbounds[0])
            return query_cpos
        elif region == 'j':
            gl_tpos = int(self.tryp_positions[gene])
            query_tpos = gl_tpos - glbounds[0]  #+ qrbounds[0]  # NOTE this coordinate is w/ respect to the j query *only*. Adjust below for coord in the whole query sequence
            query_tpos += qrbounds[0]  # add the length of the v match, and the length of query seq we'll be matching to d
            return query_tpos
            # print '%3d %s' % (tpos - glbounds[0], query_match_seq[tpos - glbounds[0]: tpos - glbounds[0] + 3])
        else:
            return -1

    # ----------------------------------------------------------------------------------------
    def check_conserved_codons(self, region, gene, query_seq, glbounds, qrbounds):  # TODO fix name conflict with fcn in utils
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
    def print_match(self, region, gene, query_seq, score, glbounds, qrbounds, codon_pos, warnings, skipping=False):
        if self.args.debug:
            buff_str = (17 - len(gene)) * ' '
            tmp_val = 0.0
            try:
                tmp_val = score / self.get_choice_prob(region, gene)
            except ZeroDivisionError:
                pass
            print '%8s%s%s%9.1e * %3.0f = %-6.1f' % (' ', utils.color_gene(gene), buff_str, self.get_choice_prob(region, gene), tmp_val, score),
            print '%4d%4d   %s' % (glbounds[0], glbounds[1], self.germline_seqs[region][gene][glbounds[0]:glbounds[1]]),
            if region == 'v':
                v_right_length = len(self.germline_seqs[region][gene]) - glbounds[0]  # germline v length minus (germline) start of v match
                print '(v right %d)' % v_right_length,
            print ''
            print '%48s  %4d%4d' % ('', qrbounds[0], qrbounds[1]),
            print '  %s (%d)' % (utils.color_mutants(self.germline_seqs[region][gene][glbounds[0]:glbounds[1]], query_seq[qrbounds[0]:qrbounds[1]]), codon_pos),
            if warnings[gene] != '':
                print 'WARNING',warnings[gene],
            if skipping:
                print 'skipping!',
            print ''                

    # ----------------------------------------------------------------------------------------
    def shift_overlapping_boundaries(self, qrbounds, glbounds, query_name, best):
        """ s-w allows d and j matches (and v and d matches) to overlap... which makes no sense, so arbitrarily give the disputed territory to the righthand region  """  # TODO do something better than using the righthand one
        for rpairs in ({'left':'v', 'right':'d'}, {'left':'d', 'right':'j'}):
            left_gene = best[rpairs['left']]
            right_gene = best[rpairs['right']]
            overlap = qrbounds[left_gene][1] - qrbounds[right_gene][0]
            if overlap > 0:
                print 'WARNING %s giving %d bases from %s match to %s match' % (query_name, overlap, rpairs['left'], rpairs['right'])
                qrbounds[left_gene] = (qrbounds[left_gene][0], qrbounds[left_gene][1] - overlap)
                glbounds[left_gene] = (glbounds[left_gene][0], glbounds[left_gene][1] - overlap)

    # ----------------------------------------------------------------------------------------
    def summarize_query(self, query_name, query_seq, raw_best, all_match_names, all_query_bounds, all_germline_bounds, perfplotter, warnings):
        # best_scores = {}
        best, match_names, n_matches = {}, {}, {}
        n_used = {'v':0, 'd':0, 'j':0}
        k_v_min, k_d_min = 999, 999
        k_v_max, k_d_max = 0, 0
        for region in utils.regions:
            all_match_names[region] = sorted(all_match_names[region], reverse=True)
            match_names[region] = []
        codon_positions = {'v':-1, 'd':-1, 'j':-1}  # conserved codon positions (v:cysteine, d:dummy, j:tryptophan)
        if self.args.debug:
            print '  smith waterman matches for %s:' % query_name
        for region in utils.regions:
            n_matches[region] = len(all_match_names[region])
            n_skipped = 0
            for score,gene in all_match_names[region]:
                glbounds = all_germline_bounds[gene]
                qrbounds = all_query_bounds[gene]
                assert qrbounds[1] <= len(query_seq)  # NOTE I'm putting these up avove as well (in process_query), so in time I should remove them from here
                assert glbounds[1] <= len(self.germline_seqs[region][gene])
                glmatchseq = self.germline_seqs[region][gene][glbounds[0]:glbounds[1]]
                if codon_positions[region] == -1:  # set to the position in the best match
                    codon_positions[region] = self.get_conserved_codon_position(region, gene, glbounds, qrbounds)  # position in the query sequence, that is

                if region == 'j' and qrbounds[1] != len(query_seq):  # TODO it's debatable whether this is a good idea. As far as I can tell, this is certainly going to be a crap match
                    print 'WARNING couldn\'t expand to right side of query seq, so skipping %s ' % utils.color_gene(gene)
                    continue
            
                # only use the best few matches
                if n_used[region] >= self.args.n_max_per_region:  # only take the top few from each region. TODO should use *lots* of d matches, but fewer vs and js
                    break

                # only use a specified set of genes
                if self.args.only_genes != None and gene not in self.args.only_genes:
                    n_skipped += 1
                    continue

                # add match to the list
                n_used[region] += 1
                match_names[region].append(gene)
                self.print_match(region, gene, query_seq, score, glbounds, qrbounds, codon_positions[region], warnings, skipping=False)

                # if the germline match and the query match aren't the same length, s-w likely added an insert, which we shouldn't get since the gap-open penalty is jacked up so high
                if len(glmatchseq) != len(query_seq[qrbounds[0]:qrbounds[1]]):
                    print glmatchseq
                    print '',query_seq[qrbounds[0]:qrbounds[1]]
                    print gene,query_name
                    print qrbounds[0],qrbounds[1]
                    print glbounds[0],glbounds[1]
                    print warnings[gene]
                assert len(glmatchseq) == len(query_seq[qrbounds[0]:qrbounds[1]])  # neurotic double check (um, I think) EDIT hey this totally saved my ass

                if region == 'v':
                    this_k_v = all_query_bounds[gene][1]
                    k_v_min = min(this_k_v, k_v_min)
                    k_v_max = max(this_k_v, k_v_max)
                if region == 'd':
                    this_k_d = all_query_bounds[gene][1] - all_query_bounds[raw_best['v']][1]  # end of d minus end of v
                    k_d_min = min(this_k_d, k_d_min)
                    k_d_max = max(this_k_d, k_d_max)

                # check consistency with best match (since the best match is excised in s-w code, and because ham is run with *one* k_v k_d set)
                if region not in best:
                    best[region] = gene
                    best[region + '_gl_seq'] = self.germline_seqs[region][gene][glbounds[0]:glbounds[1]]
                    best[region + '_qr_seq'] = query_seq[qrbounds[0]:qrbounds[1]]
                    # best_scores[region] = score

            if self.args.debug and n_skipped > 0:
                print '%8s skipped %d %s genes' % ('', n_skipped, region)
                        
        # print how many of the available matches we used
        if self.args.debug:
            print '         used',
            for region in utils.regions:
                if region != 'v':
                    print '             ',
                print ' %d / %d in %s' % (n_used[region], n_matches[region], region)

        for region in utils.regions:
            if region not in best:
                print ' no',region,'match found for',query_name  # TODO if no d match found, should just assume entire d was eroded
                if not self.args.is_data:
                    print '    true:'
                    utils.print_reco_event(self.germline_seqs, self.reco_info[query_name], 0, 0, extra_str='    ')
                return

        # s-w allows d and j matches to overlap... which makes no sense, so arbitrarily give the disputed territory to j
        self.shift_overlapping_boundaries(all_query_bounds, all_germline_bounds, query_name, best)

        # best k_v, k_d:
        k_v = all_query_bounds[best['v']][1]  # end of v match
        k_d = all_query_bounds[best['d']][1] - all_query_bounds[best['v']][1]  # end of d minus end of v

        if k_d_max < 5:  # since the s-w step matches to the longest possible j and then excises it, this sometimes gobbles up the d, resulting in a very short d alignment.
            if self.args.debug:
                print '  expanding k_d'
            k_d_max = max(8, k_d_max)
            
        v_right_length = len(self.germline_seqs['v'][best['v']]) - all_germline_bounds[best['v']][0]  # germline v length minus (germline) start of v match
        if 'IGHJ4*' in best['j'] and self.germline_seqs['d'][best['d']][-5:] == 'ACTAC':  # the end of some d versions is the same as the start of some j versions, so the s-w frequently kicks out the 'wrong' alignment
            if self.args.debug:
                print '  doubly expanding k_d'
            if k_d_max-k_d_min < 8:
                k_d_min -= 5
                k_d_max += 2

        k_v_min = max(0, k_v_min - self.args.default_v_fuzz)  # ok, so I don't *actually* want it to be zero... oh, well
        k_v_max += self.args.default_v_fuzz
        k_d_min = max(1, k_d_min - self.args.default_d_fuzz)
        k_d_max += self.args.default_d_fuzz
        assert k_v_min > 0 and k_d_min > 0 and k_v_max > 0 and k_d_max > 0 and v_right_length > 0

        assert query_name not in self.info
        self.info[query_name] = {}
        self.info[query_name]['k_v'] = {'best':k_v, 'min':k_v_min, 'max':k_v_max}
        self.info[query_name]['k_d'] = {'best':k_d, 'min':k_d_min, 'max':k_d_max}
        self.info[query_name]['v_right_length'] = v_right_length
        self.info[query_name]['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])

        assert codon_positions['v'] != -1
        assert codon_positions['j'] != -1
        self.info[query_name]['cdr3_length'] = codon_positions['j'] - codon_positions['v'] + 3  #tryp_position_in_joined_seq - self.cyst_position + 3

        # erosion, insertion, mutation info for best match
        self.info[query_name]['v_5p_del'] = all_germline_bounds[best['v']][0]
        self.info[query_name]['v_3p_del'] = len(self.germline_seqs['v'][best['v']]) - all_germline_bounds[best['v']][1]  # len(germline v) - gl_match_end
        self.info[query_name]['d_5p_del'] = all_germline_bounds[best['d']][0]
        self.info[query_name]['d_3p_del'] = len(self.germline_seqs['d'][best['d']]) - all_germline_bounds[best['d']][1]
        self.info[query_name]['j_5p_del'] = all_germline_bounds[best['j']][0]
        self.info[query_name]['j_3p_del'] = len(self.germline_seqs['j'][best['j']]) - all_germline_bounds[best['j']][1]

        self.info[query_name]['vd_insertion'] = query_seq[all_query_bounds[best['v']][1] : all_query_bounds[best['d']][0]]
        self.info[query_name]['dj_insertion'] = query_seq[all_query_bounds[best['d']][1] : all_query_bounds[best['j']][0]]

        for region in utils.regions:
            self.info[query_name][region + '_gene'] = best[region]
            self.info[query_name][region + '_gl_seq'] = best[region + '_gl_seq']
            self.info[query_name][region + '_qr_seq'] = best[region + '_qr_seq']
            self.info['all_best_matches'].add(best[region])

        self.info[query_name]['seq'] = query_seq  # only need to add this so I can pass it to print_reco_event
        if self.args.debug:
            utils.print_reco_event(self.germline_seqs, self.info[query_name], extra_str='          ')

        if self.pcounter != None:
            self.pcounter.increment(self.info[query_name])
        if perfplotter != None:
            perfplotter.evaluate(self.reco_info[query_name], self.info[query_name], query_name)
