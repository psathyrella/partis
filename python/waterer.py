import time
import sys
import json
import math
import re
import csv
import os
import itertools
import operator
import pysam
import contextlib
from subprocess import check_call, check_output, Popen, PIPE

import utils
from opener import opener
from parametercounter import ParameterCounter

# ----------------------------------------------------------------------------------------
class Waterer(object):
    """ Run smith-waterman on the query sequences in <infname> """
    def __init__(self, args, input_info, reco_info, germline_seqs, parameter_dir, write_parameters=False, plotdir=None):
        self.parameter_dir = parameter_dir
        self.plotdir = plotdir
        self.args = args
        self.input_info = input_info
        self.reco_info = reco_info
        self.germline_seqs = germline_seqs
        self.pcounter, self.true_pcounter = None, None
        if write_parameters:
            self.pcounter = ParameterCounter(self.germline_seqs)
            if not self.args.is_data:
                self.true_pcounter = ParameterCounter(self.germline_seqs)
        self.info = {}
        self.info['queries'] = []
        self.info['all_best_matches'] = set()  # set of all the matches we found (for *all* queries)
        self.info['skipped_unproductive_queries'] = []  # list of unproductive queries
        self.info['skipped_indel_queries'] = []  # list of queries that had indels
        if self.args.apply_choice_probs_in_sw:
            if self.args.debug:
                print '  reading gene choice probs from', parameter_dir
            self.gene_choice_probs = utils.read_overall_gene_probs(parameter_dir)

        with opener('r')(self.args.datadir + '/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
            self.cyst_positions = json.load(json_file)
        with opener('r')(self.args.datadir + '/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region (TGG)
            tryp_reader = csv.reader(csv_file)
            self.tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

        self.outfile = None
        if self.args.outfname != None:
            self.outfile = open(self.args.outfname, 'a')

        self.n_unproductive = 0
        self.n_total = 0

    # ----------------------------------------------------------------------------------------
    def __del__(self):
        if self.args.outfname != None:
            self.outfile.close()

    # ----------------------------------------------------------------------------------------
    def clean(self):
        if self.pcounter != None:
            self.pcounter.clean()
        if self.true_pcounter != None:
            self.true_pcounter.clean()

    # ----------------------------------------------------------------------------------------
    def run(self):
        print 'smith-waterman'
        # start = time.time()

        base_infname = 'query-seqs.fa'
        base_outfname = 'query-seqs.bam'
        sys.stdout.flush()
        n_procs = self.args.n_fewer_procs
        self.write_vdjalign_input(base_infname, n_procs=n_procs)
        if n_procs == 1:
            cmd_str = self.get_vdjalign_cmd_str(self.args.workdir, base_infname, base_outfname)
            check_call(cmd_str.split())
            if not self.args.no_clean:
                os.remove(self.args.workdir + '/' + base_infname)
        else:
            procs = []
            for iproc in range(n_procs):
                cmd_str = self.get_vdjalign_cmd_str(self.args.workdir + '/sw-' + str(iproc), base_infname, base_outfname)
                procs.append(Popen(cmd_str.split(), stdout=PIPE, stderr=PIPE))
                time.sleep(0.1)
            for iproc in range(len(procs)):
                out, err = procs[iproc].communicate()
                utils.process_out_err(out, err, extra_str=str(iproc))
            if not self.args.no_clean:
                for iproc in range(n_procs):
                    os.remove(self.args.workdir + '/sw-' + str(iproc) + '/' + base_infname)

        sys.stdout.flush()
        self.read_output(base_outfname, plot_performance=self.args.plot_performance, n_procs=n_procs)
        # print '    sw time: %.3f' % (time.time()-start)
        if self.n_unproductive > 0:
            print '      unproductive skipped %d / %d = %.2f' % (self.n_unproductive, self.n_total, float(self.n_unproductive) / self.n_total)
        if len(self.info['skipped_indel_queries']) > 0:
            print '      indels skipped %d / %d = %.2f' % (len(self.info['skipped_indel_queries']), self.n_total, float(len(self.info['skipped_indel_queries'])) / self.n_total)
        if self.pcounter != None:
            self.pcounter.write(self.parameter_dir)
            # if self.true_pcounter != None:
            #     self.true_pcounter.write(parameter_xxx_dir, plotdir=plotdir + '/true')
            if not self.args.no_plot and self.plotdir != '':
                self.pcounter.plot(self.plotdir, subset_by_gene=True, cyst_positions=self.cyst_positions, tryp_positions=self.tryp_positions)
                if self.true_pcounter != None:
                    self.true_pcounter.plot(self.plotdir + '/true', subset_by_gene=True, cyst_positions=self.cyst_positions, tryp_positions=self.tryp_positions)

    # ----------------------------------------------------------------------------------------
    def write_vdjalign_input(self, base_infname, n_procs):
        # first make a list of query names so we can iterate over an ordered collection
        ordered_info = []
        for query_name in self.input_info:
            ordered_info.append(query_name)

        queries_per_proc = float(len(self.input_info)) / n_procs
        n_queries_per_proc = int(math.ceil(queries_per_proc))
        if n_procs == 1:  # double check for rounding problems or whatnot
            assert n_queries_per_proc == len(self.input_info)
        for iproc in range(n_procs):
            workdir = self.args.workdir
            if n_procs > 1:
                workdir += '/sw-' + str(iproc)
                utils.prep_dir(workdir)
            with opener('w')(workdir + '/' + base_infname) as sub_infile:
                for iquery in range(iproc*n_queries_per_proc, (iproc + 1)*n_queries_per_proc):
                    if iquery >= len(ordered_info):
                        break
                    query_name = ordered_info[iquery]
                    sub_infile.write('>' + query_name + ' NUKES\n')
                    sub_infile.write(self.input_info[query_name]['seq'] + '\n')

    # ----------------------------------------------------------------------------------------
    def get_vdjalign_cmd_str(self, workdir, base_infname, base_outfname):
        """
        Run smith-waterman alignment (from Connor's ighutils package) on the seqs in <base_infname>, and toss all the top matches into <base_outfname>.
        """
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        os.environ['PATH'] = os.getenv('PWD') + '/packages/samtools:' + os.getenv('PATH')
        check_output(['which', 'samtools'])
        if not os.path.exists(self.args.ighutil_dir + '/bin/vdjalign'):
            raise Exception('ERROR ighutil path d.n.e: ' + self.args.ighutil_dir + '/bin/vdjalign')
        cmd_str = self.args.ighutil_dir + '/bin/vdjalign align-fastq -q'
        if self.args.slurm:
            cmd_str = 'srun ' + cmd_str
        cmd_str += ' --max-drop 50'
        cmd_str += ' --match 5 --mismatch 3'
        cmd_str += ' --gap-open ' + str(self.args.gap_open_penalty)  #1000'  #50'
        cmd_str += ' --vdj-dir ' + self.args.datadir
        cmd_str += ' ' + workdir + '/' + base_infname + ' ' + workdir + '/' + base_outfname

        return cmd_str

    # ----------------------------------------------------------------------------------------
    def read_output(self, base_outfname, plot_performance=False, n_procs=1):
        perfplotter = None
        if plot_performance:
            assert self.args.plotdir != None
            assert not self.args.is_data
            from performanceplotter import PerformancePlotter
            perfplotter = PerformancePlotter(self.germline_seqs, self.args.plotdir + '/sw/performance', 'sw')

        n_processed = 0
        for iproc in range(n_procs):
            workdir = self.args.workdir
            if n_procs > 1:
                workdir += '/sw-' + str(iproc)
            outfname = workdir + '/' + base_outfname
            with contextlib.closing(pysam.Samfile(outfname)) as bam:
                grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
                for _, reads in grouped:  # loop over query sequences
                    self.n_total += 1
                    self.process_query(bam, list(reads), perfplotter)
                    n_processed += 1

            if not self.args.no_clean:
                os.remove(outfname)
                if n_procs > 1:  # still need the top-level workdir
                    os.rmdir(workdir)

        print '    processed %d queries' % n_processed

        if perfplotter != None:
            perfplotter.plot()

    # ----------------------------------------------------------------------------------------
    def get_choice_prob(self, region, gene):
        choice_prob = 1.0
        if gene in self.gene_choice_probs[region]:
            choice_prob = self.gene_choice_probs[region][gene]
        else:
            choice_prob = 0.0  # NOTE would it make sense to use something else here?
        return choice_prob

    # ----------------------------------------------------------------------------------------
    def print_indel_info(self, query_name, cigarstr, qrseq, glseq, gene):
        cigars = re.findall('[0-9][0-9]*[A-Z]', cigarstr)  # split cigar string into its parts
        cigars = [(cstr[-1], int(cstr[:-1])) for cstr in cigars]  # split each part into the code and the length

        codestr = ''
        for code, length in cigars:
            codestr += length * code

        qrprintstr, glprintstr = '', ''
        iqr, igl = 0, 0
        for icode in range(len(codestr)):
            code = codestr[icode]
            if code == 'M':
                qrbase = qrseq[iqr]
                if qrbase != glseq[igl]:
                    qrbase = utils.color('red', qrbase)
                qrprintstr += qrbase
                glprintstr += glseq[igl]
            elif code == 'S':
                continue
            elif code == 'I':
                qrprintstr += utils.color('red', qrseq[iqr])
                glprintstr += utils.color('red', '*')
                igl -= 1
            elif code == 'D':
                qrprintstr += utils.color('red', '*')
                glprintstr += utils.color('red', glseq[igl])
                iqr -= 1
            else:
                raise Exception('unhandled code %s' % code)

            iqr += 1
            igl += 1

        print '\nWARNING in the reasoned opinion of Mssrs. Smith, Waterman (Esq.) the query %s appears to harbor a within-gene indel:' % query_name
        print '    %12s %s' % (gene, glprintstr)
        print '    %12s %s' % ('query', qrprintstr)
        print '     cigar:', cigarstr
        print 'Since the HMM does not currently support indels we\'re gonna tell it to *skip* this sequence.'
        print 'However, if you\'re pretty sure this really is an indel, you can simply add/remove the bases in question by hand and rerun.'
        print 'Note that we could do that automatically here... but that\'s really neither very transparent nor principled, so we\'re making you do it by hand for now.'

    # ----------------------------------------------------------------------------------------
    def process_query(self, bam, reads, perfplotter=None):
        primary = next((r for r in reads if not r.is_secondary), None)
        query_seq = primary.seq
        query_name = primary.qname
        raw_best = {}  # best gene for each region, as output by sw (we may decide later that we disagree on which are the best)
        first_match_query_bounds = None  # since sw excises its favorite v match, we have to know this match's boundaries in order to calculate k_d for all the other matches
        all_match_names = {}
        warnings = {}  # ick, this is a messy way to pass stuff around
        for region in utils.regions:
            all_match_names[region] = []
        all_query_bounds, all_germline_bounds = {}, {}
        n_skipped_invalid_cpos = 0
        for read in reads:  # loop over the matches found for each query sequence
            read.seq = query_seq  # only the first one has read.seq set by default, so we need to set the rest by hand
            gene = bam.references[read.tid]
            region = utils.get_region(gene)
            raw_score = read.tags[0][1]  # raw because they don't include the gene choice probs
            score = raw_score
            if self.args.apply_choice_probs_in_sw:  # NOTE I stopped applying the gene choice probs here because the smith-waterman scores don't correspond to log-probs, so throwing on the gene choice probs was dubious (and didn't seem to work that well)
                score = self.get_choice_prob(region, gene) * raw_score  # multiply by the probability to choose this gene
            qrbounds = (read.qstart, read.qend)
            glbounds = (read.pos, read.aend)
            if region == 'v' and first_match_query_bounds is None:
                first_match_query_bounds = qrbounds
            # if region not in raw_best:  # we need these 'cause *these* are the v and j that get excised by sw, so to calculate k_d for every other match we have to 
            #     raw_best[region] = gene

            if region == 'v':  # skip matches with cpos past the end of the query seq (i.e. eroded a ton on the right side of the v)
                cpos = utils.get_conserved_codon_position(self.cyst_positions, self.tryp_positions, 'v', gene, glbounds, qrbounds, assert_on_fail=False)
                if not utils.check_conserved_cysteine(self.germline_seqs['v'][gene], self.cyst_positions[gene]['cysteine-position'], assert_on_fail=False):  # some of the damn cysteine positions in the json file were wrong, so now we check
                    raise Exception('bad cysteine in %s: %d %s' % (gene, self.cyst_positions[gene]['cysteine-position'], self.germline_seqs['v'][gene]))
                if cpos < 0 or cpos >= len(query_seq):
                    n_skipped_invalid_cpos += 1
                    # if region not in raw_best:
                    #     print '    skipping best %s match with cpos %d (negative or greater than %d)' % (gene, cpos, len(query_seq))
                    continue

            if 'I' in read.cigarstring or 'D' in read.cigarstring:  # skip indels, and tell the HMM to skip indels
                if len(all_match_names[region]) == 0:  # if this is the first (best) match for this region, allow indels (otherwise skip the match)
                    print '  LEN', len(all_match_names[region])
                    self.print_indel_info(query_name, read.cigarstring, read.seq[qrbounds[0] : qrbounds[1]], self.germline_seqs[region][gene][glbounds[0] : glbounds[1]], gene)
                    self.info['skipped_indel_queries'].append(query_name)
                    self.info[query_name] = {'indels'}
                    return
                else:
                    continue

            if qrbounds[1]-qrbounds[0] != glbounds[1]-glbounds[0]:
                raise Exception('germline match (%d %d) not same length as query match (%d %d)' % (qrbounds[0], qrbounds[1], glbounds[0], glbounds[1]))

            assert qrbounds[1] <= len(query_seq)
            if glbounds[1] > len(self.germline_seqs[region][gene]):
                print '  ', gene
                print '  ', glbounds[1], len(self.germline_seqs[region][gene])
                print '  ', self.germline_seqs[region][gene]
            assert glbounds[1] <= len(self.germline_seqs[region][gene])

            assert qrbounds[1]-qrbounds[0] == glbounds[1]-glbounds[0]

            warnings[gene] = ''
            all_match_names[region].append((score, gene))  # NOTE it is important that this is ordered such that the best match is first
            all_query_bounds[gene] = qrbounds
            all_germline_bounds[gene] = glbounds

        # if n_skipped_invalid_cpos > 0:
        #     print '      skipped %d invalid cpos values for %s' % (n_skipped_invalid_cpos, query_name)
        self.summarize_query(query_name, query_seq, raw_best, all_match_names, all_query_bounds, all_germline_bounds, perfplotter, warnings, first_match_query_bounds)

    # ----------------------------------------------------------------------------------------
    def print_match(self, region, gene, query_seq, score, glbounds, qrbounds, codon_pos, warnings, skipping=False):
        if self.args.debug < 2:
            return
        out_str_list = []
        buff_str = (20 - len(gene)) * ' '
        tmp_val = score
        if self.args.apply_choice_probs_in_sw and self.get_choice_prob(region, gene) != 0.0:
            tmp_val = score / self.get_choice_prob(region, gene)
        if self.args.apply_choice_probs_in_sw:
            out_str_list.append('%8s%s%s%9.1e * %3.0f = %-6.1f' % (' ', utils.color_gene(gene), buff_str, self.get_choice_prob(region, gene), tmp_val, score))
        else:
            out_str_list.append('%8s%s%s%9s%3s %6.0f        ' % (' ', utils.color_gene(gene), '', '', buff_str, score))
        out_str_list.append('%4d%4d   %s\n' % (glbounds[0], glbounds[1], self.germline_seqs[region][gene][glbounds[0]:glbounds[1]]))
        out_str_list.append('%46s  %4d%4d' % ('', qrbounds[0], qrbounds[1]))
        out_str_list.append('   %s ' % (utils.color_mutants(self.germline_seqs[region][gene][glbounds[0]:glbounds[1]], query_seq[qrbounds[0]:qrbounds[1]])))
        if region != 'd':
            out_str_list.append('(%s %d)' % (utils.conserved_codon_names[region], codon_pos))
        if warnings[gene] != '':
            out_str_list.append('WARNING ' + warnings[gene])
        if skipping:
            out_str_list.append('skipping!')
        if self.args.outfname == None:
            print ''.join(out_str_list)
        else:
            out_str_list.append('\n')
            self.outfile.write(''.join(out_str_list))

    # ----------------------------------------------------------------------------------------
    def shift_overlapping_boundaries(self, qrbounds, glbounds, query_name, query_seq, best):
        # NOTE this does pretty much the same thing as resolve_overlapping_matches in joinparser.py
        """ s-w allows d and j matches (and v and d matches) to overlap... which makes no sense, so apportion the disputed territory between the two regions """
        for region_pairs in ({'left':'v', 'right':'d'}, {'left':'d', 'right':'j'}):
            l_reg = region_pairs['left']
            r_reg = region_pairs['right']
            l_gene = best[l_reg]
            r_gene = best[r_reg]
            overlap = qrbounds[l_gene][1] - qrbounds[r_gene][0]
            if overlap > 0:
                l_length = qrbounds[l_gene][1] - qrbounds[l_gene][0]
                r_length = qrbounds[r_gene][1] - qrbounds[r_gene][0]
                l_portion, r_portion = 0, 0
                while l_portion + r_portion < overlap:
                    if l_length <= 1 and r_length <= 1:  # don't want to erode match (in practice it'll be the d match) all the way to zero
                        print '      ERROR both lengths went to zero'
                        assert False
                    elif l_length > 1 and r_length > 1:  # if both have length left, alternate back and forth
                        if (l_portion + r_portion) % 2 == 0:
                            l_portion += 1  # give one base to the left
                            l_length -= 1
                        else:
                            r_portion += 1  # and one to the right
                            r_length -= 1
                    elif l_length > 1:
                        l_portion += 1
                        l_length -= 1
                    elif r_length > 1:
                        r_portion += 1
                        r_length -= 1

                if self.args.debug:
                    print '      WARNING %s apportioning %d bases between %s (%d) match and %s (%d) match' % (query_name, overlap, l_reg, l_portion, r_reg, r_portion)
                assert l_portion + r_portion == overlap
                qrbounds[l_gene] = (qrbounds[l_gene][0], qrbounds[l_gene][1] - l_portion)
                glbounds[l_gene] = (glbounds[l_gene][0], glbounds[l_gene][1] - l_portion)
                qrbounds[r_gene] = (qrbounds[r_gene][0] + r_portion, qrbounds[r_gene][1])
                glbounds[r_gene] = (glbounds[r_gene][0] + r_portion, glbounds[r_gene][1])

                best[l_reg + '_gl_seq'] = self.germline_seqs[l_reg][l_gene][glbounds[l_gene][0] : glbounds[l_gene][1]]
                best[l_reg + '_qr_seq'] = query_seq[qrbounds[l_gene][0]:qrbounds[l_gene][1]]
                best[r_reg + '_gl_seq'] = self.germline_seqs[r_reg][r_gene][glbounds[r_gene][0] : glbounds[r_gene][1]]
                best[r_reg + '_qr_seq'] = query_seq[qrbounds[r_gene][0]:qrbounds[r_gene][1]]

    # ----------------------------------------------------------------------------------------
    def add_to_info(self, query_name, query_seq, kvals, match_names, best, all_germline_bounds, all_query_bounds, codon_positions, perfplotter=None):
        assert query_name not in self.info
        self.info['queries'].append(query_name)
        self.info[query_name] = {}
        self.info[query_name]['unique_id'] = query_name  # redundant, but used somewhere down the line
        self.info[query_name]['k_v'] = kvals['v']
        self.info[query_name]['k_d'] = kvals['d']
        self.info[query_name]['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])

        # assert codon_positions['v'] != -1
        # assert codon_positions['j'] != -1
        self.info[query_name]['cdr3_length'] = codon_positions['j'] - codon_positions['v'] + 3  #tryp_position_in_joined_seq - self.cyst_position + 3
        self.info[query_name]['cyst_position'] = codon_positions['v']
        self.info[query_name]['tryp_position'] = codon_positions['j']

        # erosion, insertion, mutation info for best match
        self.info[query_name]['v_5p_del'] = all_germline_bounds[best['v']][0]
        self.info[query_name]['v_3p_del'] = len(self.germline_seqs['v'][best['v']]) - all_germline_bounds[best['v']][1]  # len(germline v) - gl_match_end
        self.info[query_name]['d_5p_del'] = all_germline_bounds[best['d']][0]
        self.info[query_name]['d_3p_del'] = len(self.germline_seqs['d'][best['d']]) - all_germline_bounds[best['d']][1]
        self.info[query_name]['j_5p_del'] = all_germline_bounds[best['j']][0]
        self.info[query_name]['j_3p_del'] = len(self.germline_seqs['j'][best['j']]) - all_germline_bounds[best['j']][1]

        self.info[query_name]['fv_insertion'] = query_seq[ : all_query_bounds[best['v']][0]]
        self.info[query_name]['vd_insertion'] = query_seq[all_query_bounds[best['v']][1] : all_query_bounds[best['d']][0]]
        self.info[query_name]['dj_insertion'] = query_seq[all_query_bounds[best['d']][1] : all_query_bounds[best['j']][0]]
        self.info[query_name]['jf_insertion'] = query_seq[all_query_bounds[best['j']][1] : ]

        for region in utils.regions:
            self.info[query_name][region + '_gene'] = best[region]
            self.info[query_name][region + '_gl_seq'] = best[region + '_gl_seq']
            self.info[query_name][region + '_qr_seq'] = best[region + '_qr_seq']
            self.info['all_best_matches'].add(best[region])

        self.info[query_name]['seq'] = query_seq  # only need to add this so I can pass it to print_reco_event
        if self.args.debug:
            if not self.args.is_data:
                utils.print_reco_event(self.germline_seqs, self.reco_info[query_name], extra_str='      ', label='true:')
            utils.print_reco_event(self.germline_seqs, self.info[query_name], extra_str='      ', label='inferred:')

        if self.pcounter != None:
            self.pcounter.increment_reco_params(self.info[query_name])
            self.pcounter.increment_mutation_params(self.info[query_name])
        if self.true_pcounter != None:
            self.true_pcounter.increment_reco_params(self.reco_info[query_name])
            self.true_pcounter.increment_mutation_params(self.reco_info[query_name])
        if perfplotter != None:
            perfplotter.evaluate(self.reco_info[query_name], self.info[query_name])  #, subtract_unphysical_erosions=True)

    # ----------------------------------------------------------------------------------------
    def summarize_query(self, query_name, query_seq, raw_best, all_match_names, all_query_bounds, all_germline_bounds, perfplotter, warnings, first_match_query_bounds):
        if self.args.debug:
            print '%s' % query_name

        # for region in utils.regions:
        #     if region not in raw_best:
        #         print '    skipping %s: no %s gene matches' % (query_name, region)
        #         return

        best, match_names, n_matches = {}, {}, {}
        n_used = {'v':0, 'd':0, 'j':0}
        k_v_min, k_d_min = 999, 999
        k_v_max, k_d_max = 0, 0
        for region in utils.regions:
            all_match_names[region] = sorted(all_match_names[region], reverse=True)
            match_names[region] = []
        codon_positions = {'v':-1, 'd':-1, 'j':-1}  # conserved codon positions (v:cysteine, d:dummy, j:tryptophan)
        for region in utils.regions:
            n_matches[region] = len(all_match_names[region])
            n_skipped = 0
            for score, gene in all_match_names[region]:
                glbounds = all_germline_bounds[gene]
                qrbounds = all_query_bounds[gene]
                assert qrbounds[1] <= len(query_seq)  # NOTE I'm putting these up avove as well (in process_query), so in time I should remove them from here
                assert glbounds[1] <= len(self.germline_seqs[region][gene])
                assert qrbounds[0] >= 0
                assert glbounds[0] >= 0
                glmatchseq = self.germline_seqs[region][gene][glbounds[0]:glbounds[1]]

                # TODO since I'm no longer skipping the genes after the first <args.n_max_per_region>, the OR of k-space below is overly conservative

                # only use a specified set of genes
                if self.args.only_genes != None and gene not in self.args.only_genes:
                    n_skipped += 1
                    continue

                # add match to the list
                n_used[region] += 1
                match_names[region].append(gene)

                self.print_match(region, gene, query_seq, score, glbounds, qrbounds, -1, warnings, skipping=False)

                # if the germline match and the query match aren't the same length, s-w likely added an insert, which we shouldn't get since the gap-open penalty is jacked up so high
                if len(glmatchseq) != len(query_seq[qrbounds[0]:qrbounds[1]]):  # neurotic double check (um, I think) EDIT hey this totally saved my ass
                    print 'ERROR %d not same length' % query_name
                    print glmatchseq, glbounds[0], glbounds[1]
                    print query_seq[qrbounds[0]:qrbounds[1]]
                    assert False

                if region == 'v':
                    this_k_v = all_query_bounds[gene][1]  # NOTE even if the v match doesn't start at the left hand edge of the query sequence, we still measure k_v from there.
                                                          # In other words, sw doesn't tell the hmm about it
                    k_v_min = min(this_k_v, k_v_min)
                    k_v_max = max(this_k_v, k_v_max)
                if region == 'd':
                    # this_k_d = all_query_bounds[gene][1] - all_query_bounds[raw_best['v']][1]  # end of d minus end of v
                    this_k_d = all_query_bounds[gene][1] - first_match_query_bounds[1]  # end of d minus end of v
                    k_d_min = min(this_k_d, k_d_min)
                    k_d_max = max(this_k_d, k_d_max)

                # check consistency with best match (since the best match is excised in s-w code, and because ham is run with *one* k_v k_d set)
                if region not in best:
                    best[region] = gene
                    best[region + '_gl_seq'] = self.germline_seqs[region][gene][glbounds[0]:glbounds[1]]
                    best[region + '_qr_seq'] = query_seq[qrbounds[0]:qrbounds[1]]
                    best[region + '_score'] = score

            if self.args.debug and n_skipped > 0:
                print '%8s skipped %d %s genes' % ('', n_skipped, region)

        for region in utils.regions:
            if region not in best:
                print '    no', region, 'match found for', query_name  # NOTE if no d match found, we should really should just assume entire d was eroded
                if not self.args.is_data:
                    print '    true:'
                    utils.print_reco_event(self.germline_seqs, self.reco_info[query_name], extra_str='    ')
                return

        # s-w allows d and j matches to overlap... which makes no sense, so arbitrarily give the disputed territory to j
        try:
            self.shift_overlapping_boundaries(all_query_bounds, all_germline_bounds, query_name, query_seq, best)
        except AssertionError:
            print '      ERROR %s apportionment failed' % query_name
            return

        # check for unproductive rearrangements
        for region in utils.regions:
            codon_positions[region] = utils.get_conserved_codon_position(self.cyst_positions, self.tryp_positions, region, best[region], all_germline_bounds[best[region]], all_query_bounds[best[region]], assert_on_fail=False)  # position in the query sequence, that is
        codons_ok = utils.check_both_conserved_codons(query_seq, codon_positions['v'], codon_positions['j'], debug=self.args.debug, extra_str='      ', assert_on_fail=False)
        cdr3_length = codon_positions['j'] - codon_positions['v'] + 3
        in_frame_cdr3 = (cdr3_length % 3 == 0)
        if self.args.debug and not in_frame_cdr3:
                print '      out of frame cdr3: %d %% 3 = %d' % (cdr3_length, cdr3_length % 3)
        no_stop_codon = utils.stop_codon_check(query_seq, codon_positions['v'], debug=self.args.debug)
        if not codons_ok or not in_frame_cdr3 or not no_stop_codon:
            if self.args.debug:
                print '       unproductive rearrangement in waterer codons_ok: %s   in_frame_cdr3: %s   no_stop_codon: %s' % (codons_ok, in_frame_cdr3, no_stop_codon)
            if self.args.skip_unproductive:
                if self.args.debug:
                    print '            ...skipping'
                self.n_unproductive += 1
                self.info['skipped_unproductive_queries'].append(query_name)
                return

        # best k_v, k_d:
        k_v = all_query_bounds[best['v']][1]  # end of v match
        k_d = all_query_bounds[best['d']][1] - all_query_bounds[best['v']][1]  # end of d minus end of v

        if k_d_max < 5:  # since the s-w step matches to the longest possible j and then excises it, this sometimes gobbles up the d, resulting in a very short d alignment.
            if self.args.debug:
                print '  expanding k_d'
            k_d_max = max(8, k_d_max)

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
        assert k_v_min > 0 and k_d_min > 0 and k_v_max > 0 and k_d_max > 0

        if self.args.debug:
            print '         k_v: %d [%d-%d)' % (k_v, k_v_min, k_v_max)
            print '         k_d: %d [%d-%d)' % (k_d, k_d_min, k_d_max)
            print '         used',
            for region in utils.regions:
                print ' %s: %d/%d' % (region, n_used[region], n_matches[region]),
            print ''


        kvals = {}
        kvals['v'] = {'best':k_v, 'min':k_v_min, 'max':k_v_max}
        kvals['d'] = {'best':k_d, 'min':k_d_min, 'max':k_d_max}
        self.add_to_info(query_name, query_seq, kvals, match_names, best, all_germline_bounds, all_query_bounds, codon_positions=codon_positions, perfplotter=perfplotter)
