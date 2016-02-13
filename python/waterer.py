import time
import sys
import math
import re
import os
import itertools
import operator
import pysam
import contextlib
from subprocess import check_output, Popen, PIPE
from collections import OrderedDict

import utils
from opener import opener
from parametercounter import ParameterCounter
from performanceplotter import PerformancePlotter

# ----------------------------------------------------------------------------------------
class Waterer(object):
    """ Run smith-waterman on the query sequences in <infname> """
    def __init__(self, args, input_info, reco_info, glfo, parameter_dir, write_parameters, genes_to_use):
        self.parameter_dir = parameter_dir
        self.args = args
        self.debug = self.args.debug if self.args.sw_debug is None else self.args.sw_debug

        self.max_insertion_length = 35  # if vdjalign reports an insertion longer than this, rerun the query (typically with different match/mismatch ratio)

        self.input_info = input_info
        self.remaining_queries = set()  # we remove queries from this set when we're satisfied with the current output (in general we may have to rerun some queries with different match/mismatch scores)
        for query in self.input_info.keys():
            self.remaining_queries.add(query)
        self.new_indels = 0  # number of new indels that were kicked up this time through

        self.genes_to_use = genes_to_use  # if None, we use all of 'em. NOTE do *not* use self.args.only_genes in this file (see partitiondriver)

        self.reco_info = reco_info
        self.glfo = glfo
        self.pcounter, self.true_pcounter, self.perfplotter = None, None, None
        if write_parameters:
            self.pcounter = ParameterCounter(self.glfo['seqs'])
            if not self.args.is_data:
                self.true_pcounter = ParameterCounter(self.glfo['seqs'])
        if self.args.plot_performance:
            self.perfplotter = PerformancePlotter(self.glfo['seqs'], 'sw')
        self.info = {}
        self.info['queries'] = []  # list of queries that *passed* sw, i.e. for which we have information
        self.info['all_best_matches'] = set()  # set of all the matches we found (for *all* queries)
        self.info['indels'] = {}
        if self.args.apply_choice_probs_in_sw:
            if self.debug:
                print '  reading gene choice probs from', parameter_dir
            self.gene_choice_probs = utils.read_overall_gene_probs(parameter_dir)

        self.outfile = None
        if self.args.outfname is not None:
            self.outfile = open(self.args.outfname, 'a')

        self.nth_try = 1
        self.unproductive_queries = set()

        print 'smith-waterman'

    # ----------------------------------------------------------------------------------------
    def __del__(self):
        if self.args.outfname is not None:
            self.outfile.close()

    # ----------------------------------------------------------------------------------------
    def clean(self):
        if self.pcounter is not None:
            self.pcounter.clean()
        if self.true_pcounter is not None:
            self.true_pcounter.clean()

    # ----------------------------------------------------------------------------------------
    def run(self):
        # start = time.time()
        base_infname = 'query-seqs.fa'
        base_outfname = 'query-seqs.bam'
        sys.stdout.flush()

        while len(self.remaining_queries) > 0:  # we remove queries from <self.remaining_queries> as we're satisfied with their output
            self.write_vdjalign_input(base_infname, n_procs=self.args.n_fewer_procs)
            self.execute_commands(base_infname, base_outfname, self.args.n_fewer_procs)
            self.read_output(base_outfname, n_procs=self.args.n_fewer_procs)
            if self.nth_try > 3:
                break
            self.nth_try += 1  # it's set to 1 before we begin the first try, and increases to 2 just before we start the second try

        self.finalize()

    # ----------------------------------------------------------------------------------------
    def finalize(self):
        if self.perfplotter is not None:
            self.perfplotter.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots)
        # print '    sw time: %.3f' % (time.time()-start)
        print '      info for %d' % len(self.info['queries']),
        skipped_unproductive = len(self.unproductive_queries)
        n_remaining = len(self.remaining_queries)
        if skipped_unproductive > 0 or n_remaining > 0:
            print '     (skipped',
            print '%d / %d = %.3f unproductive' % (skipped_unproductive, len(self.input_info), float(skipped_unproductive) / len(self.input_info)),
            if n_remaining > 0:
                print '   %d / %d = %.3f other' % (n_remaining, len(self.input_info), float(n_remaining) / len(self.input_info)),
            print ')',
        print ''
        if n_remaining > 0:
            printstr = '   %s %d missing annotations' % (utils.color('red', 'warning'), n_remaining)
            if n_remaining < 15:
                printstr += ' (' + ':'.join(self.remaining_queries) + ')'
            print printstr
        if self.debug and len(self.info['indels']) > 0:
            print '      indels: %s' % ':'.join(self.info['indels'].keys())
        assert len(self.info['queries']) + skipped_unproductive + n_remaining == len(self.input_info)
        if self.debug and n_remaining > 0:
            print 'true annotations for remaining events:'
            for qry in self.remaining_queries:
                utils.print_reco_event(self.glfo['seqs'], self.reco_info[qry], extra_str='      ', label='true:')  #, indelfo=self.reco_info[query_name]['indels'])
        if self.pcounter is not None:
            self.pcounter.write(self.parameter_dir)
            if self.true_pcounter is not None:
                assert self.parameter_dir[-1] != '/'
                self.true_pcounter.write(self.parameter_dir + '-true')
            if self.args.plotdir is not None:
                self.pcounter.plot(self.args.plotdir + '/sw', subset_by_gene=True, cyst_positions=self.glfo['cyst-positions'], tryp_positions=self.glfo['tryp-positions'], only_csv=self.args.only_csv_plots)
                if self.true_pcounter is not None:
                    self.true_pcounter.plot(self.args.plotdir + '/sw-true', subset_by_gene=True, cyst_positions=self.glfo['cyst-positions'], tryp_positions=self.glfo['tryp-positions'], only_csv=self.args.only_csv_plots)

        utils.pad_seqs_to_same_length(self.info['queries'], self.info, self.glfo, self.info['indels'])  # adds padded info to self.info (returns if stuff has already been padded)
        self.info['remaining_queries'] = self.remaining_queries

    # ----------------------------------------------------------------------------------------
    def execute_commands(self, base_infname, base_outfname, n_procs):
        # rewrite input germline sets
        datadir = self.args.datadir
        if self.genes_to_use is not None:
            datadir = self.args.workdir + '/germline-sets'
            rewritten_files = utils.rewrite_germline_fasta(self.args.datadir, datadir, self.genes_to_use)

        if n_procs == 1:
            cmd_str = self.get_vdjalign_cmd_str(self.args.workdir, base_infname, base_outfname, datadir)
            proc = Popen(cmd_str.split(), stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate()
            utils.process_out_err(out, err)
            if not self.args.no_clean:
                os.remove(self.args.workdir + '/' + base_infname)
        else:
            shell = True

            cmd_strs, workdirs = [], []
            for iproc in range(n_procs):
                workdirs.append(self.args.workdir + '/sw-' + str(iproc))
                cmd_strs.append(self.get_vdjalign_cmd_str(workdirs[iproc], base_infname, base_outfname, datadir, n_procs, shell))

            # start all procs for the first time
            procs, n_tries = [], []
            for iproc in range(n_procs):
                procs.append(Popen(cmd_strs[iproc], shell=shell))
                n_tries.append(1)
                time.sleep(0.1)

            # ----------------------------------------------------------------------------------------
            def get_outfname(iproc):
                return workdirs[iproc] + '/' + base_outfname

            # ----------------------------------------------------------------------------------------
            # deal with a process once it's finished (i.e. check if it failed, and restart if so)
            def finish_process(iproc):
                procs[iproc].communicate()
                utils.process_out_err('', '', extra_str=str(iproc), subworkdir=workdirs[iproc])
                if procs[iproc].returncode == 0 and os.path.exists(get_outfname(iproc)):
                    procs[iproc] = None  # job succeeded
                elif n_tries[iproc] > 5:
                    raise Exception('exceeded max number of tries for command\n    %s\nlook for output in %s' % (cmd_strs[iproc], workdirs[iproc]))
                else:
                    print '    rerunning proc %d (exited with %d' % (iproc, procs[iproc].returncode),
                    if not os.path.exists(get_outfname(iproc)):
                        print ', output %s d.n.e.' % get_outfname(iproc),
                    print ')'
                    procs[iproc] = Popen(cmd_strs[iproc], shell=shell)
                    n_tries[iproc] += 1

            # keep looping over the procs until they're all done
            while procs.count(None) != len(procs):  # we set each proc to None when it finishes
                for iproc in range(n_procs):
                    if procs[iproc] is None:  # already finished
                        continue
                    # read_progress(iproc)
                    if procs[iproc].poll() is not None:  # it's finished
                        finish_process(iproc)
                sys.stdout.flush()
                time.sleep(1)

            if not self.args.no_clean:
                for iproc in range(n_procs):
                    os.remove(workdirs[iproc] + '/' + base_infname)

        if self.genes_to_use is not None:  # this'll rewrite them for each time through, which is ok and probably preferable
            for fname in rewritten_files:
                os.remove(fname)
            os.rmdir(datadir)

        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def write_vdjalign_input(self, base_infname, n_procs):
        n_remaining = len(self.remaining_queries)
        queries_per_proc = float(n_remaining) / n_procs
        n_queries_per_proc = int(math.ceil(queries_per_proc))
        if n_procs == 1:  # double check for rounding problems or whatnot
            assert n_queries_per_proc == n_remaining
        for iproc in range(n_procs):
            workdir = self.args.workdir
            if n_procs > 1:
                workdir += '/sw-' + str(iproc)
                utils.prep_dir(workdir)
            with opener('w')(workdir + '/' + base_infname) as sub_infile:
                iquery = 0
                for query_name in self.remaining_queries:  # NOTE this is wasteful to loop of all the remaining queries for each process... but maybe not that wasteful
                    if iquery >= n_remaining:
                        break
                    if iquery < iproc*n_queries_per_proc or iquery >= (iproc + 1)*n_queries_per_proc:  # not for this process
                        iquery += 1
                        continue
                    sub_infile.write('>' + query_name + ' NUKES\n')

                    seq = self.input_info[query_name]['seq']
                    if query_name in self.info['indels']:
                        seq = self.info['indels'][query_name]['reversed_seq']  # use the query sequence with shm insertions and deletions reversed
                    sub_infile.write(seq + '\n')
                    iquery += 1

    # ----------------------------------------------------------------------------------------
    def get_vdjalign_cmd_str(self, workdir, base_infname, base_outfname, datadir, n_procs=None, shell=False):
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
        if self.args.slurm or utils.auto_slurm(n_procs):
            cmd_str = 'srun ' + cmd_str
        cmd_str += ' --max-drop 50'
        match, mismatch = self.args.match_mismatch
        cmd_str += ' --match ' + str(match) + ' --mismatch ' + str(mismatch)
        cmd_str += ' --gap-open ' + str(self.args.gap_open_penalty)  #1000'  #50'
        cmd_str += ' --vdj-dir ' + datadir  # NOTE not necessarily <self.args.datadir>
        cmd_str += ' ' + workdir + '/' + base_infname + ' ' + workdir + '/' + base_outfname

        if shell:
            cmd_str += ' 1>' + workdir + '/out' + ' 2>' + workdir + '/err'

        return cmd_str

    # ----------------------------------------------------------------------------------------
    def read_output(self, base_outfname, n_procs=1):
        queries_to_rerun = OrderedDict()  # This is to keep track of every query that we don't add to self.info (i.e. it does *not* include unproductive queries that we ignore/skip entirely because we were told to by a command line argument)
                                          # ...whereas <self.unproductive_queries> is to keep track of the queries that were definitively unproductive (i.e. we removed them from self.remaining_queries) when we were told to skip unproductives by a command line argument
        for reason in ['unproductive', 'no-match', 'weird-annot.', 'nonsense-bounds', 'invalid-codon']:
            queries_to_rerun[reason] = set()

        self.new_indels = 0
        n_processed = 0
        for iproc in range(n_procs):
            workdir = self.args.workdir
            if n_procs > 1:
                workdir += '/sw-' + str(iproc)
            outfname = workdir + '/' + base_outfname
            with contextlib.closing(pysam.Samfile(outfname)) as bam:
                grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
                for _, reads in grouped:  # loop over query sequences
                    self.process_query(bam.references, list(reads), queries_to_rerun)
                    n_processed += 1

            if not self.args.no_clean:
                os.remove(outfname)
                if n_procs > 1:  # still need the top-level workdir
                    os.rmdir(workdir)

        if self.nth_try == 1:
            print '        processed       remaining      new-indels          rerun: ' + '      '.join([reason for reason in queries_to_rerun])
        print '      %8d' % n_processed,
        if len(self.remaining_queries) > 0:
            printstr = '       %8d' % len(self.remaining_queries)
            printstr += '       %8d' % self.new_indels
            printstr += '            '
            n_to_rerun = 0
            for reason in queries_to_rerun:
                printstr += '        %8d' % len(queries_to_rerun[reason])
                n_to_rerun += len(queries_to_rerun[reason])
            print printstr,
            if n_to_rerun + self.new_indels != len(self.remaining_queries):
                print n_to_rerun, self.new_indels, len(self.remaining_queries)
                raise Exception('I\'m an exception!')
            if self.nth_try < 2 or self.new_indels == 0:  # increase the mismatch score if it's the first try, or if there's no new indels
                print '            increasing mismatch score (%d --> %d) and rerunning them' % (self.args.match_mismatch[1], self.args.match_mismatch[1] + 1)
                self.args.match_mismatch[1] += 1
            elif self.new_indels > 0:  # if there were some indels, rerun with the same parameters (but when the input is written the indel will be "reversed' in the sequences that's passed to ighutil)
                print '            rerunning for indels'
                self.new_indels = 0
            else:  # shouldn't get here
                assert False
        else:
            print '        all done'

    # ----------------------------------------------------------------------------------------
    def get_choice_prob(self, region, gene):
        choice_prob = 1.0
        if gene in self.gene_choice_probs[region]:
            choice_prob = self.gene_choice_probs[region][gene]
        else:
            choice_prob = 0.0  # NOTE would it make sense to use something else here?
        return choice_prob

    # ----------------------------------------------------------------------------------------
    def get_indel_info(self, query_name, cigarstr, qrseq, glseq, gene):
        cigars = re.findall('[0-9][0-9]*[A-Z]', cigarstr)  # split cigar string into its parts
        cigars = [(cstr[-1], int(cstr[:-1])) for cstr in cigars]  # split each part into the code and the length

        codestr = ''
        qpos = 0  # position within query sequence
        indelfo = {'reversed_seq' : '', 'indels' : []}  # replacement_seq: query seq with insertions removed and germline bases inserted at the position of deletions
        tmp_indices = []
        for code, length in cigars:
            codestr += length * code
            if code == 'I':  # advance qr seq but not gl seq
                indelfo['indels'].append({'type' : 'insertion', 'pos' : qpos, 'len' : length, 'seqstr' : ''})  # insertion begins at <pos>
                tmp_indices += [len(indelfo['indels']) - 1  for _ in range(length)]# indel index corresponding to this position in the alignment
            elif code == 'D':  # advance qr seq but not gl seq
                indelfo['indels'].append({'type' : 'deletion', 'pos' : qpos, 'len' : length, 'seqstr' : ''})  # first deleted base is <pos> (well, first base which is in the position of the first deleted base)
                tmp_indices += [len(indelfo['indels']) - 1  for _ in range(length)]# indel index corresponding to this position in the alignment
            else:
                tmp_indices += [None  for _ in range(length)]  # indel index corresponding to this position in the alignment
            qpos += length

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
                indelfo['reversed_seq'] += qrseq[iqr]  # add the base to the overall sequence with all indels reversed
            elif code == 'S':
                continue
            elif code == 'I':
                qrprintstr += utils.color('light_blue', qrseq[iqr])
                glprintstr += utils.color('light_blue', '*')
                indelfo['indels'][tmp_indices[icode]]['seqstr'] += qrseq[iqr]  # and to the sequence of just this indel
                igl -= 1
            elif code == 'D':
                qrprintstr += utils.color('light_blue', '*')
                glprintstr += utils.color('light_blue', glseq[igl])
                indelfo['reversed_seq'] += glseq[igl]  # add the base to the overall sequence with all indels reversed
                indelfo['indels'][tmp_indices[icode]]['seqstr'] += glseq[igl]  # and to the sequence of just this indel
                iqr -= 1
            else:
                raise Exception('unhandled code %s' % code)

            iqr += 1
            igl += 1

        if self.debug:
            print '\n      indels in %s' % query_name
            print '          %20s %s' % (gene, glprintstr)
            print '          %20s %s' % ('query', qrprintstr)
            for idl in indelfo['indels']:
                print '          %10s: %d bases at %d (%s)' % (idl['type'], idl['len'], idl['pos'], idl['seqstr'])
        # utils.undo_indels(indelfo)
        # print '                       %s' % self.input_info[query_name]['seq']

        return indelfo

    # ----------------------------------------------------------------------------------------
    def process_query(self, references, reads, queries_to_rerun):
        primary = next((r for r in reads if not r.is_secondary), None)
        query_seq = primary.seq
        query_name = primary.qname
        first_match_query_bounds = None  # since sw excises its favorite v match, we have to know this match's boundaries in order to calculate k_d for all the other matches
        all_match_names = {}
        warnings = {}  # ick, this is a messy way to pass stuff around
        for region in utils.regions:
            all_match_names[region] = []
        all_query_bounds, all_germline_bounds = {}, {}
        for read in reads:  # loop over the matches found for each query sequence
            # set this match's values
            read.seq = query_seq  # only the first one has read.seq set by default, so we need to set the rest by hand
            gene = references[read.tid]
            region = utils.get_region(gene)
            raw_score = read.tags[0][1]  # raw because they don't include the gene choice probs
            score = raw_score
            if self.args.apply_choice_probs_in_sw:  # NOTE I stopped applying the gene choice probs here because the smith-waterman scores don't correspond to log-probs, so throwing on the gene choice probs was dubious (and didn't seem to work that well)
                score = self.get_choice_prob(region, gene) * raw_score  # multiply by the probability to choose this gene
            qrbounds = (read.qstart, read.qend)
            glbounds = (read.pos, read.aend)
            if region == 'v' and first_match_query_bounds is None:
                first_match_query_bounds = qrbounds

            # perform a few checks and see if we want to skip this match
            # TODO I wish this wasn't here and I suspect I don't really need it (any more)
            if region == 'v':  # skip matches with cpos past the end of the query seq (i.e. eroded a ton on the right side of the v)
                cpos = self.glfo['cyst-positions'][gene] - glbounds[0] + qrbounds[0]  # position within original germline gene, minus the position in that germline gene at which the match starts, plus the position in the query sequence at which the match starts
                if cpos < 0 or cpos >= len(query_seq):
                    continue

            if 'I' in read.cigarstring or 'D' in read.cigarstring:  # skip indels, and tell the HMM to skip indels (you won't see any unless you decrease the <self.args.gap_open_penalty>)
                if self.args.no_indels:  # you can forbid indels on the command line
                    continue
                if self.nth_try < 2:  # we also forbid indels on the first try (we want to increase the mismatch score before we conclude it's "really" an indel)
                    continue
                if len(all_match_names[region]) == 0:  # if this is the first (best) match for this region, allow indels (otherwise skip the match)
                    if query_name not in self.info['indels']:
                        self.info['indels'][query_name] = self.get_indel_info(query_name, read.cigarstring, query_seq[qrbounds[0] : qrbounds[1]], self.glfo['seqs'][region][gene][glbounds[0] : glbounds[1]], gene)
                        self.info['indels'][query_name]['reversed_seq'] = query_seq[ : qrbounds[0]] + self.info['indels'][query_name]['reversed_seq'] + query_seq[qrbounds[1] : ]
                        self.new_indels += 1
                        # TODO this 'return' used to be after and indented from the else below, and that continue wasn't there. I should make sure this is how I want it
                        return  # don't process this query any further -- since it's now in the indel info it'll get run next time through
                    else:
                        if self.debug:
                            print '     ignoring subsequent indels for %s' % query_name
                        continue  # hopefully there's a later match without indels
                else:
                    continue

            if qrbounds[1]-qrbounds[0] != glbounds[1]-glbounds[0]:
                raise Exception('germline match (%d %d) not same length as query match (%d %d)' % (qrbounds[0], qrbounds[1], glbounds[0], glbounds[1]))

            assert qrbounds[1] <= len(query_seq)
            if glbounds[1] > len(self.glfo['seqs'][region][gene]):
                print '  ', gene
                print '  ', glbounds[1], len(self.glfo['seqs'][region][gene])
                print '  ', self.glfo['seqs'][region][gene]
            assert glbounds[1] <= len(self.glfo['seqs'][region][gene])
            assert qrbounds[1]-qrbounds[0] == glbounds[1]-glbounds[0]

            # and finally add this match's information
            warnings[gene] = ''
            all_match_names[region].append((score, gene))  # NOTE it is important that this is ordered such that the best match is first
            all_query_bounds[gene] = qrbounds
            all_germline_bounds[gene] = glbounds

        self.summarize_query(query_name, query_seq, all_match_names, all_query_bounds, all_germline_bounds, warnings, first_match_query_bounds, queries_to_rerun)

    # ----------------------------------------------------------------------------------------
    def print_match(self, region, gene, query_seq, score, glbounds, qrbounds, codon_pos, warnings, skipping=False):
        out_str_list = []
        buff_str = (20 - len(gene)) * ' '
        tmp_val = score
        if self.args.apply_choice_probs_in_sw and self.get_choice_prob(region, gene) != 0.0:
            tmp_val = score / self.get_choice_prob(region, gene)
        if self.args.apply_choice_probs_in_sw:
            out_str_list.append('%8s%s%s%9.1e * %3.0f = %-6.1f' % (' ', utils.color_gene(gene), buff_str, self.get_choice_prob(region, gene), tmp_val, score))
        else:
            out_str_list.append('%8s%s%s%9s%3s %6.0f        ' % (' ', utils.color_gene(gene), '', '', buff_str, score))
        out_str_list.append('%4d%4d   %s\n' % (glbounds[0], glbounds[1], self.glfo['seqs'][region][gene][glbounds[0]:glbounds[1]]))
        out_str_list.append('%46s  %4d%4d' % ('', qrbounds[0], qrbounds[1]))
        out_str_list.append('   %s ' % (utils.color_mutants(self.glfo['seqs'][region][gene][glbounds[0]:glbounds[1]], query_seq[qrbounds[0]:qrbounds[1]])))
        if region != 'd':
            out_str_list.append('(%s %d)' % (utils.conserved_codon_names[region], codon_pos))
        if warnings[gene] != '':
            out_str_list.append('WARNING ' + warnings[gene])
        if skipping:
            out_str_list.append('skipping!')
        if self.args.outfname is None:
            print ''.join(out_str_list)
        else:
            out_str_list.append('\n')
            self.outfile.write(''.join(out_str_list))

    # ----------------------------------------------------------------------------------------
    def get_overlap_and_available_space(self, rpair, best, qrbounds):
        l_reg = rpair['left']
        r_reg = rpair['right']
        l_gene = best[l_reg]
        r_gene = best[r_reg]
        overlap = qrbounds[l_gene][1] - qrbounds[r_gene][0]
        available_space = qrbounds[r_gene][1] - qrbounds[l_gene][0]
        return overlap, available_space

    # ----------------------------------------------------------------------------------------
    def check_boundaries(self, rpair, qrbounds, glbounds, query_name, query_seq, best, debug=False):
        # NOTE this duplicates code in shift_overlapping_boundaries(), which makes me cranky, but this setup avoids other things I dislike more
        l_reg = rpair['left']
        r_reg = rpair['right']
        l_gene = best[l_reg]
        r_gene = best[r_reg]

        overlap, available_space = self.get_overlap_and_available_space(rpair, best, qrbounds)

        if debug:
            print '  %s %s    overlap %d    available space %d' % (l_reg, r_reg, overlap, available_space)

        status = 'ok'
        if overlap > 0:  # positive overlap means they actually overlap
            status = 'overlap'
        if overlap > available_space:  # call it nonsense if the boundaries are really whack (i.e. there isn't enough space to resolve the overlap) -- we'll presumably either toss the query or rerun with different match/mismatch
            status = 'nonsense'

        if debug:
            print '  overlap status: %s' % status

        if status == 'nonsense' and l_reg == 'd' and self.nth_try > 2:  # on rare occasions with very high mutation, vdjalign refuses to give us a j match that's at all to the right of the d match
            assert l_reg == 'd' and r_reg == 'j'
            if debug:
                print '  %s: synthesizing d match' % query_name
            leftmost_position = min(qrbounds[l_gene][0], qrbounds[r_gene][0])
            qrbounds[l_gene] = (leftmost_position, leftmost_position + 1)  # swap whatever crummy nonsense d match we have now for a one-base match at the left end of things (things in practice should be left end of j match)
            glbounds[l_gene] = (0, 1)
            status = self.check_boundaries(rpair, qrbounds, glbounds, query_name, query_seq, best, debug)

        return status

    # ----------------------------------------------------------------------------------------
    def shift_overlapping_boundaries(self, rpair, qrbounds, glbounds, query_name, query_seq, best, debug=False):
        # NOTE this does pretty much the same thing as resolve_overlapping_matches in joinparser.py
        """
        s-w allows d and j matches (and v and d matches) to overlap... which makes no sense, so apportion the disputed territory between the two regions.
        Note that this still works if, say, v is the entire sequence, i.e. one match is entirely subsumed by another.
        """
        l_reg = rpair['left']
        r_reg = rpair['right']
        l_gene = best[l_reg]
        r_gene = best[r_reg]

        overlap, available_space = self.get_overlap_and_available_space(rpair, best, qrbounds)

        if overlap <= 0:  # nothing to do, they're already consistent
            print 'shouldn\'t get here any more if there\'s no overlap'
            return

        if overlap > available_space:
            raise Exception('overlap %d bigger than available space %d between %s and %s for %s' % (overlap, available_space, l_reg, r_reg, query_name))

        if debug:
            print '%s%s:  %d-%d overlaps with %d-%d by %d' % (l_reg, r_reg, qrbounds[l_gene][0], qrbounds[l_gene][1], qrbounds[r_gene][0], qrbounds[r_gene][1], overlap)

        l_length = qrbounds[l_gene][1] - qrbounds[l_gene][0]  # initial length of lefthand gene match
        r_length = qrbounds[r_gene][1] - qrbounds[r_gene][0]  # and same for the righthand one
        l_portion, r_portion = 0, 0  # portion of the initial overlap that we give to each side
        if debug:
            print '    lengths        portions     '
        while l_portion + r_portion < overlap:
            if debug:
                print '  %4d %4d      %4d %4d      %s %s' % (l_length, r_length, l_portion, r_portion, '', '')
            if l_length <= 1 and r_length <= 1:  # don't want to erode match (in practice it'll be the d match) all the way to zero
                raise Exception('both lengths went to zero for %s: %s %s' % (query_name, qrbounds[l_gene], qrbounds[r_gene]))
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

        if debug:
            print '  %4d %4d    %4d %4d      %s %s' % (l_length, r_length, l_portion, r_portion, '', '')
            print '      %s apportioning %d bases between %s (%d) match and %s (%d) match' % (query_name, overlap, l_reg, l_portion, r_reg, r_portion)
        assert l_portion + r_portion == overlap
        qrbounds[l_gene] = (qrbounds[l_gene][0], qrbounds[l_gene][1] - l_portion)
        glbounds[l_gene] = (glbounds[l_gene][0], glbounds[l_gene][1] - l_portion)
        qrbounds[r_gene] = (qrbounds[r_gene][0] + r_portion, qrbounds[r_gene][1])
        glbounds[r_gene] = (glbounds[r_gene][0] + r_portion, glbounds[r_gene][1])

        best[l_reg + '_gl_seq'] = self.glfo['seqs'][l_reg][l_gene][glbounds[l_gene][0] : glbounds[l_gene][1]]
        best[l_reg + '_qr_seq'] = query_seq[qrbounds[l_gene][0]:qrbounds[l_gene][1]]
        best[r_reg + '_gl_seq'] = self.glfo['seqs'][r_reg][r_gene][glbounds[r_gene][0] : glbounds[r_gene][1]]
        best[r_reg + '_qr_seq'] = query_seq[qrbounds[r_gene][0]:qrbounds[r_gene][1]]

    # ----------------------------------------------------------------------------------------
    def add_to_info(self, query_name, query_seq, kvals, match_names, best, all_germline_bounds, all_query_bounds, codon_positions):
        assert query_name not in self.info
        self.info['queries'].append(query_name)
        self.info[query_name] = {}
        self.info[query_name]['unique_id'] = query_name  # redundant, but used somewhere down the line
        self.info[query_name]['k_v'] = kvals['v']
        self.info[query_name]['k_d'] = kvals['d']
        self.info[query_name]['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])

        self.info[query_name]['cdr3_length'] = codon_positions['j'] - codon_positions['v'] + 3  #tryp_position_in_joined_seq - self.cyst_position + 3
        self.info[query_name]['cyst_position'] = codon_positions['v']
        self.info[query_name]['tryp_position'] = codon_positions['j']

        # erosion, insertion, mutation info for best match
        self.info[query_name]['v_5p_del'] = all_germline_bounds[best['v']][0]
        self.info[query_name]['v_3p_del'] = len(self.glfo['seqs']['v'][best['v']]) - all_germline_bounds[best['v']][1]  # len(germline v) - gl_match_end
        self.info[query_name]['d_5p_del'] = all_germline_bounds[best['d']][0]
        self.info[query_name]['d_3p_del'] = len(self.glfo['seqs']['d'][best['d']]) - all_germline_bounds[best['d']][1]
        self.info[query_name]['j_5p_del'] = all_germline_bounds[best['j']][0]
        self.info[query_name]['j_3p_del'] = len(self.glfo['seqs']['j'][best['j']]) - all_germline_bounds[best['j']][1]

        self.info[query_name]['fv_insertion'] = query_seq[ : all_query_bounds[best['v']][0]]
        self.info[query_name]['vd_insertion'] = query_seq[all_query_bounds[best['v']][1] : all_query_bounds[best['d']][0]]
        self.info[query_name]['dj_insertion'] = query_seq[all_query_bounds[best['d']][1] : all_query_bounds[best['j']][0]]
        self.info[query_name]['jf_insertion'] = query_seq[all_query_bounds[best['j']][1] : ]

        for region in utils.regions:
            self.info[query_name][region + '_gene'] = best[region]
            self.info[query_name][region + '_gl_seq'] = best[region + '_gl_seq']
            self.info[query_name][region + '_qr_seq'] = best[region + '_qr_seq']
            self.info['all_best_matches'].add(best[region])

        self.info[query_name]['seq'] = query_seq  # NOTE this is the seq output by vdjalign, i.e. if we reversed any indels it is the reversed sequence

        utils.add_implicit_info(self.glfo, self.info[query_name])

        if self.debug:
            if not self.args.is_data:
                utils.print_reco_event(self.glfo['seqs'], self.reco_info[query_name], extra_str='      ', label='true:', indelfo=self.reco_info[query_name]['indels'])
            utils.print_reco_event(self.glfo['seqs'], self.info[query_name], extra_str='      ', label='inferred:', indelfo=self.info['indels'].get(query_name, None))

        if self.pcounter is not None:
            self.pcounter.increment_per_family_params(self.info[query_name])
            self.pcounter.increment_per_sequence_params(self.info[query_name])
        if self.true_pcounter is not None:
            self.true_pcounter.increment_per_family_params(self.reco_info[query_name])
            self.true_pcounter.increment_per_sequence_params(self.reco_info[query_name])
        if self.perfplotter is not None:
            if query_name in self.info['indels']:
                print '    skipping performance evaluation of %s because of indels' % query_name  # I just have no idea how to handle naive hamming fraction when there's indels
            else:
                self.perfplotter.evaluate(self.reco_info[query_name], self.info[query_name])

        self.remaining_queries.remove(query_name)

    # ----------------------------------------------------------------------------------------
    def summarize_query(self, query_name, query_seq, all_match_names, all_query_bounds, all_germline_bounds, warnings, first_match_query_bounds, queries_to_rerun):
        if self.debug:
            print query_name

        best, match_names = {}, {}
        k_v_min, k_d_min = 999, 999
        k_v_max, k_d_max = 0, 0
        for region in utils.regions:
            all_match_names[region] = sorted(all_match_names[region], reverse=True)
            match_names[region] = []
        for region in utils.regions:
            for score, gene in all_match_names[region]:
                glbounds = all_germline_bounds[gene]
                qrbounds = all_query_bounds[gene]
                assert qrbounds[1] <= len(query_seq)  # NOTE I'm putting these up above as well (in process_query), so in time I should remove them from here
                assert glbounds[1] <= len(self.glfo['seqs'][region][gene])
                assert qrbounds[0] >= 0
                assert glbounds[0] >= 0
                glmatchseq = self.glfo['seqs'][region][gene][glbounds[0]:glbounds[1]]

                # TODO remove this
                if self.genes_to_use is not None and gene not in self.genes_to_use:
                    raise Exception('wtf %s not in %s' % (gene, self.genes_to_use))

                match_names[region].append(gene)

                if self.debug >= 2:
                    self.print_match(region, gene, query_seq, score, glbounds, qrbounds, -1, warnings, skipping=False)

                # if the germline match and the query match aren't the same length, s-w likely added an insert, which we shouldn't get since the gap-open penalty is jacked up so high
                if len(glmatchseq) != len(query_seq[qrbounds[0]:qrbounds[1]]):  # neurotic double check (um, I think) EDIT hey this totally saved my ass
                    print 'ERROR %d not same length' % query_name
                    print glmatchseq, glbounds[0], glbounds[1]
                    print query_seq[qrbounds[0]:qrbounds[1]]
                    assert False

                # NOTE since I'm no longer skipping the genes after the first <args.n_max_per_region>, the OR of k-space below is overly conservative. UPDATE not sure if this is still relevant, but I'll move it down here in case I feel like thinking about it later
                if region == 'v':
                    this_k_v = all_query_bounds[gene][1]  # NOTE even if the v match doesn't start at the left hand edge of the query sequence, we still measure k_v from there.
                                                          # In other words, sw doesn't tell the hmm about it
                    k_v_min = min(this_k_v, k_v_min)
                    k_v_max = max(this_k_v, k_v_max)
                if region == 'd':
                    this_k_d = all_query_bounds[gene][1] - first_match_query_bounds[1]  # end of d minus end of v
                    k_d_min = min(this_k_d, k_d_min)
                    k_d_max = max(this_k_d, k_d_max)

                # check consistency with best match (since the best match is excised in s-w code, and because ham is run with *one* k_v k_d set)
                if region not in best:
                    best[region] = gene
                    best[region + '_gl_seq'] = self.glfo['seqs'][region][gene][glbounds[0]:glbounds[1]]
                    best[region + '_qr_seq'] = query_seq[qrbounds[0]:qrbounds[1]]
                    best[region + '_score'] = score

        for region in utils.regions:
            if region not in best:
                if self.debug:
                    print '      no', region, 'match found for', query_name  # NOTE if no d match found, we should really just assume entire d was eroded
                queries_to_rerun['no-match'].add(query_name)
                return

        # s-w allows d and j matches to overlap, so we need to apportion the disputed bases
        region_pairs = ({'left':'v', 'right':'d'}, {'left':'d', 'right':'j'})
        for rpair in region_pairs:
            overlap_status = self.check_boundaries(rpair, all_query_bounds, all_germline_bounds, query_name, query_seq, best)
            if overlap_status == 'overlap':
                self.shift_overlapping_boundaries(rpair, all_query_bounds, all_germline_bounds, query_name, query_seq, best)
            elif overlap_status == 'nonsense':
                queries_to_rerun['nonsense-bounds'].add(query_name)
                return
            else:
                assert overlap_status == 'ok'

        # check for suspiciously bad annotations
        vd_insertion = query_seq[all_query_bounds[best['v']][1] : all_query_bounds[best['d']][0]]
        dj_insertion = query_seq[all_query_bounds[best['d']][1] : all_query_bounds[best['j']][0]]
        if self.nth_try < 2:
            if len(vd_insertion) > self.max_insertion_length or len(dj_insertion) > self.max_insertion_length:
                if self.debug:
                    print '      suspiciously long insertion, rerunning'
                queries_to_rerun['weird-annot.'].add(query_name)
                return

        # set and check conserved codon positions
        tmp_gl_positions = {'v' : self.glfo['cyst-positions'], 'j' : self.glfo['tryp-positions']}  # hack hack hack
        codon_positions = {}
        for region in ['v', 'j']:
            pos = tmp_gl_positions[region][best[region]] - all_germline_bounds[best[region]][0] + all_query_bounds[best[region]][0]  # position within original germline gene, minus the position in that germline gene at which the match starts, plus the position in the query sequence at which the match starts
            if pos < 0 or pos >= len(query_seq):
                if self.debug:
                    print '      invalid %s codon position (%d in seq of length %d), rerunning' % (region, pos, len(query_seq))
                queries_to_rerun['invalid-codon'].add(query_name)
                return
            codon_positions[region] = pos

        # check for unproductive rearrangements
        codons_ok = utils.check_both_conserved_codons(query_seq, codon_positions['v'], codon_positions['j'], debug=self.debug, extra_str='      ', assert_on_fail=False)
        cdr3_length = codon_positions['j'] - codon_positions['v'] + 3
        in_frame_cdr3 = (cdr3_length % 3 == 0)
        if self.debug and not in_frame_cdr3:
            print '      out of frame cdr3: %d %% 3 = %d' % (cdr3_length, cdr3_length % 3)
        no_stop_codon = utils.stop_codon_check(query_seq, codon_positions['v'], debug=self.debug)
        if not codons_ok or not in_frame_cdr3 or not no_stop_codon:
            if self.debug:
                print '       unproductive rearrangement:',
                if not codons_ok:
                    print '  bad codons',
                if not in_frame_cdr3:
                    print '  out of frame cdr3',
                if not no_stop_codon:
                    print '  stop codon'
                print ''

            if self.nth_try < 2 and (not codons_ok or not in_frame_cdr3):  # rerun with higher mismatch score (sometimes unproductiveness is the result of a really screwed up annotation rather than an actual unproductive sequence). Note that stop codons aren't really indicative of screwed up annotations, so they don't count.
                if self.debug:
                    print '            ...rerunning'
                queries_to_rerun['unproductive'].add(query_name)
                return
            elif self.args.skip_unproductive:
                if self.debug:
                    print '            ...skipping'
                self.unproductive_queries.add(query_name)
                self.remaining_queries.remove(query_name)
                return
            else:
                pass  # this is here so you don't forget that if neither of the above is true, we fall through and add the query to self.info

        # best k_v, k_d:
        k_v = all_query_bounds[best['v']][1]  # end of v match
        k_d = all_query_bounds[best['d']][1] - all_query_bounds[best['v']][1]  # end of d minus end of v

        if k_d_max < 5:  # since the s-w step matches to the longest possible j and then excises it, this sometimes gobbles up the d, resulting in a very short d alignment.
            if self.debug:
                print '  expanding k_d'
            k_d_max = max(8, k_d_max)

        if 'IGHJ4*' in best['j'] and self.glfo['seqs']['d'][best['d']][-5:] == 'ACTAC':  # the end of some d versions is the same as the start of some j versions, so the s-w frequently kicks out the 'wrong' alignment
            if self.debug:
                print '  doubly expanding k_d'
            if k_d_max-k_d_min < 8:
                k_d_min -= 5
                k_d_max += 2

        k_v_min = max(0, k_v_min - self.args.default_v_fuzz)  # ok, so I don't *actually* want it to be zero... oh, well
        k_v_max += self.args.default_v_fuzz
        k_d_min = max(1, k_d_min - self.args.default_d_fuzz)
        k_d_max += self.args.default_d_fuzz
        assert k_v_min > 0 and k_d_min > 0 and k_v_max > 0 and k_d_max > 0

        if self.debug:
            print '         k_v: %d [%d-%d)' % (k_v, k_v_min, k_v_max)
            print '         k_d: %d [%d-%d)' % (k_d, k_d_min, k_d_max)


        kvals = {}
        kvals['v'] = {'best':k_v, 'min':k_v_min, 'max':k_v_max}
        kvals['d'] = {'best':k_d, 'min':k_d_min, 'max':k_d_max}
        self.add_to_info(query_name, query_seq, kvals, match_names, best, all_germline_bounds, all_query_bounds, codon_positions=codon_positions)
