import time
import copy
import sys
import math
import re
import os
import itertools
import operator
import pysam
import contextlib
from collections import OrderedDict

import utils
import glutils
from opener import opener
from parametercounter import ParameterCounter
from performanceplotter import PerformancePlotter
from allelefinder import AlleleFinder

# ----------------------------------------------------------------------------------------
class Waterer(object):
    """ Run smith-waterman on the query sequences in <infname> """
    def __init__(self, args, input_info, reco_info, glfo, my_datadir, parameter_dir, write_parameters=False, find_new_alleles=False):
        self.parameter_dir = parameter_dir.rstrip('/')
        self.args = args
        self.debug = self.args.debug if self.args.sw_debug is None else self.args.sw_debug

        self.max_insertion_length = 35  # if vdjalign reports an insertion longer than this, rerun the query (typically with different match/mismatch ratio)
        self.absolute_max_insertion_length = 200  # just ignore them if it's longer than this

        self.input_info = input_info
        self.remaining_queries = set([q for q in self.input_info.keys()])  # we remove queries from this set when we're satisfied with the current output (in general we may have to rerun some queries with different match/mismatch scores)
        self.new_indels = 0  # number of new indels that were kicked up this time through

        self.match_mismatch = copy.deepcopy(self.args.initial_match_mismatch)  # don't want to modify it!
        self.gap_open_penalty = self.args.gap_open_penalty  # not modifying it now, but just to make sure we don't in the future

        self.reco_info = reco_info
        self.glfo = glfo
        self.info = {}
        self.info['queries'] = []  # list of queries that *passed* sw, i.e. for which we have information
        self.info['all_best_matches'] = set()  # every gene that was a best match for at least one query
        self.info['all_matches'] = {r : set() for r in utils.regions}  # every gene that was *any* match for at least one query
        self.info['indels'] = {}

        self.nth_try = 1
        self.unproductive_queries = set()

        # rewrite input germline sets (if needed)
        self.my_datadir = my_datadir

        self.alfinder, self.pcounter, self.true_pcounter, self.perfplotter = None, None, None, None
        if find_new_alleles:  # NOTE *not* the same as <self.args.find_new_alleles>
            self.alfinder = AlleleFinder(self.glfo, self.args)
        if write_parameters:  # NOTE *not* the same as <self.args.cache_parameters>
            self.pcounter = ParameterCounter(self.glfo, self.args)
            if not self.args.is_data:
                self.true_pcounter = ParameterCounter(self.glfo, self.args)
        if self.args.plot_performance:
            self.perfplotter = PerformancePlotter(self.glfo, 'sw')

        if not os.path.exists(self.args.ighutil_dir + '/bin/vdjalign'):
            raise Exception('ERROR ighutil path d.n.e: ' + self.args.ighutil_dir + '/bin/vdjalign')

    # ----------------------------------------------------------------------------------------
    def run(self):
        # start = time.time()
        base_infname = 'query-seqs.fa'
        base_outfname = 'query-seqs.bam'
        sys.stdout.flush()

        n_procs = self.args.n_fewer_procs
        initial_queries_per_proc = float(len(self.remaining_queries)) / n_procs
        while len(self.remaining_queries) > 0:  # we remove queries from <self.remaining_queries> as we're satisfied with their output
            if self.nth_try > 1 and float(len(self.remaining_queries)) / n_procs < initial_queries_per_proc:
                n_procs = int(max(1., float(len(self.remaining_queries)) / initial_queries_per_proc))
            self.write_vdjalign_input(base_infname, n_procs)
            self.execute_commands(base_infname, base_outfname, n_procs)
            self.read_output(base_outfname, n_procs)
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
            print '%d / %d = %.2f unproductive' % (skipped_unproductive, len(self.input_info), float(skipped_unproductive) / len(self.input_info)),
            if n_remaining > 0:
                print '   %d / %d = %.2f other' % (n_remaining, len(self.input_info), float(n_remaining) / len(self.input_info)),
            print ')',
        print ''
        sys.stdout.flush()
        if n_remaining > 0:
            printstr = '   %s %d missing %s' % (utils.color('red', 'warning'), n_remaining, utils.plural_str('annotation', n_remaining))
            if n_remaining < 15:
                printstr += ' (' + ':'.join(self.remaining_queries) + ')'
            print printstr
        if self.debug and len(self.info['indels']) > 0:
            print '      indels: %s' % ':'.join(self.info['indels'].keys())
        assert len(self.info['queries']) + skipped_unproductive + n_remaining == len(self.input_info)
        if self.debug and not self.args.is_data and n_remaining > 0:
            print 'true annotations for remaining events:'
            for qry in self.remaining_queries:
                utils.print_reco_event(self.glfo['seqs'], self.reco_info[qry], extra_str='      ', label='true:')
        if self.alfinder is not None:
            self.alfinder.finalize(debug=self.args.debug_new_allele_finding)
            self.info['new-alleles'] = self.alfinder.new_allele_info
            if self.args.plotdir is not None:
                self.alfinder.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots)

        # add padded info to self.info (returns if stuff has already been padded)
        self.pad_seqs_to_same_length()  # NOTE this uses *all the gene matches (not just the best ones), so it has to come before we call pcounter.write(), since that fcn rewrites the germlines removing genes that weren't best matches. But NOTE also that I'm not sure what but that the padding actually *needs* all matches (rather than just all *best* matches)

        if self.pcounter is not None:
            if self.args.plotdir is not None:
                self.pcounter.plot(self.args.plotdir + '/sw', subset_by_gene=True, codon_positions={r : self.glfo[c + '-positions'] for r, c in utils.conserved_codons.items()}, only_csv=self.args.only_csv_plots)
                if self.true_pcounter is not None:
                    self.true_pcounter.plot(self.args.plotdir + '/sw-true', subset_by_gene=True, codon_positions={r : self.glfo[c + '-positions'] for r, c in utils.conserved_codons.items()}, only_csv=self.args.only_csv_plots)
            self.pcounter.write(self.parameter_dir, self.my_datadir)
            if self.true_pcounter is not None:
                self.true_pcounter.write(self.parameter_dir + '-true')

        self.info['remaining_queries'] = self.remaining_queries

    # ----------------------------------------------------------------------------------------
    def subworkdir(self, iproc, n_procs):
        if n_procs == 1:
            return self.args.workdir
        else:
            return self.args.workdir + '/sw-' + str(iproc)

    # ----------------------------------------------------------------------------------------
    def execute_commands(self, base_infname, base_outfname, n_procs):
        # ----------------------------------------------------------------------------------------
        def get_outfname(iproc):
            return self.subworkdir(iproc, n_procs) + '/' + base_outfname
        # ----------------------------------------------------------------------------------------
        def get_cmd_str(iproc):
            return self.get_vdjalign_cmd_str(self.subworkdir(iproc, n_procs), base_infname, base_outfname, n_procs)

        # start all procs for the first time
        procs, n_tries = [], []
        for iproc in range(n_procs):
            procs.append(utils.run_cmd(get_cmd_str(iproc), self.subworkdir(iproc, n_procs)))
            n_tries.append(1)
            time.sleep(0.1)

        # keep looping over the procs until they're all done
        while procs.count(None) != len(procs):  # we set each proc to None when it finishes
            for iproc in range(n_procs):
                if procs[iproc] is None:  # already finished
                    continue
                if procs[iproc].poll() is not None:  # it's finished
                    utils.finish_process(iproc, procs, n_tries, self.subworkdir(iproc, n_procs), get_outfname(iproc), get_cmd_str(iproc))
            sys.stdout.flush()
            time.sleep(1)

        for iproc in range(n_procs):
            os.remove(self.subworkdir(iproc, n_procs) + '/' + base_infname)

        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def write_vdjalign_input(self, base_infname, n_procs):
        n_remaining = len(self.remaining_queries)
        queries_per_proc = float(n_remaining) / n_procs
        n_queries_per_proc = int(math.ceil(queries_per_proc))
        written_queries = set()  # make sure we actually write each query NOTE I should be able to remove this when I work out where they're disappearing to. But they don't seem to be disappearing any more, *sigh*
        if n_procs == 1:  # double check for rounding problems or whatnot
            assert n_queries_per_proc == n_remaining
        for iproc in range(n_procs):
            workdir = self.subworkdir(iproc, n_procs)
            if n_procs > 1:
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
                    written_queries.add(query_name)
                    iquery += 1
        not_written = self.remaining_queries - written_queries
        if len(not_written) > 0:
            raise Exception('didn\'t write %s to %s' % (':'.join(not_written), self.args.workdir))

    # ----------------------------------------------------------------------------------------
    def get_vdjalign_cmd_str(self, workdir, base_infname, base_outfname, n_procs=None):
        """
        Run smith-waterman alignment (from Connor's ighutils package) on the seqs in <base_infname>, and toss all the top matches into <base_outfname>.
        """
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        cmd_str = self.args.ighutil_dir + '/bin/vdjalign align-fastq -q'
        if self.args.slurm or utils.auto_slurm(n_procs):
            cmd_str = 'srun ' + cmd_str
        cmd_str += ' --locus ' + 'IG' + self.args.chain.upper()
        cmd_str += ' --max-drop 50'
        match, mismatch = self.match_mismatch
        cmd_str += ' --match ' + str(match) + ' --mismatch ' + str(mismatch)
        cmd_str += ' --gap-open ' + str(self.gap_open_penalty)
        cmd_str += ' --vdj-dir ' + self.my_datadir + '/' + self.args.chain
        cmd_str += ' --samtools-dir ' + self.args.partis_dir + '/packages/samtools'
        cmd_str += ' ' + workdir + '/' + base_infname + ' ' + workdir + '/' + base_outfname
        return cmd_str

    # ----------------------------------------------------------------------------------------
    def read_output(self, base_outfname, n_procs=1):
        queries_to_rerun = OrderedDict()  # This is to keep track of every query that we don't add to self.info (i.e. it does *not* include unproductive queries that we ignore/skip entirely because we were told to by a command line argument)
                                          # ...whereas <self.unproductive_queries> is to keep track of the queries that were definitively unproductive (i.e. we removed them from self.remaining_queries) when we were told to skip unproductives by a command line argument
        for reason in ['unproductive', 'no-match', 'weird-annot.', 'nonsense-bounds', 'invalid-codon']:
            queries_to_rerun[reason] = set()

        self.new_indels = 0
        queries_read_from_file = set()  # should be able to remove this, eventually
        for iproc in range(n_procs):
            outfname = self.subworkdir(iproc, n_procs) + '/' + base_outfname
            with contextlib.closing(pysam.Samfile(outfname)) as bam:
                grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
                for _, reads in grouped:  # loop over query sequences
                    qinfo = self.read_query(bam.references, list(reads))
                    self.summarize_query(qinfo, queries_to_rerun)  # returns before adding to <self.info> if it thinks we should rerun the query
                    queries_read_from_file.add(qinfo['name'])

        not_read = self.remaining_queries - queries_read_from_file
        if len(not_read) > 0:
            raise Exception('didn\'t read %s from %s' % (':'.join(not_read), self.args.workdir))

        if self.nth_try == 1:
            print '        processed       remaining      new-indels          rerun: ' + '      '.join([reason for reason in queries_to_rerun])
        print '      %8d' % len(queries_read_from_file),
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
                print ''
                raise Exception('numbers don\'t add up in sw output reader (n_to_rerun + new_indels != remaining_queries): %d + %d != %d   (look in %s)' % (n_to_rerun, self.new_indels, len(self.remaining_queries), self.args.workdir))
            if self.nth_try < 2 or self.new_indels == 0:  # increase the mismatch score if it's the first try, or if there's no new indels
                print '            increasing mismatch score (%d --> %d) and rerunning them' % (self.match_mismatch[1], self.match_mismatch[1] + 1)
                self.match_mismatch[1] += 1
            elif self.new_indels > 0:  # if there were some indels, rerun with the same parameters (but when the input is written the indel will be "reversed' in the sequences that's passed to ighutil)
                print '            rerunning for indels'
                self.new_indels = 0
            else:  # shouldn't get here
                assert False
        else:
            print '        all done'

        for iproc in range(n_procs):
            workdir = self.subworkdir(iproc, n_procs)
            os.remove(workdir + '/' + base_outfname)
            if n_procs > 1:  # still need the top-level workdir
                os.rmdir(workdir)

    # ----------------------------------------------------------------------------------------
    def get_indel_info(self, query_name, cigarstr, qrseq, glseq, gene):
        cigars = re.findall('[0-9][0-9]*[A-Z]', cigarstr)  # split cigar string into its parts
        cigars = [(cstr[-1], int(cstr[:-1])) for cstr in cigars]  # split each part into the code and the length

        codestr = ''
        qpos = 0  # position within query sequence
        indelfo = utils.get_empty_indel()  # replacement_seq: query seq with insertions removed and germline bases inserted at the position of deletions
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
    def add_dummy_d_match(self, qinfo, debug=True):
        dummy_d = glutils.dummy_d_genes[self.args.chain]
        qinfo['matches']['d'].append((1, dummy_d))
        v_end = qinfo['first_match_qrbounds'][1]
        qinfo['qrbounds'][dummy_d] = (v_end, v_end)
        qinfo['glbounds'][dummy_d] = (0, 0)

    # ----------------------------------------------------------------------------------------
    def read_query(self, references, reads):
        """ convert bam crap to python dict """
        primary = next((r for r in reads if not r.is_secondary), None)
        qinfo = {
            'name' : primary.qname,
            'seq' : primary.seq,
            'matches' : {r : [] for r in utils.regions},
            'qrbounds' : {},
            'glbounds' : {},
            'first_match_qrbounds' : None,  # since sw excises its favorite v match, we have to know this match's boundaries in order to calculate k_d for all the other matches
            'new_indel' : False
        }
        for read in reads:  # loop over the matches found for each query sequence
            # set this match's values
            read.seq = qinfo['seq']  # only the first one has read.seq set by default, so we need to set the rest by hand
            gene = references[read.tid]
            region = utils.get_region(gene)
            raw_score = read.tags[0][1]  # raw because they don't include the gene choice probs
            score = raw_score
            qrbounds = (read.qstart, read.qend)
            glbounds = (read.pos, read.aend)
            if region == 'v' and qinfo['first_match_qrbounds'] is None:
                qinfo['first_match_qrbounds'] = qrbounds

            # TODO I wish this wasn't here and I suspect I don't really need it (any more) UPDATE I dunno, this definitely eliminates some stupid (albeit rare) matches
            if region == 'v':  # skip matches with cpos past the end of the query seq (i.e. eroded a ton on the right side of the v)
                cpos = self.glfo['cyst-positions'][gene] - glbounds[0] + qrbounds[0]  # position within original germline gene, minus the position in that germline gene at which the match starts, plus the position in the query sequence at which the match starts
                if cpos < 0 or cpos >= len(qinfo['seq']):
                    continue

            if 'I' in read.cigarstring or 'D' in read.cigarstring:  # skip indels, and tell the HMM to skip indels (you won't see any unless you decrease the <self.gap_open_penalty>)
                if self.args.no_indels:  # you can forbid indels on the command line
                    continue
                if self.nth_try < 2:  # we also forbid indels on the first try (we want to increase the mismatch score before we conclude it's "really" an indel)
                    continue
                if len(qinfo['matches'][region]) == 0:  # if this is the first (best) match for this region, allow indels (otherwise skip the match)
                    if qinfo['name'] not in self.info['indels']:
                        self.info['indels'][qinfo['name']] = self.get_indel_info(qinfo['name'], read.cigarstring, qinfo['seq'][qrbounds[0] : qrbounds[1]], self.glfo['seqs'][region][gene][glbounds[0] : glbounds[1]], gene)
                        self.info['indels'][qinfo['name']]['reversed_seq'] = qinfo['seq'][ : qrbounds[0]] + self.info['indels'][qinfo['name']]['reversed_seq'] + qinfo['seq'][qrbounds[1] : ]
                        self.new_indels += 1
                        qinfo['new_indel'] = True
                        return qinfo  # don't process this query any further -- since it's now in the indel info it'll get run next time through
                    else:
                        if self.debug:
                            print '     ignoring subsequent indels for %s' % qinfo['name']
                        continue  # hopefully there's a later match without indels
                else:
                    continue

            # and finally add this match's information
            qinfo['matches'][region].append((score, gene))  # NOTE it is important that this is ordered such that the best match is first UPDATE huh, maybe I wrote this before I added the sorted() below? In any case it seems to always be sorted in the bam file, so this is really just a double-check
            qinfo['qrbounds'][gene] = qrbounds
            qinfo['glbounds'][gene] = glbounds

        if self.args.chain != 'h':
            self.add_dummy_d_match(qinfo)

        for region in utils.regions:
            qinfo['matches'][region] = sorted(qinfo['matches'][region], reverse=True)

        return qinfo

    # ----------------------------------------------------------------------------------------
    def print_match(self, region, gene, score, qseq, glbounds, qrbounds, skipping=False):
        out_str_list = []
        buff_str = (19 - len(gene)) * ' '
        out_str_list.append('%8s%s%s%9s%3s %6.0f        ' % (' ', utils.color_gene(gene), '', '', buff_str, score))
        out_str_list.append('%4d%4d   %s\n' % (glbounds[0], glbounds[1], self.glfo['seqs'][region][gene][glbounds[0]:glbounds[1]]))
        out_str_list.append('%46s  %4d%4d' % ('', qrbounds[0], qrbounds[1]))
        out_str_list.append('   %s ' % (utils.color_mutants(self.glfo['seqs'][region][gene][glbounds[0]:glbounds[1]], qseq[qrbounds[0]:qrbounds[1]])))
        if skipping:
            out_str_list.append('skipping!')

        print ''.join(out_str_list)

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
    def check_boundaries(self, rpair, qinfo, best, recursed=False, debug=False):
        # NOTE this duplicates code in shift_overlapping_boundaries(), which makes me cranky, but this setup avoids other things I dislike more
        l_reg = rpair['left']
        r_reg = rpair['right']
        l_gene = best[l_reg]
        r_gene = best[r_reg]

        overlap, available_space = self.get_overlap_and_available_space(rpair, best, qinfo['qrbounds'])

        if debug:
            print '  %s %s    overlap %d    available space %d' % (l_reg, r_reg, overlap, available_space)

        status = 'ok'
        if overlap > 0:  # positive overlap means they actually overlap
            status = 'overlap'
        if overlap > available_space or overlap == 1 and available_space == 1:  # call it nonsense if the boundaries are really whack (i.e. there isn't enough space to resolve the overlap) -- we'll presumably either toss the query or rerun with different match/mismatch
            status = 'nonsense'

        if debug:
            print '  overlap status: %s' % status

        if not recursed and status == 'nonsense' and l_reg == 'd' and self.nth_try > 2:  # on rare occasions with very high mutation, vdjalign refuses to give us a j match that's at all to the right of the d match
            assert l_reg == 'd' and r_reg == 'j'
            if debug:
                print '  %s: synthesizing d match' % qinfo['name']
            leftmost_position = min(qinfo['qrbounds'][l_gene][0], qinfo['qrbounds'][r_gene][0])
            qinfo['qrbounds'][l_gene] = (leftmost_position, leftmost_position + 1)  # swap whatever crummy nonsense d match we have now for a one-base match at the left end of things (things in practice should be left end of j match)
            qinfo['glbounds'][l_gene] = (0, 1)
            status = self.check_boundaries(rpair, qinfo, best, recursed=True, debug=debug)
            if status == 'overlap':
                if debug:
                    print '  \'overlap\' status after synthesizing d match. Setting to \'nonsense\', I can\'t deal with this bullshit'
                status = 'nonsense'

        return status

    # ----------------------------------------------------------------------------------------
    def shift_overlapping_boundaries(self, rpair, qinfo, best, debug=False):
        """
        s-w allows d and j matches (and v and d matches) to overlap... which makes no sense, so apportion the disputed territory between the two regions.
        Note that this still works if, say, v is the entire sequence, i.e. one match is entirely subsumed by another.
        """
        l_reg = rpair['left']
        r_reg = rpair['right']
        l_gene = best[l_reg]
        r_gene = best[r_reg]

        overlap, available_space = self.get_overlap_and_available_space(rpair, best, qinfo['qrbounds'])

        if overlap <= 0:  # nothing to do, they're already consistent
            print 'shouldn\'t get here any more if there\'s no overlap'
            return

        if overlap > available_space:
            raise Exception('overlap %d bigger than available space %d between %s and %s for %s' % (overlap, available_space, l_reg, r_reg, qinfo['name']))

        if debug:
            print '%s%s:  %d-%d overlaps with %d-%d by %d' % (l_reg, r_reg, qinfo['qrbounds'][l_gene][0], qinfo['qrbounds'][l_gene][1], qinfo['qrbounds'][r_gene][0], qinfo['qrbounds'][r_gene][1], overlap)

        l_length = qinfo['qrbounds'][l_gene][1] - qinfo['qrbounds'][l_gene][0]  # initial length of lefthand gene match
        r_length = qinfo['qrbounds'][r_gene][1] - qinfo['qrbounds'][r_gene][0]  # and same for the righthand one
        l_portion, r_portion = 0, 0  # portion of the initial overlap that we give to each side
        if debug:
            print '    lengths        portions     '
        while l_portion + r_portion < overlap:
            if debug:
                print '  %4d %4d      %4d %4d' % (l_length, r_length, l_portion, r_portion)
            if l_length <= 1 and r_length <= 1:  # don't want to erode match (in practice it'll be the d match) all the way to zero
                raise Exception('both lengths went to one without resolving overlap for %s: %s %s' % (qinfo['name'], qinfo['qrbounds'][l_gene], qinfo['qrbounds'][r_gene]))
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
            print '      %s apportioning %d bases between %s (%d) match and %s (%d) match' % (qinfo['name'], overlap, l_reg, l_portion, r_reg, r_portion)
        assert l_portion + r_portion == overlap
        qinfo['qrbounds'][l_gene] = (qinfo['qrbounds'][l_gene][0], qinfo['qrbounds'][l_gene][1] - l_portion)
        qinfo['glbounds'][l_gene] = (qinfo['glbounds'][l_gene][0], qinfo['glbounds'][l_gene][1] - l_portion)
        qinfo['qrbounds'][r_gene] = (qinfo['qrbounds'][r_gene][0] + r_portion, qinfo['qrbounds'][r_gene][1])
        qinfo['glbounds'][r_gene] = (qinfo['glbounds'][r_gene][0] + r_portion, qinfo['glbounds'][r_gene][1])

    # ----------------------------------------------------------------------------------------
    def add_to_info(self, qinfo, best, codon_positions):
        qname = qinfo['name']
        assert qname not in self.info

        match_names = {r : [g for _, g in qinfo['matches'][r]] for r in utils.regions}

        self.info['queries'].append(qname)
        self.info[qname] = {}
        self.info[qname]['unique_id'] = qname  # redundant, but used somewhere down the line

        kbounds = self.get_kbounds(qinfo, best)
        self.info[qname]['k_v'] = kbounds['v']
        self.info[qname]['k_d'] = kbounds['d']

        self.info[qname]['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])  # all gene matches for this query

        self.info[qname]['cdr3_length'] = codon_positions['j'] - codon_positions['v'] + 3
        self.info[qname]['codon_positions'] = copy.deepcopy(codon_positions)

        # erosion, insertion, mutation info for best match
        self.info[qname]['v_5p_del'] = qinfo['glbounds'][best['v']][0]
        self.info[qname]['v_3p_del'] = len(self.glfo['seqs']['v'][best['v']]) - qinfo['glbounds'][best['v']][1]  # len(germline v) - gl_match_end
        self.info[qname]['d_5p_del'] = qinfo['glbounds'][best['d']][0]
        self.info[qname]['d_3p_del'] = len(self.glfo['seqs']['d'][best['d']]) - qinfo['glbounds'][best['d']][1]
        self.info[qname]['j_5p_del'] = qinfo['glbounds'][best['j']][0]
        self.info[qname]['j_3p_del'] = len(self.glfo['seqs']['j'][best['j']]) - qinfo['glbounds'][best['j']][1]

        self.info[qname]['fv_insertion'] = qinfo['seq'][ : qinfo['qrbounds'][best['v']][0]]
        self.info[qname]['vd_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['v']][1] : qinfo['qrbounds'][best['d']][0]]
        self.info[qname]['dj_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['d']][1] : qinfo['qrbounds'][best['j']][0]]
        self.info[qname]['jf_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['j']][1] : ]

        self.info[qname]['indelfo'] = self.info['indels'].get(qname, utils.get_empty_indel())

        for region in utils.regions:
            self.info[qname][region + '_gene'] = best[region]
            self.info['all_best_matches'].add(best[region])
            self.info['all_matches'][region] |= set(match_names[region])

        self.info[qname]['seq'] = qinfo['seq']  # NOTE this is the seq output by vdjalign, i.e. if we reversed any indels it is the reversed sequence

        existing_implicit_keys = tuple(['cdr3_length', 'codon_positions'])
        utils.add_implicit_info(self.glfo, self.info[qname], multi_seq=False, existing_implicit_keys=existing_implicit_keys)

        if self.debug:
            if not self.args.is_data:
                utils.print_reco_event(self.glfo['seqs'], self.reco_info[qname], extra_str='      ', label='true:')
            utils.print_reco_event(self.glfo['seqs'], self.info[qname], extra_str='      ', label='inferred:')

        if self.alfinder is not None:
            self.alfinder.increment(self.info[qname])
        if self.pcounter is not None:
            self.pcounter.increment_all_params(self.info[qname])
            if self.true_pcounter is not None:
                self.true_pcounter.increment_all_params(self.reco_info[qname])
        if self.perfplotter is not None:
            if qname in self.info['indels']:
                print '    skipping performance evaluation of %s because of indels' % qname  # I just have no idea how to handle naive hamming fraction when there's indels
            else:
                self.perfplotter.evaluate(self.reco_info[qname], self.info[qname])

        self.remaining_queries.remove(qname)

    # ----------------------------------------------------------------------------------------
    def summarize_query(self, qinfo, queries_to_rerun):
        """
        Fiddle with a few things, but mostly decide whether we're satisfied with the current matches.
        If not, return without calling add_to_info().
        """
        qname = qinfo['name']
        assert qname not in self.info

        if self.debug >= 2:
            print qname
            for region in utils.regions:
                for score, gene in qinfo['matches'][region]:  # sorted by decreasing match quality
                    self.print_match(region, gene, score, qinfo['seq'], qinfo['glbounds'][gene], qinfo['qrbounds'][gene], skipping=False)

        if qinfo['new_indel']:
            if self.debug:
                print '    new indel -- rerunning'
            return

        # do we have a match for each region?
        for region in utils.getregions(self.args.chain):
            if len(qinfo['matches'][region]) == 0:
                if self.debug:
                    print '      no', region, 'match found for', qname  # TODO if no d match found, we should really just assume entire d was eroded
                queries_to_rerun['no-match'].add(qname)
                return

        best = {r : qinfo['matches'][r][0][1] for r in utils.regions}  # already made sure there's at least one match for each region

        # s-w allows d and j matches to overlap, so we need to apportion the disputed bases
        for rpair in utils.region_pairs():
            overlap_status = self.check_boundaries(rpair, qinfo, best)
            if overlap_status == 'overlap':
                self.shift_overlapping_boundaries(rpair, qinfo, best)
            elif overlap_status == 'nonsense':
                queries_to_rerun['nonsense-bounds'].add(qname)
                return
            else:
                assert overlap_status == 'ok'

        # check for suspiciously bad annotations
        for rp in utils.region_pairs():
            insertion_length = qinfo['qrbounds'][best[rp['right']]][0] - qinfo['qrbounds'][best[rp['left']]][1]  # start of right match minus end of left one
            if insertion_length > self.absolute_max_insertion_length or (self.nth_try < 2 and insertion_length > self.max_insertion_length):
                if self.debug:
                    print '      suspiciously long insertion in %s, rerunning' % qname
                queries_to_rerun['weird-annot.'].add(qname)
                return

        if self.debug == 1:
            print qname

        # set and check conserved codon positions
        codon_positions = {}
        for region, codon in utils.conserved_codons[self.args.chain].items():
            # position within original germline gene, minus the position in that germline gene at which the match starts, plus the position in the query sequence at which the match starts
            pos = self.glfo[codon + '-positions'][best[region]] - qinfo['glbounds'][best[region]][0] + qinfo['qrbounds'][best[region]][0]
            if pos < 0 or pos >= len(qinfo['seq']):
                if self.debug:
                    print '      invalid %s codon position (%d in seq of length %d), rerunning' % (codon, pos, len(qinfo['seq']))
                queries_to_rerun['invalid-codon'].add(qname)
                return
            codon_positions[region] = pos

        # check for unproductive rearrangements
        codons_ok = utils.both_codons_ok(self.args.chain, qinfo['seq'], codon_positions)
        cdr3_length = codon_positions['j'] - codon_positions['v'] + 3

        if cdr3_length < 6:  # NOTE six is also hardcoded in utils
            if self.debug:
                print '      negative cdr3 length %d' % (cdr3_length)
            queries_to_rerun['invalid-codon'].add(qname)
            return

        in_frame_cdr3 = (cdr3_length % 3 == 0)
        no_stop_codon = utils.stop_codon_check(qinfo['seq'], codon_positions['v'])
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
                queries_to_rerun['unproductive'].add(qname)
                return
            elif self.args.skip_unproductive:
                if self.debug:
                    print '            ...skipping'
                self.unproductive_queries.add(qname)
                self.remaining_queries.remove(qname)
                return
            else:
                pass  # this is here so you don't forget that if neither of the above is true, we fall through and add the query to self.info

        self.add_to_info(qinfo, best, codon_positions)

    # ----------------------------------------------------------------------------------------
    def get_kbounds(self, qinfo, best):
        # OR of k-space for all the matches
        k_v_min, k_d_min = 999, 999
        k_v_max, k_d_max = 0, 0
        for region in utils.regions:  # NOTE since here I'm not yet skipping genes beyond the first <args.n_max_per_region>, this is overly conservative. Don't really care, though... k space integration is mostly pretty cheap what with chunk caching
            for _, gene in qinfo['matches'][region]:
                if region == 'v':
                    this_k_v = qinfo['qrbounds'][gene][1]  # NOTE even if the v match doesn't start at the left hand edge of the query sequence, we still measure k_v from there. In other words, sw doesn't tell the hmm about it
                    k_v_min = min(this_k_v, k_v_min)
                    k_v_max = max(this_k_v, k_v_max)
                elif region == 'd':
                    this_k_d = qinfo['qrbounds'][gene][1] - qinfo['first_match_qrbounds'][1]  # end of d minus end of v
                    k_d_min = min(this_k_d, k_d_min)
                    k_d_max = max(this_k_d, k_d_max)

        # best k_v, k_d:
        k_v = qinfo['qrbounds'][best['v']][1]  # end of v match
        k_d = qinfo['qrbounds'][best['d']][1] - qinfo['qrbounds'][best['v']][1]  # end of d minus end of v

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

        k_v_min = max(1, k_v_min - self.args.default_v_fuzz)  # ok, so I don't *actually* want it to be zero... oh, well
        k_v_max += self.args.default_v_fuzz
        k_d_min = max(1, k_d_min - self.args.default_d_fuzz)
        k_d_max += self.args.default_d_fuzz
        assert k_v_min > 0 and k_d_min > 0 and k_v_max > 0 and k_d_max > 0

        if self.debug:
            print '         k_v: %d [%d-%d)' % (k_v, k_v_min, k_v_max)
            print '         k_d: %d [%d-%d)' % (k_d, k_d_min, k_d_max)

        kbounds = {}
        kbounds['v'] = {'best' : k_v, 'min' : k_v_min, 'max' : k_v_max}
        kbounds['d'] = {'best' : k_d, 'min' : k_d_min, 'max' : k_d_max}
        return kbounds

    # ----------------------------------------------------------------------------------------
    def get_padding_parameters(self, debug=False):
        maxima = {'gl_cpos' : None, 'gl_cpos_to_j_end' : None}
        for query in self.info['queries']:
            swfo = self.info[query]

            # find biggest cyst position among all gl matches
            fvstuff = max(0, len(swfo['fv_insertion']) - swfo['v_5p_del'])  # we always want to pad out to the entire germline sequence, so don't let this go negative
            for v_match in self.info['all_matches']['v']:  # includes matches for *other* sequences, because we want bcrham to be able to compare any sequence to any other (although, could probably use all *best* matches rather than all *all*)
                gl_cpos = self.glfo['cyst-positions'][v_match] + fvstuff
                if maxima['gl_cpos'] is None or gl_cpos > maxima['gl_cpos']:
                    maxima['gl_cpos'] = gl_cpos

            # Since we only store j_3p_del for the best match, we can't loop over all of 'em. But j stuff doesn't vary too much, so it works ok.
            cpos = swfo['codon_positions']['v']  # cyst position in query sequence (as opposed to gl_cpos, which is in germline allele)
            jfstuff = max(0, len(swfo['jf_insertion']) - swfo['j_3p_del'])
            gl_cpos_to_j_end = len(swfo['seq']) - cpos + swfo['j_3p_del'] + jfstuff
            if maxima['gl_cpos_to_j_end'] is None or gl_cpos_to_j_end > maxima['gl_cpos_to_j_end']:
                maxima['gl_cpos_to_j_end'] = gl_cpos_to_j_end

        if debug:
            print '    maxima:',
            for k, v in maxima.items():
                print '%s %d    ' % (k, v),
            print ''

        return maxima

    # ----------------------------------------------------------------------------------------
    def pad_seqs_to_same_length(self, debug=False):
        """
        Pad all sequences in <seqinfo> to the same length to the left and right of their conserved cysteine positions.
        Next, pads all sequences further out (if necessary) such as to eliminate all v_5p and j_3p deletions.
        """

        maxima = self.get_padding_parameters(debug=debug)

        for query in self.info['queries']:
            swfo = self.info[query]
            if 'padded' in swfo:  # already added padded information (we're probably partitioning, and this is not the first step)
                return
            seq = swfo['seq']
            cpos = swfo['codon_positions']['v']
            if cpos < 0 or cpos >= len(seq):
                print 'hm now what do I want to do here?'
            k_v = swfo['k_v']

            padleft = maxima['gl_cpos'] - cpos  # left padding: biggest germline cpos minus cpos in this sequence
            padright = maxima['gl_cpos_to_j_end'] - (len(seq) - cpos)
            if padleft < 0 or padright < 0:
                raise Exception('bad padding %d %d for %s' % (padleft, padright, query))

            padfo = {}
            padfo['seq'] = padleft * utils.ambiguous_bases[0] + seq + padright * utils.ambiguous_bases[0]
            if query in self.info['indels']:  # also pad the reversed sequence
                self.info['indels'][query]['reversed_seq'] = padleft * utils.ambiguous_bases[0] + self.info['indels'][query]['reversed_seq'] + padright * utils.ambiguous_bases[0]
            padfo['k_v'] = {'min' : k_v['min'] + padleft, 'max' : k_v['max'] + padleft}
            padfo['cyst_position'] = swfo['codon_positions']['v'] + padleft
            padfo['padleft'] = padleft
            padfo['padright'] = padright
            if debug:
                print '      pad %d %d   %s' % (padleft, padright, query)
                print '     %d --> %d (%d-%d --> %d-%d)' % (len(seq), len(padfo['seq']), k_v['min'], k_v['max'], padfo['k_v']['min'], padfo['k_v']['max'])
            swfo['padded'] = padfo

        if debug:
            for query in self.info['queries']:
                print '%20s %s' % (query, self.info[query]['padded']['seq'])
