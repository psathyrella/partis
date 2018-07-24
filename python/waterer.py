import time
import copy
import sys
import math
import os
import itertools
import operator
import pysam
import contextlib
from collections import OrderedDict
import csv
import numpy
import traceback

import utils
import glutils
import indelutils
from parametercounter import ParameterCounter
from performanceplotter import PerformancePlotter

# best mismatch (with a match score of 5):
# mfreq     boundaries     mutation
# 0.01         5              5
# 0.05         5              5-
# 0.15         3+             3
# 0.2         2-3             2-3
# 0.3          2-             2-
# +: next higher number isn't that much worse
# -: [...]
# mfreq was I think the sequence-wide mfreq, but was close enough to the v value that it doesn't matter

# ----------------------------------------------------------------------------------------
class Waterer(object):
    """ Run smith-waterman on the query sequences in <infname> """
    def __init__(self, args, glfo, input_info, simglfo, reco_info,
                 count_parameters=False, parameter_out_dir=None, plot_annotation_performance=False,
                 duplicates=None, pre_failed_queries=None, aligned_gl_seqs=None, vs_info=None):
        self.args = args
        self.input_info = input_info  # NOTE do *not* modify this, since it's this original input info from partitiondriver
        self.reco_info = reco_info
        self.glfo = glfo  # NOTE gets overwritten by read_cachefile()
        self.count_parameters = count_parameters  # NOTE *not* the same as <self.args.cache_parameters>
        self.plot_annotation_performance = plot_annotation_performance  # NOTE *not* the same as <self.args.plot_annotation_performance>
        self.simglfo = simglfo
        self.parameter_out_dir = parameter_out_dir
        self.duplicates = {} if duplicates is None else duplicates
        self.debug = self.args.debug if self.args.sw_debug is None else self.args.sw_debug
        self.aligned_gl_seqs = aligned_gl_seqs
        self.vs_info = vs_info

        self.absolute_max_insertion_length = 120  # but if it's longer than this, we always skip the annotation

        self.gap_open_penalty = self.args.gap_open_penalty  # not modifying it now, but just to make sure we don't in the future
        self.match_score = 5  # see commented table above ^
        if self.vs_info is None:
            self.mismatch = 4
        else:
            self.default_mfreq = 0.075  # what to use if vsearch failed on the query
            self.mfreq_mismatch_vals = [
                (0.01, 5),
                (0.05, 5),
                (0.15, 3),
                (0.25, 2),
            ]

        self.info = {}
        self.info['queries'] = []  # list of queries that *passed* sw, i.e. for which we have information
        self.info['all_best_matches'] = set()  # every gene that was a best match for at least one query
        self.info['all_matches'] = {r : set() for r in utils.regions}  # every gene that was *any* match (up to <self.args.n_max_per_region[ireg]>) for at least one query NOTE there is also an 'all_matches' in each query's info
        self.info['indels'] = {}  # NOTE if we find shm indels in a sequence, we store the indel info in here, and rerun sw with the reversed sequence (i.e. <self.info> contains the sw inference on the reversed sequence -- if you want the original sequence, get that from <self.input_info>)
        self.info['failed-queries'] = set() if pre_failed_queries is None else copy.deepcopy(pre_failed_queries)  # not really sure about the deepcopy(), but it's probably safer?
        self.info['passed-queries'] = set()
        self.info['duplicates'] = self.duplicates  # it would be really nice to rationalize this
        self.info['removed-queries'] = set()  # ...and this

        self.remaining_queries = set(self.input_info) - self.info['failed-queries']  # we remove queries from this set when we're satisfied with the current output (in general we may have to rerun some queries with different match/mismatch scores)
        self.vs_indels = set()
        self.indel_reruns = set()  # queries that either failed during indel handling, or had successful indel handling: in both cases we rerun them, with a super large gap open to prevent further indels

        self.skipped_unproductive_queries, self.kept_unproductive_queries = set(), set()

        self.my_gldir = self.args.workdir + '/sw-' + glutils.glfo_dir
        glutils.write_glfo(self.my_gldir, self.glfo)  # NOTE gets overwritten by read_cachefile()

        if not os.path.exists(self.args.ig_sw_binary):
            raise Exception('ig-sw binary d.n.e: %s' % self.args.ig_sw_binary)

    # ----------------------------------------------------------------------------------------
    def run(self, cachefname=None):
        start = time.time()
        base_infname = 'query-seqs.fa'
        base_outfname = 'query-seqs.sam'

        if self.vs_info is not None:  # if we're reading a cache file, we should make sure to read the exact same info from there
            self.add_vs_indels()

        itry = 0
        while True:  # if we're not running vsearch, we still gotta run twice to get shm indeld sequences
            mismatches, gap_opens, queries_for_each_proc = self.split_queries(self.args.n_procs)  # NOTE can tell us to run more than <self.args.n_procs> (we run at least one proc for each different mismatch score)
            self.write_input_files(base_infname, queries_for_each_proc)

            print '    running %d proc%s for %d seq%s' % (len(mismatches), utils.plural(len(mismatches)), len(self.remaining_queries), utils.plural(len(self.remaining_queries)))
            sys.stdout.flush()
            self.execute_commands(base_infname, base_outfname, mismatches, gap_opens)

            processing_start = time.time()
            self.read_output(base_outfname, len(mismatches))

            if itry > 1 or len(self.indel_reruns) == 0:
                break
            itry += 1

        self.finalize(cachefname)
        print '    water time: %.1f  (ig-sw %.1f  processing %.1f)' % (time.time() - start, time.time() - processing_start, self.ig_sw_time)

    # ----------------------------------------------------------------------------------------
    def clean_cache(self, cache_path):
        for suffix in ['.csv', '.yaml']:
            if os.path.exists(cache_path + suffix):
                print '  removing old sw cache %s%s' % (cache_path, suffix)
                os.remove(cache_path + suffix)
        if os.path.exists(cache_path + '-glfo'):
            print '  removing old sw cache glfo %s-glfo' % cache_path
            glutils.remove_glfo_files(cache_path + '-glfo', self.args.locus)

    # ----------------------------------------------------------------------------------------
    def read_cachefile(self, cachefname):
        start = time.time()
        print '        reading sw results from %s' % cachefname

        if utils.getsuffix(cachefname) == '.csv':  # old way
            cachebase = utils.getprefix(cachefname)
            if os.path.exists(cachebase + '-glfo'):  # NOTE replaces original <self.glfo>
                self.glfo = glutils.read_glfo(cachebase + '-glfo', self.args.locus)
            else:
                print '    %s didn\'t find a germline info dir along with sw cache file, but trying to read it anyway' % utils.color('red', 'warning')
            cachefile = open(cachefname)  # closes on function exit, and no this isn't a great way of doing it (but it needs to stay open for the loop over <reader>)
            reader = csv.DictReader(cachefile)
        elif utils.getsuffix(cachefname) == '.yaml':  # new way
            self.glfo, reader, _ = utils.read_yaml_output(cachefname, dont_add_implicit_info=True)  # add implicit info below, so we can skip some of 'em and use aligned gl seqs
        else:
            raise Exception('unhandled sw cache file suffix %s' % cachefname)

        for line in reader:  # NOTE failed queries are *not* written to the cache file -- they're assumed to be whatever's in input info that's missing
            if utils.getsuffix(cachefname) == '.csv':
                utils.process_input_line(line)
                for key in [k for k in [r + '_per_gene_support' for r in utils.regions] if k in line]:  # new files shouldn't have this, but I think I need to leave it for reading older files
                    del line[key]
            assert len(line['unique_ids']) == 1  # would only fail if this was not actually an sw cache file, but it's still nice to check since so many places in waterer assume it's length 1
            if line['unique_ids'][0] not in self.input_info:
                continue
            utils.add_implicit_info(self.glfo, line, aligned_gl_seqs=self.aligned_gl_seqs)
            if indelutils.has_indels(line['indelfos'][0]):
                self.info['indels'][line['unique_ids'][0]] = line['indelfos'][0]
            self.add_to_info(line)

        glutils.write_glfo(self.my_gldir, self.glfo)

        self.finalize(cachefname=None, just_read_cachefile=True)
        print '        water time: %.1f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def write_cachefile(self, cachefname):
        if self.args.write_trimmed_and_padded_seqs_to_sw_cachefname:  # hackey workaround: (in case you want to use trimmed/padded seqs for something, but shouldn't be used in general)
            self.pad_seqs_to_same_length()

        print '        writing sw results to %s' % cachefname

        if utils.getsuffix(cachefname) == '.csv':
            cachebase = utils.getprefix(cachefname)
            glutils.write_glfo(cachebase + '-glfo', self.glfo)

        # NOTE do _not_ add extra headers here since if they're in the sw cache file I'd have to deal with removing them when I read it
        # NOTE does *not* write failed queries also
        utils.write_annotations(cachefname, self.glfo, [self.info[q]for q in self.info['queries']], utils.sw_cache_headers)

    # ----------------------------------------------------------------------------------------
    def finalize(self, cachefname=None, just_read_cachefile=False):
        if self.debug:
            print '%s' % utils.color('green', 'finalizing')

        self.info['queries'] = [q for q in self.input_info if q in self.info['passed-queries']]  # it might be cleaner eventually to just make self.info['queries'] a set, but it's used in so many other places that I'm worried about changing it
        self.info['failed-queries'] |= set(self.remaining_queries)  # perhaps it doesn't make sense to keep both of 'em around after finishing?

        if len(self.info['queries']) == 0:
            print '%s no queries passing smith-waterman, exiting' % utils.color('red', 'warning')
            sys.exit(1)

        print '      info for %d / %d = %.3f   (%d failed)' % (len(self.info['queries']), len(self.input_info), float(len(self.info['queries'])) / len(self.input_info), len(self.info['failed-queries']))
        if len(self.kept_unproductive_queries) > 0:
            print '      kept %d (%.3f) unproductive' % (len(self.kept_unproductive_queries), float(len(self.kept_unproductive_queries)) / len(self.input_info))

        if not just_read_cachefile:  # it's past tense!
            if len(self.skipped_unproductive_queries) > 0:
                print '         skipped %d unproductive' % len(self.skipped_unproductive_queries)
            if self.debug:
                if len(self.info['indels']) > 0:
                    print '      indels: %s' % ':'.join(self.info['indels'].keys())
                if len(self.info['failed-queries']) > 0:
                    print '            missing annotations: ' + ' '.join(self.info['failed-queries'])
                    if not self.args.is_data:
                        print 'true annotations for remaining events:'
                        for qry in self.info['failed-queries']:
                            utils.print_reco_event(self.reco_info[qry], extra_str='      ', label='true:')

            assert len(self.info['queries']) + len(self.skipped_unproductive_queries) + len(self.info['failed-queries']) == len(self.input_info)

        if self.count_parameters:
            pcounter = ParameterCounter(self.glfo, self.args)
            for qname in self.info['queries']:
                pcounter.increment(self.info[qname])
        if self.plot_annotation_performance:
            perfplotter = PerformancePlotter('sw')
            for qname in self.info['queries']:
                perfplotter.evaluate(self.reco_info[qname], self.info[qname], simglfo=self.simglfo)

        # remove queries with cdr3 length different to the seed sequence
        if self.args.seed_unique_id is not None:  # it might be nice to get the seed cdr3 length before running anything, and then not add seqs with different cdr3 length to start with, so those other sequences' gene matches don't get mixed in? then again  aren't we pretty much always reading cached values if we're seed partitioning?
            if self.args.seed_unique_id in self.info['queries']:
                seed_cdr3_length = self.info[self.args.seed_unique_id]['cdr3_length']
            else:  # if it failed
                print '    seed unique id \'%s\' not in final s-w queries, so removing all queries' % self.args.seed_unique_id
                seed_cdr3_length = -1
            initial_n_queries = len(self.info['queries'])
            for query in copy.deepcopy(self.info['queries']):
                if self.info[query]['cdr3_length'] != seed_cdr3_length:
                    self.remove_query(query)
            n_removed = initial_n_queries - len(self.info['queries'])
            if n_removed > 0:
                print '      removed %d / %d = %.2f sequences with cdr3 length different from seed sequence (leaving %d)' % (n_removed, initial_n_queries, float(n_removed) / initial_n_queries, len(self.info['queries']))

        if not self.args.dont_remove_framework_insertions and self.args.is_data:  # don't want to do this on simulation -- it's too much trouble to keep things consistent with the simulation info (it would also screw up the purity/completeness calculation)
            self.remove_framework_insertions()
            self.remove_duplicate_sequences()

        # want to do this *before* we pad sequences, so that when we read the cache file we're reading unpadded sequences and can pad them below
        if cachefname is not None:
            self.write_cachefile(cachefname)

        self.pad_seqs_to_same_length()  # NOTE this uses all the gene matches (not just the best ones), so it has to come before we call pcounter.write(), since that fcn rewrites the germlines removing genes that weren't best matches. But NOTE also that I'm not sure what but that the padding actually *needs* all matches (rather than just all *best* matches)

        if self.plot_annotation_performance:
            perfplotter.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots)

        if self.count_parameters:
            if self.args.plotdir is not None:
                pcounter.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots, only_overall=self.args.only_overall_plots)
            if self.parameter_out_dir is not None and not self.args.dont_write_parameters:
                pcounter.write(self.parameter_out_dir)

        glutils.remove_glfo_files(self.my_gldir, self.args.locus)
        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def add_vs_indels(self):
        vsfo = self.vs_info['annotations']
        queries_with_indels = [q for q in self.remaining_queries if q in vsfo and indelutils.has_indels(vsfo[q]['indelfo'])]
        for query in queries_with_indels:
            self.vs_indels.add(query)  # this line is to tell us that this query has an indel stemming from vsearch, while the next line tells us that there's an indel (combine_indels() gets confused if we don't differentiate between the two)
            self.info['indels'][query] = copy.deepcopy(vsfo[query]['indelfo'])
            # print indelutils.get_dbg_str(vsfo[query]['indelfo'])

        if self.debug and len(queries_with_indels) > 0:
            print '    added %d vsearch indel%s%s' % (len(queries_with_indels), utils.plural(len(queries_with_indels)), (' (%s)' % ' '.join(queries_with_indels)) if len(queries_with_indels) < 100 else '')

    # ----------------------------------------------------------------------------------------
    def subworkdir(self, iproc, n_procs):
        if n_procs == 1:
            return self.args.workdir
        else:
            return self.args.workdir + '/sw-' + str(iproc)

    # ----------------------------------------------------------------------------------------
    def execute_commands(self, base_infname, base_outfname, mismatches, gap_opens):
        start = time.time()
        def get_cmd_str(iproc):
            return self.get_ig_sw_cmd_str(self.subworkdir(iproc, n_procs), base_infname, base_outfname, mismatches[iproc], gap_opens[iproc])
            # return self.get_vdjalign_cmd_str(self.subworkdir(iproc, n_procs), base_infname, base_outfname, mismatches[iproc], gap_opens[iproc] xxx update this)

        n_procs = len(mismatches)
        cmdfos = [{'cmd_str' : get_cmd_str(iproc),
                   'workdir' : self.subworkdir(iproc, n_procs),
                   'outfname' : self.subworkdir(iproc, n_procs) + '/' + base_outfname}
                  for iproc in range(n_procs)]
        utils.run_cmds(cmdfos, batch_system=self.args.batch_system, batch_options=self.args.batch_options, batch_config_fname=self.args.batch_config_fname)

        for iproc in range(n_procs):
            os.remove(self.subworkdir(iproc, n_procs) + '/' + base_infname)
        sys.stdout.flush()
        self.ig_sw_time = time.time() - start

    # ----------------------------------------------------------------------------------------
    def split_queries_by_match_mismatch(self, input_queries, n_procs, debug=False):
        def get_query_mfreq(q):
            return self.vs_info['annotations'][q]['v_mut_freq'] if q in self.vs_info['annotations'] else self.default_mfreq
        def best_mismatch(q):
            mfreq_q = self.vs_info['annotations'][q]['v_mut_freq'] if q in self.vs_info['annotations'] else self.default_mfreq
            def keyfunc(pair):
                mf, mm = pair
                return abs(mf - mfreq_q)
            nearest_mfreq, nearest_mismatch = min(self.mfreq_mismatch_vals, key=keyfunc)  # take the optimized value whose mfreq is closest to this sequence's mfreq
            return nearest_mismatch

        query_groups = utils.group_seqs_by_value(input_queries, best_mismatch)
        mismatch_vals = [best_mismatch(queries[0]) for queries in query_groups]

        # note: it'd be nice to be able to give ig-sw a different match:mismatch for each sequence (rather than running separate procs for each match:mismatch), but it initializes a matrix using the match:mismatch values before looping over sequences, so that's probably infeasible

        if debug:
            print 'start'
            for m, queries in zip(mismatch_vals, query_groups):
                print '  %d   %d    %s' % (m, len(queries), ' '.join(queries))

        # first give one proc to each mismatch (note that this takes at least as many procs as there are different mismatch values, even if --n-procs is smaller)
        while len(mismatch_vals) < n_procs:
            if debug:
                print '%d < %d' % (len(mismatch_vals), n_procs)
            largest_group = max(query_groups, key=len)
            piece_a, piece_b = largest_group[ : len(largest_group) / 2], largest_group[len(largest_group) / 2 : ]  # split up the largest group into two [almost] equal pieces
            ilargest = query_groups.index(largest_group)
            query_groups = query_groups[:ilargest] + [piece_a, piece_b] + query_groups[ilargest + 1:]
            mismatch_vals = mismatch_vals[:ilargest] + [mismatch_vals[ilargest], mismatch_vals[ilargest]] + mismatch_vals[ilargest + 1:]
            if debug:
                for m, queries in zip(mismatch_vals, query_groups):
                    print '  %d   %d    %s' % (m, len(queries), ' '.join(queries))

        return mismatch_vals, query_groups

    # ----------------------------------------------------------------------------------------
    def split_queries_evenly_among_procs(self, input_queries, n_procs):
        query_iprocs = [(input_queries[iq], iq % n_procs) for iq in range(len(input_queries))]  # loop over queries, cycling through the procs
        queries_for_each_proc = [[q for q, i in query_iprocs if i == iproc] for iproc in range(n_procs)]  # then pull out the queries for each proc
        mismatches = [self.mismatch for _ in range(n_procs)]  # all of 'em have the same mismatch
        return mismatches, queries_for_each_proc

    # ----------------------------------------------------------------------------------------
    def get_gap_opens(self, mismatches, queries_for_each_proc):
        failed_mismatches = {}
        for mismatch, queries in zip(mismatches, queries_for_each_proc):
            fqueries_in_this_mismatch = set(queries) & self.indel_reruns
            for fquery in fqueries_in_this_mismatch:  # if it's destined for this proc, keep track of the match/mismatch info, and remove it from this proc's list
                if mismatch not in failed_mismatches:  # reminder: more than one proc can have the same mismatch
                    failed_mismatches[mismatch] = []
                failed_mismatches[mismatch].append(fquery)
                queries.remove(fquery)
        assert set([q for qlist in failed_mismatches.values() for q in qlist]) == self.indel_reruns
        self.indel_reruns.clear()

        gap_opens = [self.gap_open_penalty for _ in mismatches]
        for fmismatch, fqueries in failed_mismatches.items():
            gap_opens += [self.args.no_indel_gap_open_penalty]
            mismatches += [fmismatch]
            queries_for_each_proc += [fqueries]
            # print '  adding new iproc for %s with gap open %d mismatch %d' % (' '.join(fqueries), fmismatch, self.args.no_indel_gap_open_penalty)

        return mismatches, gap_opens, queries_for_each_proc

    # ----------------------------------------------------------------------------------------
    def split_queries(self, n_procs):
        input_queries = list(self.remaining_queries)
        if self.vs_info is None:
            mismatches, queries_for_each_proc = self.split_queries_evenly_among_procs(input_queries, n_procs)
        else:
            mismatches, queries_for_each_proc = self.split_queries_by_match_mismatch(input_queries, n_procs)

        mismatches, gap_opens, queries_for_each_proc = self.get_gap_opens(mismatches, queries_for_each_proc)  # they're all the same unless we have some indel fails

        missing_queries = self.remaining_queries - set([q for proc_queries in queries_for_each_proc for q in proc_queries])
        if len(missing_queries) > 0:
            raise Exception('didn\'t write %s to %s' % (':'.join(missing_queries), self.args.workdir))

        return mismatches, gap_opens, queries_for_each_proc

    # ----------------------------------------------------------------------------------------
    def write_input_files(self, base_infname, queries_for_each_proc):
        n_procs = len(queries_for_each_proc)
        for iproc in range(n_procs):
            workdir = self.subworkdir(iproc, n_procs)
            if n_procs > 1:
                utils.prep_dir(workdir)
            with open(workdir + '/' + base_infname, 'w') as sub_infile:
                for query_name in queries_for_each_proc[iproc]:
                    if query_name in self.info['indels']:
                        seq = self.info['indels'][query_name]['reversed_seq']  # use the query sequence with shm insertions and deletions reversed
                    else:
                        assert len(self.input_info[query_name]['seqs']) == 1  # sw can't handle multiple simultaneous sequences, but it's nice to have the same headers/keys everywhere, so we use the plural versions (with lists) even here (where "it's nice" means "it used to be the other way and it fucking sucked and a fuckton of effort went into synchronizing the treatments")
                        seq = self.input_info[query_name]['seqs'][0]
                    sub_infile.write('>%s NUKES\n%s\n' % (query_name, seq))

    # # ----------------------------------------------------------------------------------------
    # def get_vdjalign_cmd_str(self, workdir, base_infname, base_outfname, mismatch):
    #     # large gap-opening penalty: we want *no* gaps in the middle of the alignments
    #     # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
    #     cmd_str = os.getenv('HOME') + '/.local/bin/vdjalign align-fastq -q'
    #     cmd_str += ' --locus ' + self.args.locus.upper()
    #     cmd_str += ' --max-drop 50'
    #     cmd_str += ' --match ' + str(self.match_score) + ' --mismatch ' + str(mismatch)
    #     cmd_str += ' --gap-open ' + str(self.gap_open_penalty)
    #     cmd_str += ' --vdj-dir ' + self.my_gldir + '/' + self.args.locus
    #     cmd_str += ' --samtools-dir ' + self.args.partis_dir + '/packages/samtools'
    #     cmd_str += ' ' + workdir + '/' + base_infname + ' ' + workdir + '/' + base_outfname
    #     return cmd_str

    # ----------------------------------------------------------------------------------------
    def get_ig_sw_cmd_str(self, workdir, base_infname, base_outfname, mismatch, gap_open):
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        cmd_str = self.args.ig_sw_binary
        cmd_str += ' -l ' + self.args.locus.upper()
        cmd_str += ' -d 50'  # max drop
        cmd_str += ' -m ' + str(self.match_score) + ' -u ' + str(mismatch)
        cmd_str += ' -o ' + str(gap_open)
        cmd_str += ' -p ' + self.my_gldir + '/' + self.args.locus + '/'  # NOTE needs the trailing slash
        cmd_str += ' ' + workdir + '/' + base_infname + ' ' + workdir + '/' + base_outfname
        return cmd_str

    # # ----------------------------------------------------------------------------------------
    # # NOTE wrote this as a way to remove the ig-sw matches with different cigar and read lengths. It works fine, except the fact we have to do all this string parsing means it's really too slow
    # def remove_length_discrepant_matches(self, samfname):
    #     start = time.time()
    #     # Remove matches for which cigar and query are different length (in all checked cases, cigar was the party at fault).
    #     # Have to rewrite the stupid file on disk because pysam doesn't seem to be able to not crash when it encounters the offending match, and ig-sw (well, klib: https://github.com/attractivechaos/klib) is too obfuscated to make fixing it there a good use of time, at least a.t.m.
    #     sam_headers = ('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM/RNEXT', 'MPOS/PNEXT', 'ISIZE/TLEN', 'SEQ', 'QUAL', 'TAGs', 'wtf')  # what kind of loser puts their headers in their tsv file, anyway?
    #     def lget(llist, sname):
    #         return llist[sam_headers.index(sname)]

    #     with open(samfname + '.tmp', 'w') as tmpfile:
    #         with open(samfname) as samfile:
    #             seqstr = None
    #             missing_primary = False
    #             for line in samfile:
    #                 if line[0] != '@':
    #                     linelist = line.strip().split()
    #                     if len(linelist) != len(sam_headers):
    #                         raise Exception('sam file line (len %d) in %s doesn\'t correspond to known headers (len %d):\n%s\n\n%s' % (len(linelist), samfname, len(sam_headers), '\n'.join(linelist), '\n'.join(sam_headers)))
    #                     cigarstr = lget(linelist, 'CIGAR')
    #                     cigarlen = sum([int(l) for l in re.findall('[0-9][0-9]*', cigarstr)])
    #                     if lget(linelist, 'FLAG') == '0':
    #                         assert seqstr != '*'
    #                         seqstr = lget(linelist, 'SEQ')
    #                     else:
    #                         assert lget(linelist, 'SEQ')
    #                     assert seqstr is not None
    #                     if cigarlen != len(seqstr):
    #                         if lget(linelist, 'FLAG') == '0':
    #                             missing_primary = True  # this is the first line, so we're removing the only one with the sequence, so we need to convert the first one that's ok to a primary
    #                         continue
    #                 if missing_primary:
    #                     iflag = sam_headers.index('FLAG')
    #                     iseq = sam_headers.index('SEQ')
    #                     linelist = linelist[:iflag] + ['0'] + linelist[iflag + 1 : iseq] + [seqstr] + linelist[iseq + 1:]
    #                     line = '\t'.join(linelist) + '\n'
    #                     missing_primary = False
    #                 tmpfile.write(line)
    #     subprocess.check_call(['mv', samfname + '.tmp', samfname])
    #     print '        time to rewrite same file: %.2f' % (time.time() - start)

    # ----------------------------------------------------------------------------------------
    def read_output(self, base_outfname, n_procs=1):
        if self.debug:
            print '%s' % utils.color('green', 'reading output')

        queries_read_from_file = set()  # should be able to remove this, eventually
        for iproc in range(n_procs):
            outfname = self.subworkdir(iproc, n_procs) + '/' + base_outfname
            # self.remove_length_discrepant_matches(outfname)
            with contextlib.closing(pysam.Samfile(outfname)) as sam:  # changed bam to sam because ig-sw outputs sam files
                grouped = itertools.groupby(iter(sam), operator.attrgetter('qname'))
                for _, reads in grouped:  # loop over query sequences
                    try:
                        readlist = list(reads)
                    except:  # should no longer happen (was a result of pysam barfing when ig-sw gave it cigar and query sequences that were different lengths, but now ig-sw should skip matches for which that's true) it would be better if ig-sw didn't make those matches to start with, but that would require understanding a lot more about ig-sw
                        raise Exception('failed to convert sam reads')
                    qinfo = self.read_query(sam.references, readlist)
                    self.summarize_query(qinfo)  # returns before adding to <self.info> if it thinks we should rerun the query
                    queries_read_from_file.add(qinfo['name'])

        not_read = self.remaining_queries - queries_read_from_file
        if len(not_read) > 0:  # ig-sw (now) doesn't write matches for cases in which cigar and read length differ, which means there are now queries for which it finds zero matches (well, it didn't seem to happen before... but not sure that it couldn't have)
            print '\n%s didn\'t read %s from %s' % (utils.color('red', 'warning'), ' '.join(not_read), self.args.workdir)

        for iproc in range(n_procs):
            workdir = self.subworkdir(iproc, n_procs)
            os.remove(workdir + '/' + base_outfname)
            if n_procs > 1:  # still need the top-level workdir
                os.rmdir(workdir)

        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def remove_query(self, query):
        # NOTE you're iterating over a deep copy of <self.info['queries']>, right? you better be!
        del self.info[query]  # this still leaves this query's gene matches in <self.info> (there may also be other traces of it)
        self.info['queries'].remove(query)
        if query in self.info['indels']:
            del self.info['indels'][query]
        self.info['removed-queries'].add(query)

    # ----------------------------------------------------------------------------------------
    def add_dummy_d_match(self, qinfo, first_v_qr_end):
        dummy_d = glutils.dummy_d_genes[self.args.locus]
        qinfo['matches']['d'].append((1, dummy_d))
        qinfo['qrbounds'][dummy_d] = (first_v_qr_end, first_v_qr_end)
        qinfo['glbounds'][dummy_d] = (1, 1)

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
            'new_indels' : {}
        }

        last_scores = {r : None for r in utils.regions}
        for read in reads:  # loop over the matches found for each query sequence
            read.seq = qinfo['seq']  # only the first one has read.seq set by default, so we need to set the rest by hand
            gene = references[read.tid]
            region = utils.get_region(gene)
            qrbounds = (read.qstart, read.qend)
            glbounds = (read.pos, read.aend)
            score = read.tags[0][1]
            if last_scores[region] is not None and score > last_scores[region]:
                raise Exception('[sb]am file from smith-waterman not ordered by match score')
            last_scores[region] = score

            # NOTE it is very important not to ever skip the best match for a region, since we need its qrbounds later to know how much was trimmed

            if len(qinfo['matches'][region]) >= self.args.n_max_per_region[utils.regions.index(region)]:
                assert len(qinfo['matches'][region]) == self.args.n_max_per_region[utils.regions.index(region)]  # there better not be a way to get more than we asked for
                continue

            indelfo = indelutils.get_indelfo_from_cigar(read.cigarstring, qinfo['seq'], qrbounds, self.glfo['seqs'][region][gene], glbounds, {region : gene}, uid=qinfo['name'])  # note that qinfo['seq'] differs from self.input_info[qinfo['name']]['seqs'][0] if we've already reversed an indel in this sequence
            if indelutils.has_indels(indelfo):
                if len(qinfo['matches'][region]) > 0:  # skip any gene matches with indels after the first one for each region (if we want to handle [i.e. reverse] an indel, we will have stored the indel info for the first match, and we'll be rerunning)
                    continue
                assert region not in qinfo['new_indels']  # only to double-check the continue just above
                qinfo['new_indels'][region] = indelfo

            # and finally add this match's information
            qinfo['matches'][region].append((score, gene))  # NOTE it is important that this is ordered such that the best match is first
            qinfo['qrbounds'][gene] = qrbounds
            qinfo['glbounds'][gene] = glbounds

        if not utils.has_d_gene(self.args.locus) and len(qinfo['matches']['v']) > 0:
            _, first_v_match = qinfo['matches']['v'][0]
            self.add_dummy_d_match(qinfo, first_v_qr_end=qinfo['qrbounds'][first_v_match][1])

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
    def super_high_mutation(self, qinfo, best):
        gl_seqs, qr_seqs = [], []
        for region in ['v', 'j']:  # insertion mutation rates are meaningless, and we expect D to frequently be really high
            gl_bounds = qinfo['glbounds'][best[region]]
            gl_seqs.append(self.glfo['seqs'][region][best[region]][gl_bounds[0] : gl_bounds[1]])
            qr_bounds = qinfo['qrbounds'][best[region]]
            qr_seqs.append(qinfo['seq'][qr_bounds[0] : qr_bounds[1]])

        vj_mute_freq = utils.hamming_fraction(''.join(gl_seqs), ''.join(qr_seqs))
        if vj_mute_freq > self.args.max_vj_mut_freq:
            if self.debug:
                print '      rerun: super high vj mutation (%.3f > %.3f)' % (vj_mute_freq, self.args.max_vj_mut_freq)
            return True
        else:
            return False

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

        if not recursed and status == 'nonsense' and l_reg == 'd':  # on rare occasions with very high mutation, vdjalign refuses to give us a j match that's at all to the right of the d match
            assert l_reg == 'd' and r_reg == 'j'
            if debug:
                print '  %s: synthesizing d match' % qinfo['name']
            leftmost_position = min(qinfo['qrbounds'][l_gene][0], qinfo['qrbounds'][r_gene][0])
            qinfo['qrbounds'][l_gene] = (leftmost_position, leftmost_position + 1)  # swap whatever crummy nonsense d match we have now for a one-base match at the left end of things (things in practice should be left end of j match)
            qinfo['glbounds'][l_gene] = (0, 1)
            if l_reg in qinfo['new_indels']:
                del qinfo['new_indels'][l_reg]
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

        # ----------------------------------------------------------------------------------------
        def dbg_print(l_length, r_length, l_portion, r_portion):
            print '  %4d %4d      %4d %4d' % (l_length, r_length, l_portion, r_portion)
            l_gene_seq = qinfo['seq'][qinfo['qrbounds'][l_gene][0] : qinfo['qrbounds'][l_gene][1] - l_portion]
            r_gene_seq = qinfo['seq'][qinfo['qrbounds'][r_gene][0] + r_portion : qinfo['qrbounds'][r_gene][1]]
            l_offset = min(qinfo['qrbounds'][r_gene][0] + r_portion, qinfo['qrbounds'][l_gene][0])
            print '                                %s%s' % ((qinfo['qrbounds'][l_gene][0] - l_offset) * ' ', l_gene_seq)
            print '                                %s%s' % ((qinfo['qrbounds'][r_gene][0] + r_portion - l_offset) * ' ', r_gene_seq)

        # ----------------------------------------------------------------------------------------
        def check_indel_interference(l_length, r_length, l_portion, r_portion):
            if l_reg in qinfo['new_indels'] and l_portion > 0:
                ifo = qinfo['new_indels'][l_reg]
                if utils.gap_len(ifo['qr_gap_seq'][-l_portion : ]) > 0 or utils.gap_len(ifo['gl_gap_seq'][-l_portion : ]) > 0:  # this should really check that l_portion is shorter than the qr gap seq... but it _should_ be impossible
                    return True
            if r_reg in qinfo['new_indels'] and r_portion > 0:
                ifo = qinfo['new_indels'][r_reg]
                if utils.gap_len(ifo['qr_gap_seq'][ : r_portion]) > 0 or utils.gap_len(ifo['gl_gap_seq'][ : r_portion]) > 0:
                    return True
            return False

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
            dbg_print(l_length, r_length, l_portion, r_portion)
        while l_portion + r_portion < overlap:
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
            if check_indel_interference(l_length, r_length, l_portion, r_portion):
                if debug:
                    print '  failed: apportionment encountered an indel'
                return True
            if debug:
                dbg_print(l_length, r_length, l_portion, r_portion)

        if debug:
            print '      %s apportioned %d bases between %s (%d) match and %s (%d) match' % (qinfo['name'], overlap, l_reg, l_portion, r_reg, r_portion)
        assert l_portion + r_portion == overlap
        qinfo['qrbounds'][l_gene] = (qinfo['qrbounds'][l_gene][0], qinfo['qrbounds'][l_gene][1] - l_portion)
        qinfo['glbounds'][l_gene] = (qinfo['glbounds'][l_gene][0], qinfo['glbounds'][l_gene][1] - l_portion)
        qinfo['qrbounds'][r_gene] = (qinfo['qrbounds'][r_gene][0] + r_portion, qinfo['qrbounds'][r_gene][1])
        qinfo['glbounds'][r_gene] = (qinfo['glbounds'][r_gene][0] + r_portion, qinfo['glbounds'][r_gene][1])
        if l_reg in qinfo['new_indels'] and l_portion > 0:  # it would be nice to check indel consistency after doing this, but my indel consistency checkers seem to only operate on <line>s, and I don't have one of those yet
            ifo = qinfo['new_indels'][l_reg]
            ifo['qr_gap_seq'] = ifo['qr_gap_seq'][ : -l_portion]  # this should really check that l_portion is shorter than the qr gap seq... but it _should_ be impossible
            ifo['gl_gap_seq'] = ifo['gl_gap_seq'][ : -l_portion]
            if debug:
                print '    removed %d base%s from right side of %s indel gap seqs' % (l_portion, utils.plural(l_portion), l_reg)
        if r_reg in qinfo['new_indels'] and r_portion > 0:
            ifo = qinfo['new_indels'][r_reg]
            ifo['qr_gap_seq'] = ifo['qr_gap_seq'][r_portion : ]
            ifo['gl_gap_seq'] = ifo['gl_gap_seq'][r_portion : ]
            if debug:
                print '    removed %d base%s from left side of %s indel gap seqs' % (r_portion, utils.plural(r_portion), r_reg)

        return False

    # ----------------------------------------------------------------------------------------
    def remove_probably_spurious_deletions(self, qinfo, best, debug=False):  # remove probably-spurious v_5p and j_3p deletions
        for erosion in utils.effective_erosions:
            region = erosion[0]
            if region == 'v':
                d_len = qinfo['glbounds'][best['v']][0]  # NOTE corresponds to *best* match (but it's what we subtract from everybody's bounds)
                insertion = qinfo['seq'][ : qinfo['qrbounds'][best['v']][0]]
            elif region == 'j':
                d_len = len(self.glfo['seqs']['j'][best['j']]) - qinfo['glbounds'][best['j']][1]  # NOTE corresponds to *best* match (but it's what we subtract from everybody's bounds)
                insertion = qinfo['seq'][qinfo['qrbounds'][best['j']][1] : ]
            else:
                assert False

            if d_len == 0 or len(insertion) == 0:
                continue
            if len(insertion) >= d_len:
                spurious_bases = insertion[len(insertion) - d_len:]  # i.e. the last <v_5p_del> bases of the fv insertion (or first <j_3p_del> bases of the jf insertion)
            else:
                spurious_bases = insertion  # the whole damn thing
            if spurious_bases.count(utils.ambiguous_bases[0]) == len(spurious_bases):  # don't do it if it's all Ns
                continue
            if debug:
                print 'EXPANDING %s %s d_len: %d   insertion: %s (len %d)   spurious: %s (len %d)' % (qinfo['name'], region, d_len, insertion, len(insertion), spurious_bases, len(spurious_bases))
            for gene in [g for g in qinfo['glbounds'] if utils.get_region(g) == region]:  # it's ok to assume glbounds and qrbounds have the same keys
                if region == 'v':
                    qinfo['glbounds'][gene] = (qinfo['glbounds'][gene][0] - len(spurious_bases), qinfo['glbounds'][gene][1])
                    qinfo['qrbounds'][gene] = (qinfo['qrbounds'][gene][0] - len(spurious_bases), qinfo['qrbounds'][gene][1])
                elif region == 'j':
                    qinfo['glbounds'][gene] = (qinfo['glbounds'][gene][0], qinfo['glbounds'][gene][1] + len(spurious_bases))
                    qinfo['qrbounds'][gene] = (qinfo['qrbounds'][gene][0], qinfo['qrbounds'][gene][1] + len(spurious_bases))
                else:
                    assert False

    # ----------------------------------------------------------------------------------------
    def convert_qinfo(self, qinfo, best, codon_positions):
        """ convert <qinfo> (which is from reading sam files) to format for <self.info> (this is so add_to_info() can be used by the cache file reader, as well) """
        qname = qinfo['name']
        assert qname not in self.info

        infoline = {}
        infoline['unique_ids'] = [qname, ]  # redundant, but used somewhere down the line
        infoline['seqs'] = [qinfo['seq'], ]  # NOTE this is the seq output by vdjalign, i.e. if we reversed any indels it is the reversed sequence, also NOTE many, many things depend on this list being of length one
        infoline['input_seqs'] = [self.input_info[qname]['seqs'][0], ]

        # erosion, insertion, mutation info for best match
        infoline['v_5p_del'] = qinfo['glbounds'][best['v']][0]
        infoline['v_3p_del'] = len(self.glfo['seqs']['v'][best['v']]) - qinfo['glbounds'][best['v']][1]  # len(germline v) - gl_match_end
        infoline['d_5p_del'] = qinfo['glbounds'][best['d']][0]
        infoline['d_3p_del'] = len(self.glfo['seqs']['d'][best['d']]) - qinfo['glbounds'][best['d']][1]
        infoline['j_5p_del'] = qinfo['glbounds'][best['j']][0]
        infoline['j_3p_del'] = len(self.glfo['seqs']['j'][best['j']]) - qinfo['glbounds'][best['j']][1]

        infoline['fv_insertion'] = qinfo['seq'][ : qinfo['qrbounds'][best['v']][0]]
        infoline['vd_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['v']][1] : qinfo['qrbounds'][best['d']][0]]
        infoline['dj_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['d']][1] : qinfo['qrbounds'][best['j']][0]]
        infoline['jf_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['j']][1] : ]

        if qname in self.info['indels']:  # NOTE at this piont indel info isn't updated for any change in the genes (the indelfo gets updated during utils.add_implicit_info())
            infoline['indelfos'] = [self.info['indels'][qname]]  # NOTE this makes it so that self.info[uid]['indelfos'] *is* self.info['indels'][uid]. It'd still be nicer to eventually do away with self.info['indels'], although I'm not sure that's really either feasible or desirable given other constraints
        else:
            infoline['indelfos'] = [indelutils.get_empty_indel()]

        infoline['duplicates'] = [self.duplicates.get(qname, []), ]  # note that <self.duplicates> doesn't handle simultaneous seqs, i.e. it's for just a single sequence

        infoline['all_matches'] = {r : [g for _, g in qinfo['matches'][r]] for r in utils.regions}  # get lists with no scores, just the names (still ordered by match quality, though)
        for region in utils.regions:
            infoline[region + '_gene'] = best[region]

        if self.args.linearham:
            sortmatches = {r : [g for _, g in qinfo['matches'][r]] for r in utils.regions}
            infoline['flexbounds'] = {}
            def span(boundlist):
                return (min(boundlist) - 2, max(boundlist) + 2)
            for region in utils.regions:
                bounds_l, bounds_r = zip(*[qinfo['qrbounds'][g] for g in sortmatches[region]])  # left- (and right-) bounds for each gene
                bounds_l = [bound - len(infoline['fv_insertion']) for bound in bounds_l]  # the left-bounds need to be adjusted for V 5' framework insertions
                bounds_r = [bound - len(infoline['fv_insertion']) for bound in bounds_r]  # the right-bounds need to be adjusted for V 5' framework insertions
                infoline['flexbounds'][region] = {'l' : span(bounds_l), 'r' : span(bounds_r)}
            infoline['relpos'] = {gene: qinfo['qrbounds'][gene][0] - glbound[0] - len(infoline['fv_insertion']) for gene, glbound in qinfo['glbounds'].items()}  # position in the query sequence of the start of each uneroded germline match (the relpos needs to be adjusted for V 5' framework insertions)

        infoline['cdr3_length'] = codon_positions['j'] - codon_positions['v'] + 3
        infoline['codon_positions'] = codon_positions

        infoline['padlefts'] = ['']
        infoline['padrights'] = ['']

        return infoline

    # ----------------------------------------------------------------------------------------
    def check_simulation_kbounds(self, line, true_line):
        true_kbounds = {'v' : {}, 'd' : {}}
        true_kbounds['v']['best'] = true_line['regional_bounds']['v'][1]
        true_kbounds['d']['best'] = true_line['regional_bounds']['d'][1] - true_line['regional_bounds']['v'][1]

        def print_kbound_warning():
            print '  %s true kset (%s) not within kbounds (%s) for %s' % (utils.color('red', 'warning'), utils.kbound_str(true_kbounds), utils.kbound_str({r : line['k_' + r] for r in true_kbounds}), ':'.join(line['unique_ids']))

        for region in true_kbounds:
            if true_kbounds[region]['best'] < line['k_' + region]['min'] or true_kbounds[region]['best'] >= line['k_' + region]['max']:
                print_kbound_warning()

    # ----------------------------------------------------------------------------------------
    def add_to_info(self, infoline):
        assert len(infoline['unique_ids'])
        qname = infoline['unique_ids'][0]

        self.info['passed-queries'].add(qname)
        self.info[qname] = infoline

        # add this query's matches into the overall gene match sets
        for region in utils.regions:
            self.info['all_best_matches'].add(self.info[qname][region + '_gene'])
            self.info['all_matches'][region] |= set(self.info[qname]['all_matches'][region])  # NOTE there's an 'all_matches' in this query's info, and also in <self.info>

        # everything this is flagging seems to be cases where the insertion (or the overlapping germline region) is the same as the germline base, i.e. the kbound getter is working
        # if not self.args.is_data:
        #     self.check_simulation_kbounds(self.info[qname], self.reco_info[qname])

        if self.debug:
            inf_label = ' ' + utils.kbound_str({r : infoline['k_' + r] for r in ['v', 'd']})
            if not self.args.is_data:
                inf_label = 'inf: ' + inf_label
                utils.print_reco_event(self.reco_info[qname], extra_str='    ', label=utils.color('green', 'true:'))
            utils.print_reco_event(self.info[qname], extra_str='    ', label=inf_label)

        if not utils.is_functional(self.info[qname]):
            self.kept_unproductive_queries.add(qname)
        self.remaining_queries.remove(qname)

    # ----------------------------------------------------------------------------------------
    def summarize_query(self, qinfo):
        """
        Fiddle with a few things, but mostly decide whether we're satisfied with the current matches.
        If not, return without calling add_to_info().
        """
        def dbgfcn(dbgstr):  # doesn't return anything, but it makes things more concise
            if self.debug:
                print '      rerun: %s' % dbgstr

        qname = qinfo['name']
        qseq = qinfo['seq']
        assert qname not in self.info
        if self.debug:
            print '  %s' % qname

        # do we have a match for each region?
        for region in utils.getregions(self.args.locus):
            if len(qinfo['matches'][region]) == 0:
                return dbgfcn('no %s match' % region)  # if no d match found, maybe we should just assume entire d was eroded?

        best = {r : qinfo['matches'][r][0][1] for r in utils.regions}  # already made sure there's at least one match for each region

        # s-w allows d and j matches to overlap, so we need to apportion the disputed bases
        for rpair in utils.region_pairs():
            overlap_status = self.check_boundaries(rpair, qinfo, best)  # I think all this overlap fixing only adjusts the bounds for the best match, but I guess that's ok? It's certainly been like that for ages
            if overlap_status == 'overlap':
                overlap_indel_fail = self.shift_overlapping_boundaries(rpair, qinfo, best)  # this is kind of a crappy way to return the information, but I can't think of anything better a.t.m.
                if overlap_indel_fail:
                    self.indel_reruns.add(qname)
                    return dbgfcn('overlap/indel fails')
            elif overlap_status == 'nonsense':
                return dbgfcn('nonsense overlap bounds')
            else:
                assert overlap_status == 'ok'

        if len(qinfo['new_indels']) > 0:  # if any of the best matches had new indels this time through
            if qname in self.info['indels']:  # no really necessary, but I'm nervous about it because I screwed it up once before
                assert qname in self.vs_indels
            indelfo = self.combine_indels(qinfo, best)  # the next time through, when we're writing ig-sw input, we look to see if each query is in <self.info['indels']>, and if it is we pass ig-sw the indel-reversed sequence, rather than the <input_info> sequence
            if indelfo is None:
                self.indel_reruns.add(qname)
                return dbgfcn('indel fails')
            else:
                self.info['indels'][qinfo['name']] = indelfo
                self.indel_reruns.add(qname)
                return dbgfcn(' new indels in %s' % ' '.join(qinfo['new_indels'].keys()))  # utils.pad_lines(indelutils.get_dbg_str(self.info['indels'][qinfo['name']]), 10)

        if self.debug >= 2:
            for region in utils.regions:
                for score, gene in qinfo['matches'][region]:  # sorted by decreasing match quality
                    self.print_match(region, gene, score, qseq, qinfo['glbounds'][gene], qinfo['qrbounds'][gene], skipping=False)

        if self.super_high_mutation(qinfo, best):
            return dbgfcn('super high mutation')

        # force v 5p and j 3p matches to (in most cases) go to the end of the sequence
        self.remove_probably_spurious_deletions(qinfo, best)

        # check for suspiciously bad annotations
        for rp in utils.region_pairs():
            insertion_length = qinfo['qrbounds'][best[rp['right']]][0] - qinfo['qrbounds'][best[rp['left']]][1]  # start of right match minus end of left one
            if insertion_length > self.absolute_max_insertion_length:
                return dbgfcn('suspiciously long insertion')

        # set conserved codon positions
        codon_positions = {}
        for region, codon in utils.conserved_codons[self.args.locus].items():
            # position within original germline gene, minus the position in that germline gene at which the match starts, plus the position in the query sequence at which the match starts
            pos = self.glfo[codon + '-positions'][best[region]] - qinfo['glbounds'][best[region]][0] + qinfo['qrbounds'][best[region]][0]
            if pos < 0 or pos >= len(qseq):
                return dbgfcn('invalid %s codon position (%d in seq of length %d), rerunning' % (codon, pos, len(qseq)))
            codon_positions[region] = pos

        # check cdr3 length
        cdr3_length = codon_positions['j'] - codon_positions['v'] + 3
        if cdr3_length < 6:  # NOTE six is also hardcoded in utils
            return dbgfcn('negative cdr3 length %d' % cdr3_length)

        # convert to regular format used elsewhere, and add implicit info
        infoline = self.convert_qinfo(qinfo, best, codon_positions)
        try:
            utils.add_implicit_info(self.glfo, infoline, aligned_gl_seqs=self.aligned_gl_seqs, reset_indel_genes=True)
        except:  # gah, I don't really like just swallowing everything... but then I *expect* it to fail here, and when I call it elsewhere, where I don't expect it to fail, shit doesn't get swallowed (and I want all the different exceptions, i.e. I kind of do want to just swallow everything here)
            exc_type, exc_value, exc_traceback = sys.exc_info()
            lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
            print utils.pad_lines(''.join(lines))
            print '      rerun: implicit info adding failed for %s (see above), rerunning' % qname  # shouldn't be able to happen, so print even if debug isn't set
            return dbgfcn('see above')

        # deal with unproductive rearrangements
        if not utils.is_functional(infoline):
            if self.args.skip_unproductive:
                if self.debug:
                    print '      skipping unproductive (%s)' % utils.is_functional_dbg_str(infoline)
                self.skipped_unproductive_queries.add(qname)
                self.remaining_queries.remove(qname)
                return
            else:
                pass  # this is here so you don't forget that if neither of the above is true, we fall through and add the query to self.info

        kbounds = self.get_kbounds(infoline, qinfo, best)  # gets the boundaries of the non-best matches from <qinfo>
        if kbounds is None:
            return dbgfcn('nonsense kbounds')
        infoline['k_v'] = kbounds['v']
        infoline['k_d'] = kbounds['d']

        self.add_to_info(infoline)

    # ----------------------------------------------------------------------------------------
    def get_kbounds(self, line, qinfo, best, debug=False):
        # NOTE
        #  - k_v is index of first d/dj insert base (i.e. length of v match)
        #  - k_v + k_d is index of first j/dj insert base (i.e. length of v + vd insert + d match)

        best_k_v = line['regional_bounds']['v'][1]  # end of v match
        best_k_d = line['regional_bounds']['d'][1] - line['regional_bounds']['v'][1]  # end of d minus end of v

        k_v_min, k_v_max = best_k_v, best_k_v
        k_d_min, k_d_max = best_k_d, best_k_d

        n_matched_to_break = 4  # once we see this many consecutive unmutated bases, assume we're actually within a germline v/d/j match

        # first make sure the hmm will check for cases in which sw over-expanded v
        if debug:
            print 'k_v_min --- %d' % k_v_min
        n_matched = 0  # if <n_matched_to_break> bases match, assume we're definitely within the germline
        qrseq = line['fv_insertion'] + line['v_qr_seqs'][0]
        glseq = line['fv_insertion'] + line['v_gl_seq']
        if k_v_min > len(qrseq):
            if debug:
                print 'k_v_min too big %d %d' % (k_v_min, len(qrseq))
            return None
        icheck = k_v_min
        while icheck > line['codon_positions']['v'] + 3:  # i.e. stop when the last v base is the last base of the cysteine
            icheck -= 1
            if debug:
                print '    check %d' % icheck
            if qrseq[icheck] == glseq[icheck]:
                n_matched += 1
                if debug:
                    print '      match number %d' % n_matched
                if n_matched >= n_matched_to_break:
                    if debug:
                        print '      break on %d matches' % n_matched_to_break
                    break
            else:  # set <k_v_min> to <icheck>, i.e. move first possible d/dj insert base to index <icheck>
                if debug:
                    print '      set k_v_min to %d' % icheck
                k_v_min = icheck
                n_matched = 0


        # then make sure the hmm will check for cases in which sw over-expanded j
        if debug:
            print 'k_d_max --- %d' % k_d_max
        k_v_min_CHK = k_v_min  # make sure we don't accidentally change it
        n_matched = 0  # if <n_matched_to_break> bases match, assume we're definitely within the germline
        qrseq, glseq = line['j_qr_seqs'][0], line['j_gl_seq']
        j_start = line['regional_bounds']['j'][0]
        icheck = k_d_max  # <icheck> is what we're considering changing k_d_max to
        # dammit, this is wrong, should be something other that len(qrseq), but oh well, it's still handled below (and it'll be much easier to really fix when I get sequence that actually causes the crash so I can dissect it)
        # while k_v_min + icheck + 1 < min(line['codon_positions']['j'], len(qrseq)):  # usually: stop when the first j/dj base is the first base of the tryp (even for the smallest k_v), but add min() with len(qrseq) to avoid extremely rare pathologies
        while k_v_min + icheck + 1 < line['codon_positions']['j']:  # usually: stop when the first j/dj base is the first base of the tryp (even for the smallest k_v), but add min() with len(qrseq) to avoid extremely rare pathologies
            icheck += 1
            if debug:
                print '    check %d' % icheck,
            if k_v_min + icheck - j_start >= len(qrseq) or k_v_min + icheck - j_start >= len(glseq):  # shouldn't happen any more (min() in while statement should avoid it), but just in case
                print '    %s for query %s: k_v_min + icheck - j_start = %d + %d - %d = %d [ < 0 or >= len(qrseq) = %d   or   len(glseq) = %d]' % (utils.color('red', 'warning'), ' '.join(line['unique_ids']), k_v_min, icheck, j_start, k_v_min + icheck - j_start, len(qrseq), len(glseq))
                return None
            if k_v_min + icheck < j_start:  # make sure we're at least as far as the start of the j (i.e. let the hmm arbitrate between d match and dj insertion)
                if debug:
                    print '      not yet to start of j (%d + %d < %d)' % (k_v_min, icheck, j_start)
                k_d_max = icheck + 1
            elif qrseq[k_v_min + icheck - j_start] == glseq[k_v_min + icheck - j_start]:
                n_matched += 1
                if debug:
                    print '      match number %d' % n_matched
                if n_matched >= n_matched_to_break:
                    if debug:
                        print '      break on %d matches' % n_matched_to_break
                    break
            else:  # set <k_d_max> to <icheck>, i.e. move first possible d/dj insert base to index <icheck>
                if debug:
                    print '      set k_d_max to %d' % (icheck + 1)
                k_d_max = icheck + 1  # i.e. if <icheck> doesn't match, we want the first d/dj base to be the *next* one (whereas for k_v, if <icheck> didn't match, we wanted <icheck> to be the first d/vd match)
                n_matched = 0

        assert k_v_min_CHK == k_v_min

        if debug:
            print 'k_d_min --- %d' % k_d_min
        k_v_max_CHK = k_v_max  # make sure we don't accidentally change it
        n_matched = 0  # if <n_matched_to_break> bases match, assume we're definitely within the germline
        qrseq, glseq = line['d_qr_seqs'][0], line['d_gl_seq']
        d_start = line['regional_bounds']['d'][0]
        icheck = k_d_min  # <icheck> is what we're considering changing k_d_min to
        while icheck > d_start - k_v_max + 1:
            icheck -= 1
            if debug:
                print '    check %d' % icheck
            if qrseq[k_v_max + icheck - d_start] == glseq[k_v_max + icheck - d_start]:  # note that it's <k_v_max> here, since even with the longest k_v we want to make sure to check as short of a d as we mean to
                n_matched += 1
                if debug:
                    print '      match number %d' % n_matched
                if n_matched >= n_matched_to_break:
                    if debug:
                        print '      break on %d matches' % n_matched_to_break
                    break
            else:  # set <k_d_min> to <icheck>, i.e. move first possible d/dj insert base to index <icheck>
                if debug:
                    print '      set k_d_min to %d' % icheck
                k_d_min = icheck
                n_matched = 0

        assert k_v_max_CHK == k_v_max

        # NOTE I think there isn't any reason to increase k_v_max or decrease k_d_min -- sw already pretty much expands the v and j matches as far as it can (i.e. these fuzzing algorithms are mostly trying to decide if the last region to be matched by sw (d) should've been given more bases)
        # UPDATE not so sure, maybe I should, but it's not so important

        # also take the OR of k-space for the best <self.args.n_max_per_region[ireg]> matches (note that this happens after the expanding above because that expanding can break if it gets bounds from other matches)

        # k_v
        tmpreg = 'v'
        n_v_genes = min(self.args.n_max_per_region[utils.regions.index(tmpreg)], len(qinfo['matches'][tmpreg]))  # NOTE don't need the n_max_per_region thing any more
        for igene in range(n_v_genes):
            _, gene = qinfo['matches'][tmpreg][igene]
            this_k_v = qinfo['qrbounds'][gene][1]  # NOTE even if the v match doesn't start at the left hand edge of the query sequence, we still measure k_v from there. In other words, sw doesn't tell the hmm about it
            k_v_min = min(this_k_v, k_v_min)
            k_v_max = max(this_k_v, k_v_max)
            if debug:
                print '    %s %d %d' % (utils.color_gene(gene, width=15), k_v_min, k_v_max)

        # k_d
        tmpreg = 'd'
        n_d_genes = min(self.args.n_max_per_region[utils.regions.index(tmpreg)], len(qinfo['matches'][tmpreg]))  # NOTE don't need the n_max_per_region thing any more
        for igene in range(n_d_genes):
            _, gene = qinfo['matches'][tmpreg][igene]
            this_k_d = qinfo['qrbounds'][gene][1] - qinfo['qrbounds'][best['v']][1]  # end of d minus end of first/best v
            k_d_min = min(max(1, this_k_d), k_d_min)
            k_d_max = max(this_k_d, k_d_max)
            if debug:
                print '    %s %d %d' % (utils.color_gene(gene, width=15), k_d_min, k_d_max)


        # ----------------------------------------------------------------------------------------
        # switch to usual indexing conventions, i.e. that start with min and do not include max NOTE would be clearer to do this more coherently
        # NOTE i.e. k_[vd]_max means different things before here and after here
        # ----------------------------------------------------------------------------------------
        k_v_max += 1
        k_d_max += 1

        if not utils.has_d_gene(self.args.locus):
            best_k_d = 1
            k_d_min = 1
            k_d_max = 2

        if best_k_v < k_v_min or best_k_v > k_v_max or best_k_d < k_d_min or best_k_d > k_d_max:
            print '  %s inconsistent best kset for %s (v: %d (%d %d)  d: %d (%d %d)' % (utils.color('red', 'error'), qinfo['name'], best_k_v, k_v_min, k_v_max, best_k_d, k_d_min, k_d_max)
            return None
        if k_v_min <= 0 or k_d_min <= 0 or k_v_min >= k_v_max or k_d_min >= k_d_max:
            print '  %s nonsense k bounds for %s (v: %d %d  d: %d %d)' % (utils.color('red', 'error'), qinfo['name'], k_v_min, k_v_max, k_d_min, k_d_max)
            return None
        if not self.args.dont_remove_framework_insertions and self.args.is_data and k_v_min - len(line['fv_insertion']) < 0:  # it happened. like, once
            print '%s trimming fwk insertion would take k_v min to less than zero for %s: %d - %d = %d' % (utils.color('yellow', 'warning'), ' '.join(line['unique_ids']), k_v_min, len(line['fv_insertion']), k_v_min - len(line['fv_insertion']))
            return None

        kbounds = {'v' : {'best' : best_k_v, 'min' : k_v_min, 'max' : k_v_max},
                   'd' : {'best' : best_k_d, 'min' : k_d_min, 'max' : k_d_max}}
        return kbounds

    # ----------------------------------------------------------------------------------------
    def remove_framework_insertions(self, debug=False):
        for query in self.info['queries']:
            swfo = self.info[query]
            assert len(swfo['seqs']) == 1

            if debug:
                print '  %-12s' % swfo['unique_ids'][0]
                print '     fv  %s' % utils.color('blue', swfo['fv_insertion'])
                print '     jf  %s' % utils.color('blue', swfo['jf_insertion'])

            fv_len = len(swfo['fv_insertion'])
            jf_len = len(swfo['jf_insertion'])
            if fv_len == 0 and jf_len == 0:
                continue

            # it would be nice to combine these shenanigans with their counterparts in pad_seqs_to_same_length() [e.g. add a fcn utils.modify_fwk_insertions(), although padding doesn't just modify the fwk insertions, so...], but I'm worried they need to be slightly different and don't want to test extensively a.t.m.
            for seqkey in ['seqs', 'input_seqs']:
                swfo[seqkey][0] = swfo[seqkey][0][fv_len : len(swfo[seqkey][0]) - jf_len]
            swfo['naive_seq'] = swfo['naive_seq'][fv_len : len(swfo['naive_seq']) - jf_len]
            if query in self.info['indels']:
                assert self.info['indels'][query] is swfo['indelfos'][0]  # it would be nice to eventually not need the dict in both these places
                indelutils.trim_indel_info(swfo, 0, swfo['fv_insertion'], swfo['jf_insertion'], 0, 0)
            for key in swfo['k_v']:
                swfo['k_v'][key] -= fv_len
            swfo['fv_insertion'] = ''
            swfo['jf_insertion'] = ''
            swfo['codon_positions']['v'] -= fv_len
            swfo['codon_positions']['j'] -= fv_len
            for region in utils.regions:
                swfo['regional_bounds'][region] = tuple([rb - fv_len for rb in swfo['regional_bounds'][region]])  # I kind of want to just use a list now, but a.t.m. don't much feel like changing it everywhere else

            if debug:
                print '    after %s' % swfo['seqs'][0]

            # *sigh* not super happy about it, but I think the best way to handle this is to also remove these bases from the simulation info
            if self.reco_info is not None:
                raise Exception('needs fixing (and maybe actually shouldn\'t be fixed) -- i.e. you probably should have turned off fwk insertion removal')
                # simfo = self.reco_info[query]
                # utils.remove_all_implicit_info(simfo)
                # simfo['seqs'][0] = simfo['seqs'][0][fv_len : len(simfo['seqs'][0]) - jf_len]
                # if indelutils.has_indels(simfo['indelfos'][0]):
                #     simfo['indelfos'][0]['reversed_seq'] = simfo['seqs'][0]
                # for indel in reversed(simfo['indelfos'][0]['indels']):
                #     indel['pos'] -= fv_len
                # simfo['fv_insertion'] = ''
                # simfo['jf_insertion'] = ''
                # utils.add_implicit_info(self.glfo, simfo)

    # ----------------------------------------------------------------------------------------
    def remove_duplicate_sequences(self, debug=False):
        # ----------------------------------------------------------------------------------------
        def getseq(uid):
            return_seq = self.info[uid]['input_seqs'][0]
            # return_seq = utils.get_seq_with_indels_reinstated(self.info[uid])
            # if return_seq not in self.input_info[uid]['seqs'][0]:  # make sure we reinstated the indels properly
            #     print '%s reinstated seq not in input sequence:\n    %s\n    %s' % (utils.color('yellow', 'warning'), reinstated_seq, self.input_info[uid]['seqs'][0])
            return return_seq
        # ----------------------------------------------------------------------------------------
        def get_key_seq(uid):  # return the sequence which will serve as the key for <uid>
            seq = getseq(uid)
            if not self.args.also_remove_duplicate_sequences_with_different_lengths:  # this is probably significantly slower, otherwise I might make it the default
                return seq
            else:
                for kseq, kuids in seqs_to_keep.items():
                    if self.info[kuids[0]]['cdr3_length'] != self.info[uid]['cdr3_length']:
                        continue
                    if seq in kseq or kseq in seq:  # note that this keeps the first one we happen to come across -- it'd be better to keep the longest one, but this is fine for now
                        if debug:
                            print '      using keyseq from %s instead of %s' % (kuids, uid)
                        return kseq
                return seq  # if we fall through, it didn't match anybody

        uids_to_pre_keep = set()  # add these uids/seqs before looping through all the queries
        if self.args.seed_unique_id is not None:
            uids_to_pre_keep.add(self.args.seed_unique_id)
        if self.args.queries is not None:
            uids_to_pre_keep |= set(self.args.queries)
        if self.args.queries_to_include is not None:  # note that the seed seq is added to queries_to_include, which is good if we're parameter caching since seed_unique_id will be set to None
            uids_to_pre_keep |= set(self.args.queries_to_include)

        if debug and len(uids_to_pre_keep) > 0:
            print 'pre-keeping %s' % ' '.join(uids_to_pre_keep)
            print '  checking pre-keepers'
        seqs_to_keep = {}  # seq : [uids that correspond to seq]
        for utpk in uids_to_pre_keep:
            if utpk not in self.info:
                print 'requested uid %s not in sw info (probably failed above)' % utpk
                continue
            keyseq = get_key_seq(utpk)
            if keyseq not in seqs_to_keep:  # NOTE it's kind of weird to have more than one uid for a sequence, but it probably just means the user specified some duplicate sequences with --queries
                seqs_to_keep[keyseq] = []
                if debug:
                    print '      new key for %s' % utpk
            seqs_to_keep[keyseq].append(utpk)
            if debug:
                if len(seqs_to_keep[keyseq]) > 1:
                    print '    add to existing key: %s' % ' '.join(seqs_to_keep[keyseq])

        if debug:
            print '  checking non-pre-kept'
        n_kept, n_removed = len(uids_to_pre_keep), 0
        for uid in copy.deepcopy(self.info['queries']):
            if uid in uids_to_pre_keep:
                continue
            keyseq = get_key_seq(uid)
            if keyseq in seqs_to_keep:
                seqs_to_keep[keyseq].append(uid)
                self.remove_query(uid)
                n_removed += 1
                if debug:
                    print '    removing %s' % uid
            else:
                seqs_to_keep[keyseq] = [uid]
                n_kept += 1
                if debug:
                    print '    new key for %s' % uid

        for seq, uids in seqs_to_keep.items():
            assert uids[0] in self.info  # the first one should've been the one we kept
            self.info[uids[0]]['duplicates'][0] = list(set(self.info[uids[0]]['duplicates'][0] + uids[1:]))  # this is wasteful (and overly verbose), but I just want to hack in a way to remove the duplicated duplicates (I think they're sneaking in from multiple waterer runs) without screwing anything up
            self.duplicates[uids[0]] = self.info[uids[0]]['duplicates'][0]  # just copies over from previous line. Note that <self.duplicates> is really just so partitiondriver can pass in previous duplicates, but then now that we have the information in two places we need to keep it synchronized

        if n_removed > 0:
            print '      removed %d / %d = %.2f duplicate sequences after trimming framework insertions (leaving %d)' % (n_removed, n_removed + n_kept, n_removed / float(n_removed + n_kept), len(self.info['queries']))

    # ----------------------------------------------------------------------------------------
    def get_padding_parameters(self, debug=False):
        padnames = ['gl_cpos', 'gl_cpos_to_j_end']

        def get_empty_maxima():
            return {pn : None for pn in padnames}

        def check_set_maxima(name, val, cdr3):
            if maxima[name] is None or val > maxima[name]:
                maxima[name] = val
            if cdr3 not in per_cdr3_maxima:
                per_cdr3_maxima[cdr3] = get_empty_maxima()
            if per_cdr3_maxima[cdr3][name] is None or val > per_cdr3_maxima[cdr3][name]:
                per_cdr3_maxima[cdr3][name] = val

        maxima = get_empty_maxima()
        per_cdr3_maxima = {}
        for query in self.info['queries']:
            swfo = self.info[query]

            # find biggest cyst position among all gl matches (NOTE this pads more than it really needs to -- it only needs to be the max cpos over genes that have this cdr3 length)
            fvstuff = max(0, len(swfo['fv_insertion']) - swfo['v_5p_del'])  # we always want to pad out to the entire germline sequence, so don't let this go negative
            # loop over all matches for all sequences (up to n_max_per_region), because we want bcrham to be able to compare any sequence to any other (although, could probably use all *best* matches rather than all *all* UPDATE no, I kinda think not)
            for v_match in self.info['all_matches']['v']:
                gl_cpos = self.glfo['cyst-positions'][v_match] + fvstuff
                check_set_maxima('gl_cpos', gl_cpos, swfo['cdr3_length'])

            # Since we only store j_3p_del for the best match, we can't loop over all of 'em. But j stuff doesn't vary too much, so it works ok.
            cpos = swfo['codon_positions']['v']  # cyst position in query sequence (as opposed to gl_cpos, which is in germline allele)
            # jfstuff = max(0, len(swfo['jf_insertion']) - swfo['j_3p_del'])  # I'm not really sure why what this was for -- maybe I needed it when fwk insertion trimming was before/after this? -- in any case I'm pretty sure it's wrong to include it now
            gl_cpos_to_j_end = len(swfo['seqs'][0]) - cpos + swfo['j_3p_del']  # + jfstuff
            check_set_maxima('gl_cpos_to_j_end', gl_cpos_to_j_end, swfo['cdr3_length'])

        if debug:
            print '  maxima:',
            for k in padnames:
                print '%s %d    ' % (k, maxima[k]),
            print ''

            print '    per-cdr3 maxima:'
            print '         %s  %s  %s' % ('cdr3', padnames[0], padnames[1])
            for cdr3 in per_cdr3_maxima:
                print '         %3d' % cdr3,
                for k in padnames:
                    print '   %d    ' % (per_cdr3_maxima[cdr3][k]),
                print ''

        return maxima, per_cdr3_maxima

    # ----------------------------------------------------------------------------------------
    def pad_seqs_to_same_length(self, debug=False):
        """
        Pad all sequences in <seqinfo> to the same length to the left and right of their conserved cysteine positions.
        Next, pads all sequences further out (if necessary) such as to eliminate all v_5p and j_3p deletions.
        """

        cluster_different_cdr3_lengths = False  # if you want glomerator.cc to try to cluster different cdr3 lengths, you need to pass it *everybody* with the same N padding... but then you're padding way more than you need to on almost every sequence, which is really wasteful and sometimes confuses bcrham

        if debug:
            print 'padding %d seqs to same length (%s cdr3 length classes)' % (len(self.info['queries']), 'within' if not cluster_different_cdr3_lengths else 'merging')

        maxima, per_cdr3_maxima = self.get_padding_parameters(debug=debug)

        if debug:
            print '    left  right    uid'
        for query in self.info['queries']:
            swfo = self.info[query]
            assert len(swfo['seqs']) == 1

            cpos = swfo['codon_positions']['v']
            if cluster_different_cdr3_lengths:
                padleft = maxima['gl_cpos'] - cpos  # left padding: biggest germline cpos minus cpos in this sequence
                padright = maxima['gl_cpos_to_j_end'] - (len(swfo['seqs'][0]) - cpos)
            else:
                padleft = per_cdr3_maxima[swfo['cdr3_length']]['gl_cpos'] - cpos  # left padding: biggest germline cpos minus cpos in this sequence
                padright = per_cdr3_maxima[swfo['cdr3_length']]['gl_cpos_to_j_end'] - (len(swfo['seqs'][0]) - cpos)
            if padleft < 0 or padright < 0:
                raise Exception('bad padding %d %d for %s' % (padleft, padright, query))

            leftstr = padleft * utils.ambiguous_bases[0]
            rightstr = padright * utils.ambiguous_bases[0]
            swfo['fv_insertion'] = leftstr + swfo['fv_insertion']
            swfo['jf_insertion'] = swfo['jf_insertion'] + rightstr
            for seqkey in ['seqs', 'input_seqs']:
                swfo[seqkey][0] = leftstr + swfo[seqkey][0] + rightstr
            swfo['naive_seq'] = leftstr + swfo['naive_seq'] + rightstr
            if query in self.info['indels']:
                assert self.info['indels'][query] is swfo['indelfos'][0]
                indelutils.pad_indelfo(swfo['indelfos'][0], leftstr, rightstr)
            for key in swfo['k_v']:
                swfo['k_v'][key] += padleft
            swfo['codon_positions']['v'] += padleft
            swfo['codon_positions']['j'] += padleft
            for region in utils.regions:
                swfo['regional_bounds'][region] = tuple([rb + padleft for rb in swfo['regional_bounds'][region]])  # I kind of want to just use a list now, but a.t.m. don't much feel like changing it everywhere else
            swfo['padlefts'] = [padleft, ]
            swfo['padrights'] = [padright, ]
            if debug:
                print '    %3d   %3d    %s' % (padleft, padright, query)

        if debug:
            print '    cdr3        uid                 padded seq'
            for query in sorted(self.info['queries'], key=lambda q: self.info[q]['cdr3_length']):
                print '    %3d   %20s    %s' % (self.info[query]['cdr3_length'], query, self.info[query]['seqs'][0])

    # ----------------------------------------------------------------------------------------
    def combine_indels(self, qinfo, best, debug=False):
        # debug = 2
        if debug:
            print '  %s: combine with %d sw indels: %s' % (qinfo['name'], len(qinfo['new_indels']), ' '.join(qinfo['new_indels'].keys()))
            for ifo in qinfo['new_indels'].values():
                print indelutils.get_dbg_str(ifo)

        regional_indelfos = {}
        qrbounds = {r : qinfo['qrbounds'][best[r]] for r in utils.regions}
        full_qrseq = qinfo['seq']
        if self.vs_info is not None and qinfo['name'] in self.vs_indels:
            if debug:
                print '    has a vsearch v indel'

            if 'v' in qinfo['new_indels']:  # if sw kicks up an additional v indel that vsearch didn't find, we rerun sw with <self.args.no_indel_gap_open_penalty>
                if debug:
                    print '      sw called a v indel, but there\'s already a vsearch v indel, sogive up (rerun sw forbidding indels)'
                return None

            vs_indelfo = self.info['indels'][qinfo['name']]
            assert 'v' in vs_indelfo['genes']  # a.t.m. vsearch is only looking for v genes, and if that changes in the future we'd need to rewrite this
            non_v_bases = len(qinfo['seq']) - qrbounds['v'][1]  # have to trim things to correspond to the new (and potentially different) sw bounds (note that qinfo['seq'] corresponds to the indel reversion from vs, but not from sw)
            vs_indelfo['qr_gap_seq'] = vs_indelfo['qr_gap_seq'][qrbounds['v'][0] : len(vs_indelfo['qr_gap_seq']) - non_v_bases]
            vs_indelfo['gl_gap_seq'] = vs_indelfo['gl_gap_seq'][qrbounds['v'][0] : len(vs_indelfo['gl_gap_seq']) - non_v_bases]
            assert len(vs_indelfo['qr_gap_seq']) == len(vs_indelfo['gl_gap_seq'])

            net_v_indel_length = indelutils.net_length(vs_indelfo)
            qrbounds['v'] = (qrbounds['v'][0], qrbounds['v'][1] + net_v_indel_length)
            for region in ['d', 'j']:
                qrbounds[region] = (qrbounds[region][0] + net_v_indel_length, qrbounds[region][1] + net_v_indel_length)
                if region in qinfo['new_indels']:
                    for ifo in qinfo['new_indels'][region]['indels']:
                        ifo['pos'] += net_v_indel_length
            full_qrseq = self.input_info[qinfo['name']]['seqs'][0]  # should in principle replace qinfo['seq'] as well, since it's the reversed seq from vsearch, but see note below
            del self.info['indels'][qinfo['name']]
            regional_indelfos['v'] = vs_indelfo
        elif 'v' in qinfo['new_indels']:
            regional_indelfos['v'] = qinfo['new_indels']['v']

        if 'd' in qinfo['new_indels']:
            regional_indelfos['d'] = qinfo['new_indels']['d']
        if 'j' in qinfo['new_indels']:
            regional_indelfos['j'] = qinfo['new_indels']['j']

        # NOTE qinfo won't be consistent with the indel reversed seq after this, but that's kind of the point, since we're just rerunning anyway

        return indelutils.combine_indels(regional_indelfos, full_qrseq, qrbounds, uid=qinfo['name'], debug=debug)
