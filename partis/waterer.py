from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
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

from . import utils
from . import glutils
from . import indelutils
from .parametercounter import ParameterCounter
from .performanceplotter import PerformancePlotter
from . import seqfileopener
from io import open

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
                 duplicates=None, pre_failed_queries=None, aligned_gl_seqs=None, vs_info=None, msa_vs_info=None, locus=None):
        self.args = args
        self.input_info = input_info if input_info is not None else OrderedDict()  # NOTE do *not* modify <input_info>, since it's this original input info from partitiondriver
        self.reco_info = reco_info
        self.glfo = copy.deepcopy(glfo)  # NOTE gets overwritten by read_cachefile()
        self.count_parameters = count_parameters  # NOTE *not* the same as <self.args.cache_parameters>
        self.plot_annotation_performance = plot_annotation_performance  # NOTE *not* the same as <self.args.plot_annotation_performance>
        self.simglfo = simglfo
        self.parameter_out_dir = parameter_out_dir
        self.duplicates = {} if duplicates is None else duplicates
        self.already_added_queries = []  # this is confusing since I'm adding it very late, but just want to keep track of queries that occur more than once in cache file (which somehow seems not to have come up before now)
        self.debug = self.args.debug if self.args.sw_debug is None else self.args.sw_debug
        self.aligned_gl_seqs = aligned_gl_seqs
        self.vs_info = vs_info
        self.msa_vs_info = msa_vs_info

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

        self.skipped_unproductive_queries, self.skipped_in_frame_queries, self.kept_unproductive_queries = set(), set(), set()

        self.my_gldir = self.args.workdir + '/sw-' + glutils.glfo_dir
        if self.glfo is None:  # reading cache file from bin/partis, rather than normal operation from python/partitiondriver.py
            self.args.locus = locus
        else:  # default, normal operation
            glutils.write_glfo(self.my_gldir, self.glfo)  # NOTE gets overwritten by read_cachefile()

        if not os.path.exists(self.args.ig_sw_binary):
            raise Exception('ig-sw binary d.n.e: %s' % self.args.ig_sw_binary)

    # ----------------------------------------------------------------------------------------
    def run(self, cachefname=None):
        start = time.time()
        base_infname = 'query-seqs.fa'
        base_outfname = 'query-seqs.sam'

        if self.vs_info is not None or self.msa_vs_info is not None:  # if we're reading a cache file, we should make sure to read the exact same info from there
            self.add_vs_indels()

        itry = 0
        while True:  # if we're not running vsearch, we still gotta run twice to get shm indeld sequences
            mismatches, gap_opens, queries_for_each_proc = self.split_queries(self.args.n_procs)  # NOTE can tell us to run more than <self.args.n_procs> (we run at least one proc for each different mismatch score)
            self.write_input_files(base_infname, queries_for_each_proc)

            print('    running %d proc%s for %d seq%s' % (len(mismatches), utils.plural(len(mismatches)), len(self.remaining_queries), utils.plural(len(self.remaining_queries))))
            sys.stdout.flush()
            self.execute_commands(base_infname, base_outfname, mismatches, gap_opens)

            processing_start = time.time()
            self.read_output(base_outfname, len(mismatches))

            if itry > 1 or len(self.indel_reruns) == 0:
                break
            itry += 1

        self.finalize(cachefname)
        print('    water time: %.1f  (ig-sw %.1f  processing %.1f)' % (time.time() - start, time.time() - processing_start, self.ig_sw_time))

    # ----------------------------------------------------------------------------------------
    def clean_cache(self, cache_path):
        for suffix in ['.csv', '.yaml']:
            if os.path.exists(cache_path + suffix):
                print('  removing old sw cache %s%s' % (cache_path, suffix))
                os.remove(cache_path + suffix)
        if os.path.exists(cache_path + '-glfo'):
            print('  removing old sw cache glfo %s-glfo' % cache_path)
            glutils.remove_glfo_files(cache_path + '-glfo', self.args.locus)

    # ----------------------------------------------------------------------------------------
    def read_cachefile(self, cachefname, ignore_seed_unique_id=False, quiet=False):
        no_input_info = len(self.input_info) == 0
        if not quiet:  # ok i'm also using this as an indicator of how much dbg i want, which maybe is dumb
            start = time.time()
            print('        reading sw results from %s' % cachefname)

        if utils.getsuffix(cachefname) == '.csv':  # old way
            cachebase = utils.getprefix(cachefname)
            if os.path.exists(cachebase + '-glfo'):  # NOTE replaces original <self.glfo>
                self.glfo = glutils.read_glfo(cachebase + '-glfo', self.args.locus)
            else:
                print('    %s didn\'t find a germline info dir along with sw cache file, but trying to read it anyway' % utils.color('red', 'warning'))
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
            uid = line['unique_ids'][0]
            if no_input_info:  # this means we're reading the cache file from bin/partis for paired clustering, so don't have input info available, so have to make it here
                self.input_info[uid] = line
            else:  # normal operation
                if uid not in self.input_info:
                    continue
            if self.args.ignore_sw_pair_info and 'paired-uids' in line:
                del line['paired-uids']
            utils.add_implicit_info(self.glfo, line, aligned_gl_seqs=self.aligned_gl_seqs)
            if indelutils.has_indels_line(line, 0):
                self.info['indels'][uid] = line['indelfos'][0]
            self.add_to_info(line)
            if len(self.remaining_queries) > 0:
                self.remaining_queries -= set(line['duplicates'][0])
            self.duplicates[uid] = line['duplicates'][0]  # Note that <self.duplicates> is really just so partitiondriver can pass in previous duplicates, but then now that we have the information in two places we need to keep it synchronized (see similar code in self.remove_duplicate_sequences())

        glutils.write_glfo(self.my_gldir, self.glfo)

        self.finalize(cachefname=None, just_read_cachefile=True, ignore_seed_unique_id=ignore_seed_unique_id, quiet=quiet)
        if not quiet:
            print('        water time: %.1f' % (time.time()-start))

    # ----------------------------------------------------------------------------------------
    def write_cachefile(self, cachefname):
        print('        writing sw results to %s' % cachefname)

        if utils.getsuffix(cachefname) == '.csv':
            cachebase = utils.getprefix(cachefname)
            glutils.write_glfo(cachebase + '-glfo', self.glfo)

        # NOTE do _not_ add extra headers here since if they're in the sw cache file I'd have to deal with removing them when I read it
        # NOTE does *not* write failed queries
        utils.write_annotations(cachefname, self.glfo, [self.info[q]for q in self.info['queries']], utils.sw_cache_headers, use_pyyaml=self.args.write_full_yaml_output, dont_write_git_info=self.args.dont_write_git_info)

    # ----------------------------------------------------------------------------------------
    def finalize(self, cachefname=None, just_read_cachefile=False, ignore_seed_unique_id=False, quiet=False):
        if self.debug:
            print('%s' % utils.color('green', 'finalizing'))

        self.info['queries'] = [q for q in self.input_info if q in self.info['passed-queries']]  # it might be cleaner eventually to just make self.info['queries'] a set, but it's used in so many other places that I'm worried about changing it
        self.info['failed-queries'] |= set(self.remaining_queries)  # perhaps it doesn't make sense to keep both of 'em around after finishing?

        if len(self.info['queries']) == 0:
            print('%s no queries passing smith-waterman, exiting' % utils.color('red', 'warning'))
            sys.exit(0)

        if not quiet:  # don't want any debug printing if we're ignoring seed id (i.e. just reading cache file for paired h/l seed clustering)
            # kind of messy, but if we just read the cache file then we need to print duplicates, but if we didn't just read the cache file then we haven't yet removed duplicates so we can't
            dupl_str, n_dupes = '', 0
            if just_read_cachefile:
                n_dupes = sum([len(dupes) for dupes in self.duplicates.values()])
                dupl_str = ', %d duplicates [removed when cache file was written]' % n_dupes
            tmp_pass_frac = float(len(self.info['queries']) + n_dupes + len(self.skipped_unproductive_queries) + len(self.skipped_in_frame_queries)) / len(self.input_info)
            print('      info for %d / %d = %.3f   (removed: %d failed%s)' % (len(self.info['queries']), len(self.input_info), tmp_pass_frac, len(self.info['failed-queries']), dupl_str))
            if tmp_pass_frac < 0.80:
                print('  %s smith-waterman step failed to find alignments for a large fraction of input sequences (see previous line)   %s'  % (utils.color('red', 'warning'), utils.reverse_complement_warning()))
            if len(self.kept_unproductive_queries) > 0:
                print('      kept %d (%.3f) unproductive' % (len(self.kept_unproductive_queries), float(len(self.kept_unproductive_queries)) / len(self.input_info)))

        if not just_read_cachefile:  # it's past tense!
            if len(self.skipped_unproductive_queries) > 0:
                print('         skipped %d unproductive' % len(self.skipped_unproductive_queries))
            if len(self.skipped_in_frame_queries) > 0:
                print('         skipped %d with in frame rearrangement (i.e. cdr3 len % 3 == 0)' % len(self.skipped_in_frame_queries))
            if self.debug:
                if len(self.info['indels']) > 0:
                    print('      indels: %s' % ':'.join(list(self.info['indels'].keys())))
                if len(self.info['failed-queries']) > 0:
                    print('            missing annotations: ' + ' '.join(sorted(self.info['failed-queries'])))
                    if not self.args.is_data:
                        print('true annotations for remaining events:')
                        for qry in sorted(self.info['failed-queries']):
                            utils.print_reco_event(self.reco_info[qry], extra_str='      ', label='true:')

            assert len(self.info['queries']) + len(self.skipped_unproductive_queries) + len(self.skipped_in_frame_queries) + len(self.info['failed-queries']) == len(self.input_info)

        if len(self.already_added_queries) > 0:
            count_str = ' '.join('%d'%len(list(ugroup)) for uid, ugroup in itertools.groupby(sorted(self.already_added_queries)))
            print('      %s %s %d queries more than once (counts: %s)' % (utils.wrnstr(), 'read' if just_read_cachefile else 'found', len(set(self.already_added_queries)), count_str))

        if self.count_parameters:
            pcounter = ParameterCounter(self.glfo, self.args, count_correlations=self.args.count_correlations)
            for qname in self.info['queries']:
                pcounter.increment(self.info[qname])
        if self.plot_annotation_performance:
            perfplotter = PerformancePlotter('sw')
            for qname in self.info['queries']:
                perfplotter.evaluate(self.reco_info[qname], self.info[qname], simglfo=self.simglfo)
        if self.args.align_constant_regions:
            utils.parse_constant_regions(self.args.species, self.args.locus, [self.info[q] for q in self.info['queries']], self.args.workdir, debug=self.args.debug)

        # remove queries with cdr3 length different to the seed sequence
        if self.args.seed_unique_id is not None and not ignore_seed_unique_id:
            assert not self.count_parameters  # informative exception is raised in bin/partis
            if self.args.seed_unique_id in self.info['queries']:
                seed_cdr3_length = self.info[self.args.seed_unique_id]['cdr3_length']
            else:  # if it failed
                print('    %s seed unique id \'%s\' not in final s-w queries, so removing all queries (if seed seq is coming from --queries-to-include-fname, maybe need to run parameter caching with this arg set?)' % (utils.color('yellow', 'note:'), self.args.seed_unique_id))
                seed_cdr3_length = -1
            initial_n_queries = len(self.info['queries'])
            for query in copy.deepcopy(self.info['queries']):
                if self.info[query]['cdr3_length'] != seed_cdr3_length:
                    self.remove_query(query)
            n_removed = initial_n_queries - len(self.info['queries'])
            if n_removed > 0 and not quiet:
                print('      removed %d / %d = %.2f sequences with cdr3 length different from seed sequence (leaving %d)' % (n_removed, initial_n_queries, float(n_removed) / initial_n_queries, len(self.info['queries'])))

        seqfileopener.add_input_metafo(self.input_info, [self.info[q] for q in self.info['queries']], keys_not_to_overwrite=['multiplicities'] if just_read_cachefile else None)  # need to do this before removing duplicates, so the duplicate info (from waterer) can get combined with multiplicities (from input metafo). And, if we just read the cache file, then we already collapsed duplicates so we don't want to overwrite multiplicity info
        if self.args.is_data:  # it's too much trouble updating reco_info on simulation, and I don't think we can add fwk insertions in simulation atm anyway
            self.remove_framework_insertions()
        if self.args.collapse_duplicate_sequences and not just_read_cachefile:
            self.remove_duplicate_sequences()

        # want to do this *before* we pad sequences, so that when we read the cache file we're reading unpadded sequences and can pad them below
        if cachefname is not None:
            self.write_cachefile(cachefname)

        self.pad_seqs_to_same_length()  # NOTE this uses all the gene matches (not just the best ones), so it has to come before we call pcounter.write(), since that fcn rewrites the germlines removing genes that weren't best matches. But NOTE also that I'm not sure what but that the padding actually *needs* all matches (rather than just all *best* matches)

        if self.plot_annotation_performance and self.args.plotdir is not None:
            perfplotter.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots)

        if self.count_parameters:
            if self.args.plotdir is not None:
                pcounter.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots, only_overall=not self.args.make_per_gene_plots)
            if self.parameter_out_dir is not None and not self.args.dont_write_parameters:
                pcounter.write(self.parameter_out_dir)

        glutils.remove_glfo_files(self.my_gldir, self.args.locus)
        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def add_vs_indels(self):
        vsfo = utils.non_none([self.msa_vs_info, self.vs_info])['annotations']  # use msa info if we have it (note that we still need the regular vs_info for the v genes)
        queries_with_indels = [q for q in sorted(self.remaining_queries) if q in vsfo and indelutils.has_indels(vsfo[q]['indelfo'])]
        for query in queries_with_indels:
            self.vs_indels.add(query)  # this line is to tell us that this query has an indel stemming from vsearch, while the next line tells us that there's an indel (combine_indels() gets confused if we don't differentiate between the two)
            self.info['indels'][query] = copy.deepcopy(vsfo[query]['indelfo'])
            if self.debug:
                print(indelutils.get_dbg_str(vsfo[query]['indelfo'], uid=query))

        if self.debug and len(queries_with_indels) > 0:
            print('    added %d %s indel%s%s' % (len(queries_with_indels), 'vsearch' if self.msa_vs_info is None else 'mafft msa', utils.plural(len(queries_with_indels)), (' (%s)' % ' '.join(queries_with_indels)) if len(queries_with_indels) < 100 else ''))

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
            print('start')
            for m, queries in zip(mismatch_vals, query_groups):
                print('  %d   %d    %s' % (m, len(queries), ' '.join(queries)))

        # first give one proc to each mismatch (note that this takes at least as many procs as there are different mismatch values, even if --n-procs is smaller)
        while len(mismatch_vals) < n_procs:
            if debug:
                print('%d < %d' % (len(mismatch_vals), n_procs))
            largest_group = max(query_groups, key=len)
            piece_a, piece_b = largest_group[ : len(largest_group) // 2], largest_group[len(largest_group) // 2 : ]  # split up the largest group into two [almost] equal pieces
            ilargest = query_groups.index(largest_group)
            query_groups = query_groups[:ilargest] + [piece_a, piece_b] + query_groups[ilargest + 1:]
            mismatch_vals = mismatch_vals[:ilargest] + [mismatch_vals[ilargest], mismatch_vals[ilargest]] + mismatch_vals[ilargest + 1:]
            if debug:
                for m, queries in zip(mismatch_vals, query_groups):
                    print('  %d   %d    %s' % (m, len(queries), ' '.join(queries)))

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
        input_queries = sorted(self.remaining_queries)
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
            print('%s' % utils.color('green', 'reading output'))

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
            print('\n%s didn\'t read %s from %s' % (utils.color('red', 'warning'), ' '.join(not_read), self.args.workdir))

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

        last_scores, bad_cigars = {r : None for r in utils.regions}, []
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

            if 'M' not in read.cigarstring:  # cigar str doesn't actually have any matches, which means the cigar parsing stuff will fail
                bad_cigars.append(read.cigarstring)
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

        if len(bad_cigars) > 0:
            print('    %s no M in %d / %d cigar strs for %s: %s' % (utils.wrnstr(), len(bad_cigars), len(reads), qinfo['name'], ' '.join(bad_cigars)))

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

        print(''.join(out_str_list))

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
                print('      rerun: super high vj mutation (%.3f > %.3f)' % (vj_mute_freq, self.args.max_vj_mut_freq))
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
            lb = qinfo['qrbounds'][l_gene]  # shorthands
            rb = qinfo['qrbounds'][r_gene]
            print('   %schecking %s/%s boundary overlap status' % ('recursively ' if recursed else '', l_reg, r_reg))
            print('       %s %s' % (' ', qinfo['seq']))
            print('       %s %s%s' % (l_reg, ' ' * lb[0], qinfo['seq'][lb[0] : lb[1]]))
            print('       %s %s%s' % (r_reg, ' ' * rb[0], qinfo['seq'][rb[0] : rb[1]]))
            print('     %s %s    overlap %d    available space %d' % (l_reg, r_reg, overlap, available_space))

        status = 'ok'
        if overlap > 0:  # positive overlap means they actually overlap
            status = 'overlap'
        if overlap > available_space or overlap == 1 and available_space == 1:  # call it nonsense if the boundaries are really whack (i.e. there isn't enough space to resolve the overlap) -- we'll presumably either toss the query or rerun with different match/mismatch
            status = 'nonsense'

        if debug:
            print('  overlap status: %s' % status)

        if not recursed and status == 'nonsense' and l_reg == 'd':  # on rare occasions with very high mutation, vdjalign refuses to give us a j match that's at all to the right of the d match
            assert l_reg == 'd' and r_reg == 'j'
            qrb, glb = [qinfo[t+'bounds'] for t in ('qr', 'gl')]  # shorthands
            if debug:
                print('  %s: synthesizing d match' % qinfo['name'])
                if qrb[r_gene][0] > qrb[l_gene][0]:  # in most (maybe all) cases we get here, it's because the d match is to the right of the j match
                    print('   huh, i should only get here if the start of the r_gene (j) is to the *left* of the start of the l_gene (d) [i.e. they\'re on the wrong side of each other], but here it looks like the opposite is true. Not really sure it\'s a problem, but it\'s probably worth looking in to')
            v_end = qrb[best['v']][1]
            lpos = v_end + (qrb[r_gene][0] - v_end) // 2  # halfway between end of v and start of j
            qrb[l_gene] = (lpos, lpos + 1)  # swap whatever crummy nonsense d match we have now for a one-base match to the left of the j match (since unlike d, the j match is probably decent)
            d_midpoint = len(self.glfo['seqs'][l_reg][l_gene]) // 2
            glb[l_gene] = (d_midpoint, d_midpoint + 1)
            if l_reg in qinfo['new_indels']:  # d is not one base long, so it can't very well have an indel (also, we just moved it to a completely different location)
                del qinfo['new_indels'][l_reg]
            if qrb[l_gene][1] > qrb[r_gene][0]:  # if left side of j needs adjusting (won't need to do this unless there was no space between the v and j, and this probably never happens any more -- it's here because I used to put the synthetic d base at the first position of the j match)
                qrb[r_gene] = (qrb[r_gene][0] + 1, qrb[r_gene][1])  # this'll give a zero-length match if the match was originally length 1, but that should just result in nonsense bounds below I think, so it should be ok
                glb[r_gene] = (glb[r_gene][0] + 1, glb[r_gene][1])
                if r_reg in qinfo['new_indels']:  # start of j match moves one base to right, so have to adjust its indelfo
                    ifo = qinfo['new_indels'][r_reg]
                    ifo['qr_gap_seq'] = ifo['qr_gap_seq'][1 : ]
                    ifo['gl_gap_seq'] = ifo['gl_gap_seq'][1 : ]
                    if debug:
                        print('    removed one base from left side of %s indel gap seqs (to match j bound moved because of new synthetic d match)' % r_reg)

            status = self.check_boundaries(rpair, qinfo, best, recursed=True, debug=debug)
            if status == 'overlap':
                if debug:
                    print('  \'overlap\' status after synthesizing d match. Setting to \'nonsense\', this is too darn hard.')
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
        qrb, glb = [qinfo[t+'bounds'] for t in ('qr', 'gl')]
        # ----------------------------------------------------------------------------------------
        def dbg_print(l_length, r_length, l_portion, r_portion):
            lb = qrb[l_gene]  # shorthands
            rb = qrb[r_gene]
            print('      %4d %4d      %4d %4d' % (l_length, r_length, l_portion, r_portion))
            l_gene_seq = qinfo['seq'][lb[0] : lb[1] - l_portion]
            r_gene_seq = qinfo['seq'][rb[0] + r_portion : rb[1]]
            l_offset = min(rb[0] + r_portion, lb[0])
            print('                                    %s%s' % ((lb[0] - l_offset) * ' ', l_gene_seq))
            print('                                    %s%s' % ((rb[0] + r_portion - l_offset) * ' ', r_gene_seq))
        # ----------------------------------------------------------------------------------------
        def shift_bounds(l_portion, r_portion):  # shift match boundaries of best left gene <l_gene> and best right gene <r_gene> according to <l_portion> and <r_portion>, then also shift other/non-best gene boundaries
            for gtmp in qrb:
                for bdict in (qrb, glb):
                    istart, istop = bdict[gtmp]
                    if utils.get_region(gtmp) == l_reg:
                        bdict[gtmp] = (istart, max(istart, istop - l_portion))  # max() is only necessary for non-best genes
                    elif utils.get_region(gtmp) == r_reg:
                        bdict[gtmp] = (min(istop, istart + r_portion), istop)  # min() is only necessary for non-best genes
                    else:
                        continue
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

        # ----------------------------------------------------------------------------------------
        overlap, available_space = self.get_overlap_and_available_space(rpair, best, qrb)

        if overlap <= 0:  # nothing to do, they're already consistent
            print('shouldn\'t get here any more if there\'s no overlap')
            return

        if overlap > available_space:
            raise Exception('overlap %d bigger than available space %d between %s and %s for %s' % (overlap, available_space, l_reg, r_reg, qinfo['name']))

        if debug:
            print('   fixing %s%s overlap:  %d-%d overlaps with %d-%d by %d' % (l_reg, r_reg, qrb[l_gene][0], qrb[l_gene][1], qrb[r_gene][0], qrb[r_gene][1], overlap))

        l_length = qrb[l_gene][1] - qrb[l_gene][0]  # initial length of lefthand gene match
        r_length = qrb[r_gene][1] - qrb[r_gene][0]  # and same for the righthand one
        l_portion, r_portion = 0, 0  # portion of the initial overlap that we give to each side
        if debug:
            print('        lengths        portions     ')
            dbg_print(l_length, r_length, l_portion, r_portion)
        while l_portion + r_portion < overlap:
            if l_length <= 1 and r_length <= 1:  # don't want to erode match (in practice it'll be the d match) all the way to zero
                raise Exception('both lengths went to one without resolving overlap for %s: %s %s' % (qinfo['name'], qrb[l_gene], qrb[r_gene]))
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
                    print('  failed: apportionment encountered an indel')
                return True
            # if debug:  # uncomment to get more verbosity
            #     dbg_print(l_length, r_length, l_portion, r_portion)
        if debug:
            dbg_print(l_length, r_length, l_portion, r_portion)

        if debug:
            print('      %s apportioned %d bases between %s (%d) match and %s (%d) match' % (qinfo['name'], overlap, l_reg, l_portion, r_reg, r_portion))
        assert l_portion + r_portion == overlap
        shift_bounds(l_portion, r_portion)
        if l_reg in qinfo['new_indels'] and l_portion > 0:  # it would be nice to check indel consistency after doing this, but my indel consistency checkers seem to only operate on <line>s, and I don't have one of those yet
            ifo = qinfo['new_indels'][l_reg]
            ifo['qr_gap_seq'] = ifo['qr_gap_seq'][ : -l_portion]  # this should really check that l_portion is shorter than the qr gap seq... but it _should_ be impossible
            ifo['gl_gap_seq'] = ifo['gl_gap_seq'][ : -l_portion]
            if debug:
                print('    removed %d base%s from right side of %s indel gap seqs' % (l_portion, utils.plural(l_portion), l_reg))
        if r_reg in qinfo['new_indels'] and r_portion > 0:
            ifo = qinfo['new_indels'][r_reg]
            ifo['qr_gap_seq'] = ifo['qr_gap_seq'][r_portion : ]
            ifo['gl_gap_seq'] = ifo['gl_gap_seq'][r_portion : ]
            if debug:
                print('    removed %d base%s from left side of %s indel gap seqs' % (r_portion, utils.plural(r_portion), r_reg))

        return False

    # ----------------------------------------------------------------------------------------
    def fix_conserved_codon_deletions(self, qinfo, best, debug=False):
        assert False  # needs to fix indel stuff if i want to turn it on
        for region, side in [('v', '3p'), ('j', '5p')]:
            for _, gene in qinfo['matches'][region]:
                qrb, glb = [qinfo[t+'bounds'] for t in ('qr', 'gl')]
                if region == 'v' and glb[gene][1] < utils.cdn_pos(self.glfo, region, gene) + 3:
                    n_to_expand = len(self.glfo['seqs'][region][gene]) - glb[gene][1]  # expand all the way to the end of the gl seq, since if it deleted way past the conserved codon it probably has no idea what's actually right, so we want the hmm to check everything
                    if debug:
                        print('      %s_%s %s match deletes conserved codon, so increasing righthand bounds: qr %d --> %d, gl %d --> %d' % (region, side, utils.color_gene(gene), qrb[gene][1], qrb[gene][1] + n_to_expand, glb[gene][1], glb[gene][1] + n_to_expand))
                    qrb[gene] = (qrb[gene][0], qrb[gene][1] + n_to_expand)
                    glb[gene] = (glb[gene][0], glb[gene][1] + n_to_expand)
                elif region == 'j' and glb[gene][0] > utils.cdn_pos(self.glfo, region, gene):
                    n_to_expand = glb[gene][0]  # i guess also extend the j all the way? this is starting to make less sense
                    if debug:
                        print('      %s_%s %s match deletes conserved codon, so increasing lefthand bounds: qr %d --> %d, gl %d --> %d' % (region, side, utils.color_gene(gene), qrb[gene][1], qrb[gene][1] + n_to_expand, glb[gene][1], glb[gene][1] + n_to_expand))
                    qrb[gene] = (qrb[gene][0] - n_to_expand, qrb[gene][1])
                    glb[gene] = (glb[gene][0] - n_to_expand, glb[gene][1])  # this just takes the start bound it to 0

    # ----------------------------------------------------------------------------------------
    def remove_probably_spurious_deletions(self, qinfo, best, debug=False):  # remove probably-spurious v_5p and j_3p deletions
        if debug:
            print('  looking for spurious v_5p and j_3p deletions')
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

            if debug:
                print('      %s   deletion len: %3d   insertion: %s' % (erosion, d_len, insertion))

            if d_len == 0 or len(insertion) == 0:
                continue
            if len(insertion) >= d_len:
                spurious_bases = insertion[len(insertion) - d_len:]  # i.e. the last <v_5p_del> bases of the fv insertion (or first <j_3p_del> bases of the jf insertion)
            else:
                spurious_bases = insertion  # the whole damn thing
            if spurious_bases.count(utils.ambig_base) == len(spurious_bases):  # don't do it if it's all Ns
                continue
            if debug:
                print('EXPANDING %s %s d_len: %d   insertion: %s (len %d)   spurious: %s (len %d)' % (qinfo['name'], region, d_len, insertion, len(insertion), spurious_bases, len(spurious_bases)))
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

        line = {}
        line['unique_ids'] = [qname, ]  # redundant, but used somewhere down the line
        line['seqs'] = [qinfo['seq'], ]  # NOTE this is the seq output by vdjalign, i.e. if we reversed any indels it is the reversed sequence, also NOTE many, many things depend on this list being of length one
        line['input_seqs'] = [self.input_info[qname]['seqs'][0], ]

        for region in utils.regions:
            line[region + '_gene'] = best[region]

        # erosion, insertion, mutation info for best match
        line['v_5p_del'] = qinfo['glbounds'][best['v']][0]
        line['v_3p_del'] = len(self.glfo['seqs']['v'][best['v']]) - qinfo['glbounds'][best['v']][1]  # len(germline v) - gl_match_end
        line['d_5p_del'] = qinfo['glbounds'][best['d']][0]
        line['d_3p_del'] = len(self.glfo['seqs']['d'][best['d']]) - qinfo['glbounds'][best['d']][1]
        line['j_5p_del'] = qinfo['glbounds'][best['j']][0]
        line['j_3p_del'] = len(self.glfo['seqs']['j'][best['j']]) - qinfo['glbounds'][best['j']][1]

        line['fv_insertion'] = qinfo['seq'][ : qinfo['qrbounds'][best['v']][0]]
        line['vd_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['v']][1] : qinfo['qrbounds'][best['d']][0]]
        line['dj_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['d']][1] : qinfo['qrbounds'][best['j']][0]]
        line['jf_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['j']][1] : ]
        line['leader_seqs'] = [line['fv_insertion']]  # yes it sucks to have both these and fv/jf insertions, but see notes elsewhere -- it'd be a mess to try to remove fv/jf insertions
        line['c_gene_seqs'] = [line['jf_insertion']]

        if qname in self.info['indels']:  # NOTE at this piont indel info isn't updated for any change in the genes (the indelfo gets updated during utils.add_implicit_info())
            line['indelfos'] = [self.info['indels'][qname]]  # NOTE this makes it so that self.info[uid]['indelfos'] *is* self.info['indels'][uid]. It'd still be nicer to eventually do away with self.info['indels'], although I'm not sure that's really either feasible or desirable given other constraints
        else:
            line['indelfos'] = [indelutils.get_empty_indel()]

        line['duplicates'] = [self.duplicates.get(qname, []), ]  # note that <self.duplicates> doesn't handle simultaneous seqs, i.e. it's for just a single sequence

        matchfo = {}
        for region in utils.regions:
            matchfo[region] = {  # NOTE could use an ordered dict, but json dump() then writes it as a regular dict, so it's no longer ordered if read from sw cache file, so it's better to just have it alwyas be unsorted to avoid confusion
                gene : {'score' : score, 'glbounds' : qinfo['glbounds'][gene], 'qrbounds' : qinfo['qrbounds'][gene]} for score, gene in qinfo['matches'][region]
            }
        line['all_matches'] = [matchfo]

        line['cdr3_length'] = codon_positions['j'] - codon_positions['v'] + 3
        line['codon_positions'] = codon_positions

        line['padlefts'] = ['']
        line['padrights'] = ['']

        return line

    # ----------------------------------------------------------------------------------------
    def check_simulation_kbounds(self, line, true_line):
        true_kbounds = {'v' : {}, 'd' : {}}
        true_kbounds['v']['best'] = true_line['regional_bounds']['v'][1]
        true_kbounds['d']['best'] = true_line['regional_bounds']['d'][1] - true_line['regional_bounds']['v'][1]

        def print_kbound_warning():
            print('  %s true kset (%s) not within kbounds (%s) for %s' % (utils.color('red', 'warning'), utils.kbound_str(true_kbounds), utils.kbound_str({r : line['k_' + r] for r in true_kbounds}), ':'.join(line['unique_ids'])))

        for region in true_kbounds:
            if true_kbounds[region]['best'] < line['k_' + region]['min'] or true_kbounds[region]['best'] >= line['k_' + region]['max']:
                print_kbound_warning()

    # ----------------------------------------------------------------------------------------
    def add_to_info(self, line):
        assert len(line['unique_ids'])
        qname = line['unique_ids'][0]
        if qname in self.info:  # i didn't end up needing this for the case I impelemented it, but it still seems worth having
            self.already_added_queries.append(qname)
            return

        self.info['passed-queries'].add(qname)
        self.info[qname] = line

        # add this query's matches into the overall gene match sets
        for region in utils.regions:
            self.info['all_best_matches'].add(self.info[qname][region + '_gene'])
            self.info['all_matches'][region] |= set(self.info[qname]['all_matches'][0][region])  # NOTE there's an 'all_matches' in this query's info, and also in <self.info>

        # everything this is flagging seems to be cases where the insertion (or the overlapping germline region) is the same as the germline base, i.e. the kbound getter is working
        # if not self.args.is_data:
        #     self.check_simulation_kbounds(self.info[qname], self.reco_info[qname])

        if self.debug:
            inf_label = ' ' + utils.kbound_str({r : line['k_' + r] for r in ['v', 'd']})
            if not self.args.is_data:
                inf_label = 'inf: ' + inf_label
                utils.print_reco_event(self.reco_info[qname], extra_str='    ', label=utils.color('green', 'true:'), extra_print_keys=self.args.extra_print_keys)
            utils.print_reco_event(self.info[qname], extra_str='    ', label=inf_label, extra_print_keys=self.args.extra_print_keys)

        if not utils.is_functional(self.info[qname], iseq=0):
            self.kept_unproductive_queries.add(qname)
        if len(self.remaining_queries) > 0:
            self.remaining_queries.remove(qname)

    # ----------------------------------------------------------------------------------------
    def summarize_query(self, qinfo):
        """
        Fiddle with a few things, but mostly decide whether we're satisfied with the current matches.
        If not, return without calling add_to_info().
        """
        def dbgfcn(dbgstr):  # doesn't return anything, but it makes things more concise
            if self.debug:
                print('      rerun: %s' % dbgstr)

        qname = qinfo['name']
        qseq = qinfo['seq']
        assert qname not in self.info
        if self.debug:
            print('  %s' % qname)

        # in very rare cases the j match gets extended so far left it touches the v, in which case we get no d match, but we don't really want to throw these out
        tmpbest = {r : qinfo['matches'][r][0][1] for r in utils.regions if len(qinfo['matches'][r]) > 0}  # this is hackey and duplicates too much code from just below, but I don't want to change anything that's down there since it's been there for so long and is so thoroughly tested
        if 'd' in utils.getregions(self.args.locus) and set(tmpbest.keys()) == set(['v', 'j']):  # missing d (not sure I actually need the first clause)
            v_end, j_start = qinfo['qrbounds'][tmpbest['v']][1], qinfo['qrbounds'][tmpbest['j']][0]
            if v_end >= j_start:  # set all the d genes to zero-length matches at that point, then the overlap resolution should make things at least consistent, and then the hmm should sort stuf out
                for dg in self.glfo['seqs']['d']:  # NOTE duplicates a bit of code in check_boundaries()
                    qinfo['matches']['d'].append((0, dg))  # add the d gene as a match with score 0
                    qinfo['qrbounds'][dg] = (v_end, v_end + 1)
                    d_midpoint = len(self.glfo['seqs']['d'][dg]) // 2
                    qinfo['glbounds'][dg] = (d_midpoint, d_midpoint + 1)

        # do we have a match for each region?
        for region in utils.getregions(self.args.locus):
            if len(qinfo['matches'][region]) == 0:
                return dbgfcn('no %s match' % region)  # if no d match found, maybe we should just assume entire d was eroded?

        best = {r : qinfo['matches'][r][0][1] for r in utils.regions}  # already made sure there's at least one match for each region

        # NOTE this works ok, but if i want to use it, it would need to also fix the indel lengths (and that doesn't seem worthwhile, since a) this is already done when building hmms and b) same goal is accomplished by the new stuff in get_kbounds()
        # if not self.args.allow_conserved_codon_deletion:
        #     self.fix_conserved_codon_deletions(qinfo, best)

        # s-w allows d and j matches to overlap, so we need to apportion the disputed bases UPDATE ok it turns out ig-sw only allowed overlaps because it had a bug, and the bug is fixed (https://github.com/psathyrella/partis/commit/471e5eac6d2b0fbdbb2b6024c81af14cdc3d9399), so maybe this stuff will never get used now.
        for rpair in utils.region_pairs(self.glfo['locus']):
            overlap_status = self.check_boundaries(rpair, qinfo, best)  #, debug=self.debug>1)  # I think all this overlap fixing only adjusts the bounds for the best match, but I guess that's ok? It's certainly been like that for ages
            if overlap_status == 'overlap':
                overlap_indel_fail = self.shift_overlapping_boundaries(rpair, qinfo, best)  #, debug=self.debug>1)  # this is kind of a crappy way to return the information, but I can't think of anything better a.t.m.
                if overlap_indel_fail:
                    self.indel_reruns.add(qname)
                    return dbgfcn('overlap/indel fails')
            elif overlap_status == 'nonsense':
                return dbgfcn('nonsense overlap bounds')
            else:
                assert overlap_status == 'ok'

        if len(qinfo['new_indels']) > 0:  # if any of the best matches had new indels this time through
            if qname in self.info['indels'] and qname not in self.vs_indels:  # if it already has a non-vsearch indel (i.e. from a previous sw round), just toss it (this used to assert that this wasn't the case, but it failed once in issue #309)
                self.indel_reruns.add(qname)
                return dbgfcn('indel fails')
            indelfo = self.combine_indels(qinfo, best)  # the next time through, when we're writing ig-sw input, we look to see if each query is in <self.info['indels']>, and if it is we pass ig-sw the indel-reversed sequence, rather than the <input_info> sequence
            if indelfo is None:
                self.indel_reruns.add(qname)
                return dbgfcn('indel fails')
            else:
                self.info['indels'][qinfo['name']] = indelfo
                self.indel_reruns.add(qname)
                return dbgfcn(' new indels in %s' % ' '.join(list(qinfo['new_indels'].keys())))  # utils.pad_lines(indelutils.get_dbg_str(self.info['indels'][qinfo['name']]), 10)

        if self.debug >= 2:  # keep this after the indel stuff, since it will crash (from different-length strings) if there's indels
            for region in utils.regions:
                for score, gene in qinfo['matches'][region]:  # sorted by decreasing match quality
                    self.print_match(region, gene, score, qseq, qinfo['glbounds'][gene], qinfo['qrbounds'][gene], skipping=False)

        if self.super_high_mutation(qinfo, best):
            return dbgfcn('super high mutation')

        # force v 5p and j 3p matches to (in most cases) go to the end of the sequence
        self.remove_probably_spurious_deletions(qinfo, best)

        # check for suspiciously bad annotations
        for rp in utils.region_pairs(self.glfo['locus']):
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
            return dbgfcn('negative cdr3 length %d (well, start and end codons overlap)' % cdr3_length)

        # convert to regular format used elsewhere, and add implicit info
        line = self.convert_qinfo(qinfo, best, codon_positions)
        try:
            utils.add_implicit_info(self.glfo, line, aligned_gl_seqs=self.aligned_gl_seqs, reset_indel_genes=True)
        except:
            elines = traceback.format_exception(*sys.exc_info())
            print(utils.pad_lines(''.join(elines)))
            print('      rerun: implicit info adding failed for %s (see above), rerunning' % qname)  # shouldn't be able to happen, so print even if debug isn't set
            return dbgfcn('see above')

        # deal with unproductive rearrangements
        if not utils.is_functional(line, iseq=0):
            if self.args.skip_unproductive:
                if self.debug:
                    print('      skipping unproductive (%s)' % utils.is_functional_dbg_str(line, iseq=0))
                self.skipped_unproductive_queries.add(qname)
                self.remaining_queries.remove(qname)
                return
        if self.args.skip_in_frame_rearrangements and line['cdr3_length'] % 3 == 0:  # NOTE *not* the same as line['in_frames'][0] (since here we're caring if the original rearrangement was productive, whereas 'in_frames' depends also on shm indels)
            if self.debug:
                print('      skipping in frame rearrangement')
            self.skipped_in_frame_queries.add(qname)
            self.remaining_queries.remove(qname)
            return

        kbounds = self.get_kbounds(line, qinfo, best)  # gets the boundaries of the non-best matches from <qinfo>
        if kbounds is None:
            return dbgfcn('nonsense kbounds')
        line['k_v'] = kbounds['v']
        line['k_d'] = kbounds['d']

        self.add_to_info(line)

    # ----------------------------------------------------------------------------------------
    def get_kbounds(self, line, qinfo, best, debug=False):
        # NOTE
        #  - k_v is index of first d/dj insert base (i.e. length of v match)
        #  - k_v + k_d is index of first j/dj insert base (i.e. length of v + vd insert + d match)
        if debug:
            print('%s for %s' % (utils.color('blue', 'get kbounds'), line['unique_ids'][0]))

        best_k_v = line['regional_bounds']['v'][1]  # end of v match
        best_k_d = line['regional_bounds']['d'][1] - line['regional_bounds']['v'][1]  # end of d minus end of v

        k_v_min, k_v_max = best_k_v, best_k_v
        k_d_min, k_d_max = best_k_d, best_k_d

        n_matched_to_break = 4  # once we see this many consecutive unmutated bases, assume we're actually within a germline v/d/j match
        typical_dj_insert_len = 5

        # first make sure the hmm will check for cases in which sw over-expanded v
        if debug:
            print('k_v_min --- %d' % k_v_min)
        n_matched = 0  # if <n_matched_to_break> bases match, assume we're definitely within the germline
        qrseq = line['fv_insertion'] + line['v_qr_seqs'][0]
        glseq = line['fv_insertion'] + line['v_gl_seq']
        if k_v_min > len(qrseq):
            if debug:
                print('k_v_min too big %d %d' % (k_v_min, len(qrseq)))
            return None
        icheck = k_v_min
        while icheck > line['codon_positions']['v'] + 3:  # i.e. stop when the last v base is the last base of the cysteine
            icheck -= 1
            if debug:
                print('    check %d' % icheck)
            if qrseq[icheck] == glseq[icheck]:
                n_matched += 1
                if debug:
                    print('      match number %d' % n_matched)
                if n_matched >= n_matched_to_break:
                    if debug:
                        print('      break on %d matches' % n_matched_to_break)
                    break
            else:  # set <k_v_min> to <icheck>, i.e. move first possible d/dj insert base to index <icheck>
                if debug:
                    print('      set k_v_min to %d' % icheck)
                k_v_min = icheck
                n_matched = 0

        # then check if sw might've eroded too much v or j
        if k_v_max < line['codon_positions']['v'] + 3:  # i.e. if first d/dj insert base (k_v_max) is within the three bases in the codon
            n_to_right_of_codon = len(self.glfo['seqs']['v'][line['v_gene']]) - utils.cdn_pos(self.glfo, 'v', line['v_gene'])  # number of bases to right of first position in codon
            if debug:
                print('k_v_max --- %d' % k_v_max)
                print('    set to length of v: %d' % (line['codon_positions']['v'] + n_to_right_of_codon))
            k_v_max = line['codon_positions']['v'] + n_to_right_of_codon  # extend all the way to end of v gene
        if k_v_min + k_d_min > line['codon_positions']['j']:  # i.e. if the first possible j base (k_v_min + k_d_min) is within the conserved codon
            n_to_left_of_codon = utils.cdn_pos(self.glfo, 'j', line['j_gene'])
            new_k_d_min = max(1, line['codon_positions']['j'] - n_to_left_of_codon - k_v_min - typical_dj_insert_len)
            if debug:
                print('k_d_min --- %d' % k_d_min)
                print('    k_v_min + k_d_min = %d + %d = %d > j cpos = %d, so set k_d_min to %d (so sum is at start of j, minus a bit)' % (k_v_min, k_d_min, k_v_min + k_d_min, line['codon_positions']['j'], new_k_d_min))
            k_d_min = new_k_d_min

        # then make sure the hmm will check for cases in which sw over-expanded j
        if debug:
            print('k_d_max --- %d' % k_d_max)
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
                print('    check %d' % icheck, end=' ')
            if k_v_min + icheck - j_start >= len(qrseq) or k_v_min + icheck - j_start >= len(glseq):  # shouldn't happen any more (min() in while statement should avoid it), but just in case
                print('    %s for query %s: k_v_min + icheck - j_start = %d + %d - %d = %d [ < 0 or >= len(qrseq) = %d   or   len(glseq) = %d]' % (utils.color('red', 'warning'), ' '.join(line['unique_ids']), k_v_min, icheck, j_start, k_v_min + icheck - j_start, len(qrseq), len(glseq)))
                return None
            if k_v_min + icheck < j_start:  # make sure we're at least as far as the start of the j (i.e. let the hmm arbitrate between d match and dj insertion)
                if debug:
                    print('      not yet to start of j (%d + %d < %d)' % (k_v_min, icheck, j_start))
                k_d_max = icheck + 1
            elif qrseq[k_v_min + icheck - j_start] == glseq[k_v_min + icheck - j_start]:
                n_matched += 1
                if debug:
                    print('      match number %d' % n_matched)
                if n_matched >= n_matched_to_break:
                    if debug:
                        print('      break on %d matches' % n_matched_to_break)
                    break
            else:  # set <k_d_max> to <icheck>, i.e. move first possible d/dj insert base to index <icheck>
                if debug:
                    print('      set k_d_max to %d' % (icheck + 1))
                k_d_max = icheck + 1  # i.e. if <icheck> doesn't match, we want the first d/dj base to be the *next* one (whereas for k_v, if <icheck> didn't match, we wanted <icheck> to be the first d/vd match)
                n_matched = 0

        assert k_v_min_CHK == k_v_min

        if debug:
            print('k_d_min --- %d' % k_d_min)
        k_v_max_CHK = k_v_max  # make sure we don't accidentally change it
        n_matched = 0  # if <n_matched_to_break> bases match, assume we're definitely within the germline
        qrseq, glseq = line['d_qr_seqs'][0], line['d_gl_seq']
        d_start = line['regional_bounds']['d'][0]
        icheck = k_d_min  # <icheck> is what we're considering changing k_d_min to
        while k_v_max < d_start and icheck > d_start - k_v_max + 1:  # i.e. icheck/k_d_min can't be smaller than the distance between d start and the furthest right possible v end (if k_v_max is larger than d_start, it's probably because we increased k_v_max above [a change which isn't reflected in d regional bounds in <line>])
            icheck -= 1
            if debug:
                print('    check %d' % icheck)
            if k_v_max + icheck - d_start < len(qrseq) and qrseq[k_v_max + icheck - d_start] == glseq[k_v_max + icheck - d_start]:  # note that it's <k_v_max> here, since even with the longest k_v we want to make sure to check as short of a d as we mean to
                n_matched += 1
                if debug:
                    print('      match number %d' % n_matched)
                if n_matched >= n_matched_to_break:
                    if debug:
                        print('      break on %d matches' % n_matched_to_break)
                    break
            else:  # set <k_d_min> to <icheck>, i.e. move first possible d/dj insert base to index <icheck>
                if debug:
                    print('      set k_d_min to %d' % icheck)
                k_d_min = icheck
                n_matched = 0

        assert k_v_max_CHK == k_v_max

        # NOTE I think there isn't any reason to increase k_v_max or decrease k_d_min -- sw already pretty much expands the v and j matches as far as it can (i.e. these fuzzing algorithms are mostly trying to decide if the last region to be matched by sw (d) should've been given more bases)
        # UPDATE not so sure, maybe I should, but it's not so important
        # SECOND UPDATE yep definitely need to expand k_v_max (and k_d_min), partly/mostly because we're now forbidding conserved codon deletion when building hmms (so now adding it above): added these above

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
                print('    %s %d %d' % (utils.color_gene(gene, width=15), k_v_min, k_v_max))

        # k_d
        tmpreg = 'd'
        n_d_genes = min(self.args.n_max_per_region[utils.regions.index(tmpreg)], len(qinfo['matches'][tmpreg]))  # NOTE don't need the n_max_per_region thing any more
        for igene in range(n_d_genes):
            _, gene = qinfo['matches'][tmpreg][igene]
            this_k_d = qinfo['qrbounds'][gene][1] - qinfo['qrbounds'][best['v']][1]  # end of d minus end of first/best v
            k_d_min = min(max(1, this_k_d), k_d_min)
            k_d_max = max(this_k_d, k_d_max)
            if debug:
                print('    %s %d %d' % (utils.color_gene(gene, width=15), k_d_min, k_d_max))

        # ----------------------------------------------------------------------------------------
        # switch to usual indexing conventions, i.e. that start with min and do not include max NOTE would be clearer to do this more coherently
        # NOTE i.e. k_[vd]_max means different things before here and after here
        # ----------------------------------------------------------------------------------------
        k_v_max += 1
        k_d_max += 1

        if not utils.has_d_gene(self.args.locus):
            k_v_min -= 3  # i'm not really sure why this can't be -1, but i think it has to do with coming up against the conserved tryp in j (whose erosion we now absolutely forbid) [this is getting fixed long after the k d modifications, i'm not sure why it didn't cause problems before]
            best_k_d = 1
            k_d_min = 1
            k_d_max = 2

        if best_k_v < k_v_min or best_k_v > k_v_max or best_k_d < k_d_min or best_k_d > k_d_max:
            if self.debug:
                print('  %s inconsistent best kset for %s (v: %d (%d %d)  d: %d (%d %d)' % (utils.color('red', 'error'), qinfo['name'], best_k_v, k_v_min, k_v_max, best_k_d, k_d_min, k_d_max))
            return None
        if k_v_min <= 0 or k_d_min <= 0 or k_v_min >= k_v_max or k_d_min >= k_d_max:
            if self.debug:
                print('  %s nonsense k bounds for %s (v: %d %d  d: %d %d)' % (utils.color('red', 'error'), qinfo['name'], k_v_min, k_v_max, k_d_min, k_d_max))
            return None
        if self.args.is_data and k_v_min - len(line['fv_insertion']) < 0:
            if self.debug:
                print('%s trimming fwk insertion would take k_v min to less than zero for %s: %d - %d = %d   %s' % (utils.color('yellow', 'warning'), ' '.join(line['unique_ids']), k_v_min, len(line['fv_insertion']), k_v_min - len(line['fv_insertion']), utils.reverse_complement_warning()))
            return None

        kbounds = {'v' : {'best' : best_k_v, 'min' : k_v_min, 'max' : k_v_max},
                   'd' : {'best' : best_k_d, 'min' : k_d_min, 'max' : k_d_max}}
        return kbounds

    # ----------------------------------------------------------------------------------------
    def remove_framework_insertions(self, debug=False):  # note that we modify non-implicit info here because it's faster than removing and re-adding it (although it makes this fcn more complicated)
        # NOTE duplicates code in utils.trim_fwk_insertions() (but they're enough different I don't want to deal with merging them, e.g. this has extra sw-specific keys, plus i care way more about speed here so don't really want to have to re-add implicit info)
        for query in self.info['queries']:
            swfo = self.info[query]
            assert len(swfo['seqs']) == 1

            if debug:
                print('  %-12s' % swfo['unique_ids'][0])
                print('     fv  %s' % utils.color('blue', swfo['fv_insertion']))
                print('     jf  %s' % utils.color('blue', swfo['jf_insertion']))

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
            swfo['fv_insertion'] = ''  # maybe don't really need these any more now that we've got leader_seqs and c_gene_seqs, but otoh they get used in lots of places (e.g. padding to same length) so it'd be really hard to remove them
            swfo['jf_insertion'] = ''
            swfo['codon_positions']['v'] -= fv_len
            swfo['codon_positions']['j'] -= fv_len
            for region in utils.regions:
                swfo['regional_bounds'][region] = tuple([rb - fv_len for rb in swfo['regional_bounds'][region]])  # I kind of want to just use a list now, but a.t.m. don't much feel like changing it everywhere else

            for region in utils.regions:
                if isinstance(swfo['all_matches'][0][region], list):  # if we just read an old sw cache file, it'll be a list of genes sorted by score, rather than a dict keyed by gene that includes scores and bounds, so there's no bounds to adjust
                    continue
                for gene, gfo in swfo['all_matches'][0][region].items():
                    gfo['qrbounds'] = tuple(b - fv_len for b in gfo['qrbounds'])

            if debug:
                print('    after %s' % swfo['seqs'][0])

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
        def get_pre_kept_queries():
            pre_kept_uids = set()
            if self.args.seed_unique_id is not None:
                pre_kept_uids.add(self.args.seed_unique_id)
            if self.args.queries is not None:
                pre_kept_uids |= set(self.args.queries)
            if self.args.queries_to_include is not None:  # note that the seed seq is added to queries_to_include, which is good if we're parameter caching since seed_unique_id will be set to None
                pre_kept_uids |= set(self.args.queries_to_include)
            info_queries = set(self.info['queries'])
            if len(pre_kept_uids - info_queries) > 0:
                print('  %s %d requested uid%s not in sw info: %s' % (utils.color('yellow', 'warning'), len(pre_kept_uids - info_queries), utils.plural(len(pre_kept_uids - info_queries)), ' '.join(pre_kept_uids - info_queries)))
            pre_kept_uids &= info_queries
            return pre_kept_uids
        # ----------------------------------------------------------------------------------------
        def process_seqs(uid_list, remove=False, lkseqs=None):
            for uid in uid_list:
                keyseq = getseq(uid) if lkseqs is None else lkseqs[uid]
                if keyseq in seqs_to_keep:
                    seqs_to_keep[keyseq].append(uid)
                    if remove:
                        self.remove_query(uid)
                        removed_queries.add(uid)
                else:
                    seqs_to_keep[keyseq] = [uid]
        # ----------------------------------------------------------------------------------------
        def get_long_seqs(tdbg=False):
            long_seqs, seq_classes = {}, {}  # <long_seqs>: map from each long/kept seq to its uid, <seq_classes>: map from each long/kept seq to the list of uids in its class
            for uid in self.info['queries']:
                useq = getseq(uid)
                found, switch = False, False
                if uid in pre_kept_uids:  # NOTE that if two pre-kept queries have the same seq, we'll just keep whichever one is last, which isn't really right but oh well
                    switch = True
                for lseq, lid in list(long_seqs.items()):
                    if self.info[lid]['cdr3_length'] != self.info[uid]['cdr3_length']:
                        continue
                    if useq in lseq:  # if lseq is longer (or they're the same), keep the one that's in there (lseq)
                        found = True
                        if uid in pre_kept_uids and len(useq) < len(lseq) and lid not in pre_kept_uids:
                            print('  %s pre-included query \'%s\' is being kept, but has shorter sequence than \'%s\', which we\'re marking as duplicate:\n    %s %s\n    %s %s' % (utils.color('yellow', 'warning'), uid, lid, useq, uid, lseq, lid))
                        break
                    elif lseq in useq:  # but useq is longer, we need to switch to useq
                        found = True
                        switch = True
                        break
                if found:
                    if switch:
                        long_seqs[useq] = uid
                        seq_classes[useq] = seq_classes[lseq] + [uid]
                        if tdbg:
                            print('  %s --> %s (%s)' % (lid, uid, ' '.join(seq_classes[useq])))
                        if lseq != useq:  # if they're the same this must be a pre-kept query
                            del long_seqs[lseq]
                            del seq_classes[lseq]
                    else:
                        seq_classes[lseq].append(uid)
                        if tdbg:
                            print('    add %s: %s' % (uid, ' '.join(seq_classes[lseq])))
                else:
                    if tdbg:
                        print('  new: %s' % uid)
                    assert useq not in long_seqs
                    long_seqs[useq] = uid
                    seq_classes[useq] = [uid]
            lkseqs = {u : lseq for lseq, uids in seq_classes.items() for u in uids}  # map from each uid to its 'keyseq', i.e. the longest seq that contains its seq
            return long_seqs, lkseqs

        # ----------------------------------------------------------------------------------------
        seqs_to_keep = {}  # seq : [uids that correspond to seq]

        # handle any pre-kept queries, just adding any duplicates to the appropriate list in <seqs_to_keep> *without* actually removing the duplicates
        pre_kept_uids = get_pre_kept_queries()  # set of seqs that we definitely keep, all as separate seqs (even if some are duplicates of each other, which occurs e.g. if there's duplicate sequences in --queries or --queries-to-include)
        process_seqs(pre_kept_uids)
        if len(seqs_to_keep) < len(pre_kept_uids):
            print('  %s keeping duplicate sequences in pre-kept queries: %s' % (utils.color('yellow', 'warning'), ',   '.join(' '.join(uids) for uids in seqs_to_keep.values() if len(uids) > 1)))

        lkseqs, long_seqs = None, {}
        if self.args.also_remove_duplicate_sequences_with_different_lengths:
            long_seqs, lkseqs = get_long_seqs()
            process_seqs(list(long_seqs.values()))  # I'm not setting remove=True, but there shouldn't/can't be any duplicates anyway

        removed_queries = set()
        process_seqs(set(self.info['queries']) - pre_kept_uids - set(long_seqs.values()), remove=True, lkseqs=lkseqs)

        any_have_mtpy = any('multiplicities' in self.info[u] for u in self.info['queries'])  # when we read input meta info, we only fill it for seqs for which a key is present, which is fine, *except* for multiplicities, since later on (in utils.get_multiplicity()) we look in a multi-sequence annotation to see if 'multiplicities' is there, so if it's there for one it has to also be correct for all of them
        for seq, uids in seqs_to_keep.items():
            kept_uid = uids[0]
            previous_duplicates = set(self.info[kept_uid]['duplicates'][0])  # probably from previous waterer run
            new_duplicates = set(uids[1:]) - pre_kept_uids  # don't actually add as duplicates uids that are in <pre_kept_uids>, since these will not have been removed from <self.info>
            all_duplicates = list(previous_duplicates | new_duplicates)
            self.info[kept_uid]['duplicates'][0] = all_duplicates
            self.duplicates[kept_uid] = all_duplicates  # copy info from previous line to <self.duplicates>, which is just so partitiondriver can pass in previous duplicates, and yes having the info in two places is dumb
            if any_have_mtpy or len(all_duplicates) > 0:  # we *also* have to add it if there's any duplicates, since some of those duplicate sequences might have multiplicities, which would then get lost if no other sequences have multiplicities (this is admittedly very rare)
                def getmult(ltmp): return ltmp['multiplicities'][0] if 'multiplicities' in ltmp else utils.input_metafile_defaults('multiplicities')
                self.info[kept_uid]['multiplicities'] = [sum(getmult(self.input_info[u]) for u in [kept_uid] + all_duplicates)]

        if len(removed_queries) > 0:
            print('      removed %d / %d = %.2f duplicate sequences after trimming framework insertions (leaving %d)' % (len(removed_queries), len(removed_queries) + len(self.info['queries']), len(removed_queries) / float(len(removed_queries) + len(self.info['queries'])), len(self.info['queries'])))

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
            print('  maxima:', end=' ')
            for k in padnames:
                print('%s %d    ' % (k, maxima[k]), end=' ')
            print('')

            print('    per-cdr3 maxima:')
            print('         %s  %s  %s' % ('cdr3', padnames[0], padnames[1]))
            for cdr3 in per_cdr3_maxima:
                print('         %3d' % cdr3, end=' ')
                for k in padnames:
                    print('   %d    ' % (per_cdr3_maxima[cdr3][k]), end=' ')
                print('')

        return maxima, per_cdr3_maxima

    # ----------------------------------------------------------------------------------------
    def pad_seqs_to_same_length(self, debug=False):
        """
        Pad all sequences in <seqinfo> to the same length to the left and right of their conserved cysteine positions.
        Next, pads all sequences further out (if necessary) such as to eliminate all v_5p and j_3p deletions.
        NOTE that this gets run *after* reading or writing cache file, i.e. the cache file has *un*padded seqs, while the padded sw annotations/seqs only show up in-memory en route to hmm input
        """

        cluster_different_cdr3_lengths = False  # if you want glomerator.cc to try to cluster different cdr3 lengths, you need to pass it *everybody* with the same N padding... but then you're padding way more than you need to on almost every sequence, which is really wasteful and sometimes confuses bcrham

        if debug:
            print('padding %d seqs to same length (%s cdr3 length classes)' % (len(self.info['queries']), 'within' if not cluster_different_cdr3_lengths else 'merging'))

        maxima, per_cdr3_maxima = self.get_padding_parameters(debug=debug)

        if debug:
            print('    left  right    uid')
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

            leftstr = padleft * utils.ambig_base
            rightstr = padright * utils.ambig_base
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

            for region in utils.regions:
                if isinstance(swfo['all_matches'][0][region], list):  # if we just read an old sw cache file, it'll be a list of genes sorted by score, rather than a dict keyed by gene that includes scores and bounds, so there's no bounds to adjust
                    continue
                for gene, gfo in swfo['all_matches'][0][region].items():
                    gfo['qrbounds'] = tuple(b + padleft for b in gfo['qrbounds'])

            swfo['padlefts'] = [padleft, ]
            swfo['padrights'] = [padright, ]
            if debug:
                print('    %3d   %3d    %s' % (padleft, padright, query))

        if debug:
            print('    cdr3        uid                 padded seq')
            for query in sorted(self.info['queries'], key=lambda q: self.info[q]['cdr3_length']):
                print('    %3d   %20s    %s' % (self.info[query]['cdr3_length'], query, self.info[query]['seqs'][0]))

    # ----------------------------------------------------------------------------------------
    def combine_indels(self, qinfo, best, debug=False):
        # debug = 2
        if debug:
            print('  %s: combine with %d sw indels: %s' % (qinfo['name'], len(qinfo['new_indels']), ' '.join(list(qinfo['new_indels'].keys()))))
            for ifo in qinfo['new_indels'].values():
                print(indelutils.get_dbg_str(ifo))

        regional_indelfos = {}
        qrbounds = {r : qinfo['qrbounds'][best[r]] for r in utils.regions}
        full_qrseq = qinfo['seq']
        if self.vs_info is not None and qinfo['name'] in self.vs_indels:
            if debug:
                print('    has a vsearch v indel')

            if 'v' in qinfo['new_indels']:  # if sw kicks up an additional v indel that vsearch didn't find, we rerun sw with <self.args.no_indel_gap_open_penalty>
                if debug:
                    print('      sw called a v indel, but there\'s already a vsearch v indel, so give up (delete vsearch indel, ignore new sw indel, and rerun sw forbidding indels)')
                if qinfo['name'] in self.info['indels']:
                    del self.info['indels'][qinfo['name']]
                return None

            vs_indelfo = self.info['indels'][qinfo['name']]
            assert 'v' in vs_indelfo['genes']  # a.t.m. vsearch is only looking for v genes, and if that changes in the future we'd need to rewrite this

            if any(ifo['pos'] >= qrbounds['v'][1] for ifo in vs_indelfo['indels']):  # if any of the vsearch indels are to the right of the sw v bounds, we want to ignore them all (well, really I guess we only want to ignore the ones that're to the right of the sw bounds, but this is so incredibly rare that I think this is ok)
                if self.debug:
                    print('    ignoring vsearch v indel that\'s to the right of the sw v bounds')
            else:
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
                regional_indelfos['v'] = vs_indelfo

            full_qrseq = self.input_info[qinfo['name']]['seqs'][0]  # should in principle replace qinfo['seq'] as well, since it's the reversed seq from vsearch, but see note below
            del self.info['indels'][qinfo['name']]
        elif 'v' in qinfo['new_indels']:
            regional_indelfos['v'] = qinfo['new_indels']['v']

        if 'd' in qinfo['new_indels']:
            regional_indelfos['d'] = qinfo['new_indels']['d']
        if 'j' in qinfo['new_indels']:
            regional_indelfos['j'] = qinfo['new_indels']['j']

        # NOTE qinfo won't be consistent with the indel reversed seq after this, but that's kind of the point, since we're just rerunning anyway

        return indelutils.combine_indels(regional_indelfos, full_qrseq, qrbounds, uid=qinfo['name'], debug=debug)
