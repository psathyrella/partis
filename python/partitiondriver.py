import json
import numpy
import time
import sys
import itertools
import math
import os
import glob
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import random
from collections import OrderedDict
from subprocess import Popen, check_call, PIPE, CalledProcessError, check_output
import copy
import multiprocessing
import operator

import utils
import glutils
import seqfileopener
from glomerator import Glomerator
from clusterpath import ClusterPath
from waterer import Waterer
from parametercounter import ParameterCounter
from performanceplotter import PerformancePlotter
from partitionplotter import PartitionPlotter
from hist import Hist

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    """ Class to parse input files, start bcrham jobs, and parse/interpret bcrham output for annotation and partitioning """
    def __init__(self, args, action, initial_gldir):  # NOTE <initial_gldir> is not, in general, the same as <args.initial_germline_dir>
        self.args = args
        self.current_action = action  # *not* necessarily the same as <self.args.action>
        utils.prep_dir(self.args.workdir)
        self.my_gldir = self.args.workdir + '/' + glutils.glfo_dir
        self.glfo = glutils.read_glfo(initial_gldir, locus=self.args.locus, only_genes=self.args.only_genes)
        self.simglfo = self.glfo
        if self.args.simulation_germline_dir is not None:
            self.simglfo = glutils.read_glfo(self.args.simulation_germline_dir, locus=self.args.locus)  # NOTE uh, I think I don't want to apply <self.args.only_genes>
        glutils.write_glfo(self.my_gldir, self.glfo)  # need a copy on disk for vdjalign and bcrham (note that what we write to <self.my_gldir> in general differs from what's in <initial_gldir>)

        self.input_info, self.reco_info = None, None
        if self.args.infname is not None:
            self.input_info, self.reco_info = seqfileopener.get_seqfile_info(self.args.infname, self.args.is_data, n_max_queries=self.args.n_max_queries, args=self.args, glfo=self.glfo, simglfo=self.simglfo)
            if len(self.input_info) > 1000:
                if self.args.n_procs == 1:
                    print '  note:! running on %d sequences spread over %d processes. This will be kinda slow, so it might be a good idea to set --n-procs N to the number of processors on your local machine, or look into non-local parallelization with --batch-system.\n' % (len(self.input_info), self.args.n_procs)
                if self.args.outfname is None and self.current_action != 'cache-parameters':
                    print '  note: running on a lot of sequences without setting --outfname. Which is ok! But there\'ll be no persistent record of the results'
            self.default_sw_cachefname = self.args.parameter_dir + '/sw-cache-' + repr(abs(hash(''.join(self.input_info.keys())))) + '.csv'  # maybe I shouldn't abs it? collisions are probably still unlikely, and I don't like the extra dash in my file name
        elif self.current_action != 'view-annotations' and self.current_action != 'view-partitions' and self.current_action != 'view-alternative-naive-seqs':
            raise Exception('--infname is required for action \'%s\'' % args.action)

        self.sw_info = None
        self.duplicates = {}
        self.bcrham_proc_info = None
        self.timing_info = []  # TODO clean up this and bcrham_proc_info

        self.unseeded_seqs = None  # all the queries that we *didn't* cluster with the seed uid
        self.small_cluster_seqs = None  # all the queries that we removed after a few partition steps 'cause they were in small clusters

        self.all_new_allele_info = []

        self.sw_param_dir = self.args.parameter_dir + '/sw'
        self.hmm_param_dir = self.args.parameter_dir + '/hmm'
        self.sub_param_dir = self.args.parameter_dir + '/' + self.args.parameter_type

        self.hmm_infname = self.args.workdir + '/hmm_input.csv'
        self.hmm_cachefname = self.args.workdir + '/hmm_cached_info.csv'
        self.hmm_outfname = self.args.workdir + '/hmm_output.csv'
        self.annotation_fname = self.hmm_outfname.replace('.csv', '_annotations.csv')

        if self.args.outfname is not None:
            utils.prep_dir(dirname=None, fname=self.args.outfname, allow_other_files=True)

        self.deal_with_persistent_cachefile()

        self.cached_naive_hamming_bounds = self.args.naive_hamming_bounds  # just so we don't get them every iteration through the clustering loop

        self.aligned_gl_seqs = None
        if self.args.aligned_germline_fname is not None:
            self.aligned_gl_seqs = glutils.read_aligned_gl_seqs(self.args.aligned_germline_fname, self.glfo)

    # ----------------------------------------------------------------------------------------
    def clean(self):
        glutils.remove_glfo_files(self.my_gldir, self.args.locus)

        # merge persistent and current cache files into the persistent cache file
        if self.args.persistent_cachefname is not None:
            lockfname = self.args.persistent_cachefname + '.lock'
            while os.path.exists(lockfname):
                print '  waiting for lock on %s' % lockfname
                time.sleep(0.5)
            lockfile = open(lockfname, 'w')
            if not os.path.exists(self.args.persistent_cachefname):
                open(self.args.persistent_cachefname, 'w').close()
            self.merge_files(infnames=[self.args.persistent_cachefname, self.hmm_cachefname], outfname=self.args.persistent_cachefname, dereplicate=True)
            lockfile.close()
            os.remove(lockfname)
        if os.path.exists(self.hmm_cachefname):
            os.remove(self.hmm_cachefname)

        try:
            os.rmdir(self.args.workdir)
        except OSError:
            raise Exception('workdir (%s) not empty: %s' % (self.args.workdir, ' '.join(os.listdir(self.args.workdir))))  # hm... you get weird recursive exceptions if you get here. Oh, well, it still works

    # ----------------------------------------------------------------------------------------
    def deal_with_persistent_cachefile(self):
        if self.args.persistent_cachefname is None or not os.path.exists(self.args.persistent_cachefname):  # nothin' to do (ham'll initialize it)
            return

        with open(self.args.persistent_cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            if set(reader.fieldnames) == set(utils.annotation_headers):
                raise Exception('doesn\'t work yet')
                print '  parsing annotation output file %s to partition cache file %s' % (self.args.persistent_cachefname, self.hmm_cachefname)
                with open(self.hmm_cachefname, 'w') as outcachefile:
                    writer = csv.DictWriter(outcachefile, utils.partition_cachefile_headers)
                    writer.writeheader()
                    for line in reader:
                        if line['v_gene'] == '':  # failed
                            continue
                        utils.process_input_line(line)
                        outrow = {'unique_ids' : line['unique_ids'], 'naive_seq' : line['padlefts'][0] * utils.ambiguous_bases[0] + line['naive_seq'] + line['padrights'][0] * utils.ambiguous_bases[0]}
                        writer.writerow(outrow)
            elif set(reader.fieldnames) == set(utils.partition_cachefile_headers):  # headers are ok, so can just copy straight over
                check_call(['cp', self.args.persistent_cachefname, self.hmm_cachefname])
            else:
                raise Exception('--persistent-cachefname %s has unexpected header list %s' % (self.args.persistent_cachefname, reader.fieldnames))

    # ----------------------------------------------------------------------------------------
    def run_waterer(self, count_parameters=False, write_parameters=False, remove_less_likely_alleles=False, find_new_alleles=False, itry=None, write_cachefile=False, look_for_cachefile=False):
        print 'smith-waterman',
        if write_parameters:
            print '  (writing parameters)',
        if find_new_alleles:
            print '  (looking for new alleles)',
        if remove_less_likely_alleles:
            print '  (removing less-likely alleles)',
        print ''
        sys.stdout.flush()

        # can probably remove this... I just kind of want to know if it happens EDIT hell no don't remove this
        if not count_parameters and not find_new_alleles and not remove_less_likely_alleles and os.path.exists(self.sub_param_dir + '/hmms'):
            genes_with_hmms = set(utils.find_genes_that_have_hmms(self.sub_param_dir))
            expected_genes = set([g for r in utils.regions for g in self.glfo['seqs'][r].keys()])  # this'll be the & of the gldir (maybe rewritten, maybe not)
            if self.args.only_genes is None and len(genes_with_hmms - expected_genes) > 0:
                print '  %s yamels in %s for %d genes that aren\'t in glfo' % (utils.color('red', 'warning'), self.sub_param_dir, len(genes_with_hmms - expected_genes))
            if len(expected_genes - genes_with_hmms) > 0:
                print '  %s %d genes in glfo that don\'t have yamels in %s' % (utils.color('red', 'warning'), len(expected_genes - genes_with_hmms), self.sub_param_dir)

        pre_failed_queries = self.sw_info['failed-queries'] if self.sw_info is not None else None  # don't re-run on failed queries if this isn't the first sw run (i.e., if we're parameter caching)
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.glfo,
                          count_parameters=count_parameters,  # (remove_less_likely_alleles or parameter_out_dir is not None),
                          parameter_out_dir=self.sw_param_dir if write_parameters else None,
                          remove_less_likely_alleles=remove_less_likely_alleles,
                          find_new_alleles=find_new_alleles,
                          plot_performance=(self.args.plot_performance and not remove_less_likely_alleles and not find_new_alleles),
                          simglfo=self.simglfo, itry=itry, duplicates=self.duplicates, pre_failed_queries=pre_failed_queries)
        cachefname = self.default_sw_cachefname if self.args.sw_cachefname is None else self.args.sw_cachefname
        if not look_for_cachefile and os.path.exists(cachefname):  # i.e. if we're not explicitly told to look for it, and it's there, then it's probably out of date
            print '  removing old sw cache %s' % cachefname.replace('.csv', '')
            os.remove(cachefname)
            if os.path.exists(cachefname.replace('.csv', '-glfo')):  # it should always be there now, but there could be some old sw cache files lying around from before they had their own glfo dirs
                glutils.remove_glfo_files(cachefname.replace('.csv', '-glfo'), self.args.locus)
        if look_for_cachefile and os.path.exists(cachefname):  # run sw if we either don't want to do any caching (None) or if we are planning on writing the results after we run
            waterer.read_cachefile(cachefname)
        else:
            waterer.run(cachefname if write_cachefile else None)

        self.sw_info = waterer.info
        for uid, dupes in waterer.duplicates.items():  # <waterer.duplicates> is <self.duplicates> OR'd into any new duplicates from this run
            self.duplicates[uid] = dupes

    # ----------------------------------------------------------------------------------------
    def find_new_alleles(self):
        """ look for new alleles with sw, write any that you find to the germline set directory in <self.workdir>, add them to <self.glfo>, and repeat until you don't find any. """
        assert len(self.all_new_allele_info) == 0
        # alleles_with_evidence = set()
        itry = 0
        while True:
            self.run_waterer(find_new_alleles=True, itry=itry)
            if len(self.sw_info['new-alleles']) == 0:
                break
            if os.path.exists(self.default_sw_cachefname):
                print '    removing sw cache file %s (it has outdated germline info)' % self.default_sw_cachefname
                os.remove(self.default_sw_cachefname)

            self.all_new_allele_info += [afo for afo in self.sw_info['new-alleles'] if '+' in afo['gene']]  # i.e. skip new/inferred alleles that turned out to be previously known (previous: in original glfo)
            # alleles_with_evidence |= self.sw_info['alleles-with-evidence']
            glutils.restrict_to_genes(self.glfo, list(self.sw_info['all_best_matches']))
            glutils.add_new_alleles(self.glfo, self.sw_info['new-alleles'])
            # if self.args.generate_germline_set:
            #     for alfo in self.sw_info['new-alleles']:
            #         if alfo['template-gene'] not in alleles_with_evidence:  # XXX [update comment] if the new allele is actually new (i.e. not in imgt), and if we never had explicit evidence for the template gene (i.e. it was just the best match we had) then remove the template gene
            #             print '    removing template gene %s' % utils.color_gene(alfo['template-gene'])
            #             glutils.remove_gene(self.glfo, alfo['template-gene'])
            glutils.write_glfo(self.my_gldir, self.glfo)  # write glfo modifications to disk

            itry += 1
            if itry >= self.args.n_max_allele_finding_iterations:
                break

        if self.args.new_allele_fname is not None:
            n_new_alleles = len(self.all_new_allele_info)
            print '  writing %d new %s to %s' % (n_new_alleles, utils.plural_str('allele', n_new_alleles), self.args.new_allele_fname)
            with open(self.args.new_allele_fname, 'w') as outfile:
                for allele_info in self.all_new_allele_info:
                    outfile.write('>%s\n' % allele_info['gene'])
                    outfile.write('%s\n' % allele_info['seq'])

    # ----------------------------------------------------------------------------------------
    def restrict_to_observed_alleles(self, subpdir):
        # TODO do I still need this now I'm using alleleremover?
        """ Restrict <self.glfo> to genes observed in <subpdir>, and write the changes to <self.my_gldir>. """
        if self.args.debug:
            print '  restricting self.glfo (and %s) to alleles observed in %s' % (self.my_gldir, subpdir)
        only_genes = set()
        for region in utils.regions:
            with open(subpdir + '/' + region + '_gene-probs.csv', 'r') as pfile:
                reader = csv.DictReader(pfile)
                for line in reader:
                    only_genes.add(line[region + '_gene'])
        glutils.restrict_to_genes(self.glfo, only_genes, debug=False)
        glutils.write_glfo(self.my_gldir, self.glfo, debug=False)  # write glfo modifications to disk

    # ----------------------------------------------------------------------------------------
    def cache_parameters(self):
        """ Infer full parameter sets and write hmm files for sequences from <self.input_info>, first with Smith-Waterman, then using the SW output as seed for the HMM """
        print 'caching parameters'
        if len(self.input_info) < 10 * self.args.min_observations_to_write:
            print """
            %s: number of input sequences (%d) isn\'t very large compared to --min-observations-to-write (%d), i.e. when we write hmm files we\'re going to be doing a lot of interpolation and smoothing.
            This is not necessarily terrible -- if you really only have %d input sequences, you will, in general, get sensible answers.
            But:
              - if this is a subset of a larger file, you should cache parameters using the entire file (well, a random subset of, say, ~50k sequences is typically sufficient)
              - if you have another data set that you believe is similar to this small one (e.g. same human, so the germlines are the same), you will get more accurate results if you cache parameters on that data set, and then run inference using those parameters (i.e. use the --parameter-dir argument) on this small data set.

            For now, we assume the first case (you actually want to infer parameters on this small data set), so we reset --min-observations-to-write to 1 and charge ahead.
            It would also be sensible to compare the hmm output results (--outfname) to the smith-waterman results (which are cached in a file whose path should be printed just below).
            """ % (utils.color('red', 'warning'), len(self.input_info), self.args.min_observations_to_write, len(self.input_info), )

            self.args.min_observations_to_write = 1

        if not self.args.dont_remove_unlikely_alleles:
            self.run_waterer(remove_less_likely_alleles=True, count_parameters=True)
            glutils.remove_genes(self.glfo, self.sw_info['genes-to-remove'])
            glutils.write_glfo(self.my_gldir, self.glfo)
        if self.args.find_new_alleles:
            self.find_new_alleles()
        self.run_waterer(count_parameters=True, write_parameters=True, write_cachefile=True)
        self.restrict_to_observed_alleles(self.sw_param_dir)
        self.write_hmms(self.sw_param_dir)
        if self.args.only_smith_waterman:
            return

        print 'hmm'
        self.run_hmm('viterbi', parameter_in_dir=self.sw_param_dir, parameter_out_dir=self.hmm_param_dir, count_parameters=True)
        self.restrict_to_observed_alleles(self.hmm_param_dir)
        self.write_hmms(self.hmm_param_dir)

        if len(self.all_new_allele_info) > 0:
            print '  %d new allele%s written to {sw,hmm}/%s subdirs of parameter dir %s' % (len(self.all_new_allele_info), utils.plural_str('', len(self.all_new_allele_info)), glutils.glfo_dir, self.args.parameter_dir)

    # ----------------------------------------------------------------------------------------
    def run_algorithm(self, algorithm):
        """ Just run <algorithm> (either 'forward' or 'viterbi') on sequences in <self.input_info> and exit. You've got to already have parameters cached in <self.args.parameter_dir> """
        print 'running %s' % algorithm
        self.run_waterer(look_for_cachefile=True)
        if self.args.only_smith_waterman:
            return
        print 'hmm'
        self.run_hmm(algorithm, parameter_in_dir=self.sub_param_dir)

    # ----------------------------------------------------------------------------------------
    def view_existing_annotations(self):
        with open(self.args.outfname) as csvfile:
            failed_queries = set()
            reader = csv.DictReader(csvfile)
            n_queries_read = 0
            for line in reader:
                if line['v_gene'] == '':
                    failed_queries.add(line['unique_ids'])
                    continue
                if self.args.queries is not None and self.args.queries[0] not in line['unique_ids']:  # check that at least the first requested query is present in the unconverted unique_ids string, since conversion is slow
                    continue
                utils.process_input_line(line)
                if self.args.queries is not None and len(set(self.args.queries) & set(line['unique_ids'])) == 0:  # actually make sure this is the precise set of queries we want (note that --queries and line['unique_ids'] are both ordered, and this ignores that... oh, well, sigh.)
                    continue
                if self.args.infname is not None and self.reco_info is not None:
                    utils.print_true_events(self.glfo, self.reco_info, line, extra_str='')
                utils.add_implicit_info(self.glfo, line)
                print 'inferred:',
                if len(line['unique_ids']) > 1:
                    print '   %s' % ':'.join(line['unique_ids'])
                else:
                    print ''
                utils.print_reco_event(line, extra_str='  ')
                n_queries_read += 1
                if self.args.n_max_queries > 0 and n_queries_read >= self.args.n_max_queries:
                    break
        if len(failed_queries) > 0:
            print '\n%d failed queries' % len(failed_queries)

    # ----------------------------------------------------------------------------------------
    def view_existing_partitions(self):
        cp = ClusterPath()
        cp.readfile(self.args.outfname)
        cp.print_partitions(abbreviate=self.args.abbreviate, reco_info=self.reco_info)

    # ----------------------------------------------------------------------------------------
    def view_alternative_naive_seqs(self):
        if self.args.queries is None:
            print '%s in order to view alternative naive sequences, you have to specify a set of uids in which you\'re interested. Choose something from here:' % utils.color('red', 'error')
            cp = ClusterPath()
            cp.readfile(self.args.outfname)
            cp.print_partitions(abbreviate=self.args.abbreviate, reco_info=self.reco_info)
            raise Exception('see above')
        self.print_subcluster_naive_seqs(self.args.queries)

    # ----------------------------------------------------------------------------------------
    def partition(self):
        """ Partition sequences in <self.input_info> into clonally related lineages """
        print 'partitioning'
        self.run_waterer(look_for_cachefile=True)  # run smith-waterman
        if len(self.sw_info['queries']) == 0:
            if self.args.outfname is not None:
                check_call(['touch', self.args.outfname])
            return
        if self.args.only_smith_waterman:
            return

        print 'hmm'
        # cache hmm naive seq for each single query NOTE <self.current_action> is (and needs to be) still set to partition for this
        if not self.args.dont_precache_naive_seqs and (len(self.sw_info['queries']) > 50 or self.args.naive_vsearch or self.args.naive_swarm):
            print '--> caching all %d naive sequences' % len(self.sw_info['queries'])
            self.run_hmm('viterbi', self.sub_param_dir, n_procs=self.auto_nprocs(len(self.sw_info['queries'])), precache_all_naive_seqs=True)

        if self.args.naive_vsearch or self.args.naive_swarm:
            cpath = self.cluster_with_naive_vsearch_or_swarm(parameter_dir=self.sub_param_dir)
        else:
            cpath = self.cluster_with_bcrham()

        cluster_annotations = self.get_cluster_annotations(cpath.partitions[cpath.i_best])

        if self.args.plotdir is not None:
            partplotter = PartitionPlotter()
            partplotter.plot(self.args.plotdir + '/partitions', partition=cpath.partitions[cpath.i_best], annotations=cluster_annotations, only_csv=self.args.only_csv_plots)

        if self.args.debug:
            print 'final'
            cpath.print_partitions(self.reco_info, print_header=True, calc_missing_values='all' if (len(self.input_info) < 500) else 'best')
            if not self.args.is_data:
                true_cp = ClusterPath(seed_unique_id=self.args.seed_unique_id)
                true_cp.add_partition(utils.get_true_partition(self.reco_info), -1., 1)
                print 'true:'
                true_cp.print_partitions(self.reco_info, print_header=False, calc_missing_values='best')

        self.check_partition(cpath.partitions[cpath.i_best])
        if self.args.outfname is not None:
            self.write_clusterpaths(self.args.outfname, cpath)  # [last agglomeration step]

    # ----------------------------------------------------------------------------------------
    def split_seeded_clusters(self, old_cpath):
        seeded_clusters, unseeded_clusters = utils.split_partition_with_criterion(old_cpath.partitions[old_cpath.i_best_minus_x], lambda cluster: self.args.seed_unique_id in cluster)
        self.unseeded_seqs = [uid for uclust in unseeded_clusters for uid in uclust]  # they should actually all be singletons, since we shouldn't have merged anything that wasn't seeded
        if len(unseeded_clusters) != len(self.unseeded_seqs):
            print '%s unseeded clusters not all singletons' % utils.color('red', 'warning')
        seeded_cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
        seeded_singleton_set = set([uid for sclust in seeded_clusters for uid in sclust])  # in case there's duplicates
        seeded_singletons = [[uid, ] for uid in seeded_singleton_set]
        seeded_cpath.add_partition(seeded_singletons, -1., 1)
        print '      removing %d sequences in unseeded clusters, and splitting %d seeded clusters into %d singletons' % (len(self.unseeded_seqs), len(seeded_clusters), len(seeded_cpath.partitions[seeded_cpath.i_best_minus_x]))

        return seeded_cpath

    # ----------------------------------------------------------------------------------------
    def remove_small_clusters(self, old_cpath):
        big_clusters, small_clusters = utils.split_partition_with_criterion(old_cpath.partitions[old_cpath.i_best_minus_x], lambda cluster: len(cluster) not in self.args.small_clusters_to_ignore)
        self.small_cluster_seqs = [sid for sclust in small_clusters for sid in sclust]
        new_cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
        new_cpath.add_partition(big_clusters, -1., 1)
        print '      removing %d sequences in %d small clusters (with sizes among %s)' % (len(self.small_cluster_seqs), len(small_clusters), ' '.join([str(sz) for sz in self.args.small_clusters_to_ignore]))
        return new_cpath

    # ----------------------------------------------------------------------------------------
    def scale_n_procs_for_new_n_clusters(self, initial_nseqs, initial_nprocs, cpath):
        new_n_clusters = len(cpath.partitions[cpath.i_best_minus_x])  # when removing small clusters, this is the number of clusters, not the number of sequences, but it's maybe still ok
        int_initial_seqs_per_proc = max(1, int(float(initial_nseqs) / initial_nprocs))
        new_n_procs = max(1, int(float(new_n_clusters) / int_initial_seqs_per_proc))
        if new_n_clusters > 20:
            new_n_procs *= 3  # multiply by something 'cause we're turning off the seed uid for the last few times through
        if self.args.batch_system is None:
            new_n_procs = min(new_n_procs, multiprocessing.cpu_count())
        new_n_procs = min(new_n_procs, self.args.n_procs)  # don't let it be bigger than whatever was initially specified
        print '        new n_procs %d (initial seqs/proc: %.2f   new seqs/proc: %.2f' % (new_n_procs, float(initial_nseqs) / initial_nprocs, float(new_n_clusters) / new_n_procs)
        return new_n_procs

    # ----------------------------------------------------------------------------------------
    def shall_we_reduce_n_procs(self, last_n_procs, n_proc_list):
        if self.timing_info[-1]['total'] < self.args.min_hmm_step_time:  # mostly for when you're running on really small samples
            return True
        n_calcd_per_process = self.get_n_calculated_per_process()
        if n_calcd_per_process < self.args.n_max_to_calc_per_process and last_n_procs > 2:  # should be replaced by time requirement, since especially in later iterations, the larger clusters make this a crappy metric (2 is kind of a special case, becase, well, small integers and all)
            return True
        times_to_try_this_n_procs = max(3, last_n_procs)  # if we've already milked this number of procs for most of what it's worth (once you get down to 2 or 3, you don't want to go lower)
        if n_proc_list.count(last_n_procs) >= times_to_try_this_n_procs:
            return True

        return False

    # ----------------------------------------------------------------------------------------
    def prepare_next_iteration(self, n_proc_list, cpath, initial_nseqs):
        last_n_procs = n_proc_list[-1]
        next_n_procs = last_n_procs

        factor = 1.3
        if self.shall_we_reduce_n_procs(last_n_procs, n_proc_list):
            next_n_procs = int(next_n_procs / float(factor))

        def time_to_remove_some_seqs(n_proc_threshold):
            return len(n_proc_list) >= n_proc_threshold or next_n_procs == 1

        if self.args.small_clusters_to_ignore is not None and self.small_cluster_seqs is None and time_to_remove_some_seqs(self.args.n_steps_after_which_to_ignore_small_clusters):
            cpath = self.remove_small_clusters(cpath)
            next_n_procs = self.scale_n_procs_for_new_n_clusters(initial_nseqs, n_proc_list[0], cpath)
        if self.args.seed_unique_id is not None and self.unseeded_seqs is None and time_to_remove_some_seqs(3):  # if we didn't already remove the unseeded clusters in a partition previous step
            cpath = self.split_seeded_clusters(cpath)
            next_n_procs = self.scale_n_procs_for_new_n_clusters(initial_nseqs, n_proc_list[0], cpath)

        return next_n_procs, cpath

    # ----------------------------------------------------------------------------------------
    def get_n_calculated_per_process(self):
        assert self.bcrham_proc_info is not None

        total = 0.  # sum over each process
        for procinfo in self.bcrham_proc_info:
            if 'vtb' not in procinfo['calcd'] or 'fwd' not in procinfo['calcd']:
                print '%s couldn\'t find vtb/fwd in:\n%s' % (utils.color('red', 'warning'), procinfo['calcd'])  # may as well not fail, it probably just means we lost some stdout somewhere. Which, ok, is bad, but let's say it shouldn't be fatal.
                return 1.  # er, or something?
            if self.args.naive_hamming:  # make sure we didn't accidentally calculate some fwds
                assert procinfo['calcd']['fwd'] == 0.
            total += procinfo['calcd']['vtb'] + procinfo['calcd']['fwd']
        if self.args.debug:
            print '          vtb + fwd calcd: %d (%.1f per proc)' % (total, float(total) / len(self.bcrham_proc_info))
        return float(total) / len(self.bcrham_proc_info)

    # ----------------------------------------------------------------------------------------
    def are_we_finished_clustering(self, n_procs, cpath):
        if n_procs == 1:
            return True
        elif self.args.n_final_clusters is not None and len(cpath.partitions[cpath.i_best]) <= self.args.n_final_clusters:  # NOTE I *think* I want the best, not best-minus-x here (hardish to be sure a.t.m., since I'm not really using the minus-x part right now)
            print '  stopping with %d (<= %d) clusters' % (len(cpath.partitions[cpath.i_best]), self.args.n_final_clusters)
            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------
    def cluster_with_bcrham(self):
        n_procs = self.args.n_procs
        cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
        initial_nsets = [[q, ] for q in self.sw_info['queries']]
        cpath.add_partition(initial_nsets, logprob=0., n_procs=n_procs)  # NOTE sw info excludes failed sequences (and maybe also sequences with different cdr3 length)
        n_proc_list = []
        start = time.time()
        while n_procs > 0:
            print '--> %d clusters with %d proc%s' % (len(cpath.partitions[cpath.i_best_minus_x]), n_procs, utils.plural(n_procs))
            cpath = self.run_hmm('forward', self.sub_param_dir, n_procs=n_procs, partition=cpath.partitions[cpath.i_best_minus_x], shuffle_input=True)  # NOTE that a.t.m. i_best and i_best_minus_x are usually the same, since we're usually not calculating log probs of partitions (well, we're trying to avoid calculating any extra log probs, which means we usually don't know the log prob of the entire partition)
            n_proc_list.append(n_procs)
            if self.are_we_finished_clustering(n_procs, cpath):
                break
            n_procs, cpath = self.prepare_next_iteration(n_proc_list, cpath, len(initial_nsets))

        print '      loop time: %.1f' % (time.time()-start)
        return cpath

    # ----------------------------------------------------------------------------------------
    def check_partition(self, partition):
        uids = set([uid for cluster in partition for uid in cluster])
        input_ids = set(self.sw_info['queries'])  # note that this does not include queries that were removed in sw
        missing_ids = input_ids - uids
        if self.unseeded_seqs is not None:
            missing_ids -= set(self.unseeded_seqs)
        if self.small_cluster_seqs is not None:
            missing_ids -= set(self.small_cluster_seqs)
        if len(missing_ids) > 0:
            warnstr = 'queries missing from partition: ' + str(len(missing_ids))
            print '  ' + utils.color('red', 'warning') + ' ' + warnstr

    # ----------------------------------------------------------------------------------------
    def auto_nprocs(self, nseqs):
        if self.args.n_precache_procs is not None:  # command line override
            return self.args.n_precache_procs

        seqs_per_proc = 500  # 2.5 mins (at something like 0.3 sec/seq)
        if nseqs > 3000:
            seqs_per_proc *= 2
        if nseqs > 10000:
            seqs_per_proc *= 1.5
        n_precache_procs = int(math.ceil(float(nseqs) / seqs_per_proc))
        n_precache_procs = min(n_precache_procs, self.args.n_max_procs)  # I can't get more'n a few hundred slots at a time, so it isn't worth using too much more than that
        if self.args.batch_system is None:  # if we're not on a batch system, make sure it's less than the number of cpus
            n_precache_procs = min(n_precache_procs, multiprocessing.cpu_count())

        return n_precache_procs

    # ----------------------------------------------------------------------------------------
    def write_clusterpaths(self, outfname, cpath):
        outfile, writer = cpath.init_outfile(outfname, self.args.is_data)
        true_partition = None
        if not self.args.is_data:
            true_partition = utils.get_true_partition(self.reco_info)
        cpath.write_partitions(writer=writer, reco_info=self.reco_info, true_partition=true_partition, is_data=self.args.is_data, n_to_write=self.args.n_partitions_to_write, calc_missing_values='best')
        outfile.close()

        if self.args.presto_output:
            outstr = check_output(['mv', '-v', self.args.outfname, self.args.outfname + '.partis'])
            print '    backing up partis output before converting to presto: %s' % outstr.strip()
            cpath.write_presto_partitions(self.args.outfname, self.input_info)

    # ----------------------------------------------------------------------------------------
    def get_cluster_annotations(self, partition):
        if len(partition) == 0:
            return
        action_cache = self.current_action
        self.current_action = 'run-viterbi'
        partition = sorted(partition, key=len, reverse=True)  # as opposed to in clusterpath, where we *don't* want to sort by size, it's nicer to have them sorted by size here, since then as you're scanning down a long list of cluster annotations you know once you get to the singletons you won't be missing something big
        n_procs = min(self.args.n_procs, len(partition))  # we want as many procs as possible, since the large clusters can take a long time (depending on if we're translating...), but in general we treat <self.args.n_procs> as the maximum allowable number of processes
        print '--> getting annotations for final partition'
        self.run_hmm('viterbi', self.sub_param_dir, n_procs=n_procs, partition=partition, read_output=False)  # it would be nice to rearrange <self.read_hmm_output()> so I could remove this option
        if n_procs > 1:
            _ = self.merge_all_hmm_outputs(n_procs, precache_all_naive_seqs=False)
        annotations = self.read_annotation_output(self.hmm_outfname, outfname=self.args.cluster_annotation_fname, print_annotations=self.args.print_cluster_annotations)
        if os.path.exists(self.hmm_infname):
            os.remove(self.hmm_infname)
        self.current_action = action_cache

        return annotations

    # ----------------------------------------------------------------------------------------
    def run_swarm(self, naive_seqs, threshold, outfname):
        raise Exception('needs updating/fixing')
        # # ----------------------------------------------------------------------------------------
        # if self.args.naive_swarm:
        #     print '    NOTE: replacing N with A for input to swarm'
        # with open(fastafname, 'w') as fastafile:
        #     for query, naive_seq in naive_seqs.items():
        #         if self.args.naive_swarm:
        #             query += '_1'
        #             naive_seq = utils.remove_ambiguous_ends(naive_seq)
        #             naive_seq = naive_seq.replace('N', 'A')
        #         fastafile.write('>' + query + '\n' + naive_seq + '\n')
        # # ----------------------------------------------------------------------------------------
        #     clusterfname = self.args.workdir + '/swarm-clusters.txt'
        #     cmd = './bin/swarm-2.1.1-linux-x86_64 ' + fastafname
        #     cmd += ' -t 5'  # five threads TODO set this more intelligently
        #     # cmd += ' -f'
        #     cmd += ' --match-reward ' + str(self.args.match_mismatch[0])
        #     cmd += ' --mismatch-penalty ' + str(self.args.match_mismatch[1])
        #     cmd += ' --gap-opening-penalty ' + str(self.args.gap_open_penalty)
        #     # cmd += ' --gap-extension-penalty'
        #     tmpstart = time.time()
        #     total = 0.
        #     for key in self.sw_info['queries']:
        #         seq = self.input_info[key]['seqs'][0]
        #         total += float(len(seq))
        #     mean_length = total / len(self.sw_info['queries'])
        #     raise Exception('update for new thresholds')
        #     bound = self.get_naive_hamming_threshold(parameter_dir, 'tight') /  2.  # yay for heuristics! (I did actually optimize this...)
        #     differences = int(round(mean_length * bound))
        #     print '        d = mean len * mut freq bound = %f * %f = %f --> %d' % (mean_length, bound, mean_length * bound, differences)
        #     print '      swarm average time: %.1f' % (time.time()-tmpstart)
        #     cmd += ' --differences ' + str(differences)
        #     cmd += ' --uclust-file ' + clusterfname
        #     check_call(cmd.split())
        # # ----------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------
    def run_vsearch(self, naive_seqs, threshold):
        # write input
        infname = self.args.workdir + '/naive-seqs.fasta'
        outfname = self.args.workdir + '/vsearch-clusters.txt'
        with open(infname, 'w') as fastafile:
            for query, naive_seq in naive_seqs.items():
                fastafile.write('>' + query + '\n' + naive_seq + '\n')

        # run
        id_fraction = 1. - threshold
        cmd = self.args.partis_dir + '/bin/vsearch-1.1.3-linux-x86_64 --threads ' + str(self.args.n_procs) + ' --uc ' + outfname + ' --cluster_fast ' + infname + ' --id ' + str(id_fraction) + ' --maxaccept 0 --maxreject 0'
        cmdfos = [{'cmd_str' : cmd, 'outfname' : outfname, 'workdir' : self.args.workdir, 'threads' : self.args.n_procs}, ]
        utils.run_cmds(cmdfos, batch_system=self.args.batch_system, batch_options=self.args.batch_options, batch_config_fname=self.args.batch_config_fname)

        # read output
        id_clusters = {}
        with open(outfname) as clusterfile:
            reader = csv.DictReader(clusterfile, fieldnames=['type', 'cluster_id', '3', '4', '5', '6', '7', 'crap', 'query', 'morecrap'], delimiter='\t')
            for line in reader:
                if line['type'] == 'C':  # batshit output format: some lines are a cluster, and some are a query sequence. Skip the cluster ones.
                    continue
                cluster_id = int(line['cluster_id'])
                if cluster_id not in id_clusters:
                    id_clusters[cluster_id] = []
                uid = line['query']
                if self.args.naive_swarm and uid[-2:] == '_1':  # remove (dummy) abundance information
                    uid = uid[:-2]
                id_clusters[cluster_id].append(uid)
        partition = id_clusters.values()

        os.remove(infname)
        os.remove(outfname)
        return partition

    # ----------------------------------------------------------------------------------------
    def cluster_with_naive_vsearch_or_swarm(self, parameter_dir=None, read_hmm_cachefile=True):
        start = time.time()

        naive_seq_list = []
        if read_hmm_cachefile:
            assert parameter_dir is not None
            threshold = self.get_naive_hamming_bounds(parameter_dir)[0]  # lo and hi are the same
            cached_naive_seqs = {}
            with open(self.hmm_cachefname) as cachefile:
                reader = csv.DictReader(cachefile)
                for line in reader:
                    unique_ids = line['unique_ids'].split(':')
                    if len(unique_ids) == 1:  # if it's a cache file left over from a previous partitioning, there'll be clusters in it, too
                        cached_naive_seqs[unique_ids[0]] = line['naive_seq']
            for uid in self.sw_info['queries']:
                if uid not in cached_naive_seqs:
                    raise Exception('naive sequence for %s not found in %s' % (uid, self.hmm_cachefname))
                naive_seq_list.append((uid, cached_naive_seqs[uid]))
        else:
            assert parameter_dir is None  # i.e. get mut freq from sw info
            threshold = self.get_naive_hamming_bounds(parameter_dir=None, overall_mute_freq=self.sw_info['mute-freqs']['all'])[0]  # lo and hi are the same
            naive_seq_list = [(q, self.sw_info[q]['naive_seq']) for q in self.sw_info['queries']]

        all_naive_seqs, naive_seq_hashes = utils.collapse_naive_seqs(naive_seq_list, self.sw_info)

        print '    using hfrac bound for vsearch %.3f' % threshold

        partition = []
        print '    running vsearch %d times (once for each cdr3 length class):' % len(all_naive_seqs),
        for cdr3_length, sub_naive_seqs in all_naive_seqs.items():
            sub_hash_partition = self.run_vsearch(sub_naive_seqs, threshold)
            sub_uid_partition = [[uid for hashstr in hashcluster for uid in naive_seq_hashes[hashstr]] for hashcluster in sub_hash_partition]
            partition += sub_uid_partition
            print '.',
            sys.stdout.flush()
        print ''

        ccfs = [None, None]
        if not self.args.is_data:  # it's ok to always calculate this since it's only ever for one partition
            queries_without_annotations = set(self.input_info) - set(self.sw_info['queries'])
            tmp_partition = copy.deepcopy(partition) + [[q, ] for q in queries_without_annotations]  # just add the missing ones as singletons
            self.check_partition(tmp_partition)
            true_partition = utils.get_true_partition(self.reco_info)
            ccfs = utils.new_ccfs_that_need_better_names(tmp_partition, true_partition, self.reco_info)
        cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
        cpath.add_partition(partition, logprob=0.0, n_procs=1, ccfs=ccfs)

        print '      vsearch time: %.1f' % (time.time()-start)
        return cpath

    # ----------------------------------------------------------------------------------------
    def get_naive_hamming_bounds(self, parameter_dir=None, overall_mute_freq=None):
        if self.cached_naive_hamming_bounds is not None:  # only run the stuff below once
            return self.cached_naive_hamming_bounds

        if parameter_dir is not None:
            assert overall_mute_freq is None
            mutehist = Hist(fname=parameter_dir + '/all-mean-mute-freqs.csv')
            mute_freq = mutehist.get_mean(ignore_overflows=True)  # TODO should I not ignore overflows here? or should I ignore them when I set waterer.info['mute-freqs']?
        else:
            assert overall_mute_freq is not None
            mute_freq = overall_mute_freq

        # just use a line based on two points (mute_freq, threshold)
        x1, x2 = 0.05, 0.2  # 0.5x, 3x (for 10 leaves)

        if self.args.naive_hamming:  # set lo and hi to the same thing, so we don't use log prob ratios, i.e. merge if less than this, don't merge if greater than this
            y1, y2 = 0.035, 0.06
            lo = utils.intexterpolate(x1, y1, x2, y2, mute_freq)
            hi = lo
        elif self.args.naive_vsearch:  # set lo and hi to the same thing, so we don't use log prob ratios, i.e. merge if less than this, don't merge if greater than this
            y1, y2 = 0.02, 0.05
            lo = utils.intexterpolate(x1, y1, x2, y2, mute_freq)
            hi = lo
        else:  # these are a bit larger than the tight ones and should almost never merge non-clonal sequences, i.e. they're appropriate for naive hamming preclustering if you're going to run the full likelihood on nearby sequences
            y1, y2 = 0.015, 0.015  # TODO get better numbers for this
            lo = utils.intexterpolate(x1, y1, x2, y2, mute_freq)  # ...and never merge 'em if it's bigger than this
            y1, y2 = 0.08, 0.15
            hi = utils.intexterpolate(x1, y1, x2, y2, mute_freq)  # ...and never merge 'em if it's bigger than this

        print '       naive hfrac bounds: %.3f %.3f   (%.3f mean mutation in parameter dir %s)' % (lo, hi, mute_freq, parameter_dir)
        self.cached_naive_hamming_bounds = (lo, hi)
        return self.cached_naive_hamming_bounds

    # ----------------------------------------------------------------------------------------
    def get_hmm_cmd_str(self, algorithm, csv_infname, csv_outfname, parameter_dir, precache_all_naive_seqs, n_procs):
        """ Return the appropriate bcrham command string """
        cmd_str = self.args.partis_dir + '/packages/ham/bcrham'
        cmd_str += ' --algorithm ' + algorithm
        if self.args.debug > 0:
            cmd_str += ' --debug ' + str(self.args.debug)
        cmd_str += ' --hmmdir ' + os.path.abspath(parameter_dir) + '/hmms'
        cmd_str += ' --datadir ' + self.my_gldir
        cmd_str += ' --infile ' + csv_infname
        cmd_str += ' --outfile ' + csv_outfname
        cmd_str += ' --locus ' + self.args.locus
        cmd_str += ' --random-seed ' + str(self.args.seed)
        if self.args.cache_naive_hfracs:
            cmd_str += ' --cache-naive-hfracs'
        if n_procs > 1:  # only cache vals for sequence sets with newly-calculated vals (initial cache file is copied to each subdir)
            cmd_str += ' --only-cache-new-vals'

        if self.args.dont_rescale_emissions:
            cmd_str += ' --dont-rescale-emissions'
        if self.current_action == 'partition':
            if os.path.exists(self.hmm_cachefname):
                cmd_str += ' --input-cachefname ' + self.hmm_cachefname
            cmd_str += ' --output-cachefname ' + self.hmm_cachefname
            if precache_all_naive_seqs:
                cmd_str += ' --cache-naive-seqs'
            else:  # actually partitioning
                cmd_str += ' --partition'
                cmd_str += ' --max-logprob-drop ' + str(self.args.max_logprob_drop)

                hfrac_bounds = self.get_naive_hamming_bounds(parameter_dir)
                if self.args.naive_hamming:  # shouldn't be able to happen, but...
                    assert hfrac_bounds[0] == hfrac_bounds[1]
                cmd_str += ' --hamming-fraction-bound-lo ' + str(hfrac_bounds[0])
                cmd_str += ' --hamming-fraction-bound-hi ' + str(hfrac_bounds[1])
                cmd_str += ' --logprob-ratio-threshold ' + str(self.args.logprob_ratio_threshold)
                cmd_str += ' --biggest-naive-seq-cluster-to-calculate ' + str(self.args.biggest_naive_seq_cluster_to_calculate)
                cmd_str += ' --biggest-logprob-cluster-to-calculate ' + str(self.args.biggest_logprob_cluster_to_calculate)
                cmd_str += ' --n-partitions-to-write ' + str(self.args.n_partitions_to_write)  # don't write too many, since calculating the extra logprobs is kind of expensive
                if n_procs == 1:  # if this is the last time through, with one process, we want glomerator.cc to calculate the total logprob of each partition NOTE this is quite expensive, since we have to turn off translation entirely
                    cmd_str += '  --write-logprob-for-each-partition'

                if self.args.seed_unique_id is not None and self.unseeded_seqs is None:  # if we're in the last few cycles (i.e. we've removed unseeded clusters) we want bcrham to not know about the seed (this gives more accurate clustering 'cause we're really doing hierarchical agglomeration)
                    cmd_str += ' --seed-unique-id ' + self.args.seed_unique_id

                if n_procs == 1 and self.args.n_final_clusters is not None:
                    cmd_str += ' --n-final-clusters ' + str(self.args.n_final_clusters)

        assert len(utils.ambiguous_bases) == 1  # could allow more than one, but it's not implemented a.t.m.
        cmd_str += ' --ambig-base ' + utils.ambiguous_bases[0]

        return cmd_str

    # ----------------------------------------------------------------------------------------
    def subworkdir(self, iproc, n_procs):
        if n_procs == 1:
            return self.args.workdir
        else:
            return self.args.workdir + '/hmm-' + str(iproc)

    # ----------------------------------------------------------------------------------------
    def print_partition_dbgfo(self):
        if self.bcrham_proc_info is None:
            return
        summaryfo = utils.summarize_bcrham_dbgstrs(self.bcrham_proc_info)

        pwidth = str(len(str(len(self.input_info))))  # close enough
        if sum(summaryfo['read-cache'].values()) == 0:
            print '                no/empty cache file'
        else:
            print ('          read from cache:  naive-seqs %' + pwidth + 'd   logprobs %' + pwidth + 'd') % (summaryfo['read-cache']['naive-seqs'], summaryfo['read-cache']['logprobs'])
        print ('                    calcd:         vtb %' + pwidth + 'd        fwd %' + pwidth + 'd') % (summaryfo['calcd']['vtb'], summaryfo['calcd']['fwd'])
        print ('                   merged:       hfrac %' + pwidth + 'd     lratio %' + pwidth + 'd') % (summaryfo['merged']['hfrac'], summaryfo['merged']['lratio'])

        if len(self.bcrham_proc_info) == 1:
            print '                     time:  %.1f sec' % summaryfo['time']['bcrham'][0]
        else:
            print '             min-max time:  %.1f - %.1f sec' % (summaryfo['time']['bcrham'][0], summaryfo['time']['bcrham'][1])

    # ----------------------------------------------------------------------------------------
    def check_wait_times(self, wait_time):
        max_bcrham_time = max([procinfo['time']['bcrham'] for procinfo in self.bcrham_proc_info])
        if max_bcrham_time > 0. and wait_time / max_bcrham_time > 1.5 and wait_time > 30.:  # if we were waiting for a lot longer than the slowest process took, and if it took long enough for us to care
            print '    spent much longer waiting for bcrham (%.1fs) than bcrham reported taking (max per-proc time %.1fs)' % (wait_time, max_bcrham_time)

    # ----------------------------------------------------------------------------------------
    def execute(self, cmd_str, n_procs):
        # ----------------------------------------------------------------------------------------
        def get_outfname(iproc):
            return self.hmm_outfname.replace(self.args.workdir, self.subworkdir(iproc, n_procs))
        # ----------------------------------------------------------------------------------------
        def get_cmd_str(iproc):
            strlist = cmd_str.split()
            for istr in range(len(strlist)):
                if strlist[istr] == self.hmm_infname or strlist[istr] == self.hmm_cachefname or strlist[istr] == self.hmm_outfname:
                    strlist[istr] = strlist[istr].replace(self.args.workdir, self.subworkdir(iproc, n_procs))
            return ' '.join(strlist)

        print '    running %d proc%s' % (n_procs, utils.plural(n_procs))
        sys.stdout.flush()
        start = time.time()

        self.bcrham_proc_info = [{} for _ in range(n_procs)]
        cmdfos = [{'cmd_str' : get_cmd_str(iproc),
                   'workdir' : self.subworkdir(iproc, n_procs),
                   'outfname' : get_outfname(iproc),
                   'dbgfo' : self.bcrham_proc_info[iproc]}
                  for iproc in range(n_procs)]
        utils.run_cmds(cmdfos, batch_system=self.args.batch_system, batch_options=self.args.batch_options, batch_config_fname=self.args.batch_config_fname, debug='print' if self.args.debug else None)
        if self.current_action == 'partition':
            self.print_partition_dbgfo()

        self.check_wait_times(time.time()-start)
        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, algorithm, parameter_in_dir, parameter_out_dir='', count_parameters=False, n_procs=None, precache_all_naive_seqs=False, partition=None, shuffle_input=False, read_output=True):
        """ 
        Run bcrham, possibly with many processes, and parse and interpret the output.
        NOTE the local <n_procs>, which overrides the one from <self.args>
        """
        start = time.time()
        if len(self.sw_info['queries']) == 0:
            print '  %s no input queries for hmm' % utils.color('red', 'warning')
            return

        if n_procs is None:
            n_procs = self.args.n_procs

        self.prepare_for_hmm(algorithm, parameter_in_dir, partition, shuffle_input=shuffle_input)

        cmd_str = self.get_hmm_cmd_str(algorithm, self.hmm_infname, self.hmm_outfname, parameter_dir=parameter_in_dir, precache_all_naive_seqs=precache_all_naive_seqs, n_procs=n_procs)

        if n_procs > 1:
            self.split_input(n_procs, self.hmm_infname)

        exec_start = time.time()
        self.execute(cmd_str, n_procs)
        exec_time = time.time() - exec_start

        new_cpath = None
        if read_output:
            new_cpath = self.read_hmm_output(algorithm, n_procs, count_parameters, parameter_out_dir, precache_all_naive_seqs)

        step_time = time.time() - start
        if step_time - exec_time > 0.1:
            print '         infra time: %.1f' % (step_time - exec_time)  # i.e. time for non-executing, infrastructure time
        print '      hmm step time: %.1f' % step_time
        self.timing_info.append({'exec' : exec_time, 'total' : step_time})  # NOTE in general, includes pre-cache step

        return new_cpath

    # ----------------------------------------------------------------------------------------
    def read_hmm_cachefile(self):
        cachefo = {}
        if not os.path.exists(self.hmm_cachefname):
            return cachefo
        with open(self.hmm_cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                utils.process_input_line(line, hmm_cachefile=True)
                cachefo[':'.join(line['unique_ids'])] = line
        return cachefo

    # ----------------------------------------------------------------------------------------
    def print_subcluster_naive_seqs(self, uids_of_interest):
        uids_of_interest = set(uids_of_interest)
        uidstr_of_interest, ref_naive_seq = None, None  # we don't know what order they're in the cache file yet

        cachefo = self.read_hmm_cachefile()
        sub_uidstrs = []
        for uidstr, info in cachefo.items():
            if info['naive_seq'] == '':
                continue
            uids = set(info['unique_ids'])
            if len(uids & uids_of_interest) > 0:  # first see if there's some overlap with what we're interested in
                sub_uidstrs.append(uidstr)
            if len(uids - uids_of_interest) == 0 and len(uids_of_interest - uids) == 0:  # then see if it's the actual cluster we're interested in
                uidstr_of_interest = uidstr

        sub_uidstrs = sorted(sub_uidstrs, key=lambda x: x.count(':'))
        if uidstr_of_interest is None:
            uidstr_of_interest = sub_uidstrs[-1]
            uids_of_interest = None  # make sure we don't use it later
            print '%s couldn\'t find the exact requested cluster, so using the biggest cluster that has some overlap with the requested cluster for the reference sequence' % utils.color('yellow', 'warning')
        ref_naive_seq = cachefo[uidstr_of_interest]['naive_seq']
        print '  subcluster naive sequences for %s (in %s below)' % (uidstr_of_interest, utils.color('blue', 'blue'))
        # print '      %s  %s' % (ref_naive_seq, utils.color('blue', uidstr_of_interest))
        for uidstr in sub_uidstrs:
            post_str = uidstr
            if uidstr == uidstr_of_interest:
                post_str = utils.color('blue', post_str)
            print utils.color_mutants(ref_naive_seq, cachefo[uidstr]['naive_seq'], extra_str='      ', print_isnps=True, post_str='  ' + post_str)

    # ----------------------------------------------------------------------------------------
    def get_padded_true_naive_seq(self, qry):
        assert len(self.sw_info[qry]['padlefts']) == 1
        return self.sw_info[qry]['padlefts'][0] * utils.ambiguous_bases[0] + self.reco_info[qry]['naive_seq'] + self.sw_info[qry]['padrights'][0] * utils.ambiguous_bases[0]

    # ----------------------------------------------------------------------------------------
    def get_sw_naive_seqs(self, info, namekey):

        naive_seqs = {}
        for line in info:
            query = line[namekey]
            if len(query.split(':')) == 1:  # ...but if we don't have them, use smith-waterman (should only be for single queries)
               naive_seqs[query] = self.sw_info[query]['naive_seq']
            elif len(query.split(':')) > 1:
                naive_seqs[query] = self.sw_info[query.split(':')[0]]['naive_seq']  # just arbitrarily use the naive seq from the first one. This is ok partly because if we cache the logprob but not the naive seq, that's because we thought about merging two clusters but did not -- so they're naive seqs should be similar. Also, this is just for divvying queries.
            else:
                raise Exception('no naive sequence found for ' + str(query))
            if naive_seqs[query] == '':
                raise Exception('zero-length naive sequence found for ' + str(query))
        return naive_seqs

    # ----------------------------------------------------------------------------------------
    def split_input(self, n_procs, infname):

        # should we pull out the seeded clusters, and carefully re-inject them into each process?
        separate_seeded_clusters = self.args.seed_unique_id is not None and self.unseeded_seqs is None

        # read single input file
        info = []
        seeded_clusters = {}
        with open(infname, 'r') as infile:
            reader = csv.DictReader(infile, delimiter=' ')
            for line in reader:
                if separate_seeded_clusters and self.args.seed_unique_id in set(line['names'].split(':')):
                    if len(seeded_clusters) > 0 and ':' not in line['names']:  # the first time through, we add the seed uid to *every* process. So, when we read those results back in, the procs that didn't merge the seed with anybody will have it as a singleton still, and we only need the singleton once
                        continue
                    seeded_clusters[line['names']] = line
                    continue  # don't want the seeded clusters mixed in with the non-seeded clusters just yet (see below)
                info.append(line)

        # find the smallest seeded cluster
        if separate_seeded_clusters:
            if len(seeded_clusters) == 0:
                raise Exception('couldn\'t find info for seed query %s in %s' % (self.args.seed_unique_id, infname))
            smallest_seed_cluster_str = None
            for unique_id_str in seeded_clusters:
                if smallest_seed_cluster_str is None or len(unique_id_str.split(':')) < len(smallest_seed_cluster_str.split(':')):
                    smallest_seed_cluster_str = unique_id_str

        # ----------------------------------------------------------------------------------------
        def get_sub_outfile(siproc, mode):
            return open(self.subworkdir(siproc, n_procs) + '/' + os.path.basename(infname), mode)
        def get_writer(sub_outfile):
            return csv.DictWriter(sub_outfile, reader.fieldnames, delimiter=' ')
        def copy_cache_files(n_procs):
            cmdfos = [{'cmd_str' : 'cp ' + self.hmm_cachefname + ' ' + self.subworkdir(iproc, n_procs) + '/',
                       'workdir' : self.subworkdir(iproc, n_procs),
                       'outfname' : self.subworkdir(iproc, n_procs) + '/' + os.path.basename(self.hmm_cachefname)}
                      for iproc in range(n_procs)]
            utils.run_cmds(cmdfos)

        # initialize output/cache files
        for iproc in range(n_procs):
            utils.prep_dir(self.subworkdir(iproc, n_procs))
            sub_outfile = get_sub_outfile(iproc, 'w')
            get_writer(sub_outfile).writeheader()
            sub_outfile.close()  # can't leave 'em all open the whole time 'cause python has the thoroughly unreasonable idea that one oughtn't to have thousands of files open at once
        if self.current_action == 'partition' and os.path.exists(self.hmm_cachefname):  # copy cachefile to this subdir (first clause is just for so when we're getting cluster annotations we don't copy over the cache files)
            copy_cache_files(n_procs)

        seed_clusters_to_write = seeded_clusters.keys()  # the keys in <seeded_clusters> that we still need to write
        for iproc in range(n_procs):
            sub_outfile = get_sub_outfile(iproc, 'a')
            writer = get_writer(sub_outfile)

            # first deal with the seeded clusters
            if separate_seeded_clusters:  # write the seed info line to each file
                if len(seed_clusters_to_write) > 0:
                    if iproc < n_procs - 1:  # if we're not on the last proc, pop off and write the first one
                        writer.writerow(seeded_clusters[seed_clusters_to_write.pop(0)])
                    else:
                        while len(seed_clusters_to_write) > 0:  # keep adding 'em until we run out
                            writer.writerow(seeded_clusters[seed_clusters_to_write.pop(0)])
                else:  # if we don't have any more that we *need* to write (i.e. that have other seqs in them), just write the shortest one (which will frequently be a singleton)
                    writer.writerow(seeded_clusters[smallest_seed_cluster_str])

            # then loop over the non-seeded clusters
            for iquery in range(len(info)):
                if iquery % n_procs != iproc:
                    continue
                writer.writerow(info[iquery])
            sub_outfile.close()

    # ----------------------------------------------------------------------------------------
    def merge_subprocess_files(self, fname, n_procs, include_outfile=False):
        subfnames = []
        for iproc in range(n_procs):
            subfnames.append(self.subworkdir(iproc, n_procs) + '/' + os.path.basename(fname))
        if include_outfile:  # also merge the output file <fname> (i.e. for the cache files, the sub files only include *new* information, so we need to also merge them with the original file)
            subfnames.append(fname)
        self.merge_files(subfnames, fname, dereplicate=False)

    # ----------------------------------------------------------------------------------------
    def merge_files(self, infnames, outfname, dereplicate):
        """ 
        Merge <infnames> into <outfname>.
        NOTE that <outfname> is overwritten with the zero-length file if it exists, otherwise it is created.
        Some of <infnames> may not exist.
        """
        # check_call(['wc', ] + [fn for fn in infnames if fn != outfname])
        # if os.path.exists(outfname):
        #     check_call(['wc', outfname])
        # else:
        #     print '  outfname d.n.e.'

        header = ''
        outfile = None
        one_real_file = False
        if outfname not in infnames or not os.path.exists(outfname):  # if it *is* in <infnames> we assume we can just tack the other infnames onto the end of it and use <outfname>'s header
            outfile = open(outfname, 'w')
        for fname in infnames:
            if not os.path.exists(fname) or os.stat(fname).st_size == 0:
                continue
            one_real_file = True
            with open(fname) as headfile:
                reader = csv.DictReader(headfile)
                header = ','.join(reader.fieldnames)
                if outfile is not None:
                    writer = csv.DictWriter(outfile, reader.fieldnames)
                    writer.writeheader()
            break  # kinda weird to do it this way, but we just need one of the infiles to get the header info (and some may be zero length)
        if outfile is not None:
            outfile.close()
        if not one_real_file:
            print '    nothing to merge into %s' % outfname
            return

        assert header != ''

        non_out_infnames = [fn for fn in infnames if fn != outfname]
        if len(non_out_infnames) == 0:
            raise Exception('merge_files() called with <infnames> consisting only of <outfname>')
        cmd = 'cat ' + ' '.join(non_out_infnames) + ' | grep -v \'' + header + '\''
        cmd += ' >>' + outfname
        try:
            check_call(cmd, shell=True)
        except CalledProcessError:
            print '    nothing to merge into %s' % outfname
            # raise Exception('only read headers from %s', ' '.join([fn for fn in infnames if fn != outfname]))

        if dereplicate:
            tmpfname = outfname + '.tmp'
            check_call('echo ' + header + ' >' + tmpfname, shell=True)
            check_call('grep -v \'' + header + '\' ' + outfname + ' | sort | uniq >>' + tmpfname, shell=True)  # NOTE there can be multiple lines with the same uid string, but this is ok -- the c++ handles it
            check_call(['mv', tmpfname, outfname])

        for infname in infnames:
            if infname != outfname:
                os.remove(infname)

    # ----------------------------------------------------------------------------------------
    def merge_all_hmm_outputs(self, n_procs, precache_all_naive_seqs):
        """ Merge any/all output files from subsidiary bcrham processes """
        cpath = None  # it would be nice to figure out a cleaner way to do this
        if self.current_action == 'partition':  # merge partitions from several files
            if n_procs > 1:
                self.merge_subprocess_files(self.hmm_cachefname, n_procs, include_outfile=True)  # sub cache files only have new info

            if not precache_all_naive_seqs:
                if n_procs == 1:
                    infnames = [self.hmm_outfname, ]
                else:
                    infnames = [self.subworkdir(iproc, n_procs) + '/' + os.path.basename(self.hmm_outfname) for iproc in range(n_procs)]
                glomerer = Glomerator(self.reco_info, seed_unique_id=self.args.seed_unique_id)
                glomerer.read_cached_agglomeration(infnames, debug=self.args.debug)  #, outfname=self.hmm_outfname)
                assert len(glomerer.paths) == 1
                cpath = glomerer.paths[0]
        else:
            self.merge_subprocess_files(self.hmm_outfname, n_procs)

        if n_procs == 1:
            os.remove(self.hmm_outfname)
        else:
            for iproc in range(n_procs):
                subworkdir = self.subworkdir(iproc, n_procs)
                os.remove(subworkdir + '/' + os.path.basename(self.hmm_infname))
                if os.path.exists(subworkdir + '/' + os.path.basename(self.hmm_outfname)):
                    os.remove(subworkdir + '/' + os.path.basename(self.hmm_outfname))
                os.rmdir(subworkdir)

        return cpath

    # ----------------------------------------------------------------------------------------
    def write_hmms(self, parameter_dir):
        """ Write hmm model files to <parameter_dir>/hmms, using information from <parameter_dir> """
        print '  writing hmms',
        sys.stdout.flush()
        start = time.time()

        from hmmwriter import HmmWriter
        hmm_dir = parameter_dir + '/hmms'
        utils.prep_dir(hmm_dir, '*.yaml')

        if self.args.debug:
            print '    to %s' % parameter_dir + '/hmms'

        for region in utils.regions:
            for gene in self.glfo['seqs'][region]:
                writer = HmmWriter(parameter_dir, hmm_dir, gene, self.glfo, self.args)
                writer.write()

        print '(%.1f sec)' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def get_existing_hmm_files(self, parameter_dir):
        fnames = [os.path.basename(fn) for fn in glob.glob(parameter_dir + '/hmms/*.yaml')]
        genes = set([utils.unsanitize_name(os.path.splitext(fn)[0]) for fn in fnames])
        if len(genes) == 0:
            raise Exception('no yamels in %s' % parameter_dir + '/hmms')
        return genes

    # ----------------------------------------------------------------------------------------
    def all_regions_present(self, gene_list, skipped_gene_matches, query_name, second_query_name=None):
        """ Check that we have at least one gene for each region """
        for region in utils.regions:
            found = False
            for gene in gene_list:
                if utils.get_region(gene) == region:
                    found = True
                    break
            if not found:
                print '       no %s genes in %s for %s %s' % (region, ':'.join(gene_list), query_name, '' if (second_query_name == None) else second_query_name)
                print '          skipped %s' % (':'.join(skipped_gene_matches))
                print 'giving up on query'
                return False

        return True

    # ----------------------------------------------------------------------------------------
    def combine_queries(self, query_names, genes_with_hmm_files, skipped_gene_matches=None):
        """ 
        Return the 'logical OR' of the queries in <query_names>, i.e. the maximal extent in k_v/k_d space and OR of only_gene sets.
        """

        # Note that this whole thing probably ought to use cached hmm info if it's available.
        # Also, this just always uses the SW mutation rate, but I should really update it with the (multi-)hmm-derived ones (same goes for k space boundaries)

        combo = {}
        combo['seqs'] = [self.sw_info[name]['seqs'][0] for name in query_names]
        combo['mut_freq'] = numpy.mean([utils.hamming_fraction(self.sw_info[name]['naive_seq'], self.sw_info[name]['seqs'][0]) for name in query_names])
        cdr3_lengths = [self.sw_info[name]['cdr3_length'] for name in query_names]
        if cdr3_lengths.count(cdr3_lengths[0]) != len(cdr3_lengths):
            uids_and_lengths = {q : self.sw_info[q]['cdr3_length'] for q in query_names}
            uids_and_lengths = sorted(uids_and_lengths.items(), key=operator.itemgetter(1))
            uids, lengths = zip(*uids_and_lengths)
            raise Exception('%s cdr3 lengths not all the same for %s (%s)' % (utils.color('red', 'warning'), ' '.join(uids), ' '.join([str(c) for c in lengths])))
        combo['cdr3_length'] = cdr3_lengths[0]

        combo['k_v'] = {'min' : 99999, 'max' : -1}
        combo['k_d'] = {'min' : 99999, 'max' : -1}
        combo['only_genes'] = []
        if self.args.linearham:
            combo['relpos'] = {}
        for name in query_names:
            swfo = self.sw_info[name]
            k_v = swfo['k_v']
            k_d = swfo['k_d']
            assert len(swfo['seqs']) == 1  # checking that when we filled in 'seqs' all was well
            combo['k_v']['min'] = min(k_v['min'], combo['k_v']['min'])
            combo['k_v']['max'] = max(k_v['max'], combo['k_v']['max'])
            combo['k_d']['min'] = min(k_d['min'], combo['k_d']['min'])
            combo['k_d']['max'] = max(k_d['max'], combo['k_d']['max'])
            if self.args.linearham:
                for gene in set(swfo['relpos']) - set(combo['relpos']):  # loop over genes that haven't come up yet in previous queries (note that this takes <pos> from the first <name> that happens to have a match to <gene>)
                    combo['relpos'][gene] = swfo['relpos'][gene]

            genes_to_use = set()  # genes from this query that'll get ORd into the ones from the previous queries
            for region in utils.regions:
                regmatches = swfo['all_matches'][region]  # the best <n_max_per_region> matches for this query, ordered by sw match quality
                skipped_gene_matches |= set([g for g in regmatches if g not in genes_with_hmm_files])
                genes_to_use |= set([g for g in regmatches if g in genes_with_hmm_files])

            # OR this query's genes into the ones from previous queries
            combo['only_genes'] = list(set(genes_to_use) | set(combo['only_genes']))  # NOTE using the OR of all sets of genes (from all query seqs) like this *really* helps,

        if self.args.linearham:
            boundfcns = {'l' : min, 'r' : max}  # take the widest set of flexbounds over all the queries
            combo['flexbounds'] = {region : {side : boundfcns[side]([self.sw_info[q]['flexbounds'][region][side] for q in query_names])
                                             for side in ['l', 'r']}
                                   for region in utils.regions}

        if not self.all_regions_present(combo['only_genes'], skipped_gene_matches, query_names):
            return {}

        for kb in ['k_v', 'k_d']:
            if combo[kb]['min'] <= 0 or combo[kb]['min'] >= combo[kb]['max']:
                raise Exception('nonsense k bounds for %s (v: %d %d  d: %d %d)' % (':'.join(query_names), combo['k_v']['min'], combo['k_v']['max'], combo['k_d']['min'], combo['k_d']['max']))

        return combo

    # ----------------------------------------------------------------------------------------
    def write_fake_cache_file(self, nsets):
        """ Write a fake cache file which, instead of the inferred naive sequences, has the *true* naive sequences. Used to generate synthetic partitions. """

        if self.reco_info is None:
            raise Exception('can\'t write fake cache file for --synthetic-distance-based-partition unless --is-simu is specified (and there\'s sim info in the input csv)')

        if os.path.exists(self.hmm_cachefname):
            print '      cache file exists, not writing fake true naive seqs'
            return

        print '      caching fake true naive seqs'
        with open(self.hmm_cachefname, 'w') as fakecachefile:
            writer = csv.DictWriter(fakecachefile, utils.partition_cachefile_headers)
            writer.writeheader()
            for query_name_list in nsets:
                writer.writerow({
                    'unique_ids' : ':'.join([qn for qn in query_name_list]),
                    'naive_seq' : self.get_padded_true_naive_seq(query_name_list[0])  # NOTE just using the first one... but a.t.m. I think I'll only run this fcn the first time through when they're all singletons, anyway
                })

    # ----------------------------------------------------------------------------------------
    def write_to_single_input_file(self, fname, nsets, parameter_dir, skipped_gene_matches, shuffle_input=False):
        csvfile = open(fname, 'w')
        header = ['names', 'k_v_min', 'k_v_max', 'k_d_min', 'k_d_max', 'mut_freq', 'cdr3_length', 'only_genes', 'seqs']
        if self.args.linearham:
            header += ['flexbounds', 'relpos']
        writer = csv.DictWriter(csvfile, header, delimiter=' ')
        writer.writeheader()

        if shuffle_input:  # shuffle nset order (this is absolutely critical when clustering with more than one process, in order to redistribute sequences among the several processes)
            random.shuffle(nsets)

        if self.args.synthetic_distance_based_partition:
            self.write_fake_cache_file(nsets)

        genes_with_hmm_files = self.get_existing_hmm_files(parameter_dir)

        glfo_genes = set([g for r in utils.regions for g in self.glfo['seqs'][r]])
        if self.args.only_genes is None and len(genes_with_hmm_files - glfo_genes) > 0:
            print '  %s hmm files for %s that aren\'t in glfo' % (utils.color('red', 'warning'), ' '.join(genes_with_hmm_files - glfo_genes))
        if len(glfo_genes - genes_with_hmm_files) > 0:
            print '  %s no hmm files for glfo genes %s' % (utils.color('red', 'warning'), ' '.join(glfo_genes - genes_with_hmm_files))

        for query_name_list in nsets:  # NOTE in principle I think I should remove duplicate singleton <seed_unique_id>s here. But I think they in effect get removed 'cause in bcrham everything's stored as hash maps, so any duplicates just overwites the original upon reading its input
            combined_query = self.combine_queries(query_name_list, genes_with_hmm_files, skipped_gene_matches=skipped_gene_matches)
            if len(combined_query) == 0:  # didn't find all regions
                continue
            row = {'names' : ':'.join([qn for qn in query_name_list]),
                   'k_v_min' : combined_query['k_v']['min'],
                   'k_v_max' : combined_query['k_v']['max'],
                   'k_d_min' : combined_query['k_d']['min'],
                   'k_d_max' : combined_query['k_d']['max'],
                   'mut_freq' : combined_query['mut_freq'],
                   'cdr3_length' : combined_query['cdr3_length'],
                   'only_genes' : ':'.join(combined_query['only_genes']),
                   'seqs' : ':'.join(combined_query['seqs'])}
            if self.args.linearham:
                row['flexbounds'] = json.dumps(combined_query['flexbounds'], separators=(',', ':'))
                row['relpos'] = json.dumps(combined_query['relpos'], separators=(',', ':'))
            writer.writerow(row)

        csvfile.close()

    # ----------------------------------------------------------------------------------------
    def prepare_for_hmm(self, algorithm, parameter_dir, partition, shuffle_input=False):
        """ Write input file for bcrham """
        if self.args.debug:
            print '    writing input'

        if partition is not None:
            nsets = copy.deepcopy(partition)  # needs to be a deep copy so we can shuffle the order
        else:
            qlist = [q for q in self.input_info if q in self.sw_info['queries']]  # we want the queries from sw (to skip failures), but the order from input_info
            if self.args.simultaneous_true_clonal_seqs:
                assert self.args.n_simultaneous_seqs is None and not self.args.is_data  # are both already checked in ./bin/partis
                nsets = utils.get_true_partition(self.reco_info, ids=qlist)
            elif self.args.n_simultaneous_seqs is None:  # plain ol' singletons
                nsets = [[q] for q in qlist]
            else:
                nlen = self.args.n_simultaneous_seqs  # shorthand
                # nsets = [qlist[iq : min(iq + nlen, len(qlist))] for iq in range(0, len(qlist), nlen)]  # this way works fine, but it's hard to get right 'cause it's hard to understand 
                nsets = [list(group) for _, group in itertools.groupby(qlist, key=lambda q: qlist.index(q) / nlen)]  # integer division

        skipped_gene_matches = set()
        self.write_to_single_input_file(self.hmm_infname, nsets, parameter_dir, skipped_gene_matches, shuffle_input=shuffle_input)  # single file gets split up later if we've got more than one process
        if self.args.debug and len(skipped_gene_matches) > 0:
            print '    not found in %s, so removing from consideration for hmm (i.e. were only the nth best, but never the best sw match for any query):' % (parameter_dir),
            for region in utils.regions:
                # print '  %s: %d' % (region, len([gene for gene in skipped_gene_matches if utils.get_region(gene) == region])),
                print '\n      %s: %s' % (region, ' '.join([utils.color_gene(gene) for gene in sorted(skipped_gene_matches) if utils.get_region(gene) == region]))
            print ''

    # ----------------------------------------------------------------------------------------
    def read_hmm_output(self, algorithm, n_procs, count_parameters, parameter_out_dir, precache_all_naive_seqs):
        cpath = None  # would be nice to figure out a cleaner way to do this
        if self.current_action == 'partition' or n_procs > 1:
            cpath = self.merge_all_hmm_outputs(n_procs, precache_all_naive_seqs)

        if self.current_action != 'partition' or count_parameters:
            if algorithm == 'viterbi':
                self.read_annotation_output(self.hmm_outfname, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir, outfname=self.args.outfname)
            elif algorithm == 'forward':
                self.read_forward_output(self.hmm_outfname)

        if os.path.exists(self.hmm_infname):
            os.remove(self.hmm_infname)

        return cpath

    # ----------------------------------------------------------------------------------------
    def check_did_bcrham_fail(self, line, errorfo):
        if line['errors'] == '':  # no problems
            return False

        failed = False
        for ecode in line['errors'].split(':'):  # I don't think I've ever seen one with more than one, but it could happen
            if ecode not in errorfo:
                errorfo[ecode] = set()
            for uid in line['unique_ids']:
                errorfo[ecode].add(uid)
            if ecode == 'no_path':  # boundary errors aren't failures, they're just telling us we had to expand the boundaries EDIT oh, wait, or does it mean it couldn't expand them enough? in any case, we still get an answer
                failed = True

        return failed

    # ----------------------------------------------------------------------------------------
    def read_forward_output(self, annotation_fname):
        probs = OrderedDict()
        with open(annotation_fname, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for line in reader:
                if line['errors'] != '':
                    print '  bcrham errors (%s) for %s' % (line['errors'], line['unique_ids'])
                probs[line['unique_ids']] = float(line['logprob'])

        if self.args.outfname is not None:
            with open(self.args.outfname, 'w') as outfile:
                writer = csv.DictWriter(outfile, ('unique_ids', 'logprob'))
                writer.writeheader()
                for uids, prob in probs.items():
                    writer.writerow({'unique_ids' : uids, 'logprob' : prob})

        os.remove(annotation_fname)

    # ----------------------------------------------------------------------------------------
    def process_dummy_d_hack(self, line, debug=False):
        """
        a.t.m. we force bcrham to give us D of length one for loci with no D genes.
        Here, we delete the dummy D base, and give it to either V, J, or the insertion.
        """
        tmpline = copy.deepcopy(line)
        utils.add_implicit_info(self.glfo, tmpline)
        if debug:
            print ''
            print '  dummy d hack for %s' % ' '.join(line['unique_ids'])
            utils.print_reco_event(tmpline, extra_str='    ', label='before')

        gl_v_base = None
        if tmpline['v_3p_del'] > 0:
            full_v_gl_seq = self.glfo['seqs']['v'][tmpline['v_gene']]
            gl_v_base = full_v_gl_seq[tmpline['lengths']['v']]

        gl_j_base = None
        if tmpline['j_5p_del'] > 0 and len(tmpline['dj_insertion']) == 0:
            full_j_gl_seq = self.glfo['seqs']['j'][tmpline['j_gene']]
            gl_j_base = full_j_gl_seq[tmpline['j_5p_del'] - 1]
        if debug:
            print '    gl_j_base', gl_j_base
            print '    gl_v_base', gl_v_base

        # take a majority vote as to whom we should give the base
        votes = {'v' : 0, 'j' : 0, 'dj_insertion' : 0}
        qr_base_votes = {n : 0 for n in utils.expected_characters}
        for iseq in range(len(tmpline['seqs'])):
            d_qr_base = tmpline['d_qr_seqs'][iseq]
            qr_base_votes[d_qr_base] += 1
            if d_qr_base == gl_v_base:
                votes['v'] += 1
            if d_qr_base == gl_j_base:  # yeah, that's right, you can vote twice
                votes['j'] += 1
            if d_qr_base != gl_v_base and d_qr_base != gl_j_base:
                votes['dj_insertion'] += 1

        sorted_votes = sorted(votes.items(), key=operator.itemgetter(1), reverse=True)
        winner = sorted_votes[0][0]
        sorted_qr_base_votes = sorted(qr_base_votes.items(), key=operator.itemgetter(1), reverse=True)
        qr_base_winner = sorted_qr_base_votes[0][0]
        if debug:
            print '   ', sorted_votes
            print '   ', sorted_qr_base_votes
            print '    winner', winner, qr_base_winner

        line['d_5p_del'] = 1  # NOTE we don't modify tmpline, i.e. we modify the line *without* implicit info, 'cause it's simpler
        if winner == 'v':
            assert line['v_3p_del'] > 0
            line['v_3p_del'] -= 1
        elif winner == 'j':
            line['j_5p_del'] -= 1
        else:
            assert winner == 'dj_insertion'
            line['dj_insertion'] = qr_base_winner + line['dj_insertion']

        after_line = copy.deepcopy(line)
        utils.add_implicit_info(self.glfo, after_line)
        if debug:
            utils.print_reco_event(after_line, extra_str='    ', label='after')

    # ----------------------------------------------------------------------------------------
    def check_for_unexpectedly_missing_keys(self, annotations, hmm_failure_ids):
        missing_input_keys = set(self.input_info)
        missing_input_keys -= set([uid for line in  annotations.values() for uid in line['unique_ids']])  # set(self.sw_info['queries'])  # all the queries for which we had decent sw annotations (sw failures are accounted for below)
        missing_input_keys -= self.sw_info['failed-queries']
        missing_input_keys -= set([d for dlist in self.sw_info['duplicates'].values() for d in dlist])
        missing_input_keys -= hmm_failure_ids
        if self.unseeded_seqs is not None:
            missing_input_keys -= set(self.unseeded_seqs)
        if self.small_cluster_seqs is not None:
            missing_input_keys -= set(self.small_cluster_seqs)
        if len(missing_input_keys) > 0:
            print '  %s couldn\'t account for %d missing input uid%s%s' % (utils.color('red', 'warning'), len(missing_input_keys), utils.plural(len(missing_input_keys)), ': %s' % ' '.join(missing_input_keys) if len(missing_input_keys) < 15 else '')

    # ----------------------------------------------------------------------------------------
    def read_annotation_output(self, annotation_fname, outfname=None, count_parameters=False, parameter_out_dir=None, print_annotations=False):
        """ Read bcrham annotation output """
        print '    read output'
        sys.stdout.flush()

        pcounter = ParameterCounter(self.glfo, self.args) if count_parameters else None
        true_pcounter = ParameterCounter(self.simglfo, self.args) if (count_parameters and not self.args.is_data) else None
        perfplotter = PerformancePlotter('hmm') if self.args.plot_performance else None

        n_lines_read, n_seqs_processed, n_events_processed, n_invalid_events = 0, 0, 0, 0
        at_least_one_mult_hmm_line = False
        eroded_annotations, padded_annotations = OrderedDict(), OrderedDict()
        hmm_failures = set()  # hm, does this duplicate info I'm already keeping track of in one of these other variables?
        errorfo = {}
        with open(annotation_fname, 'r') as hmm_csv_outfile:
            reader = csv.DictReader(hmm_csv_outfile)
            for padded_line in reader:  # line coming from hmm output is N-padded such that all the seqs are the same length

                utils.process_input_line(padded_line)
                n_lines_read += 1

                failed = self.check_did_bcrham_fail(padded_line, errorfo)
                if failed:
                    hmm_failures |= set(padded_line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                    continue

                uids = padded_line['unique_ids']
                uidstr = ':'.join(uids)
                padded_line['indelfos'] = [self.sw_info['indels'].get(uid, utils.get_empty_indel()) for uid in uids]  # reminder: hmm was given a sequence with any indels reversed (i.e. <self.sw_info['indels'][uid]['reverersed_seq']>)
                padded_line['duplicates'] = [self.duplicates.get(uid, []) for uid in uids]

                if not utils.has_d_gene(self.args.locus):
                    self.process_dummy_d_hack(padded_line)

                utils.add_implicit_info(self.glfo, padded_line, aligned_gl_seqs=self.aligned_gl_seqs)
                utils.process_per_gene_support(padded_line)  # switch per-gene support from log space to normalized probabilities
                if padded_line['invalid']:
                    n_invalid_events += 1
                    if self.args.debug:
                        print '      %s padded line invalid' % uidstr
                        utils.print_reco_event(padded_line, extra_str='    ', label='invalid:')
                    hmm_failures |= set(padded_line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                    continue

                assert uidstr not in padded_annotations
                padded_annotations[uidstr] = padded_line

                if len(uids) > 1:  # if there's more than one sequence, we need to use the padded line
                    at_least_one_mult_hmm_line = True
                    line_to_use = padded_line
                else:  # otherwise, the eroded line is kind of simpler to look at
                    # get a new dict in which we have edited the sequences to swap Ns on either end (after removing fv and jf insertions) for v_5p and j_3p deletions
                    eroded_line = utils.reset_effective_erosions_and_effective_insertions(self.glfo, padded_line, aligned_gl_seqs=self.aligned_gl_seqs)  #, padfo=self.sw_info)
                    if eroded_line['invalid']:  # not really sure why the eroded line is sometimes invalid when the padded line is not, but it's very rare and I don't really care, either
                        n_invalid_events += 1
                        hmm_failures |= set(eroded_line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                        continue
                    line_to_use = eroded_line
                    eroded_annotations[uidstr] = eroded_line  # these only get used if there aren't any multi-seq lines, so it's ok that they don't all get added if there is a multi seq line

                if self.args.debug or print_annotations:
                    self.print_hmm_output(line_to_use, print_true=True)

                n_events_processed += 1
                n_seqs_processed += len(uids)

                if pcounter is not None:
                    pcounter.increment(line_to_use)
                if true_pcounter is not None:
                    true_pcounter.increment(self.reco_info[uids[0]])  # NOTE doesn't matter which id you pass it, since they all have the same reco parameters

                if perfplotter is not None:
                    for iseq in range(len(uids)):  # TODO get perfplotter handling multi-seq lines
                        singlefo = utils.synthesize_single_seq_line(line_to_use, iseq)
                        if len(singlefo['naive_seq']) != len(self.reco_info[uids[iseq]]['naive_seq']):  # this seems to mostly depend on if it's the multi-hmm or not (but the real problem is that different js are different lengths)
                            singlefo = utils.reset_effective_erosions_and_effective_insertions(self.glfo, singlefo)
                        if uids[iseq] in self.sw_info['indels']:
                            print '    skipping performance evaluation of %s because of indels' % uids[iseq]  # I just have no idea how to handle naive hamming fraction when there's indels
                        else:
                            perfplotter.evaluate(self.reco_info[uids[iseq]], singlefo)

        # parameter and performance writing/plotting
        if pcounter is not None:
            if self.args.plotdir is not None:
                pcounter.plot(self.args.plotdir + '/hmm', only_csv=self.args.only_csv_plots, only_overall=self.args.only_overall_plots)
            pcounter.write(parameter_out_dir)
        if true_pcounter is not None:
            if self.args.plotdir is not None:
                true_pcounter.plot(self.args.plotdir + '/hmm-true', only_csv=self.args.only_csv_plots, only_overall=self.args.only_overall_plots)
            true_pcounter.write(parameter_out_dir + '-true')
        if perfplotter is not None:
            perfplotter.plot(self.args.plotdir + '/hmm', only_csv=self.args.only_csv_plots)

        print '        processed %d hmm output lines with %d sequences in %d events  (%d failures)' % (n_lines_read, n_seqs_processed, n_events_processed, len(hmm_failures))
        if n_invalid_events > 0:
            print '            %s skipped %d invalid events' % (utils.color('red', 'warning'), n_invalid_events)
        for ecode in errorfo:
            if ecode == 'no_path':
                print '          %s no valid paths: %s' % (utils.color('red', 'warning'), ' '.join(errorfo[ecode]))
            elif ecode == 'boundary':
                print '          %d boundary warnings' % len(errorfo[ecode])
                if self.args.debug:
                    print '                %s' % ' '.join(errorfo[ecode])
            else:
                print '          %s unknown ecode \'%s\': %s' % (utils.color('red', 'warning'), ecode, ' '.join(errorfo[ecode]))

        annotations_to_use = padded_annotations if at_least_one_mult_hmm_line else eroded_annotations  # if every query is a single-sequence query, then the output will be less confusing to people if the N padding isn't there. But you kinda need the padding in order to make the multi-seq stuff work

        self.check_for_unexpectedly_missing_keys(annotations_to_use, hmm_failures)  # NOTE not sure if it's really correct to use <annotations_to_use>, [maybe since <hmm_failures> has ones that failed the conversion to eroded line (and maybe other reasons)]

        # write output file
        if outfname is not None:
            self.write_annotations(annotations_to_use, outfname, hmm_failures)  # [0] takes the best annotation... if people want other ones later it's easy to change

        # annotation (VJ CDR3) clustering
        if self.args.annotation_clustering is not None:
            self.deal_with_annotation_clustering(annotations_to_use, outfname)

        os.remove(annotation_fname)
        return annotations_to_use

    # ----------------------------------------------------------------------------------------
    def deal_with_annotation_clustering(self, annotations, outfname):
        if self.args.annotation_clustering != 'vollmers':
            raise Exception('we only handle \'vollmers\' (vj cdr3 0.x) annotation clustering at the moment')

        # initialize output file
        if outfname is not None:
            outfile = open(outfname, 'w')  # NOTE overwrites annotation info that's already been written to <outfname>
            headers = ['n_clusters', 'threshold', 'partition']
            if not self.args.is_data:
                headers += ['ccf_under', 'ccf_over']
            writer = csv.DictWriter(outfile, headers)
            writer.writeheader()

        # perform annotation clustering for each threshold and write to file
        import annotationclustering
        for thresh in self.args.annotation_clustering_thresholds:
            partition = annotationclustering.vollmers(annotations, threshold=thresh, reco_info=self.reco_info)
            n_clusters = len(partition)
            if outfname is not None:
                row = {'n_clusters' : n_clusters, 'threshold' : thresh, 'partition' : utils.get_str_from_partition(partition)}
                if not self.args.is_data:
                    true_partition = utils.get_true_partition(self.reco_info)
                    ccfs = utils.new_ccfs_that_need_better_names(partition, true_partition, self.reco_info)
                    row['ccf_under'] = ccfs[0]
                    row['ccf_over'] = ccfs[1]
                writer.writerow(row)

        if outfname is not None:
            outfile.close()

    # ----------------------------------------------------------------------------------------
    def print_hmm_output(self, line, print_true=False):
        if print_true and not self.args.is_data:  # first print true event (if this is simulation)
            utils.print_true_events(self.glfo, self.reco_info, line)

        if len(line['unique_ids']) > 1:  # make it easier to cut and paste for --queries (it'd be nice if this could go on the same line as 'inferred:', but then for really big clusters it pushes the 'inserts' label out of alignment
            print '          ' + ':'.join(line['unique_ids'])
        label = 'inferred:'
        if self.args.seed_unique_id is not None and self.args.seed_unique_id in line['unique_ids']:
            label += '   (found %d sequences clonal to seed %s)' % (len(line['unique_ids']), self.args.seed_unique_id)
        utils.print_reco_event(line, extra_str='    ', label=label, seed_uid=self.args.seed_unique_id)

    # ----------------------------------------------------------------------------------------
    def write_annotations(self, annotations, outfname, hmm_failures):
        outpath = outfname
        if outpath[0] != '/':  # if full output path wasn't specified on the command line, write to current directory
            outpath = os.getcwd() + '/' + outpath

        with open(outpath, 'w') as outfile:
            writer = csv.DictWriter(outfile, utils.annotation_headers)
            writer.writeheader()
            # hmm annotations
            for line in annotations.values():
                outline = utils.get_line_for_output(line)  # convert lists to colon-separated strings and whatnot (doesn't modify <line>
                outline = {k : v for k, v in outline.items() if k in utils.annotation_headers}  # remove the columns we don't want to output
                writer.writerow(outline)

            # write empty lines for seqs that failed either in sw or the hmm
            for fid in self.sw_info['failed-queries'] | hmm_failures:  # both use single-seq ids a.t.m.
                writer.writerow({'unique_ids' : fid, 'input_seqs' : ':'.join([self.input_info[f]['seqs'][0] for f in fid.split(':')])})  # .split() stuff is to handle in the future multi-seq ids

        # presto!
        if self.args.presto_output:
            raise Exception('needs updating')
            outstr = check_output(['mv', '-v', self.args.outfname, self.args.outfname + '.partis'])
            print '    backing up partis output before converting to presto: %s' % outstr.strip()

            prestoheader = utils.presto_headers.values()
            with open(outpath, 'w') as outfile:
                writer = csv.DictWriter(outfile, prestoheader)
                writer.writeheader()

                for full_line in annotations.values():
                    outline = copy.deepcopy(full_line)  # in case we modify it

                    utils.remove_all_implicit_info(outline)
                    utils.add_implicit_info(self.glfo, outline, aligned_gl_seqs=self.aligned_gl_seqs)

                    outline = utils.convert_to_presto_headers(outline)

                    outline = utils.get_line_for_output(outline)  # convert lists to colon-separated strings and whatnot
                    outline = {k : v for k, v in outline.items() if k in prestoheader}  # remove the columns we don't want to output

                    writer.writerow(outline)

                # and write empty lines for seqs that failed either in sw or the hmm
                if len(missing_input_keys) > 0:  # NOTE assumes it's already been set by the first loop
                    print 'missing %d input keys' % len(missing_input_keys)
                    for uid in missing_input_keys:
                        col = utils.presto_headers['unique_id']
                        writer.writerow({col : uid})

