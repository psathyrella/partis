import yaml
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
import traceback

import utils
import glutils
import indelutils
import treeutils
from glomerator import Glomerator
from clusterpath import ClusterPath
from waterer import Waterer
from parametercounter import ParameterCounter
from alleleclusterer import AlleleClusterer
from alleleremover import AlleleRemover
from allelefinder import AlleleFinder
from performanceplotter import PerformancePlotter
from partitionplotter import PartitionPlotter
from hist import Hist
import seqfileopener

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    """ Class to parse input files, start bcrham jobs, and parse/interpret bcrham output for annotation and partitioning """
    def __init__(self, args, glfo, input_info, simglfo, reco_info):
        self.args = args
        self.glfo = glfo
        self.input_info = input_info
        self.simglfo = simglfo
        self.reco_info = reco_info

        utils.prep_dir(self.args.workdir)
        self.my_gldir = self.args.workdir + '/' + glutils.glfo_dir
        if args.infname is not None:
            if self.args.sw_cachefname is None:
                self.sw_cache_path = self.args.parameter_dir + '/sw-cache-' + repr(abs(hash(''.join(self.input_info.keys()))))  # remain suffix-agnostic
            else:
                self.sw_cache_path = utils.getprefix(self.args.sw_cachefname)

        self.vs_info, self.sw_info = None, None
        self.duplicates = {}
        self.bcrham_proc_info = None
        self.timing_info = []  # it would be really nice to clean up both this and bcrham_proc_info
        self.istep = None  # stupid hack to get around network file system issues (see self.subworkidr()
        self.subworkdirs = []  # arg. same stupid hack

        self.unseeded_seqs = None  # all the queries that we *didn't* cluster with the seed uid
        self.small_cluster_seqs = None  # all the queries that we removed after a few partition steps 'cause they were in small clusters

        self.sw_param_dir, self.hmm_param_dir, self.multi_hmm_param_dir = ['%s/%s' % (self.args.parameter_dir, s) for s in ['sw', 'hmm', 'multi-hmm']]
        self.sub_param_dir = self.args.parameter_dir + '/' + self.args.parameter_type

        self.hmm_infname = self.args.workdir + '/hmm_input.csv'
        self.hmm_cachefname = self.args.workdir + '/hmm_cached_info.csv'
        self.hmm_outfname = self.args.workdir + '/hmm_output.csv'
        self.cpath_progress_dir = '%s/cluster-path-progress' % self.args.workdir  # write the cluster paths for each clustering step to separate files in this dir

        if self.args.outfname is not None:
            utils.prep_dir(dirname=None, fname=self.args.outfname, allow_other_files=True)

        self.deal_with_persistent_cachefile()

        self.cached_naive_hamming_bounds = self.args.naive_hamming_bounds  # just so we don't get them every iteration through the clustering loop

        self.aligned_gl_seqs = None
        if self.args.aligned_germline_fname is not None:
            self.aligned_gl_seqs = glutils.read_aligned_gl_seqs(self.args.aligned_germline_fname, self.glfo, self.args.locus, dont_warn_about_duplicates=True)

        self.action_fcns = {
            'cache-parameters'             : self.cache_parameters,
            'annotate'                     : self.annotate,
            'partition'                    : self.partition,
            'view-output'                  : self.read_existing_output,
            'view-annotations'             : self.read_existing_output,
            'view-partitions'              : self.read_existing_output,
            'plot-partitions'              : self.read_existing_output,
            'get-tree-metrics'             : self.read_existing_output,
            'view-alternative-annotations'  : self.view_alternative_annotations,
        }

    # ----------------------------------------------------------------------------------------
    def get_cpath_progress_fname(self, istep):
        return '%s/istep-%d.csv' % (self.cpath_progress_dir, istep)

    # ----------------------------------------------------------------------------------------
    def get_all_cpath_progress_fnames(self):
        assert len(os.listdir(self.cpath_progress_dir)) == self.istep + 1  # not really checking for anything that has a decent chance of happening, it's more because I get the file names in a slightly different loop in self.clean()
        return [self.get_cpath_progress_fname(istep) for istep in range(self.istep + 1)]

    # ----------------------------------------------------------------------------------------
    def run(self, actions):
        for tmpaction in actions:
            self.current_action = tmpaction  # NOTE gets changed on the fly below, I think just in self.get_cluster_annotations() (which is kind of hackey, but I can't figure out a way to improve on it that wouldn't involve wasting a foolish amount of time rewriting things. Bottom line is that the control flow for different actions is really complicated, and that complexity is going to show up somewhere)
            self.action_fcns[tmpaction]()

    # ----------------------------------------------------------------------------------------
    def clean(self):
        if self.args.new_allele_fname is not None:
            new_allele_region = 'v'
            new_alleles = [(g, seq) for g, seq in self.glfo['seqs'][new_allele_region].items() if glutils.is_snpd(g)]
            print '  writing %d new %s to %s' % (len(new_alleles), utils.plural_str('allele', len(new_alleles)), self.args.new_allele_fname)
            with open(self.args.new_allele_fname, 'w') as outfile:
                for name, seq in new_alleles:
                    outfile.write('>%s\n' % name)
                    outfile.write('%s\n' % seq)

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

        for subd in self.subworkdirs:
            if os.path.exists(subd):  # if there was only one proc for this step, it'll have already been removed
                os.rmdir(subd)

        if os.path.exists(self.cpath_progress_dir):  # only exists for partitioning
            for cpfname in self.get_all_cpath_progress_fnames():
                os.remove(cpfname)
            os.rmdir(self.cpath_progress_dir)

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
    def run_waterer(self, count_parameters=False, write_parameters=False, write_cachefile=False, look_for_cachefile=False, dbg_str=''):
        print 'smith-waterman%s' % (('  (%s)' % dbg_str) if dbg_str != '' else '')
        sys.stdout.flush()

        self.vs_info = None  # should already be None, but we want to make sure (if --no-sw-vsearch is set we need it to be None, and if we just removed unlikely alleles we need to rerun vsearch with the likely alleles)
        if not self.args.no_sw_vsearch:
            self.set_vsearch_info(get_annotations=True)

        pre_failed_queries = self.sw_info['failed-queries'] if self.sw_info is not None else None  # don't re-run on failed queries if this isn't the first sw run (i.e., if we're parameter caching)
        waterer = Waterer(self.args, self.glfo, self.input_info, self.simglfo, self.reco_info,
                          count_parameters=count_parameters,
                          parameter_out_dir=self.sw_param_dir if write_parameters else None,
                          plot_annotation_performance=self.args.plot_annotation_performance,
                          duplicates=self.duplicates, pre_failed_queries=pre_failed_queries, aligned_gl_seqs=self.aligned_gl_seqs, vs_info=self.vs_info)

        cachefname = self.sw_cache_path + ('.yaml' if self.args.sw_cachefname is None else utils.getsuffix(self.args.sw_cachefname))  # use yaml, unless csv was explicitly set on the command line
        if look_for_cachefile:
            if os.path.exists(self.sw_cache_path + '.csv'):  # ...but if there's already an old csv, use that
                cachefname = self.sw_cache_path + '.csv'
        else:  # i.e. if we're not explicitly told to look for it (and it exists) then it should be out of date
            waterer.clean_cache(self.sw_cache_path)  # hm, should this be <cachefname> instead of <self.sw_cache_path>? i mean they're the same, but still
        if look_for_cachefile and os.path.exists(cachefname):  # run sw if we either don't want to do any caching (None) or if we are planning on writing the results after we run
            waterer.read_cachefile(cachefname)
        else:
            if look_for_cachefile:
                print '    couldn\'t find sw cache file %s, so running sw%s' % (cachefname, ' (this is probably because --seed-seq/--seed-unique-id is set to a sequence that wasn\'t in the input file on which we cached parameters [if it\'s inconvenient to put your seed sequences in your input file, you can avoid this by putting them instead in separate file and set --queries-to-include-fname])' if self.args.seed_unique_id is not None else '')
            waterer.run(cachefname if write_cachefile else None)

        self.sw_info = waterer.info
        for uid, dupes in waterer.duplicates.items():  # <waterer.duplicates> is <self.duplicates> OR'd into any new duplicates from this run
            self.duplicates[uid] = dupes

        # utils.compare_vsearch_to_sw(self.sw_info, self.vs_info)  # only compares indels a.t.m.

        # NOTE neither of these really works right here, since any time we change the germline set we in really need to go back and rerun sw. But we can't do them betweeen allele finding and running parameter counting sw, since then the counts are wrong
        # # d j allele removal based on snps/counts (just printing for now)
        # if count_parameters and self.current_action == 'cache-parameters':  # I'm not sure this is precisely the criterion I want, but it does the job of not running dj removal printing when we're just annotating with existing parameters (which was causing a crash [key error] with inconsistent glfo)
        #     print ' testing d+j snp-based allele removal'
        #     alremover = AlleleRemover(self.glfo, self.args, simglfo=self.simglfo, reco_info=self.reco_info)
        #     alremover.finalize(gene_counts=None, annotations={q : self.sw_info[q] for q in self.sw_info['queries']}, regions=['d', 'j'], debug=self.args.debug_allele_finding)
        #     print '  (not actually removing d and j alleleremover genes)'
        #     # glutils.remove_genes(self.glfo, alremover.genes_to_remove, debug=True)
        #     # glutils.write_glfo('xxx _output/glfo-test', self.glfo)
        # gene name-based allele removal:
        # if self.args.n_max_alleles_per_gene is not None:  # it would be nice to use AlleleRemover for this, but that's really set up more as a precursor to allele finding, so it ends up being pretty messy to implement
        #     gene_counts = utils.get_gene_counts_from_annotations({q : self.sw_info[q] for q in self.sw_info['queries']})
        #     glutils.remove_extra_alleles_per_gene(self.glfo, self.args.n_max_alleles_per_gene, gene_counts)
        #     # glutils.write_glfo(self.sw_param_dir + '/' + glutils.glfo_dir, self.glfo)  # don't need to rewrite glfo above, since we haven't yet written parameters, but now we do/have

    # ----------------------------------------------------------------------------------------
    def set_vsearch_info(self, get_annotations=False):  # NOTE setting match:mismatch to optimized values from sw (i.e. 5:-4) results in much worse shm indel performance, so we leave it at the vsearch defaults ('2:-4')
        seqs = {sfo['unique_ids'][0] : sfo['seqs'][0] for sfo in self.input_info.values()}
        self.vs_info = utils.run_vsearch('search', seqs, self.args.workdir + '/vsearch', threshold=0.3, glfo=self.glfo, print_time=True, vsearch_binary=self.args.vsearch_binary, get_annotations=get_annotations, no_indels=self.args.no_indels)

    # ----------------------------------------------------------------------------------------
    def cache_parameters(self):
        print 'caching parameters'

        # remove unlikely alleles (can only remove v alleles here, since we're using vsearch annotations, but that's ok since it's mostly a speed optimization)
        if not self.args.dont_remove_unlikely_alleles:
            self.set_vsearch_info(get_annotations=(self.args.debug_allele_finding and self.args.is_simu))  # we only use the annotations to print some debug info in alleleremover
            alremover = AlleleRemover(self.glfo, self.args, simglfo=self.simglfo, reco_info=self.reco_info)
            alremover.finalize({'v' : self.vs_info['gene-counts']}, annotations=(None if len(self.vs_info['annotations']) == 0 else self.vs_info['annotations']), debug=self.args.debug_allele_finding)
            glutils.remove_genes(self.glfo, alremover.genes_to_remove)
            self.vs_info = None  # don't want to keep this around, since it has alignments against all the genes we removed (also maybe memory control)
            alremover = None  # memory control (not tested)

        # (re-)add [new] alleles
        if self.args.allele_cluster:
            self.run_waterer(dbg_str='new-allele clustering')
            alclusterer = AlleleClusterer(self.args, glfo=self.glfo, reco_info=self.reco_info, simglfo=self.simglfo)
            alcluster_alleles = alclusterer.get_alleles(self.sw_info, debug=self.args.debug_allele_finding, plotdir=None if self.args.plotdir is None else self.args.plotdir + '/sw/alcluster')
            if len(alcluster_alleles) > 0:
                glutils.add_new_alleles(self.glfo, alcluster_alleles.values(), use_template_for_codon_info=False, simglfo=self.simglfo, debug=True)
                if self.aligned_gl_seqs is not None:
                    glutils.add_missing_alignments(self.glfo, self.aligned_gl_seqs, debug=True)
            alclusterer = None

        if not self.args.dont_find_new_alleles:
            self.run_waterer(dbg_str='new-allele fitting')
            alfinder = AlleleFinder(self.glfo, self.args)
            new_allele_info = alfinder.increment_and_finalize(self.sw_info, debug=self.args.debug_allele_finding)  # incrementing and finalizing are intertwined since it needs to know the distribution of 5p and 3p deletions before it can increment
            if self.args.plotdir is not None:
                alfinder.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots)
            if len(new_allele_info) > 0:
                glutils.restrict_to_genes(self.glfo, list(self.sw_info['all_best_matches']))
                glutils.add_new_alleles(self.glfo, new_allele_info, debug=True, simglfo=self.simglfo, use_template_for_codon_info=False)  # <remove_template_genes> stuff is handled in <new_allele_info> (also note, can't use template for codon info since we may have already removed it)
                if self.aligned_gl_seqs is not None:
                    glutils.add_missing_alignments(self.glfo, self.aligned_gl_seqs, debug=True)

        # get and write sw parameters
        self.run_waterer(count_parameters=True, write_parameters=True, write_cachefile=True, dbg_str='writing parameters')
        self.write_hmms(self.sw_param_dir)  # note that this modifies <self.glfo>
        if self.args.only_smith_waterman:
            if self.args.outfname is not None:  # NOTE this is _not_ identical to the sw cache file (e.g. padding, failed query writing, plus probably other stuff)
                self.write_output(None, set(), write_sw=True)
            return

        # get and write hmm parameters
        print 'hmm'
        sys.stdout.flush()
        _, annotations, hmm_failures = self.run_hmm('viterbi', parameter_in_dir=self.sw_param_dir, parameter_out_dir=self.hmm_param_dir, count_parameters=True)
        if self.args.outfname is not None:
            self.write_output(annotations.values(), hmm_failures)
        self.write_hmms(self.hmm_param_dir)  # note that this modifies <self.glfo>

    # ----------------------------------------------------------------------------------------
    def annotate(self):
        print 'annotating'
        if self.sw_info is None:
            self.run_waterer(look_for_cachefile=not self.args.write_sw_cachefile, write_cachefile=self.args.write_sw_cachefile, count_parameters=self.args.count_parameters)
        if self.args.only_smith_waterman:
            if self.args.outfname is not None:  # NOTE this is _not_ identical to the sw cache file (e.g. padding, failed query writing, plus probably other stuff)
                self.write_output(None, set(), write_sw=True)  # note that if you're auto-parameter caching, this will just be rewriting an sw output file that's already there from parameter caching, but oh, well. If you're setting --only-smith-waterman and not using cache-parameters, you have only yourself to blame
            return
        print 'hmm'
        _, annotations, hmm_failures = self.run_hmm('viterbi', parameter_in_dir=self.sub_param_dir, count_parameters=self.args.count_parameters, parameter_out_dir=self.multi_hmm_param_dir if self.args.parameter_out_dir is None else self.args.parameter_out_dir)
        if self.args.get_tree_metrics:
            self.calculate_tree_metrics(annotations, cpath=None)  # adds tree metrics to <annotations>
        if self.args.outfname is not None:
            self.write_output(annotations.values(), hmm_failures)

    # ----------------------------------------------------------------------------------------
    def calculate_tree_metrics(self, annotations, cpath=None):
        if self.current_action == 'get-tree-metrics' and self.args.input_metafname is not None:  # presumably if you're running 'get-tree-metrics' with --input-metafname set, that means you didn't add the affinities (+ other metafo) when you partitioned, so we need to add it now
            seqfileopener.read_input_metafo(self.args.input_metafname, annotations.values(), debug=True)
        treeutils.calculate_tree_metrics(annotations, self.args.min_tree_metric_cluster_size, self.args.lb_tau, cpath=cpath, reco_info=self.reco_info, treefname=self.args.treefname,
                                         use_true_clusters=self.reco_info is not None, base_plotdir=self.args.plotdir, ete_path=self.args.ete_path, workdir=self.args.workdir, debug=self.args.debug)

    # ----------------------------------------------------------------------------------------
    def parse_existing_annotations(self, annotation_lines, ignore_args_dot_queries=False, process_csv=False):
        n_queries_read = 0
        failed_query_strs = set()
        annotations = OrderedDict()
        for line in annotation_lines:
            if process_csv:
                utils.process_input_line(line)
            uidstr = ':'.join(line['unique_ids'])
            if ('invalid' in line and line['invalid']) or line['v_gene'] == '':  # first way is the new way, but we have to check the empty-v-gene way too for old files
                failed_query_strs.add(uidstr)
                continue
            if self.args.queries is not None and not ignore_args_dot_queries:  # second bit is because when printing subcluster naive seqs, we want to read all the ones that have any overlap with self.args.queries, not just the exact cluster of self.args.queries
                if len(set(self.args.queries) & set(line['unique_ids'])) == 0:  # actually make sure this is the precise set of queries we want (note that --queries and line['unique_ids'] are both ordered, and this ignores that... oh, well, sigh.)
                    continue
            if self.args.reco_ids is not None and line['reco_id'] not in self.args.reco_ids:
                continue
            utils.add_implicit_info(self.glfo, line)
            annotations[uidstr] = line

            n_queries_read += 1
            if self.args.n_max_queries > 0 and n_queries_read >= self.args.n_max_queries:
                break

        if len(failed_query_strs) > 0:
            print '\n%d failed queries' % len(failed_query_strs)

        return annotations

    # ----------------------------------------------------------------------------------------
    def view_alternative_annotations(self):
        print '  %s getting alternative annotation information from existing output file. These results will only be meaningful if you had --calculate-alternative-annotations set when writing the output file (so that all subcluster annotations were stored). We can\'t check for that here directly, so instead we print this warning to make sure you had it set ;-)' % utils.color('yellow', 'note')

        # we used to require that you set --queries to tell us which to get, but I think now it makes sense to by default just get all of them (but not sure enough to delete this yet)
        # if self.args.queries is None:
        #     _, cpath = self.read_existing_output(read_partitions=True)
        #     clusterstrs = []
        #     for cluster in sorted(cpath.partitions[cpath.i_best], key=len, reverse=True):
        #         clusterstrs.append('      %s' % ':'.join(cluster))
        #     raise Exception('in order to view alternative annotations, you have to specify (with --queries) a cluster from the final partition. Choose from the following:\n%s' % '\n'.join(clusterstrs))

        cluster_annotations, cpath = self.read_existing_output(ignore_args_dot_queries=True, read_partitions=True, read_annotations=True)  # note that even if we don't need the cpath to re-write output below, we need to set read_annotations=True, since the fcn gets confused otherwise and doesn't read the right cluster annotation file (for deprecated csv files)

        # just for backwards compatibility with old style csv output/cache file (we used to keep the hmm cache file around to read here, but now we just write everything we need to the regular output file)
        cachefo = None
        if os.path.exists(self.hmm_cachefname):
            cachefo = self.read_hmm_cachefile()
            print '  %s persistent hmm cache file %s exists (probably just copied it from %s), so we\'re assuming that all the sub cluster annotations *weren\'t* written to the output file, i.e. this is old-style output. This\'ll still work, it will just be missing a lot of the v/d/j gene info)' % (utils.color('yellow', 'note'), self.hmm_cachefname, self.args.persistent_cachefname)
            # fuck it, this isn't worth it (but not yet quite willing to delete the effort)
            # print '  %s persistent hmm cache file %s exists (probably just copied it from %s), so we\'re assuming that the sub cluster annotations *weren\'t* written to the output file, i.e. this is old-style output (if you set --infname (and --parameter-dir is either set or the default corresponds to an existing parameter directory), we\'ll just run the cluster annotations for cache file uidstrs right now. Otherwise this\'ll still work, it just won\'t have the v/d/j gene info)' % (utils.color('yellow', 'note'), self.hmm_cachefname, self.args.persistent_cachefname)
            # if self.args.parameter_dir is not None and self.args.infname is not None:
            #     self.run_waterer(look_for_cachefile=not self.args.write_sw_cachefile, write_cachefile=self.args.write_sw_cachefile)  # need sw info to run the hmm (this is not really a good idea to be running new stuff during an action that's supposed to be running on existing output, but it's only for backwards compatibility, so oh well)
            #     cache_clusters = [set(uidstr.split(':')) for uidstr in cachefo]  # all of 'em
            #     cache_clusters = [sc for sc in cache_clusters if sc <= uids_of_interest]  # just the ones that are composed entire of queries in uids_of_interest
            #     # DAMMIT neither of these ways of handling missing sw queries work. oh, well
            #     # cache_clusters = [sc for sc in cache_clusters if len(sc & set(self.sw_info['failed-queries'])) == 0]  # and remove any sw failues (arg, this shouldn't be needed, why are they failing now when they didn't before?)
            #     # cache_clusters = [sc - set(self.sw_info['failed-queries']) for sc in cache_clusters if len(sc & set(self.sw_info['failed-queries'])) == 0]  # and remove any sw failues (arg, this shouldn't be needed, why are they failing now when they didn't before?)
            #     cache_clusters = set([tuple(sc) for sc in cache_clusters])
            #     cluster_annotations = self.get_cluster_annotations(cpath=None, clusters_to_annotate=cache_clusters, dont_write_anything=True, extra_dbg_str=' (using clusters from old hmm cache file, not final partition)')  # replaces <cluster_annotations> that was just read the existing output file:
            #     cachefo = cluster_annotations  # account for the clusters we removed because they had queries that weren't in uids_of_interest

        clusters_to_use = cpath.partitions[cpath.i_best] if self.args.queries is None else [self.args.queries]
        clusters_in_partitions = set([':'.join(c) for partition in cpath.partitions for c in partition])
        for cluster in sorted(clusters_to_use, key=len, reverse=True):
            self.get_subcluster_naive_seqs(cluster, cluster_annotations, clusters_in_partitions, cachefo=cachefo, debug=True)

        print '  note: rewriting output file with newly-calculated alternative annotation info'
        self.write_output(cluster_annotations.values(), set(), cpath=cpath, dont_write_failed_queries=True)  # I *think* we want <dont_write_failed_queries> set, because the failed queries should already have been written, so now they'll just be mixed in with the others in <annotations>

    # ----------------------------------------------------------------------------------------
    def print_results(self, cpath, annotations):
        seed_uid = self.args.seed_unique_id
        if cpath is not None:
            # it's expected that sometimes you'll write a seed partition cpath, but then when you read the file you don't bother to seed the seed id on the command line. The reverse, however, shouldn't happen
            if seed_uid is not None and cpath.seed_unique_id != seed_uid:
                print '  %s seed uids from args and cpath don\'t match %s %s ' % (utils.color('red', 'error'), self.args.seed_unique_id, cpath.seed_unique_id)
            seed_uid = cpath.seed_unique_id
            print utils.color('green', 'partitions:')
            cpath.print_partitions(abbreviate=self.args.abbreviate, reco_info=self.reco_info, highlight_cluster_indices=self.args.cluster_indices, calc_missing_values=('all' if cpath.n_seqs() < 500 else 'best'))
            if not self.args.is_data and self.reco_info is not None:  # if we're reading existing output, it's pretty common to not have the reco info even when it's simulation, since you have to also pass in the simulation input file on the command line
                true_cp = ClusterPath(seed_unique_id=self.args.seed_unique_id)
                true_cp.add_partition(utils.get_true_partition(self.reco_info), -1., 1)
                print 'true:'
                true_cp.print_partitions(self.reco_info, print_header=False, calc_missing_values='best')

        if len(annotations) > 0:
            print utils.color('green', 'annotations:')
            sorted_annotations = sorted(annotations.values(), key=lambda l: len(l['unique_ids']), reverse=True)
            if self.args.cluster_indices is not None:
                sorted_annotations = [sorted_annotations[iclust] for iclust in self.args.cluster_indices]
            for line in sorted_annotations:
                if self.args.only_print_best_partition and cpath is not None and cpath.i_best is not None and line['unique_ids'] not in cpath.partitions[cpath.i_best]:
                    continue
                if (self.args.only_print_seed_clusters or self.args.seed_unique_id is not None) and seed_uid not in line['unique_ids']:  # we only use the seed id from the command line here, so you can print all the clusters even if you ran seed partitioning UPDATE wait did I change my mind? need to check
                    continue
                label, post_label = '', ''
                if self.args.infname is not None and self.reco_info is not None:
                    utils.print_true_events(self.simglfo, self.reco_info, line, extra_str='  ')
                    label = 'inferred:'
                if cpath is not None and cpath.i_best is not None and line['unique_ids'] in cpath.partitions[cpath.i_best]:
                    post_label += ' (from best partition)'
                utils.print_reco_event(line, extra_str='  ', label=label, post_label=post_label, seed_uid=seed_uid)

    # ----------------------------------------------------------------------------------------
    def read_existing_output(self, outfname=None, ignore_args_dot_queries=False, read_partitions=False, read_annotations=False):
        if outfname is None:
            outfname = self.args.outfname

        annotation_lines = []
        cpath = None
        tmpact = self.current_action  # just a shorthand for brevity
        if utils.getsuffix(outfname) == '.csv':  # old way
            if tmpact == 'view-partitions' or tmpact == 'plot-partitions' or tmpact == 'view-output' or tmpact == 'get-tree-metrics' or read_partitions:
                cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id, fname=outfname)
            if tmpact == 'view-annotations' or tmpact == 'plot-partitions' or tmpact == 'view-output' or tmpact == 'get-tree-metrics' or read_annotations:
                csvfile = open(outfname if cpath is None else self.args.cluster_annotation_fname)  # closes on function exit, and no this isn't a great way of doing it (but it needs to stay open for the loop)
                annotation_lines = csv.DictReader(csvfile)
                if 'unique_ids' not in annotation_lines.fieldnames:
                    raise Exception('not an annotation file: %s' % outfname)
        elif utils.getsuffix(outfname) == '.yaml':  # new way
            # NOTE replaces <self.glfo>, which is definitely what we want (that's the point of putting glfo in the yaml file), but it's still different behavior than if reading a csv
            assert self.glfo is None  # make sure bin/partis successfully figured out that we would be reading the glfo from the yaml output file
            self.glfo, annotation_lines, cpath = utils.read_yaml_output(outfname, n_max_queries=self.args.n_max_queries, dont_add_implicit_info=True, seed_unique_id=self.args.seed_unique_id)  # add implicit info below, so we can skip some of 'em
        else:
            raise Exception('unhandled annotation file suffix %s' % outfname)

        annotations = self.parse_existing_annotations(annotation_lines, ignore_args_dot_queries=ignore_args_dot_queries, process_csv=utils.getsuffix(outfname) == '.csv')

        if tmpact == 'get-tree-metrics':
            self.calculate_tree_metrics(annotations, cpath=cpath)  # adds tree metrics to <annotations>
            print '  note: rewriting output file %s with newly-calculated tree metrics' % outfname
            self.write_output(annotations.values(), set(), cpath=cpath, dont_write_failed_queries=True)  # I *think* we want <dont_write_failed_queries> set, because the failed queries should already have been written, so now they'll just be mixed in with the others in <annotations>

        if tmpact == 'plot-partitions':
            partplotter = PartitionPlotter(self.args)
            partplotter.plot(self.args.plotdir + '/partitions', partition=cpath.partitions[cpath.i_best], annotations=annotations, reco_info=self.reco_info, cpath=cpath)

        if tmpact in ['view-output', 'view-annotations', 'view-partitions']:
            self.print_results(cpath, annotations)

        return annotations, cpath

    # ----------------------------------------------------------------------------------------
    def partition(self):
        """ Partition sequences in <self.input_info> into clonally related lineages """
        print 'partitioning'
        if self.sw_info is None:
            self.run_waterer(look_for_cachefile=not self.args.write_sw_cachefile, write_cachefile=self.args.write_sw_cachefile, count_parameters=self.args.count_parameters)  # run smith-waterman
        if len(self.sw_info['queries']) == 0:
            if self.args.outfname is not None:
                ClusterPath().write(self.args.outfname, False)
            return
        if self.args.only_smith_waterman:
            return

        print 'hmm'

        # pre-cache hmm naive seq for each single query NOTE <self.current_action> is still 'partition' for this (so that we build the correct bcrham command line)
        if self.args.persistent_cachefname is None or not os.path.exists(self.hmm_cachefname):  # if the default (no persistent cache file), or if a not-yet-existing persistent cache file was specified
            print 'caching all %d naive sequences' % len(self.sw_info['queries'])  # this used to be a speed optimization, but now it's so we have better naive sequences for the pre-bcrham collapse
            self.run_hmm('viterbi', self.sub_param_dir, n_procs=self.auto_nprocs(len(self.sw_info['queries'])), precache_all_naive_seqs=True)

        if self.args.naive_vsearch or self.args.naive_swarm:
            cpath = self.cluster_with_naive_vsearch_or_swarm(parameter_dir=self.sub_param_dir)
        else:
            cpath = self.cluster_with_bcrham()

        best_cluster_annotations = self.get_cluster_annotations(cpath)

        if self.args.plotdir is not None:
            partplotter = PartitionPlotter(self.args)
            partplotter.plot(self.args.plotdir + '/partitions', partition=cpath.partitions[cpath.i_best], annotations=best_cluster_annotations, reco_info=self.reco_info, cpath=cpath)

        if self.args.debug:
            print 'final'
            self.print_results(cpath, best_cluster_annotations)

        self.check_partition(cpath.partitions[cpath.i_best])

    # ----------------------------------------------------------------------------------------
    def split_seeded_clusters(self, old_cpath):
        seeded_clusters, unseeded_clusters = utils.split_partition_with_criterion(old_cpath.partitions[old_cpath.i_best_minus_x], lambda cluster: self.args.seed_unique_id in cluster)
        self.unseeded_seqs = [uid for uclust in unseeded_clusters for uid in uclust]  # NOTE we no longer expect them to all be singletons, since we're merging queries with identical naive seqs before passing to glomerator.cc
        seeded_singleton_set = set([uid for sclust in seeded_clusters for uid in sclust])  # in case there's duplicates
        seeded_partition = utils.collapse_naive_seqs(self.synth_sw_info(seeded_singleton_set), split_by_cdr3=True)
        seeded_cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
        seeded_cpath.add_partition(seeded_partition, -1., 1)
        print '      removed %d sequences in unseeded clusters,' % len(self.unseeded_seqs),
        print 'split %d seeded clusters into %d singletons, and merged these into %d clusters with identical naive seqs' % (len(seeded_clusters), len(seeded_singleton_set), len(seeded_cpath.partitions[seeded_cpath.i_best_minus_x]))

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
        times_to_try_this_n_procs = max(4, last_n_procs)  # if we've already milked this number of procs for most of what it's worth (once you get down to 2 or 3, you don't want to go lower)
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
    def merge_shared_clusters(self, cpath, debug=False):  # replace the most likely partition with a new partition in which any clusters that share a sequence have been merged
        # cpath.partitions[cpath.i_best] = [['a', 'b', 'c', 'e'], ['d'], ['f', 'a'], ['g'], ['h'], ['i'], ['j', 'a'], ['x', 'y', 'z', 'd'], ['xx', 'x']]
        partition = cpath.partitions[cpath.i_best]

        if debug:
            print 'merging shared clusters'
            cpath.print_partitions()

        # find every pair of clusters that has some overlap
        cluster_groups = []
        if debug:
            print ' making cluster_groups'
        for iclust in range(len(partition)):
            for jclust in range(iclust + 1, len(partition)):
                if len(set(partition[iclust]) & set(partition[jclust])) > 0:
                    if debug:
                        print '  %d %d' % (iclust, jclust)
                    cluster_groups.append(set([iclust, jclust]))

        # merge these pairs of clusters into groups
        while True:
            no_more_merges = True
            for cp1, cp2 in itertools.combinations(cluster_groups, 2):
                if len(cp1 & cp2) > 0:
                    if debug:
                        print '  merging %s and %s' % (cp1, cp2)
                    cluster_groups.append(cp1 | cp2)
                    cluster_groups.remove(cp1)
                    cluster_groups.remove(cp2)
                    no_more_merges = False
                    break  # we've modified it now, so we have to go back and remake the iterator
            if no_more_merges:
                break

        # actually merge the groups of clusters
        new_clusters = []
        for cgroup in cluster_groups:
            new_clusters.append(list(set([uid for iclust in cgroup for uid in partition[iclust]])))
        if debug:
            print ' removing'
        for iclust in sorted([i for cgroup in cluster_groups for i in cgroup], reverse=True):
            if debug:
                print '    %d' % iclust
            partition.pop(iclust)
        for nclust in new_clusters:
            partition.append(nclust)

        if debug:
            cpath.print_partitions()

    # ----------------------------------------------------------------------------------------
    def are_we_finished_clustering(self, n_procs, cpath):
        if n_procs == 1:
            return True
        elif self.args.n_final_clusters is not None and len(cpath.partitions[cpath.i_best]) <= self.args.n_final_clusters:  # NOTE I *think* I want the best, not best-minus-x here (hardish to be sure a.t.m., since I'm not really using the minus-x part right now)
            print '  stopping with %d (<= %d) clusters' % (len(cpath.partitions[cpath.i_best]), self.args.n_final_clusters)
            return True
        elif self.args.max_cluster_size is not None and max([len(c) for c in cpath.partitions[cpath.i_best]]) > self.args.max_cluster_size:  # NOTE I *think* I want the best, not best-minus-x here (hardish to be sure a.t.m., since I'm not really using the minus-x part right now)
            print '   --max-cluster-size (partitiondriver): stopping with a cluster of size %d (> %d)' % (max([len(c) for c in cpath.partitions[cpath.i_best]]), self.args.max_cluster_size)
            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------
    def synth_sw_info(self, queries):  # only used for passing info to utils.collapse_naive_seqs()
        # this uses the cached hmm naive seqs (since we have them and they're better) but then later we pass the hmm the sw annotations, so we have to make sure the sw cdr3 length is the same within each cluster (it's very rare that it isn't)
        synth_sw_info = {q : {'naive_seq' : s, 'cdr3_length' : self.sw_info[q]['cdr3_length']} for q, s in self.get_cached_hmm_naive_seqs(queries).items()}  # NOTE code duplication in cluster_with_bcrham()
        synth_sw_info['queries'] = synth_sw_info.keys()
        return synth_sw_info

    # ----------------------------------------------------------------------------------------
    def init_cpath(self, n_procs):
        initial_nseqs = len(self.sw_info['queries'])  # NOTE um, maybe I should change this to the number of clusters, now that we're doing some preclustering here?
        initial_nsets = utils.collapse_naive_seqs(self.synth_sw_info(self.sw_info['queries']), split_by_cdr3=True, debug=True)
        cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
        cpath.add_partition(initial_nsets, logprob=0., n_procs=n_procs)  # NOTE sw info excludes failed sequences (and maybe also sequences with different cdr3 length)
        os.makedirs(self.cpath_progress_dir)
        if self.args.debug:
            cpath.print_partitions(abbreviate=self.args.abbreviate, reco_info=self.reco_info)
        return cpath, initial_nseqs

    # ----------------------------------------------------------------------------------------
    def merge_cpaths_from_previous_steps(self, final_cpath, debug=False):
        if debug:
            print 'final (unmerged) cpath:'
            final_cpath.print_partitions(abbreviate=self.args.abbreviate)
            print ''

        n_before, n_after = self.args.n_partitions_to_write, self.args.n_partitions_to_write  # this takes more than we need, since --n-partitions-to-write is the *full* width, not half-width, but oh, well
        if self.args.debug or self.args.calculate_alternative_annotations or self.args.get_tree_metrics:  # take all of 'em
            n_before, n_after = sys.maxint, sys.maxint
        elif self.args.write_additional_cluster_annotations is not None:
            n_before, n_after = [max(waca, n_) for waca, n_ in zip(self.args.write_additional_cluster_annotations, (n_before, n_after))]
        # NOTE we don't actually do anything with <n_after>, since we can't add any extra partitions here (well, we don't want to)

        cpfnames = self.get_all_cpath_progress_fnames()  # last one corresponds to <final_cpath>

        if final_cpath.i_best >= n_before or len(cpfnames) < 2:  # if we already have enough partitions, or if there was only one step, there's nothing to do
            if debug:
                print '   nothing to merge'
            return final_cpath

        icpfn = len(cpfnames) - 1
        merged_cp = ClusterPath(fname=cpfnames[icpfn], seed_unique_id=self.args.seed_unique_id)
        assert merged_cp.partitions[merged_cp.i_best] == final_cpath.partitions[final_cpath.i_best]  # shouldn't really be necessary, and is probably kind of slow
        while merged_cp.i_best < n_before and icpfn > 0:
            icpfn -= 1
            previous_cp = ClusterPath(fname=cpfnames[icpfn], seed_unique_id=self.args.seed_unique_id)
            for ip in range(len(merged_cp.partitions)):
                if len(merged_cp.partitions[ip]) == len(previous_cp.partitions[-1]):  # skip identical partitions (don't bother checking unless they're at least the same length
                    if set([tuple(c) for c in merged_cp.partitions[ip]]) == set([tuple(c) for c in previous_cp.partitions[-1]]):  # no, they're not always in the same order (I think because they get parcelled out to different processes, and then read back in random order)
                        if math.isinf(merged_cp.logprobs[ip]) and not math.isinf(previous_cp.logprobs[-1]):  # it should only be possible for the *later* partition to have non-infinite logprob, since we usually only calculate full logprobs in the last clustering step (which is why we're taking the later partition, from merged_cp), so print an error if the earlier one, that we're about to throw away, is the one that's non-infinite
                            print '%s earlier partition (that we\'re discarding) has non-infinite logprob %f, while later partition\'s is infinite %f' % (utils.color('red', 'error'), previous_cp.logprobs[-1], merged_cp.logprobs[ip])
                        previous_cp.remove_partition(len(previous_cp.partitions) - 1)  # remove it from previous_cp, since we want the one that may have a logprob set (and which has smaller n_procs, although I don't think we care about that)
                previous_cp.add_partition(list(merged_cp.partitions[ip]), merged_cp.logprobs[ip], merged_cp.n_procs[ip])
            merged_cp = previous_cp
            assert merged_cp.partitions[merged_cp.i_best] == final_cpath.partitions[final_cpath.i_best]  # shouldn't really be necessary, and is probably kind of slow
            if debug:
                print '%s' % utils.color('red', str(icpfn))
                merged_cp.print_partitions()

        if icpfn > 0:
            print '  %s not merging entire cpath history' % utils.color('yellow', 'note')

        return merged_cp

    # ----------------------------------------------------------------------------------------
    def cluster_with_bcrham(self):
        tmpstart = time.time()
        n_procs = self.args.n_procs
        cpath, initial_nseqs = self.init_cpath(n_procs)
        n_proc_list = []
        self.istep = 0
        start = time.time()
        while n_procs > 0:
            print '%d clusters with %d proc%s' % (len(cpath.partitions[cpath.i_best_minus_x]), n_procs, utils.plural(n_procs))  # NOTE that a.t.m. i_best and i_best_minus_x are usually the same, since we're usually not calculating log probs of partitions (well, we're trying to avoid calculating any extra log probs, which means we usually don't know the log prob of the entire partition)
            cpath, _, _ = self.run_hmm('forward', self.sub_param_dir, n_procs=n_procs, partition=cpath.partitions[cpath.i_best_minus_x], shuffle_input=True)  # note that this annihilates the old <cpath>, which is a memory optimization (but we write all of them to the cpath progress dir)
            n_proc_list.append(n_procs)
            if self.are_we_finished_clustering(n_procs, cpath):
                break
            n_procs, cpath = self.prepare_next_iteration(n_proc_list, cpath, initial_nseqs)
            self.istep += 1

        if self.args.max_cluster_size is not None:
            print '   --max-cluster-size (partitiondriver): merging shared clusters'
            self.merge_shared_clusters(cpath)

        cpath = self.merge_cpaths_from_previous_steps(cpath)

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

        n_max_procs = 100

        if nseqs < 1000:
            seqs_per_proc = 250
        elif nseqs < 3000:
            seqs_per_proc = 500
        else:
            seqs_per_proc = 1000
        if self.args.batch_system is not None:  # if we're using a batch system, all the overhead (and priority issues) means it makes sense to have fewer processes
            seqs_per_proc *= 4
        n_precache_procs = int(math.ceil(float(nseqs) / seqs_per_proc))
        n_precache_procs = min(n_precache_procs, n_max_procs)  # I can't get more'n a few hundred slots at a time, so it isn't worth using too much more than that
        if self.args.batch_system is None:  # if we're not on a batch system, make sure it's less than the number of cpus
            n_precache_procs = min(n_precache_procs, multiprocessing.cpu_count())
        else:
            n_precache_procs = min(n_precache_procs, self.args.n_procs)  # aw, screw it, just limit it to --n-procs

        return n_precache_procs

    # ----------------------------------------------------------------------------------------
    def get_cluster_annotations(self, cpath, extra_dbg_str=''):  # we used to try to have glomerator.cc try to guess what annoations to write (using ::WriteAnnotations()), but now we go back and rerun a separate bcrham process just to get the annotations we want, partly for that control, but also partly because we typically want all these annotations calculated without any uid translation.
        def add_additional_clusters(clusters_to_annotate):
            cluster_set = set([tuple(c) for c in clusters_to_annotate])
            if self.args.write_additional_cluster_annotations is not None:
                istart = max(0, cpath.i_best - self.args.write_additional_cluster_annotations[0])
                istop = min(len(cpath.partitions), cpath.i_best + 1 + self.args.write_additional_cluster_annotations[1])
                for ip in range(istart, istop):
                    cluster_set |= set([tuple(c) for c in cpath.partitions[ip]])
            if self.args.calculate_alternative_annotations:  # add every cluster from the entire clustering history NOTE stuff that's in self.hmm_cachefname (which we used to use for this), but not in the cpath progress files: translated clusters, clusters that we calculated but didn't merge (maybe? not sure)
                cluster_set |= set([(uid,) for cluster in cpath.partitions[cpath.i_best] for uid in cluster])  # add the singletons separately, since we don't write a singleton partition before collapsing naive sequences before the first clustering step
                cluster_set |= set([tuple(cluster) for partition in cpath.partitions for cluster in partition])  # kind of wasteful to re-add clusters from the best partition here, but oh well
            print '    added %d clusters (in addition to the %d from the best partition) before running cluster annotations' % (len(cluster_set) - len(clusters_to_annotate), len(clusters_to_annotate))
            if self.args.debug:
                print '       %s these additional clusters will also be printed below, since --debug is greater than 0' % utils.color('yellow', 'note:')
            return [list(c) for c in cluster_set]
        def remove_additional_clusters(best_annotations):
            keys_to_remove = [uidstr for uidstr in best_annotations if uidstr.split(':') not in cpath.partitions[cpath.i_best]]
            for uidstr in keys_to_remove:
                del best_annotations[uidstr]
            if len(best_annotations) != len(cpath.partitions[cpath.i_best]):
                if len(best_annotations) < len(cpath.partitions[cpath.i_best]):  # if <best_annotations> is too short, it should be because there was a failed annotation
                    print '    %s read fewer cluster annotations than there are clusters in the best partition (should be accounted for above)' % utils.color('yellow', 'warning')
                else:
                    raise Exception('something went wrong when removing extra clusters from best_annotations (%d vs %d)' % (len(best_annotations), len(cpath.partitions[cpath.i_best])))

        # ----------------------------------------------------------------------------------------
        clusters_to_annotate = cpath.partitions[cpath.i_best]
        if self.args.write_additional_cluster_annotations is not None or self.args.calculate_alternative_annotations:
            clusters_to_annotate = add_additional_clusters(clusters_to_annotate)
            extra_dbg_str += ' (including additional clusters)'

        if len(clusters_to_annotate) == 0:
            return
        action_cache = self.current_action  # hackey, but probably not worth trying (more) to improve
        self.current_action = 'annotate'
        clusters_to_annotate = sorted(clusters_to_annotate, key=len, reverse=True)  # as opposed to in clusterpath, where we *don't* want to sort by size, it's nicer to have them sorted by size here, since then as you're scanning down a long list of cluster annotations you know once you get to the singletons you won't be missing something big
        n_procs = min(self.args.n_procs, len(clusters_to_annotate))  # we want as many procs as possible, since the large clusters can take a long time (depending on if we're translating...), but in general we treat <self.args.n_procs> as the maximum allowable number of processes
        print 'getting annotations for final partition%s' % extra_dbg_str
        self.run_hmm('viterbi', self.sub_param_dir, n_procs=n_procs, partition=clusters_to_annotate, read_output=False)
        if n_procs > 1:
            self.merge_all_hmm_outputs(n_procs, precache_all_naive_seqs=False)
        best_annotations, hmm_failures = self.read_annotation_output(self.hmm_outfname, count_parameters=self.args.count_parameters, parameter_out_dir=self.multi_hmm_param_dir if self.args.parameter_out_dir is None else self.args.parameter_out_dir)
        if self.args.get_tree_metrics:
            self.calculate_tree_metrics(best_annotations, cpath=cpath)  # adds tree metrics to <annotations>

        if self.args.calculate_alternative_annotations:
            clusters_in_partitions = set([':'.join(c) for partition in cpath.partitions for c in partition])
            for cluster in sorted(cpath.partitions[cpath.i_best], key=len, reverse=True):
                self.get_subcluster_naive_seqs(cluster, best_annotations, clusters_in_partitions, debug=self.args.debug)  # NOTE modifies the annotations (adds alternative annotation info)

        if self.args.outfname is not None:  # NOTE need to write *before* removing any clusters from the non-best partition
            self.write_output(best_annotations.values(), hmm_failures, cpath=cpath, dont_write_failed_queries=True)

        if self.args.write_additional_cluster_annotations is not None:  # remove the clusters that aren't actually in the best partition (we need them for partition plotting)
            remove_additional_clusters(best_annotations)

        if os.path.exists(self.hmm_infname):
            os.remove(self.hmm_infname)
        self.current_action = action_cache

        return best_annotations

    # ----------------------------------------------------------------------------------------
    def get_cached_hmm_naive_seqs(self, queries=None):
        # would be nice to merge this with self.read_hmm_cachefile()
        expected_queries = self.sw_info['queries'] if queries is None else queries
        cached_naive_seqs = {}
        with open(self.hmm_cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                if ':' in line['unique_ids']:  # if it's a cache file left over from a previous partitioning, there'll be clusters in it, too
                    continue
                if line['unique_ids'] not in expected_queries:  # probably can only happen if self.args.persistent_cachefname is set
                    continue
                cached_naive_seqs[line['unique_ids']] = line['naive_seq']
                if len(cached_naive_seqs) == len(expected_queries):  # already got everybody
                    break

        if set(cached_naive_seqs) != set(expected_queries):  # can happen if hmm can't find a path for a sequence for which sw *did* have an annotation (but in that case the annotation is almost certainly garbage)
            extra = set(cached_naive_seqs) - set(expected_queries)
            missing = set(expected_queries) - set(cached_naive_seqs)
            if len(extra) > 0:
                print '    %s read %d extra queries from hmm cache file %s' % (utils.color('yellow', 'warning:'), len(extra), ' '.join(extra))
            if len(missing) > 0:
                print '    %s missing %d queries from hmm cache file (using sw naive sequence instead): %s' % (utils.color('yellow', 'warning:'), len(missing), ' '.join(missing))
                for uid in missing:
                    cached_naive_seqs[uid] = self.sw_info[uid]['naive_seq']

        return cached_naive_seqs

    # ----------------------------------------------------------------------------------------
    def cluster_with_naive_vsearch_or_swarm(self, parameter_dir=None):
        start = time.time()

        naive_seq_list = []
        assert parameter_dir is not None
        threshold = self.get_naive_hamming_bounds(parameter_dir)[0]  # lo and hi are the same
        cached_naive_seqs = self.get_cached_hmm_naive_seqs()
        for uid in self.sw_info['queries']:
            if uid not in cached_naive_seqs:
                raise Exception('naive sequence for %s not found in %s' % (uid, self.hmm_cachefname))
            naive_seq_list.append((uid, cached_naive_seqs[uid]))

        all_naive_seqs, naive_seq_hashes = utils.collapse_naive_seqs_with_hashes(naive_seq_list, self.sw_info)

        print '    using hfrac bound for vsearch %.3f' % threshold

        partition = []
        print '    running vsearch %d times (once for each cdr3 length class):' % len(all_naive_seqs),
        for cdr3_length, sub_naive_seqs in all_naive_seqs.items():
            sub_hash_partition = utils.run_vsearch('cluster', sub_naive_seqs, self.args.workdir + '/vsearch', threshold, vsearch_binary=self.args.vsearch_binary)
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
    def get_naive_hamming_bounds(self, parameter_dir=None, overall_mute_freq=None):  # parameterize the relationship between mutation frequency and naive sequence inaccuracy
        if self.cached_naive_hamming_bounds is not None:  # only run the stuff below once
            return self.cached_naive_hamming_bounds

        if parameter_dir is not None:
            assert overall_mute_freq is None
            mutehist = Hist(fname=parameter_dir + '/all-mean-mute-freqs.csv')
            mute_freq = mutehist.get_mean(ignore_overflows=True)  # should I not ignore overflows here?
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
            y1, y2 = 0.015, 0.015  # would be nice to get better numbers for this
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
        if n_procs > 1:  # only cache vals for sequence sets with newly-calculated vals (initial cache file is copied to each subdir)
            cmd_str += ' --only-cache-new-vals'

        if self.args.dont_rescale_emissions:
            cmd_str += ' --dont-rescale-emissions'
        if self.current_action == 'partition':
            if self.args.cache_naive_hfracs:
                cmd_str += ' --cache-naive-hfracs'
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
                cmd_str += ' --n-partitions-to-write ' + str(self.args.n_partitions_to_write)  # don't write too many, since calculating the extra logprobs is kind of expensive (note that in practice bcrham rarely has the default of 10 to work with any more, mostly because it's only looking at merges from one clustering step. Also note that we typically won't have the logprobs for all of them, for the same reason)
                if n_procs == 1:  # if this is the last time through, with one process, we want glomerator.cc to calculate the total logprob of each partition NOTE this is quite expensive, since we have to turn off translation entirely
                    cmd_str += '  --write-logprob-for-each-partition'

                if self.args.seed_unique_id is not None and self.unseeded_seqs is None:  # if we're in the last few cycles (i.e. we've removed unseeded clusters) we want bcrham to not know about the seed (this gives more accurate clustering 'cause we're really doing hierarchical agglomeration)
                    cmd_str += ' --seed-unique-id ' + self.args.seed_unique_id

                if n_procs == 1:
                    if self.args.n_final_clusters is not None:
                        cmd_str += ' --n-final-clusters ' + str(self.args.n_final_clusters)
                    if self.args.min_largest_cluster_size is not None:
                        cmd_str += ' --min-largest-cluster-size ' + str(self.args.min_largest_cluster_size)

                if self.args.max_cluster_size is not None:
                    cmd_str += ' --max-cluster-size ' + str(self.args.max_cluster_size)

        assert len(utils.ambiguous_bases) == 1  # could allow more than one, but it's not implemented a.t.m.
        cmd_str += ' --ambig-base ' + utils.ambiguous_bases[0]

        return cmd_str

    # ----------------------------------------------------------------------------------------
    def subworkdir(self, iproc, n_procs):
        if n_procs == 1:  # NOTE it kind of sucks not using the <self.istep> stuff for n_procs=1 (and also not using it for either naive sequence precaching or cluster annotation), but there's too many different places where stuff would need to change and I don't want to fix it now.
            return self.args.workdir
        else:
            subworkdir = self.args.workdir
            if self.istep is not None:  # have to use a separate darn subdir for each iteration so the network filesystem doesn't screw everything up (same thing in waterer)
                subworkdir += '/istep-%d' % self.istep
                self.subworkdirs.append(subworkdir)
            return subworkdir + '/hmm-' + str(iproc)

    # ----------------------------------------------------------------------------------------
    def print_partition_dbgfo(self):
        if self.bcrham_proc_info is None:
            return
        actionstr = self.current_action if self.current_action != 'cache-parameters' else 'annotate'
        summaryfo = utils.summarize_bcrham_dbgstrs(self.bcrham_proc_info, action=actionstr)

        pwidth = str(len(str(len(self.input_info))))  # close enough
        if 'read-cache' in utils.bcrham_dbgstrs[actionstr]:
            if sum(summaryfo['read-cache'].values()) == 0:
                print '                no/empty cache file'
            else:
                print ('          read from cache:  naive-seqs %' + pwidth + 'd   logprobs %' + pwidth + 'd') % (summaryfo['read-cache']['naive-seqs'], summaryfo['read-cache']['logprobs'])
        print ('                    calcd:         vtb %' + pwidth + 'd        fwd %' + pwidth + 'd') % (summaryfo['calcd']['vtb'], summaryfo['calcd']['fwd'])
        if 'merged' in utils.bcrham_dbgstrs[actionstr]:
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
        def get_cmd_str(iproc):  # all this does at this point is replace workdir with sub-workdir in hmm input, output, and cache file arguments
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
        glutils.write_glfo(self.my_gldir, self.glfo)

        cmd_str = self.get_hmm_cmd_str(algorithm, self.hmm_infname, self.hmm_outfname, parameter_dir=parameter_in_dir, precache_all_naive_seqs=precache_all_naive_seqs, n_procs=n_procs)

        if n_procs > 1:
            self.split_input(n_procs, self.hmm_infname)

        exec_start = time.time()
        self.execute(cmd_str, n_procs)
        exec_time = time.time() - exec_start

        glutils.remove_glfo_files(self.my_gldir, self.args.locus)

        cpath, annotations, hmm_failures = None, None, None
        if read_output:
            if self.current_action == 'partition' or n_procs > 1:
                cpath = self.merge_all_hmm_outputs(n_procs, precache_all_naive_seqs)
                if cpath is not None:
                    cpath.write(self.get_cpath_progress_fname(self.istep), self.args.is_data, reco_info=self.reco_info, true_partition=utils.get_true_partition(self.reco_info) if not self.args.is_data else None)

            if algorithm == 'viterbi' and not precache_all_naive_seqs:
                annotations, hmm_failures = self.read_annotation_output(self.hmm_outfname, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir, print_annotations=self.args.debug)

            if os.path.exists(self.hmm_infname):
                os.remove(self.hmm_infname)

        step_time = time.time() - start
        if step_time - exec_time > 0.1:
            print '         infra time: %.1f' % (step_time - exec_time)  # i.e. time for non-executing, infrastructure time
        print '      hmm step time: %.1f' % step_time
        self.timing_info.append({'exec' : exec_time, 'total' : step_time})  # NOTE in general, includes pre-cache step

        return cpath, annotations, hmm_failures

    # ----------------------------------------------------------------------------------------
    def read_hmm_cachefile(self):
        # would be nice to merge this with self.get_cached_hmm_naive_seqs()
        cachefo = {}
        if not os.path.exists(self.hmm_cachefname):
            return cachefo
        with open(self.hmm_cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                utils.process_input_line(line)
                cachefo[':'.join(line['unique_ids'])] = line
        return cachefo

    # ----------------------------------------------------------------------------------------
    def get_subcluster_naive_seqs(self, uids_of_interest, cluster_annotations, clusters_in_partitions, cachefo=None, debug=False):
        # <clusters_in_partitions> is a list of the clusters that actually occur in one of the partitions in the clusterpath. It's so we can ignore annotations for clusters that were never actually formed, i.e. where we got the annotation for a potential cluster, but decided not to make the merge (at least, I think this is why some annotations don't correspond to clusters in the cluster path)
        # NOTE that when seed partitioning, in the steps before we throw out non-seeded clusters, there's *tons* of very overlapping clusters, since I think we take the the biggest seeded cluster, and pass that to all the subprocs, so then if a different sequence gets added to it in each subproc you end up with lots of different versions (including potentially lots of different orderings of the exact same cluster)

        if cachefo is None:  # <cachefo> and <cluster_annotations> are only different (i.e. cachefo is passed in) if we're reading old-style (csv) output with a separate cache file
            cachefo = cluster_annotations

        if self.args.seed_unique_id is not None and debug and self.args.seed_unique_id not in uids_of_interest:
            print '  note: only printing alternative annotations for seed clusters (the rest are still written to disk, but if you want to print the others, don\'t set --seed-unique-id)'
            debug = False

        uids_of_interest = set(uids_of_interest)
        uidstr_of_interest = None  # we don't yet know in what order the uids in <uids_of_interest> appear in the file
        sub_info = {}  # map from naive seq : sub uid strs
        final_info = {'naive-seqs' : OrderedDict(), 'gene-calls' : {r : OrderedDict() for r in utils.regions}}  # this is the only info that's persistent, it get's added to the cluster annotation corresponding to <uids_of_interest>

        # pull all the clusters out of the cache info that have any overlap with <uids_of_interest>
        for uidstr, info in cachefo.items():
            if info['naive_seq'] == '':  # hmm cache file lines that only have logprobs
                continue
            uid_set = set(info['unique_ids'])
            if len(uid_set & uids_of_interest) == 0:
                continue
            if ':' in uidstr and uidstr not in clusters_in_partitions:  # first bit is because singletons are not in general all in a partition, since the first partition is from after initial naive seq merging (at least seems to be? I kind of thought I added the singleton partition, but I guess not)
                continue
            naive_seq = cachefo[uidstr]['naive_seq']
            if naive_seq not in sub_info:
                sub_info[naive_seq] = []
            sub_info[naive_seq].append(uidstr)
            if uid_set == uids_of_interest:
                uidstr_of_interest = uidstr

        if uidstr_of_interest is None:
            raise Exception('hmm cache file doesn\'t have a cluster with the exact requested uids (run without setting --queries to get a list of available clusters): %s' % uids_of_interest)

        cache_file_naive_seq = cachefo[uidstr_of_interest]['naive_seq']  # ok not actually from the hmm cache file with new-style output, but I don't feel like renaming it

        if uidstr_of_interest in cluster_annotations:
            line_of_interest = cluster_annotations[uidstr_of_interest]
        else:  # only possible if we're reading the hmm cache file, i.e. old-style output
            n_max_overlap, max_cluster = None, None  # see if it's in there in a different order (not sure how, but this happened once)
            for clusterstr in cluster_annotations:
                cluster = clusterstr.split(':')
                overlap = set(cluster) & uids_of_interest
                if n_max_overlap is None or len(overlap) > n_max_overlap:
                    n_max_overlap = len(overlap)
                    max_cluster = cluster
            if n_max_overlap == len(uids_of_interest):
                print '  %s uids of interest in cluster annotations are in a different order to those in cache file (not really sure how this would happen, but it should only be possible on deprecated csv files, so oh well)' % utils.color('yellow', 'warning')
                line_of_interest = cluster_annotations[':'.join(max_cluster)]
            else:
                raise Exception('annotations don\'t have cluster corresponding to that found in hmm cache file')

        genes_of_interest = {r : line_of_interest[r + '_gene'] for r in utils.regions}
        gene_strs_of_interest = {r : utils.color_gene(genes_of_interest[r]) for r in utils.regions}  # have to do this up here so we know the width before we start looping
        gswidth = str(utils.len_excluding_colors('  '.join(gene_strs_of_interest.values())))  # out of order, but doesn't matter since we only want the length

        naive_seq_of_interest = line_of_interest['naive_seq']
        if naive_seq_of_interest != cache_file_naive_seq:  # not possible with new-style output, since we don't read the hmm cache file
            print '%s naive sequences from cluster annotation and cache file aren\'t the same:' % utils.color('yellow', 'warning')
            utils.color_mutants(naive_seq_of_interest, cache_file_naive_seq, print_result=True, ref_label='cluster annotation  ', seq_label='cache file  ', extra_str='     ')
            print ''

        if debug:
            print '  subcluster naive sequences for:\n    %s\n   (in %s below)' % (uidstr_of_interest, utils.color('blue', 'blue'))
            print ''
            utils.print_reco_event(utils.synthesize_single_seq_line(line_of_interest, iseq=0), extra_str='      ', label='annotation for a single (arbitrary) sequence from the cluster:')
            print ''
            print ''
            print '%s          unique       inconsistent (%s: missing info)     cluster' % (' ' * len(cache_file_naive_seq), utils.color('blue', 'x'))  # 'missing info' means that we didn\'t have an annotation for at least one of the clusters in that line. This is probably because we're using an old hmm cache file, and at best passed in the sw cache file annotations, which are only single-sequence.
            print '%s        seqs  frac     regions         genes               sizes (+singletons)' % (' ' * len(cache_file_naive_seq))

            # ----------------------------------------------------------------------------------------
            def print_naive_line(other_genes, uid_str_list, naive_seq, n_independent_seqs, max_len_other_gene_str=20):
                gene_strs, other_gene_strs = [], []
                for region in utils.regions:
                    if len(other_genes[region]) > 0:
                        gene_strs += [utils.color('red', '%s' % region)]  # , width=utils.len_excluding_colors(gene_strs_of_interest[region]))]
                        for gene in other_genes[region]:
                            other_gene_strs += [utils.color_gene(gene, width=utils.len_excluding_colors(gene_strs_of_interest[region]))]
                    else:
                        gene_strs += [' ']  #  * utils.len_excluding_colors(gene_strs_of_interest[region])]
                gene_str = '  '.join(gene_strs)  # ('%' + gswidth + 's') % '  '.join(gene_strs)
                other_gene_str = ' '.join(other_gene_strs)  # have to split this out into a spearate string since sometimes a single line will have, e.g., three separate inconsistent d genes
                if utils.len_excluding_colors(other_gene_str) > max_len_other_gene_str:
                    max_len_other_gene_str = utils.len_excluding_colors(other_gene_str)
                other_gene_str += ' ' * (max_len_other_gene_str - utils.len_excluding_colors(other_gene_str)) # - 3)  # this will be too narrow until we get to whichever line happens to have the max number of other genes, but, oh, well (and no, I don't know why I need the -3)
                if no_info:  # at least one cluster didn't have an annotation
                    other_gene_str = '%s %s' % (utils.color('blue', 'x'), other_gene_str)

                cluster_size_strs = []
                n_singletons = 0
                for uidstr in uid_str_list:
                    size_str = str(uidstr.count(':') + 1)
                    if size_str == '1':
                        n_singletons += 1
                        continue
                    if uidstr_of_interest == uidstr:
                        size_str = utils.color('blue', size_str)
                    cluster_size_strs.append(size_str)
                if n_singletons > 0:
                    cluster_size_strs.append('(+%d)' % n_singletons)

                pre_str, post_str = '', ''
                if uidstr_of_interest in uid_str_list:
                    pre_str = utils.color('blue', '-->', width=5)
                    post_str = utils.color('blue', ' <-- requested uids', width=5)
                if self.args.print_all_annotations:
                    print '  printing annotations for all clusters supporting naive seq with fraction %.3f:' % final_info['naive-seqs'][naive_seq]
                print ('%5s %s  %4d  %4.2f     %s     %s    %s%s') % (pre_str, utils.color_mutants(cache_file_naive_seq, naive_seq), n_independent_seqs, final_info['naive-seqs'][naive_seq], gene_str, other_gene_str, ' '.join(cluster_size_strs), post_str)
                if self.args.print_all_annotations:
                    for uidstr in sorted(uid_str_list, key=lambda x: x.count(':'), reverse=True):
                        if uidstr in cluster_annotations:
                            utils.print_reco_event(cluster_annotations[uidstr], extra_str='      ', label='')
                        else:
                            print '    %s missing from cluster annotations' % uidstr
                    print '\n'

        independent_seq_info = {naive_seq : set([uid for uidstr in uid_str_list for uid in uidstr.split(':')]) for naive_seq, uid_str_list in sub_info.items()}
        total_independent_seqs = sum(len(uids) for uids in independent_seq_info.values())
        unique_seqs_for_each_gene = {r : {} for r in utils.regions}  # for each gene, keeps track of the number of unique sequences that contributed to an annotation that used that gene
        cluster_sizes_for_each_gene = {r : {} for r in utils.regions}  # same, but number/sizes of different clusters (so we can tell if a gene is only supported by like one really large cluster, but all the smaller clusters point to another gene)
        uidstrs_for_each_gene = {r : {} for r in utils.regions}  # need the actual uid strs, but only use them for --print-all-annotations
        def init_dicts(region, gene_call):  # ick
            unique_seqs_for_each_gene[region][gene_call] = set()
            cluster_sizes_for_each_gene[region][gene_call] = []
            uidstrs_for_each_gene[region][gene_call] = []
        for naive_seq in sorted(independent_seq_info, key=lambda ns: len(independent_seq_info[ns]), reverse=True):
            final_info['naive-seqs'][naive_seq] = len(independent_seq_info[naive_seq]) / float(total_independent_seqs)

            uid_str_list = sorted(sub_info[naive_seq], key=lambda uidstr: uidstr.count(':') + 1, reverse=True)

            # print the v gene along side the first naive sequence, as well as for any subsequent ones that have a different v
            other_genes = {r : set() for r in utils.regions}
            no_info = False  # gets set to true if *any* of the <uidstr>s are missing info (so <no_info> can be set to true for a line for which we also have inconsistent genes set)
            for uidstr in uid_str_list:
                for region in utils.regions:
                    if uidstr in cluster_annotations:  # we should have all of 'em now, at least for new-style output
                        gene_call = cluster_annotations[uidstr][region + '_gene']
                        if gene_call != genes_of_interest[region]:
                            other_genes[region].add(gene_call)
                        if gene_call not in unique_seqs_for_each_gene[region]:
                            init_dicts(region, gene_call)
                        unique_seqs_for_each_gene[region][gene_call] |= set(cluster_annotations[uidstr]['unique_ids'])
                        cluster_sizes_for_each_gene[region][gene_call].append(len(cluster_annotations[uidstr]['unique_ids']))
                        uidstrs_for_each_gene[region][gene_call].append(uidstr)
                    else:
                        no_info = True
            if debug:
                print_naive_line(other_genes, uid_str_list, naive_seq, len(independent_seq_info[naive_seq]))

        for region in utils.regions:
            total_unique_seqs_this_region = sum(len(uids) for uids in unique_seqs_for_each_gene[region].values())  # yes, it is in general different for each region
            for gene_call, uids in sorted(unique_seqs_for_each_gene[region].items(), key=lambda s: len(s[1]), reverse=True):
                final_info['gene-calls'][region][gene_call] = len(uids) / float(total_unique_seqs_this_region)
        if debug:
            print ''
            print '                        unique'
            print '                      seqs  frac          cluster sizes (+singletons)'
            for region in utils.regions:
                print '    %s' % utils.color('green', region)
                for gene_call, uids in sorted(unique_seqs_for_each_gene[region].items(), key=lambda s: len(s[1]), reverse=True):
                    csizes = cluster_sizes_for_each_gene[region][gene_call]
                    if self.args.print_all_annotations:
                        print '  printing annotations for all clusters supprting gene call with fraction %.3f:' % final_info['gene-calls'][region][gene_call]
                    cluster_size_strs = ['%d' % cs for cs in sorted(csizes, reverse=True) if cs > 1]
                    n_singletons = csizes.count(1)
                    if n_singletons > 0:
                        cluster_size_strs.append('(+%d)' % n_singletons)
                    print '      %s %4d  %4.2f             %s' % (utils.color_gene(gene_call, width=15), len(uids), final_info['gene-calls'][region][gene_call], ' '.join(cluster_size_strs))
                    if self.args.print_all_annotations:
                        for uidstr in sorted(uidstrs_for_each_gene[region][gene_call], key=lambda x: x.count(':'), reverse=True):
                            if uidstr in cluster_annotations:
                                utils.print_reco_event(cluster_annotations[uidstr], extra_str='      ', label='')
                            else:
                                print '    %s missing from cluster annotations' % uidstr
                        print '\n'

        # this is messy because we want to acces them as dicts in this fcn, for dbg printing, but want them in the output file as simple combinations of lists/tuples
        line_of_interest['alternative-annotations'] = {'naive-seqs' : final_info['naive-seqs'].items(), 'gene-calls' : {r : final_info['gene-calls'][r].items() for r in utils.regions}}

    # ----------------------------------------------------------------------------------------
    def get_padded_true_naive_seq(self, qry):
        assert len(self.sw_info[qry]['padlefts']) == 1
        return self.sw_info[qry]['padlefts'][0] * utils.ambiguous_bases[0] + self.reco_info[qry]['naive_seq'] + self.sw_info[qry]['padrights'][0] * utils.ambiguous_bases[0]

    # ----------------------------------------------------------------------------------------
    def split_input(self, n_procs, infname):

        # should we pull out the seeded clusters, and carefully re-inject them into each process?
        separate_seeded_clusters = self.current_action == 'partition' and self.args.seed_unique_id is not None and self.unseeded_seqs is None  # I think it's no longer possible to have seed_unique_id set if we're not partitioning, but I'll leave it just to be safe (otherwise we get the seed seq sent to every process)

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
                cpath = glomerer.read_cached_agglomeration(infnames, debug=self.args.debug)  #, outfname=self.hmm_outfname)
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
        if self.args.dont_write_parameters:
            return
        print '  writing hmms',
        sys.stdout.flush()
        start = time.time()

        from hmmwriter import HmmWriter
        hmm_dir = parameter_dir + '/hmms'
        utils.prep_dir(hmm_dir, '*.yaml')
        glutils.restrict_to_observed_genes(self.glfo, parameter_dir)  # this is kind of a weird place to put this... it would make more sense to read the glfo from the parameter dir, but I don't want to mess around with changing that a.t.m.

        if self.args.debug:
            print 'to %s' % parameter_dir + '/hmms',

        if multiprocessing.cpu_count() * utils.memory_usage_fraction() > 0.8:  # already using a lot of memory, so don't to call multiprocessing, which will duplicate all the memory for each process
            for region in utils.regions:
                for gene in self.glfo['seqs'][region]:
                    writer = HmmWriter(parameter_dir, hmm_dir, gene, self.glfo, self.args)
                    writer.write()
        else:
            def write_single_hmm(gene):
                writer = HmmWriter(parameter_dir, hmm_dir, gene, self.glfo, self.args)
                writer.write()
            procs = [multiprocessing.Process(target=write_single_hmm, args=(gene,))
                     for region in utils.regions for gene in self.glfo['seqs'][region]]
            utils.run_proc_functions(procs)  # uses all the cores (should only be for a little bit, though)

        print '(%.1f sec)' % (time.time()-start)
        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def get_existing_hmm_files(self, parameter_dir):
        fnames = [os.path.basename(fn) for fn in glob.glob(parameter_dir + '/hmms/*.yaml')]
        genes = set([utils.unsanitize_name(utils.getprefix(fn)) for fn in fnames])
        if len(genes) == 0:
            raise Exception('no yamels in %s' % parameter_dir + '/hmms')
        return genes

    # ----------------------------------------------------------------------------------------
    def combine_queries(self, query_names, available_genes):
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
            raise Exception('cdr3 lengths not all the same for %s (%s)' % (' '.join(uids), ' '.join([str(c) for c in lengths])))
        combo['cdr3_length'] = cdr3_lengths[0]

        combo['k_v'] = {'min' : 99999, 'max' : -1}
        combo['k_d'] = {'min' : 99999, 'max' : -1}
        combo['only_genes'] = set()
        for name in query_names:
            swfo = self.sw_info[name]
            k_v = swfo['k_v']
            k_d = swfo['k_d']
            assert len(swfo['seqs']) == 1  # checking that when we filled in 'seqs' all was well
            combo['k_v']['min'] = min(k_v['min'], combo['k_v']['min'])
            combo['k_v']['max'] = max(k_v['max'], combo['k_v']['max'])
            combo['k_d']['min'] = min(k_d['min'], combo['k_d']['min'])
            combo['k_d']['max'] = max(k_d['max'], combo['k_d']['max'])

            genes_to_use = set()  # genes from this query that'll get ORd into the ones from the previous queries
            for region in utils.regions:
                regmatches = set(swfo['all_matches'][region])  # the best <n_max_per_region> matches for this query, ordered by sw match quality (well, ordered before we call set())
                genes_to_use |= regmatches & available_genes

            # OR this query's genes into the ones from previous queries
            combo['only_genes'] |= genes_to_use  # NOTE using the OR of all sets of genes (from all query seqs) like this *really* helps,

        combo['only_genes'] = list(combo['only_genes'])  # maybe I don't need to convert it to a list, but it used to be a list, and I don't want to make sure there wasn't a reason for that

        # if we don't have at least one gene for each region, add all available genes from that regions
        if set(utils.regions) > set(utils.get_region(g) for g in combo['only_genes']):  # this is *very* rare
            missing_regions = set(utils.regions) - set(utils.get_region(g) for g in combo['only_genes'])  # this is wasteful, but it hardly *ever* happens
            combo['only_genes'] += [g for g in available_genes if utils.get_region(g) in missing_regions]
            if self.args.debug:
                print '    %s no %s genes for query %s, so added in all the genes from these regions' % ('note:', ' or '.join(missing_regions), query_names)  # these other genes may not work (the kbounds are probably wrong), but they might work, and it's better than just giving up on the query

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
    def write_to_single_input_file(self, fname, nsets, parameter_dir, shuffle_input=False):
        csvfile = open(fname, 'w')
        header = ['names', 'k_v_min', 'k_v_max', 'k_d_min', 'k_d_max', 'mut_freq', 'cdr3_length', 'only_genes', 'seqs']
        writer = csv.DictWriter(csvfile, header, delimiter=' ')
        writer.writeheader()

        if shuffle_input:  # shuffle nset order (this is absolutely critical when clustering with more than one process, in order to redistribute sequences among the several processes)
            random.shuffle(nsets)

        if self.current_action == 'partition' and self.args.synthetic_distance_based_partition:
            self.write_fake_cache_file(nsets)

        genes_with_hmm_files = self.get_existing_hmm_files(parameter_dir)
        genes_with_enough_counts = utils.get_genes_with_enough_counts(parameter_dir, self.args.min_allele_prevalence_fractions)  # it would be nice to do this at some earlier step, but then we have to rerun sw (note that there usually won't be any to remove for V, since the same threshold was already applied in alleleremover, but there we can't do d and j since we don't yet have annotations for them)
        glfo_genes = set([g for r in utils.regions for g in self.glfo['seqs'][r]])
        if self.args.only_genes is None and len(genes_with_hmm_files - glfo_genes) > 0:
            print '  %s hmm files for genes that aren\'t in glfo: %s' % (utils.color('red', 'warning'), utils.color_genes(genes_with_hmm_files - glfo_genes))
        if len(glfo_genes - genes_with_hmm_files) > 0:
            print '    skipping matches from %d genes that don\'t have hmm files: %s' % (len(glfo_genes - genes_with_hmm_files), utils.color_genes(glfo_genes - genes_with_hmm_files))
        if self.current_action == 'cache-parameters' and len(glfo_genes - genes_with_enough_counts) > 0:
            print '    skipping matches from %d genes without enough counts: %s' % (len(glfo_genes - genes_with_enough_counts), utils.color_genes(glfo_genes - genes_with_enough_counts))
        available_genes = genes_with_hmm_files & genes_with_enough_counts

        for query_name_list in nsets:  # NOTE in principle I think I should remove duplicate singleton <seed_unique_id>s here. But I think they in effect get removed 'cause in bcrham everything's stored as hash maps, so any duplicates just overwites the original upon reading its input
            combined_query = self.combine_queries(query_name_list, available_genes)
            if len(combined_query) == 0:  # didn't find all regions
                continue
            writer.writerow({
                'names' : ':'.join([qn for qn in query_name_list]),
                'k_v_min' : combined_query['k_v']['min'],
                'k_v_max' : combined_query['k_v']['max'],
                'k_d_min' : combined_query['k_d']['min'],
                'k_d_max' : combined_query['k_d']['max'],
                'mut_freq' : combined_query['mut_freq'],
                'cdr3_length' : combined_query['cdr3_length'],
                'only_genes' : ':'.join(combined_query['only_genes']),
                'seqs' : ':'.join(combined_query['seqs'])
            })

        csvfile.close()

    # ----------------------------------------------------------------------------------------
    @utils.timeprinter
    def prepare_for_hmm(self, algorithm, parameter_dir, partition, shuffle_input=False):
        """ Write input file for bcrham """

        if partition is not None:
            nsets = copy.deepcopy(partition)  # needs to be a deep copy so we can shuffle the order
        else:
            qlist = self.sw_info['queries']  # shorthand

            if self.args.simultaneous_true_clonal_seqs:
                assert self.args.n_simultaneous_seqs is None and not self.args.is_data  # are both already checked in ./bin/partis
                nsets = utils.get_true_partition(self.reco_info, ids=qlist)
                nsets = utils.split_clusters_by_cdr3(nsets, self.sw_info, warn=True)  # arg, have to split some clusters apart by cdr3, for rare cases where we call an shm indel in j within the cdr3
            elif self.args.all_seqs_simultaneous:  # everybody together
                nsets = [qlist]
                nsets = utils.split_clusters_by_cdr3(nsets, self.sw_info, warn=True)  # arg, have to split some clusters apart by cdr3, for rare cases where we call an shm indel in j within the cdr3
            elif self.args.n_simultaneous_seqs is not None:  # set number of simultaneous seqs
                nlen = self.args.n_simultaneous_seqs  # shorthand
                # nsets = [qlist[iq : min(iq + nlen, len(qlist))] for iq in range(0, len(qlist), nlen)]  # this way works fine, but it's hard to get right 'cause it's hard to understand 
                nsets = [list(group) for _, group in itertools.groupby(qlist, key=lambda q: qlist.index(q) / nlen)]  # integer division
            else:  # plain ol' singletons
                nsets = [[q] for q in qlist]

        self.write_to_single_input_file(self.hmm_infname, nsets, parameter_dir, shuffle_input=shuffle_input)  # single file gets split up later if we've got more than one process

    # ----------------------------------------------------------------------------------------
    def check_did_bcrham_fail(self, line, errorfo):
        if line['errors'] == '':  # no problems
            return False

        if 'v_gene' in line and line['v_gene'] == '':  # utils.process_input_line() just returns without doing anything in this case
            line['unique_ids'] = line['unique_ids'].split(':')

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
    def correct_multi_hmm_boundaries(self, line, debug=False):
        factor = 2.
        assert 'regional_bounds' not in line  # need to make sure implicit info isn't in there

        for boundary in utils.boundaries:
            delname = boundary[1] + '_5p_del'  # could just as well use the insertion length, but this is at least a reminder that the insertion is part of the righthand region
            sw_deletion_lengths = [self.sw_info[q][delname] for q in line['unique_ids']]
            median_single_length = numpy.median(sw_deletion_lengths)
            single_length_err = numpy.std(sw_deletion_lengths)

            if line[delname] < median_single_length + factor * single_length_err:  # skip it if the multi-hmm deletion isn't much longer than the single sw deletions
                continue

            if debug:
                true_del = None
                if self.reco_info is not None:
                    true_line = utils.synthesize_multi_seq_line_from_reco_info(line['unique_ids'], self.reco_info)
                    true_del = true_line[delname]
                print boundary
                print '  multi %d%s' % (line[delname], ('   true %d' % true_del) if true_del is not None else '')
                print '  %.2f +/- %.3f  (%s)' % (median_single_length, single_length_err, sw_deletion_lengths)

            while line[delname] > median_single_length and len(line[boundary + '_insertion']) > 0:
                line[delname] -= 1
                line[boundary + '_insertion'] = line[boundary + '_insertion'][:-1]

    # ----------------------------------------------------------------------------------------
    def process_dummy_d_hack(self, line, debug=False):
        """
        a.t.m. we force bcrham to give us D of length one for loci with no D genes.
        Here, we delete the dummy D base, and give it to either V, J, or the insertion.
        """
        tmpline = copy.deepcopy(line)
        utils.add_implicit_info(self.glfo, tmpline, reset_indel_genes=True)
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
        utils.add_implicit_info(self.glfo, after_line, reset_indel_genes=True)
        if debug:
            utils.print_reco_event(after_line, extra_str='    ', label='after')

    # ----------------------------------------------------------------------------------------
    def check_for_unexpectedly_missing_keys(self, annotations, hmm_failure_ids):
        missing_input_keys = set(self.input_info)
        missing_input_keys -= set([uid for line in  annotations.values() for uid in line['unique_ids']])  # set(self.sw_info['queries'])  # all the queries for which we had decent sw annotations (sw failures are accounted for below)
        missing_input_keys -= self.sw_info['failed-queries']
        missing_input_keys -= self.sw_info['removed-queries']
        missing_input_keys -= set([d for dlist in self.sw_info['duplicates'].values() for d in dlist])
        missing_input_keys -= hmm_failure_ids
        if self.unseeded_seqs is not None:
            missing_input_keys -= set(self.unseeded_seqs)
        if self.small_cluster_seqs is not None:
            missing_input_keys -= set(self.small_cluster_seqs)
        if self.unseeded_seqs is not None:
            missing_input_keys -= set(self.unseeded_seqs)
        if len(missing_input_keys) > 0:
            print '  %s couldn\'t account for %d missing input uid%s%s' % (utils.color('red', 'warning'), len(missing_input_keys), utils.plural(len(missing_input_keys)), ': %s' % ' '.join(missing_input_keys) if len(missing_input_keys) < 15 else '')

    # ----------------------------------------------------------------------------------------
    def read_annotation_output(self, annotation_fname, count_parameters=False, parameter_out_dir=None, print_annotations=False):
        """ Read bcrham annotation output """
        print '    read output'
        sys.stdout.flush()

        pcounter = ParameterCounter(self.glfo, self.args) if count_parameters else None
        true_pcounter = ParameterCounter(self.simglfo, self.args) if (count_parameters and not self.args.is_data) else None
        perfplotter = PerformancePlotter('hmm') if self.args.plot_annotation_performance else None

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

                padded_line['indelfos'] = [self.sw_info['indels'].get(uid, indelutils.get_empty_indel()) for uid in uids]  # reminder: hmm was given a sequence with any indels reversed (i.e. <self.sw_info['indels'][uid]['reverersed_seq']>)
                padded_line['input_seqs'] = [self.sw_info[uid]['input_seqs'][0] for uid in uids]  # not in <padded_line>, since the hmm doesn't know anything about the input (i.e. non-indel-reversed) sequences
                padded_line['duplicates'] = [self.duplicates.get(uid, []) for uid in uids]
                for lkey in [lk for lk in utils.input_metafile_keys.values() if lk in self.sw_info[uids[0]]]:  # if it's in one, it should be in all of them
                    padded_line[lkey] = [self.sw_info[uid][lkey][0] for uid in uids]

                if not utils.has_d_gene(self.args.locus):
                    self.process_dummy_d_hack(padded_line)

                # if self.args.correct_boundaries and len(padded_line['unique_ids']) > 1:  # this does a decent job of correcting the multi-hmm's tendency to overestimate insertion and deletion lengths, but it also removes a significant portion of the multi-hmm's advantage in naive hamming distance
                #     self.correct_multi_hmm_boundaries(padded_line)

                try:
                    utils.add_implicit_info(self.glfo, padded_line, aligned_gl_seqs=self.aligned_gl_seqs, reset_indel_genes=True)
                except:  # I really don't like just swallowing it, but it's crashing deep in the new[ish] indel code on an extraordinarily rare and I think super screwed up sequence, and I can't replicate it without running on the entire stupid huge sample
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
                    print '      %s implicit info adding failed for %s when reading hmm output (so adding to failed queries):' % (utils.color('red', 'warning'), uidstr)
                    print utils.pad_lines(''.join(lines))
                    hmm_failures |= set(padded_line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                    continue

                utils.process_per_gene_support(padded_line)  # switch per-gene support from log space to normalized probabilities

                if self.args.linearham and self.current_action == 'annotate':
                    # add flexbounds/relpos to padded line
                    status = utils.add_linearham_info(self.sw_info, padded_line)
                    if status == 'nonsense':
                        hmm_failures |= set(padded_line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                        continue

                if padded_line['invalid']:
                    n_invalid_events += 1
                    if self.args.debug:
                        print '      %s padded line invalid' % uidstr
                        utils.print_reco_event(padded_line, extra_str='    ', label='invalid:')
                    hmm_failures |= set(padded_line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                    continue

                if uidstr in padded_annotations:  # this shouldn't happen, but it's more an indicator that something else has gone wrong than that in and of itself it's catastrophic
                    print '%s uidstr %s already read from file %s' % (utils.color('yellow', 'warning'), uidstr, annotation_fname)
                padded_annotations[uidstr] = padded_line

                if len(uids) > 1 or (self.args.linearham and self.current_action == 'annotate'):  # if there's more than one sequence (or we're in linearham mode), we need to use the padded line
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

                n_events_processed += 1
                n_seqs_processed += len(uids)

                if pcounter is not None:
                    pcounter.increment(line_to_use)

                if perfplotter is not None:
                    messed_up = False  # this would be really nice to clean up
                    for iseq in range(len(line_to_use['unique_ids'])):
                        if indelutils.has_indels(self.reco_info[uids[iseq]]['indelfos'][0]) or indelutils.has_indels(line_to_use['indelfos'][iseq]):
                            simlen = indelutils.net_length(self.reco_info[uids[iseq]]['indelfos'][0])
                            inflen = indelutils.net_length(line_to_use['indelfos'][iseq])
                            if simlen != inflen:  # see similar code in performanceplotter.py
                                messed_up = True
                                break
                    if messed_up:
                        continue
                    for iseq in range(len(uids)):  # NOTE this counts rearrangement-level parameters once for every mature sequence, which is inconsistent with the pcounters... but I think might make more sense here?
                        perfplotter.evaluate(self.reco_info[uids[iseq]], utils.synthesize_single_seq_line(line_to_use, iseq), simglfo=self.simglfo)

        if true_pcounter is not None:
            for uids in utils.get_true_partition(self.reco_info, ids=self.sw_info['queries']):  # NOTE this'll include queries that passed sw but failed the hmm... there aren't usually really any of those
                true_pcounter.increment(utils.synthesize_multi_seq_line_from_reco_info(uids, self.reco_info))

        # parameter and performance writing/plotting
        if pcounter is not None:
            path_str = 'hmm' if parameter_out_dir is None else os.path.basename(parameter_out_dir)  # at the moment, this is just 'hmm' for regular parameter caching, and 'multi-hmm' for parameter caching after partitioning
            if self.args.plotdir is not None:
                pcounter.plot('%s/%s' % (self.args.plotdir, path_str), only_csv=self.args.only_csv_plots, only_overall=not self.args.make_per_gene_plots)
                if true_pcounter is not None:
                    true_pcounter.plot('%s/%s' % (self.args.plotdir, path_str.replace('hmm', 'true')), only_csv=self.args.only_csv_plots, only_overall=not self.args.make_per_gene_plots)
            if not self.args.dont_write_parameters:
                pcounter.write(parameter_out_dir)
                if true_pcounter is not None:
                    true_pcounter.write('%s/%s' % (os.path.dirname(parameter_out_dir), path_str.replace('hmm', 'true')))

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
        seqfileopener.add_input_metafo(self.input_info, annotations_to_use.values())
        if print_annotations:
            self.print_results(None, annotations_to_use)

        self.check_for_unexpectedly_missing_keys(annotations_to_use, hmm_failures)  # NOTE not sure if it's really correct to use <annotations_to_use>, [maybe since <hmm_failures> has ones that failed the conversion to eroded line (and maybe other reasons)]

        # annotation (VJ CDR3) clustering
        if self.args.annotation_clustering is not None:
            self.deal_with_annotation_clustering(annotations_to_use, outfname)

        os.remove(annotation_fname)
        return annotations_to_use, hmm_failures

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
    def write_output(self, annotation_list, hmm_failures, cpath=None, dont_write_failed_queries=False, write_sw=False):
        if write_sw:
            assert annotation_list is None
            annotation_list = [self.sw_info[q] for q in self.input_info if q in self.sw_info['queries']]

        failed_queries = None
        if not dont_write_failed_queries:  # write empty lines for seqs that failed either in sw or the hmm
            failed_queries = [{'unique_ids' : [uid], 'invalid' : True, 'input_seqs' : self.input_info[uid]['seqs']} for uid in self.sw_info['failed-queries'] | hmm_failures]  # <uid> *needs* to be single-sequence (but there shouldn't really be any way for it to not be)

        if self.args.presto_output:
            presto_annotation_fname = self.args.outfname
            if cpath is not None:  # note that if we're partitioning, we get here with <self.current_action> set to 'annotate'
                presto_annotation_fname = utils.getprefix(self.args.outfname) + '.tsv'
                cpath.write_presto_partitions(self.args.outfname, self.input_info)
            utils.write_presto_annotations(presto_annotation_fname, annotation_list, failed_queries=failed_queries)
            return
        elif self.args.airr_output:
            utils.write_airr_output(self.args.outfname, annotation_list, cpath=cpath, failed_queries=failed_queries)
            return

        partition_lines = None
        if cpath is not None:
            true_partition = utils.get_true_partition(self.reco_info) if not self.args.is_data else None
            partition_lines = cpath.get_partition_lines(self.args.is_data, reco_info=self.reco_info, true_partition=true_partition, n_to_write=self.args.n_partitions_to_write, calc_missing_values=('all' if (len(annotation_list) < 500) else 'best'))

        headers = utils.add_lists(utils.annotation_headers if not write_sw else utils.sw_cache_headers, self.args.extra_annotation_columns)
        if self.args.linearham and self.current_action == 'annotate':
            headers = utils.add_lists(headers, utils.linearham_headers)
            utils.write_linearham_seqs(self.args.outfname, annotation_list)
        if utils.getsuffix(self.args.outfname) == '.csv':
            if cpath is not None:
                cpath.write(self.args.outfname, self.args.is_data, partition_lines=partition_lines)  # don't need to pass in reco_info/true_partition since we passed them when we got the partition lines
            annotation_fname = self.args.outfname if cpath is None else self.args.cluster_annotation_fname
            utils.write_annotations(annotation_fname, self.glfo, annotation_list, headers, failed_queries=failed_queries)
        elif utils.getsuffix(self.args.outfname) == '.yaml':
            utils.write_annotations(self.args.outfname, self.glfo, annotation_list, headers, failed_queries=failed_queries, partition_lines=partition_lines, use_pyyaml=self.args.write_full_yaml_output)
        else:
            raise Exception('unhandled annotation file suffix %s' % self.args.outfname)
