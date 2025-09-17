from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import json
import numpy
import time
import sys
import itertools
import math
import os
import glob
import csv
from io import open
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import random
from collections import OrderedDict
from subprocess import check_call, CalledProcessError
import copy
import multiprocessing
import operator
import traceback

from . import utils
from . import glutils
from . import indelutils
from . import treeutils
from . import lbplotting
from .glomerator import Glomerator
from .clusterpath import ClusterPath, ptnprint
from .waterer import Waterer
from .parametercounter import ParameterCounter
from .alleleclusterer import AlleleClusterer
from .alleleremover import AlleleRemover
from .allelefinder import AlleleFinder
from .performanceplotter import PerformancePlotter
from .partitionplotter import PartitionPlotter
from .hist import Hist
from . import seqfileopener

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

        self.vs_info, self.sw_info, self.msa_vs_info = None, None, None
        self.duplicates = {}
        self.bcrham_proc_info = None
        self.timing_info = []  # it would be really nice to clean up both this and bcrham_proc_info
        self.istep = None  # stupid hack to get around network file system issues (see self.subworkidr()
        self.subworkdirs = []  # arg. same stupid hack

        self.unseeded_seqs = None  # all the queries that we *didn't* cluster with the seed uid
        self.small_cluster_seqs = None  # all the queries that we removed after a few partition steps 'cause they were in small clusters

        self.sw_param_dir, self.hmm_param_dir, self.multi_hmm_param_dir = ['%s/%s' % (self.args.parameter_dir, s) for s in ['sw', 'hmm', 'multi-hmm']]
        self.sub_param_dir = utils.parameter_type_subdir(self.args, self.args.parameter_dir)
        self.final_multi_paramdir = utils.non_none([self.args.parameter_out_dir, self.multi_hmm_param_dir])  # ick

        self.hmm_infname = self.args.workdir + '/hmm_input.csv'
        self.hmm_cachefname = self.args.workdir + '/hmm_cached_info.csv'
        self.hmm_outfname = self.args.workdir + '/hmm_output.csv'
        self.cpath_progress_dir = '%s/cluster-path-progress' % self.args.workdir  # write the cluster paths for each clustering step to separate files in this dir

        self.print_status = self.args.debug  # if set, print some extra info (e.g. hmm calculation stats) that we used to print by default, but don't want to any more

        if self.args.outfname is not None:
            utils.prep_dir(dirname=None, fname=self.args.outfname, allow_other_files=True)

        self.input_partition, self.input_cpath = None, None
        if self.args.input_partition_fname is not None:
            self.input_glfo, self.input_antn_list, self.input_cpath = utils.read_yaml_output(self.args.input_partition_fname, skip_annotations=not self.args.continue_from_input_partition)
            if self.args.continue_from_input_partition:
                self.input_antn_dict = utils.get_annotation_dict(self.input_antn_list, ignore_duplicates=True)  # NOTE not really sure ignore_duplicates should be set, but when we're reading merged subset partitions it's nice to avoid the warnings, and in general duplicates doesn't seem like a big problem
            self.input_partition = self.input_cpath.partitions[self.input_cpath.i_best if self.args.input_partition_index is None else self.args.input_partition_index]
            # print('    %s input partition has duplicates: sum of cluster sizes %d vs %d unique queries' % (utils.wrnstr(), sum(len(c) for c in self.input_partition), len(set(u for c in self.input_partition for u in c))))
            print('  --input-partition-fname: read %s partition with %d sequences in %d clusters from %s' % ('best' if self.args.input_partition_index is None else 'index-%d'%self.args.input_partition_index, sum(len(c) for c in self.input_partition), len(self.input_partition), self.args.input_partition_fname))
            input_partition_queries = set(u for c in self.input_partition for u in c)
            ids_to_rm = set(self.input_info) - input_partition_queries
            for uid in ids_to_rm:
                del self.input_info[uid]
            if len(ids_to_rm) > 0:
                print('    removed %d/%d queries from input info that were absent from input partition' % (len(ids_to_rm), len(self.input_info) + len(ids_to_rm)))

        self.deal_with_persistent_cachefile()

        self.cached_naive_hamming_bounds = self.args.naive_hamming_bounds  # this exists so we don't get the bounds every iteration through the clustering loop (and here is just set to the value from the args, but is set for real below)

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
            'get-selection-metrics'        : self.read_existing_output,
            'get-linearham-info'           : self.read_existing_output,
            'update-meta-info'             : self.read_existing_output,
            'view-alternative-annotations'  : self.view_alternative_annotations,
        }

    # ----------------------------------------------------------------------------------------
    def sw_cache_path(self, find_any=False):
        if self.args.sw_cachefname is not None:
            return utils.getprefix(self.args.sw_cachefname)
        elif None not in [self.args.parameter_dir, self.input_info]:
            if find_any:
                fnames = glob.glob(self.args.parameter_dir + '/sw-cache*')  # remain suffix-agnostic
                if len(fnames) == 0:
                    raise Exception('couldn\'t find any sw cache files in %s, despite setting <find_any>' % self.args.parameter_dir)
                return utils.getprefix(fnames[0])
            else:
                return self.args.parameter_dir + '/sw-cache-' + utils.uidhashstr(''.join(self.input_info.keys()))  # remain suffix-agnostic
        else:
            return None

    # ----------------------------------------------------------------------------------------
    def get_cpath_progress_fname(self, istep):
        return '%s/istep-%d.csv' % (self.cpath_progress_dir, istep)

    # ----------------------------------------------------------------------------------------
    def get_all_cpath_progress_fnames(self):
        assert len(os.listdir(self.cpath_progress_dir)) == self.istep + 1  # not really checking for anything that has a decent chance of happening, it's more because I get the file names in a slightly different loop in self.clean()
        return [self.get_cpath_progress_fname(istep) for istep in range(self.istep + 1)]

    # ----------------------------------------------------------------------------------------
    def run(self, actions):
        self.all_actions = actions
        for tmpaction in actions:
            self.current_action = tmpaction  # NOTE gets changed on the fly below, I think just in self.get_annotations_for_partitions() (which is kind of hackey, but I can't figure out a way to improve on it that wouldn't involve wasting a foolish amount of time rewriting things. Bottom line is that the control flow for different actions is really complicated, and that complexity is going to show up somewhere)
            self.action_fcns[tmpaction]()

    # ----------------------------------------------------------------------------------------
    def clean(self):
        if self.args.new_allele_fname is not None:
            new_allele_region = 'v'
            new_alleles = [(g, seq) for g, seq in self.glfo['seqs'][new_allele_region].items() if glutils.is_snpd(g)]
            print('  writing %d new %s to %s' % (len(new_alleles), utils.plural_str('allele', len(new_alleles)), self.args.new_allele_fname))
            with open(self.args.new_allele_fname, 'w') as outfile:
                for name, seq in new_alleles:
                    outfile.write('>%s\n' % name)
                    outfile.write('%s\n' % seq)

        # merge persistent and current cache files into the persistent cache file
        if self.args.persistent_cachefname is not None:
            lockfname = self.args.persistent_cachefname + '.lock'
            while os.path.exists(lockfname):
                print('  waiting for lock on %s' % lockfname)
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
                print('  parsing annotation output file %s to partition cache file %s' % (self.args.persistent_cachefname, self.hmm_cachefname))
                with open(self.hmm_cachefname, utils.csv_wmode()) as outcachefile:
                    writer = csv.DictWriter(outcachefile, utils.partition_cachefile_headers)
                    writer.writeheader()
                    for line in reader:
                        if line['v_gene'] == '':  # failed
                            continue
                        utils.process_input_line(line)
                        outrow = {'unique_ids' : line['unique_ids'], 'naive_seq' : line['padlefts'][0] * utils.ambig_base + line['naive_seq'] + line['padrights'][0] * utils.ambig_base}
                        writer.writerow(outrow)
            elif set(reader.fieldnames) == set(utils.partition_cachefile_headers):  # headers are ok, so can just copy straight over
                check_call(['cp', self.args.persistent_cachefname, self.hmm_cachefname])
            else:
                raise Exception('--persistent-cachefname %s has unexpected header list %s' % (self.args.persistent_cachefname, reader.fieldnames))

    # ----------------------------------------------------------------------------------------
    def run_waterer(self, count_parameters=False, write_parameters=False, write_cachefile=False, look_for_cachefile=False, require_cachefile=False, dbg_str=''):
        print('smith-waterman%s' % (('  (%s)' % dbg_str) if dbg_str != '' else ''))
        sys.stdout.flush()

        self.vs_info = None  # should already be None, but we want to make sure (if --no-sw-vsearch is set we need it to be None, and if we just removed unlikely alleles we need to rerun vsearch with the likely alleles)
        if not self.args.no_sw_vsearch:
            self.set_vsearch_info(get_annotations=True)
        if self.args.simultaneous_true_clonal_seqs:  # it might be better to just copy over the true indel info in this case? it depends what you're trying to test, and honestly really if you're using this option you just shouldn't be putting indels in your simulation to start with
            print('  note: not running msa indel stuff for --simultaneous-true-clonal-seqs, so any families with shm indels within cdr3 will be split up before running the hmm. To fix this you\'ll either need to run set_msa_info() (which is fine and easy, but slow, and requires deciding whether to make sure to run parameter caching with the arg, or else rerun smith waterman with the msa indels')
        if self.args.all_seqs_simultaneous and self.msa_vs_info is None:  # only run the first time we run sw
            self.set_msa_info(debug=self.args.debug)
            look_for_cachefile, require_cachefile = False, False
            print('  note: ignoring any existing sw cache file to ensure we\'re getting msa indel info')  # the main use case for this is with 'annotate' or 'partition' on existing parameters that were run on the whole repertoire, so a) it shouldn't be a big deal to rerun and b) you probably don't want to run msa indel info when parameter caching. Also, it's not easy to figure out if msa indel info is in the sw cached file without first reading it

        pre_failed_queries = self.sw_info['failed-queries'] if self.sw_info is not None else None  # don't re-run on failed queries if this isn't the first sw run (i.e., if we're parameter caching)
        waterer = Waterer(self.args, self.glfo, self.input_info, self.simglfo, self.reco_info,  # NOTE if we're reading a cache file, this glfo gets replaced with the glfo from the file
                          count_parameters=count_parameters,
                          parameter_out_dir=self.sw_param_dir if write_parameters else None,
                          plot_annotation_performance=self.args.plot_annotation_performance,
                          duplicates=self.duplicates, pre_failed_queries=pre_failed_queries, aligned_gl_seqs=self.aligned_gl_seqs, vs_info=self.vs_info, msa_vs_info=self.msa_vs_info)

        cache_path = self.sw_cache_path(find_any=require_cachefile)
        cachefname = cache_path + ('.yaml' if self.args.sw_cachefname is None else utils.getsuffix(self.args.sw_cachefname))  # use yaml, unless csv was explicitly set on the command line
        if look_for_cachefile or require_cachefile:
            if os.path.exists(cache_path + '.csv'):  # ...but if there's already an old csv, use that
                cachefname = cache_path + '.csv'
        else:  # i.e. if we're not explicitly told to look for it (and it exists) then it should be out of date
            waterer.clean_cache(cache_path)  # hm, should this be <cachefname> instead of <cache_path>? i mean they're the same, but still
        if (look_for_cachefile or require_cachefile) and os.path.exists(cachefname):
            waterer.read_cachefile(cachefname)
        else:
            if require_cachefile:
                raise Exception('sw cache file %s not found' % cachefname)
            if look_for_cachefile:
                print('    couldn\'t find sw cache file %s, so running sw%s' % (cachefname, ' (this is probably because --seed-unique-id is set to a sequence that wasn\'t in the input file on which we cached parameters [if it\'s inconvenient to put your seed sequences in your input file, you can avoid this by putting them instead in separate file and set --queries-to-include-fname])' if self.args.seed_unique_id is not None else ''))
            waterer.run(cachefname if write_cachefile else None)

        self.sw_info = waterer.info
        self.sw_glfo = waterer.glfo  # ick
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
    def set_msa_info(self, debug=False):  # NOTE not running this for args.simultaneous_true_clonal_seqs any more, but i'm leaving the stuff in here for that arg in case I change my mind later
        # ----------------------------------------------------------------------------------------
        def run_msa(cluster):  # NOTE that this is really slow, and could probably be sped up? But i don't really care, the only time you'd run on a lot of families is simulation with tons of indels, which just isn't an important use case
            unln_seqfos = [{'name' : q, 'seq' : self.input_info[q]['seqs'][0]} for q in cluster]  # ignore the indels that already cam from vsearch, combining them would be hard (and we want the rest of the vsearch info for other purposes)
            if self.args.simultaneous_true_clonal_seqs and len(set(len(s['seq']) for s in unln_seqfos)) == 1:  # if all the seqs are the same length, they almost certainly don't have shm indels
                if debug:
                    print('    all %d seqs the same length, skipping' % len(unln_seqfos))
                return {'gene-counts' : None, 'annotations' : OrderedDict(), 'failures' : []}
            aln_seqfos = utils.align_many_seqs(unln_seqfos, extra_str='        ', debug=debug)
            cseq = utils.cons_seq(aligned_seqfos=aln_seqfos, extra_str='        ', debug=debug)
            indeld_cseq = []  # cons seq where we remove any "indels" (well, gaps in the msa) that are present in less than half the seqs
            for ich, cons_char in enumerate(cseq):
                msa_chars = [s['seq'][ich] for s in aln_seqfos]
                n_gap_chars = len([c for c in msa_chars if c in utils.gap_chars])
                # print ''.join(msa_chars), n_gap_chars, len(msa_chars), cons_char if n_gap_chars <= len(msa_chars) // 2 else ''
                if n_gap_chars <= len(msa_chars) // 2:  # if it's less than half gap chars, we want the cons char in indeld_cseq
                    indeld_cseq.append(cons_char)
            indeld_cseq = ''.join(indeld_cseq)
            if debug:
                print('    indeld cons seq: %s' % indeld_cseq)
            fglfo = glutils.get_empty_glfo(self.args.locus)
            fglfo['seqs']['v'] = {'IGHVx-x*x' : indeld_cseq}  # it's not a real v gene, it extends through the whole (vdj) sequence, but i have to put something here, and i think this won't cause problems
            return utils.run_vsearch('search', {s['name'] : s['seq'] for s in unln_seqfos}, self.args.workdir + '/vsearch', threshold=0.3, glfo=fglfo, vsearch_binary=self.args.vsearch_binary, get_annotations=True)  # don't really need to align again, but this gets us the cigar seqs automatically, and i REALLY don't want to write anything more to do with cigars (i.e. converting aln_seqfos to cigars)
        # ----------------------------------------------------------------------------------------
        print('  running mafft+vsearch for msa indel info for --all-seqs-simultaneous/--simultaneous-true-clonal-seqs')
        if self.args.all_seqs_simultaneous:  # if you set both of these, that's your problem, it doesn't make sense anyway
            nsets = [[q for q in self.input_info]]  # maybe i should exclude any that failed sw, but otoh if you set all simultaneous, that means you want *all* simultaneous
        elif self.args.simultaneous_true_clonal_seqs:
            nsets = utils.get_partition_from_reco_info(self.reco_info)
        else:
            assert False
        all_antns, all_failed_queries = OrderedDict(), []
        for cluster in nsets:
            cfo = run_msa(cluster)
            all_antns.update(cfo['annotations'])
            all_failed_queries += cfo['failures']
        self.msa_vs_info = {'gene-counts' : None, 'annotations' : all_antns, 'failures' : all_failed_queries}

    # ----------------------------------------------------------------------------------------
    def cache_parameters(self):
        print('caching parameters')

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
                glutils.add_new_alleles(self.glfo, list(alcluster_alleles.values()), use_template_for_codon_info=False, simglfo=self.simglfo, debug=True)
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
        print('hmm')
        sys.stdout.flush()
        _, annotations, hmm_failures = self.run_hmm('viterbi', self.sw_param_dir, parameter_out_dir=self.hmm_param_dir, count_parameters=True, partition=self.input_partition)
        if self.args.outfname is not None and self.current_action == self.all_actions[-1]:
            self.write_output(list(annotations.values()), hmm_failures, cpath=self.input_cpath)
        self.write_hmms(self.hmm_param_dir)  # note that this modifies <self.glfo>

    # ----------------------------------------------------------------------------------------
    def annotate(self):
        print('annotating    (with %s)%s' % (self.sub_param_dir, ' (and star tree annotation, since --subcluster-annotation-size is None)' if self.args.subcluster_annotation_size is None else ''))
        if self.sw_info is None:
            self.run_waterer(look_for_cachefile=not self.args.write_sw_cachefile, write_cachefile=self.args.write_sw_cachefile, count_parameters=self.args.count_parameters)
        if self.args.only_smith_waterman:
            if self.args.outfname is not None:  # NOTE this is _not_ identical to the sw cache file (e.g. padding, failed query writing, plus probably other stuff)
                self.write_output(None, set(), write_sw=True)  # note that if you're auto-parameter caching, this will just be rewriting an sw output file that's already there from parameter caching, but oh, well. If you're setting --only-smith-waterman and not using cache-parameters, you have only yourself to blame
            return
        print('hmm')
        self.added_extra_clusters_to_annotate = False  # ugh (see other places this gets set, and fcn in next line gets called)
        annotations, hmm_failures = self.actually_get_annotations_for_clusters(clusters_to_annotate=self.input_partition)
        if self.args.get_selection_metrics:
            self.calc_tree_metrics(annotations)  # adds tree metrics to <annotations>
        if self.args.annotation_clustering:  # VJ CDR3 clustering (NOTE it would probably be better to have this under 'partition' action, but it's historical and also not very important)
            from . import annotationclustering
            antn_ptn = annotationclustering.vollmers(annotations, self.args.annotation_clustering_threshold)
            antn_cpath = ClusterPath(partition=antn_ptn)
            self.get_annotations_for_partitions(antn_cpath)  # get new annotations corresponding to <antn_ptn>
        else:
            if self.args.outfname is not None:
                self.write_output(list(annotations.values()), hmm_failures, cpath=antn_cpath if self.args.annotation_clustering else self.input_cpath)
            if self.args.plot_partitions or self.input_partition is not None and self.args.plotdir is not None:
                assert self.input_partition is not None
                partplotter = PartitionPlotter(self.args, glfo=self.glfo)
                partplotter.plot(self.args.plotdir + '/partitions', self.input_partition, annotations, reco_info=self.reco_info, args=self.args)
            if self.args.count_parameters and not self.args.dont_write_parameters:
                self.write_hmms(self.final_multi_paramdir)  # note that this modifies <self.glfo>

    # ----------------------------------------------------------------------------------------
    def calc_tree_metrics(self, annotation_dict, annotation_list=None, cpath=None):
        if annotation_list is None:
            annotation_list = list(annotation_dict.values())
        if self.current_action == 'get-selection-metrics' and self.args.input_metafnames is not None:  # presumably if you're running 'get-selection-metrics' with --input-metafnames set, that means you didn't add the affinities (+ other metafo) when you partitioned, so we need to add it now
            seqfileopener.read_input_metafo(self.args.input_metafnames, annotation_list)
        if self.args.seed_unique_id is not None:  # restrict to seed cluster in the best partition (clusters from non-best partition have duplicate uids, which then make fasttree barf, and it doesn't seem worth the trouble to fix it now)
            print('    --seed-unique-id: restricting selection metric calculation to seed cluster in best partition (mostly to avoid fasttree crash on duplicate uids)')
            annotation_dict = OrderedDict([(uidstr, line) for uidstr, line in annotation_dict.items() if self.args.seed_unique_id in line['unique_ids'] and line['unique_ids'] in cpath.partitions[cpath.i_best]])
        treeutils.add_smetrics(self.args, self.args.selection_metrics_to_calculate, annotation_dict, self.args.lb_tau, reco_info=self.reco_info,  # NOTE keys in <annotation_dict> may be out of sync with 'unique_ids' if we add inferred ancestral seqs here
                               use_true_clusters=self.reco_info is not None, base_plotdir=self.args.plotdir, workdir=self.args.workdir,
                               outfname=self.args.selection_metric_fname, glfo=self.glfo, tree_inference_outdir=self.args.tree_inference_outdir, debug=self.args.debug)

    # ----------------------------------------------------------------------------------------
    def parse_existing_annotations(self, annotation_list, ignore_args_dot_queries=False, process_csv=False):
        n_queries_read = 0
        failed_query_strs, fake_paired_strs = set(), set()
        new_annotation_list = []
        for line in annotation_list:
            if process_csv:
                utils.process_input_line(line)
            uidstr = ':'.join(line['unique_ids'])
            if ('invalid' in line and line['invalid']) or line['v_gene'] == '':  # first way is the new way, but we have to check the empty-v-gene way too for old files
                failed_query_strs.add(uidstr)
                if line.get('is_fake_paired', False):
                    fake_paired_strs.add(uidstr)
                continue
            if self.args.queries is not None and not ignore_args_dot_queries:  # second bit is because when printing subcluster naive seqs, we want to read all the ones that have any overlap with self.args.queries, not just the exact cluster of self.args.queries
                if len(set(self.args.queries) & set(line['unique_ids'])) == 0:  # actually make sure this is the precise set of queries we want (note that --queries and line['unique_ids'] are both ordered, and this ignores that... oh, well, sigh.)
                    continue
            if self.args.reco_ids is not None and line['reco_id'] not in self.args.reco_ids:
                continue
            utils.add_implicit_info(self.glfo, line)
            new_annotation_list.append(line)

            n_queries_read += 1
            if self.args.n_max_queries > 0 and n_queries_read >= self.args.n_max_queries:
                break

        if len(failed_query_strs) > 0:
            print('\n%d failed queries%s' % (len(failed_query_strs), '' if len(fake_paired_strs) == 0 else ' (%d were fake paired annotations)' % len(fake_paired_strs)))
        return new_annotation_list, len(fake_paired_strs)

    # ----------------------------------------------------------------------------------------
    def view_alternative_annotations(self):
        print('  %s getting alternative annotation information from existing output file. These results will only be meaningful if you had --calculate-alternative-annotations set when writing the output file (so that all subcluster annotations were stored). We can\'t check for that here directly, so instead we print this warning to make sure you had it set ;-)' % utils.color('yellow', 'note'))

        # we used to require that you set --queries to tell us which to get, but I think now it makes sense to by default just get all of them (but not sure enough to delete this yet)
        # if self.args.queries is None:
        #     _, cpath = self.read_existing_output(read_partitions=True)
        #     clusterstrs = []
        #     for cluster in sorted(cpath.partitions[cpath.i_best], key=len, reverse=True):
        #         clusterstrs.append('      %s' % ':'.join(cluster))
        #     raise Exception('in order to view alternative annotations, you have to specify (with --queries) a cluster from the final partition. Choose from the following:\n%s' % '\n'.join(clusterstrs))

        cluster_annotations, cpath = self.read_existing_output(ignore_args_dot_queries=True, read_partitions=True, read_annotations=True)  # note that even if we don't need the cpath to re-write output below, we need to set read_annotations=True, since the fcn gets confused otherwise and doesn't read the right cluster annotation file (for deprecated csv files)

        clusters_to_use = cpath.partitions[cpath.i_best] if self.args.queries is None else [self.args.queries]
        n_skipped = 0
        for cluster in sorted(clusters_to_use, key=len, reverse=True):
            if len(cluster) < self.args.min_selection_metric_cluster_size:
                n_skipped += 1
                continue
            self.process_alternative_annotations(cluster, cluster_annotations, cpath=cpath, debug=True)
        if n_skipped > 0:
            print('  skipped %d clusters smaller than --min-selection-metric-cluster-size %d' % (n_skipped, self.args.min_selection_metric_cluster_size))

        print('  note: rewriting output file with newly-calculated alternative annotation info')
        self.write_output(list(cluster_annotations.values()), set(), cpath=cpath, dont_write_failed_queries=True)  # I *think* we want <dont_write_failed_queries> set, because the failed queries should already have been written, so now they'll just be mixed in with the others in <annotations>

    # ----------------------------------------------------------------------------------------
    def get_index_restricted_clusters(self, cpath):
        tptn = cpath.best()
        if self.args.partition_index_to_print is not None:
            if self.args.partition_index_to_print > len(cpath.partitions) - 1:
                cpath.print_partitions()
                print('  %s --partition-index-to-print %d too large for cpath with length %d, so ignoring it and using best partition' % (utils.wrnstr(), self.args.partition_index_to_print, len(cpath.partitions)))
            else:
                tptn = cpath.partitions[self.args.partition_index_to_print]
        tptn = sorted(tptn, key=len, reverse=True)  # NOTE we always want this sorted, i.e. dont_sort doesn't apply to this, since here we're only doing --cluster-indices, which says in its help message that we sort
        if self.args.cluster_indices is None:
            return tptn
        else:
            return [tptn[i] for i in self.args.cluster_indices]

    # ----------------------------------------------------------------------------------------
    def print_results(self, cpath, annotation_list, dont_sort=False, label_list=None, extra_str=''):
        if label_list is not None:
            assert len(label_list) == len(annotation_list)
        seed_uid = self.args.seed_unique_id
        true_partition = None
        restricted_clusters = None
        if cpath is not None and len(cpath.partitions) > 0:
            if len(annotation_list) > 0:  # this is here just so we get a warning if any of the clusters in the best partition are missing from the annotations
                _ = utils.get_annotation_dict(annotation_list, cpath=cpath)
            # it's expected that sometimes you'll write a seed partition cpath, but then when you read the file you don't bother to seed the seed id on the command line. The reverse, however, shouldn't happen
            if seed_uid is not None and cpath.seed_unique_id != seed_uid:
                print('  %s seed uids from args and cpath don\'t match %s %s ' % (utils.color('red', 'error'), self.args.seed_unique_id, cpath.seed_unique_id))
            if self.args.cluster_indices is not None:
                restricted_clusters = self.get_index_restricted_clusters(cpath)
            seed_uid = cpath.seed_unique_id
            n_to_print, ipart_center = None, None
            if self.args.partition_index_to_print is not None:
                print('  --partition-index-to-print: using non-default partition with index %d' % self.args.partition_index_to_print)
                n_to_print, ipart_center = 1, self.args.partition_index_to_print
            print('%s%s' % (extra_str, utils.color('green', 'partitions:')))
            cpath.print_partitions(abbreviate=self.args.abbreviate, reco_info=self.reco_info, highlight_cluster_indices=self.args.cluster_indices,
                                   calc_missing_values=('all' if cpath.n_seqs() < 500 else 'best'), print_partition_indices=True, n_to_print=n_to_print, ipart_center=ipart_center)
            if not self.args.is_data and self.reco_info is not None:  # if we're reading existing output, it's pretty common to not have the reco info even when it's simulation, since you have to also pass in the simulation input file on the command line
                true_partition = utils.get_partition_from_reco_info(self.reco_info)
                true_cp = ClusterPath(seed_unique_id=self.args.seed_unique_id)
                true_cp.add_partition(true_partition, -1., 1)
                print('%strue:' % extra_str)
                # print utils.per_seq_correct_cluster_fractions(cpath.partitions[cpath.i_best], true_partition, reco_info=self.reco_info, seed_unique_id=self.args.seed_unique_id)
                true_cp.print_partitions(self.reco_info, print_header=False, calc_missing_values='best', extrastr=extra_str, print_partition_indices=True)

        if len(annotation_list) > 0:
            print('%s%s' % (extra_str, utils.color('green', 'annotations:')))
            if dont_sort:
                sorted_annotations = annotation_list
            else:
                sorted_annotations = sorted(annotation_list, key=lambda l: len(l['unique_ids']), reverse=True)
            if self.args.cluster_indices is not None:
                print('    --cluster-indices: restricting to %d cluster%s with indices: %s' % (len(self.args.cluster_indices), utils.plural(len(self.args.cluster_indices)), ' '.join(str(i) for i in self.args.cluster_indices)))
                # sorted_annotations = [sorted_annotations[iclust] for iclust in self.args.cluster_indices]  # this is what it used to be, but this is wrong
                antn_dict = utils.get_annotation_dict(sorted_annotations)
                sorted_annotations = [antn_dict.get(':'.join(rc)) for rc in restricted_clusters]
                if None in sorted_annotations:
                    print('    %s missing %d requested annotations' % (utils.color('yellow', 'warning'), sorted_annotations.count(None)))
                    sorted_annotations = [l for l in sorted_annotations if l is not None]

            for iline, line in enumerate(sorted_annotations):
                if self.args.only_print_best_partition and cpath is not None and cpath.i_best is not None and line['unique_ids'] not in cpath.partitions[cpath.i_best]:
                    continue
                if (self.args.only_print_seed_clusters or self.args.seed_unique_id is not None) and seed_uid not in line['unique_ids']:  # we only use the seed id from the command line here, so you can print all the clusters even if you ran seed partitioning UPDATE wait did I change my mind? need to check
                    continue
                if self.args.only_print_queries_to_include_clusters and len(set(self.args.queries_to_include) & set(line['unique_ids'])) == 0:  # will barf if you don't tell us what queries to include, but then that's your fault isn't it
                    continue
                if self.args.print_trees:
                    treestr = line.get('tree', lbplotting.get_tree_in_line(line, self.args.is_simu))  # ok this weird, but i want to be able to print the tree on simulation files even without setting --is-simu, sincei if --is-simu is set i may not be able to print the annotations (since for that, reco_info has to be set, but that depends how the simulation file was read. Anyway...
                    if treestr is None:
                        print('  --print-trees: no tree found in line')
                    else:
                        dtree = treeutils.get_dendro_tree(treestr=treestr)
                        print(utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=dtree)))
                    continue  # eh, maybe just continue so it doesn't crash if it sees the multi-seq true annotation
                label, post_label = [], []
                if self.args.infname is not None and self.reco_info is not None:
                    utils.print_true_events(self.simglfo, self.reco_info, line, full_true_partition=true_partition, extra_str=extra_str+'  ')
                    label += ['inferred:']
                if cpath is not None and cpath.i_best is not None:  # maybe I could do the iparts stuff even if i_best isn't set, but whatever, I think it's only really not set if the cpath is null anyway
                    iparts = cpath.find_iparts_for_cluster(line['unique_ids'])
                    ipartstr = 'none' if len(iparts) == 0 else ' '.join([str(i) for i in iparts])
                    post_label += ['   partition%s: %s' % (utils.plural(len(iparts)), ipartstr)]
                    if cpath.i_best in iparts:
                        post_label += [', %s' % utils.color('yellow', 'best')]
                queries_to_emphasize = []
                if seed_uid is not None and seed_uid in line['unique_ids']:
                    post_label += [', %s' % utils.color('red', 'seed')]
                    queries_to_emphasize += [seed_uid]
                if self.args.queries_to_include is not None and len(set(self.args.queries_to_include) & set(line['unique_ids'])) > 0:  # will barf if you don't tell us what queries to include, but then that's your fault isn't it
                    post_label += [', %s' % utils.color('red', 'queries-to-include')]
                    queries_to_emphasize += self.args.queries_to_include
                if label_list is not None:
                    post_label += [label_list[iline]]
                utils.print_reco_event(line, extra_str=extra_str+'  ', label=''.join(label), post_label=''.join(post_label), queries_to_emphasize=queries_to_emphasize, extra_print_keys=self.args.extra_print_keys)

    # ----------------------------------------------------------------------------------------
    def restrict_ex_out_clusters(self, cpath, annotation_list):  # NOTE would be nice to use [bits of] this also for result printing fcn above, but there we need the indices to line up, so have to do the 'continue' thing
        n_before = len(annotation_list)
        ptn_to_use, dbg_str = cpath.best() , []
        if self.args.only_print_best_partition and cpath is not None and cpath.i_best is not None:
            annotation_list = [l for l in annotation_list if l['unique_ids'] in cpath.partitions[cpath.i_best]]
        if self.args.only_print_seed_clusters or self.args.seed_unique_id is not None:
            annotation_list = [l for l in annotation_list if self.args.seed_unique_id in l['unique_ids']]
            ptn_to_use = [c for c in ptn_to_use if self.args.seed_unique_id in c]
        if self.args.only_print_queries_to_include_clusters:
            annotation_list = [l for l in annotation_list if len(set(self.args.queries_to_include) & set(l['unique_ids'])) > 0]  # will barf if you don't tell us what queries to include, but then that's your fault isn't it
            ptn_to_use = [c for c in ptn_to_use if len(set(self.args.queries_to_include) & set(c)) > 0]
        if self.args.n_final_clusters is not None or self.args.min_largest_cluster_size is not None:
            tptns = cpath.partitions
            if self.args.n_final_clusters is not None:
                tptns = [p for p in tptns if len(p) == self.args.n_final_clusters]
                dbg_str.append('--n-final-clusters')
            if self.args.min_largest_cluster_size is not None:
                tptns = [p for p in cpath.partitions if any(len(c) >= self.args.min_largest_cluster_size for c in p)]
                dbg_str.append('--min-largest-cluster-size')
            if len(tptns) > 1:
                print('    %s multiple partitions satisfy --n-final-clusters/--min-largest-cluster-size criteria, just picking first one' % utils.wrnstr())
            ptn_to_use = cpath.partitions[-1] if len(tptns)==0 else tptns[0]
            annotation_list = [l for l in annotation_list if l['unique_ids'] in ptn_to_use]
        if self.args.partition_index_to_print is not None or self.args.cluster_indices is not None:
            ptn_to_use = self.get_index_restricted_clusters(cpath)
            annotation_list = [l for l in annotation_list if l['unique_ids'] in ptn_to_use]
            if self.args.partition_index_to_print is not None:
                dbg_str.append('--partition-index-to-print')
            if self.args.cluster_indices is not None:
                dbg_str.append('--cluster-indices')
        if self.args.only_print_best_partition or self.args.only_print_seed_clusters or self.args.only_print_queries_to_include_clusters or len(dbg_str) > 0:
            astr = ', '.join(['--only-print-'+s for s in ['best-partition', 'seed-clusters', 'queries-to-include-clusters'] if getattr(self.args, ('only-print-'+s).replace('-', '_'))] + dbg_str)
            print('  %s: restricting to %d/%d annotations' % (astr, len(annotation_list), n_before))
        else:
            print('  note: By default we print/operate on *all* annotations in the output file, which in general can include annotations from non-best partititons and non-seed clusters (e.g. if --n-final-clusters was set).\n        If you want to restrict to particular annotations, use one of --only-print-best-partition, --only-print-seed-clusters, or --only-print-queries-to-include-clusters (or, if set during partitioning, --n-final-clusters or --min-largest-cluster-size).')

        return ptn_to_use, annotation_list

    # ----------------------------------------------------------------------------------------
    def read_existing_output(self, outfname=None, ignore_args_dot_queries=False, read_partitions=False, read_annotations=False):
        if outfname is None:
            outfname = self.args.outfname

        annotation_list = []
        cpath = None
        tmpact = self.current_action  # just a shorthand for brevity
        if utils.getsuffix(outfname) == '.csv':  # old way
            if tmpact == 'view-partitions' or tmpact == 'plot-partitions' or tmpact == 'view-output' or tmpact == 'get-selection-metrics' or read_partitions:
                cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id, fname=outfname)
            if tmpact == 'view-annotations' or tmpact == 'plot-partitions' or tmpact == 'view-output' or tmpact == 'get-selection-metrics' or read_annotations:
                csvfile = open(outfname if cpath is None else self.args.cluster_annotation_fname)  # closes on function exit, and no this isn't a great way of doing it (but it needs to stay open for the loop)
                reader = csv.DictReader(csvfile)
                if 'unique_ids' not in reader.fieldnames:
                    raise Exception('not an annotation file: %s' % outfname)
                annotation_list = list(reader)
        elif utils.getsuffix(outfname) == '.yaml':  # new way
            # NOTE replaces <self.glfo>, which is definitely what we want (that's the point of putting glfo in the yaml file), but it's still different behavior than if reading a csv
            assert self.glfo is None  # make sure bin/partis successfully figured out that we would be reading the glfo from the yaml output file
            self.glfo, annotation_list, cpath = utils.read_yaml_output(outfname, n_max_queries=self.args.n_max_queries, dont_add_implicit_info=True, seed_unique_id=self.args.seed_unique_id)  # add implicit info below, so we can skip some of 'em
        else:
            raise Exception('unhandled annotation file suffix %s' % outfname)

        annotation_list, n_fake_paired = self.parse_existing_annotations(annotation_list, ignore_args_dot_queries=ignore_args_dot_queries, process_csv=utils.getsuffix(outfname) == '.csv')  # NOTE modifies <annotation_list>
        ptn_to_use, annnotation_list = self.restrict_ex_out_clusters(cpath, annotation_list)
        if len(annotation_list) == 0:
            if cpath is not None and tmpact in ['view-output', 'view-annotations', 'view-partitions']:
                self.print_results(cpath, [])  # used to just return, but now i want to at least see the cpath
            print('zero annotations to print, exiting%s' % ('' if n_fake_paired==0 else ' (%s %d were fake paired annotations)'%(utils.color('yellow', 'note'), n_fake_paired)))
            return
        annotation_dict = utils.get_annotation_dict(annotation_list)  # returns none type if there's duplicate annotations
        extra_headers = list(set([h for l in annotation_list for h in l.keys() if h not in utils.annotation_headers]))  # note that this basically has to be hackey/wrong, since we're trying to guess what the headers were when the file was written

        if tmpact == 'get-linearham-info':
            self.input_info = OrderedDict([(u, {'unique_ids' : [u], 'seqs' : [s]}) for l in annotation_list for u, s in zip(l['unique_ids'], l['input_seqs'])])  # this is hackey, but I think is ok (note that the order won't be the same as it would've been before)
            self.run_waterer(require_cachefile=True)
            uids_missing_sw_info = [u for l in annotation_list for u in l['unique_ids'] if u not in self.sw_info]  # make sure <annotation_list> and <self.sw_info> have the same uids (they can get out of sync because we're re-running sw here with potentially different options (or versions) to whatever command was run to make <annotation_list>, e.g. if --is-simu was turned on duplicates won't have been removed before)
            dup_dict = {d : u for u in self.sw_info['queries'] for d in self.sw_info[u]['duplicates'][0]}
            for uid in uids_missing_sw_info:  # NOTE <self.sw_info> is somewhat inconsistent after we do this, but this code should only get run when we're just adding linearham info so idgaf
                if uid in dup_dict:
                    self.sw_info[uid] = self.sw_info[dup_dict[uid]]  # it would be proper to fix the duplicates in here, and probably some other things
                else:
                    pass  # switching this to pass: even though I think it can ony happen if people are running with options that don't really make sense, I don't think there's really a pressing need to crash here raise Exception('no sw info for query %s' % uid)  # I can't really do anything else, it makes no sense to go remove it from <annotation_list> when the underlying problem is (probably) that sw info was just run with different options
            self.write_output(annotation_list, set(), cpath=cpath, outfname=self.args.linearham_info_fname, dont_write_failed_queries=True, extra_headers=extra_headers)  # I *think* we want <dont_write_failed_queries> set, because the failed queries should already have been written, so now they'll just be mixed in with the others in <annotation_list>

        if tmpact == 'get-selection-metrics':
            self.calc_tree_metrics(annotation_dict, annotation_list=annotation_list, cpath=cpath)  # adds tree metrics to <annotations>
        if tmpact == 'update-meta-info':
            inids, alist_ids = set(self.input_info), set(u for l in annotation_list for u in l['unique_ids'])
            if inids != alist_ids:
                print('  %s input info (len %d) has different uids to annotation list (len %d), and meta info will only be set/correct for the ones in input info (%d missing from input info, %d missing from annotation list, %d in common)' % (utils.wrnstr(), len(inids), len(alist_ids), len(alist_ids - inids), len(inids - alist_ids), len(inids  & alist_ids)))
            # may want to add this? not sure: overwrite_all=True
            seqfileopener.add_input_metafo(self.input_info, annotation_list, keys_not_to_overwrite=['multiplicities', 'paired-uids'])  # these keys are modified by sw (multiplicities) or paired clustering (paired-uids), so if you want to update them with this action here you're out of luck
        if tmpact == 'update-meta-info' or (tmpact == 'get-selection-metrics' and self.args.add_selection_metrics_to_outfname):
            print('  rewriting output file with %s: %s' % ('newly-calculated selection metrics' if tmpact=='get-selection-metrics' else 'updated input meta info', outfname))
            if self.args.add_selection_metrics_to_outfname and 'gctree' in self.args.tree_inference_method:
                print('  %s writing gctree annotations (with inferred ancestral sequences added) to original output file, which means that if you rerun gctree things may crash/be messed up since the inferred ancestral sequences are already in the annotation' % utils.wrnstr())
            self.write_output(annotation_list, set(), cpath=cpath, dont_write_failed_queries=True, extra_headers=extra_headers)  # I *think* we want <dont_write_failed_queries> set, because the failed queries should already have been written, so now they'll just be mixed in with the others in <annotation_list>

        if self.args.align_constant_regions:
            utils.parse_constant_regions(self.args.species, self.args.locus, annotation_list, self.args.workdir, csv_outdir=os.path.realpath(os.path.dirname(self.args.outfname)) if self.args.outfname is not None else None, debug=self.args.debug)

        if tmpact == 'plot-partitions':
            partplotter = PartitionPlotter(self.args, glfo=self.glfo)
            partplotter.plot(self.args.plotdir + '/partitions', ptn_to_use, annotation_dict, reco_info=self.reco_info, args=self.args)

        if tmpact in ['view-output', 'view-annotations', 'view-partitions']:
            self.print_results(cpath, annotation_list)

        return annotation_dict, cpath

    # ----------------------------------------------------------------------------------------
    def partition(self):
        """ Partition sequences in <self.input_info> into clonally related lineages """
        print('partitioning     (with %s)' % self.sub_param_dir)
        if self.sw_info is None:
            self.run_waterer(look_for_cachefile=not self.args.write_sw_cachefile, write_cachefile=self.args.write_sw_cachefile, count_parameters=False)  # self.args.count_parameters)  # run smith-waterman
        if len(self.sw_info['queries']) == 0:
            if self.args.outfname is not None:
                self.write_output([], set())
            return
        if self.args.only_smith_waterman:
            return

        print('hmm')

        # pre-cache hmm naive seq for each single query NOTE <self.current_action> is still 'partition' for this (so that we build the correct bcrham command line)
        if self.args.persistent_cachefname is not None:
            print('  --persistent-cachefname: using existing hmm cache file %s' % self.args.persistent_cachefname)
        if self.args.persistent_cachefname is None or not os.path.exists(self.hmm_cachefname):  # if the default (no persistent cache file), or if a not-yet-existing persistent cache file was specified
            print('%scaching all %d naive sequences' % ('' if self.print_status else '  ', len(self.sw_info['queries'])), end='\n' if self.input_partition is not None and self.args.continue_from_input_partition else ' ')
            if self.args.synthetic_distance_based_partition:
                self.write_bcrham_cache_file([[q] for q in self.sw_info['queries']])
            elif self.input_partition is not None and self.args.continue_from_input_partition:
                self.write_bcrham_cache_file(self.input_partition, ctype='input')
            else:
                self.run_hmm('viterbi', self.sub_param_dir, precache_all_naive_seqs=True)  # , n_procs=self.auto_nprocs(len(self.sw_info['queries']))

        if self.args.simultaneous_true_clonal_seqs:
            print('  --simultaneous-true-clonal-seqs: using true clusters instead of partitioning')
            true_partition = [[uid for uid in cluster if uid in self.sw_info] for cluster in utils.get_partition_from_reco_info(self.reco_info)]  # mostly just to remove duplicates, although I think there might be other reasons why a uid would be missing
            cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id, partition=true_partition)
        elif self.args.all_seqs_simultaneous:
            print('  --all-seqs-simultaneous: using single cluster instead of partitioning')
            one_clust_ptn = [[u for u in self.sw_info['queries']]]
            cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id, partition=one_clust_ptn)
        elif self.input_partition is not None and not self.args.continue_from_input_partition:
            print('  --input-partition-fname: using input cpath instead of running partitioning')
            cpath = self.input_cpath
        elif self.args.naive_vsearch:  # or self.args.naive_swarm:
            cpath = self.cluster_with_naive_vsearch_or_swarm(parameter_dir=self.sub_param_dir)
        else:
            cpath = self.cluster_with_bcrham()

        self.get_annotations_for_partitions(cpath)

        self.check_partition(cpath.partitions[cpath.i_best])

    # ----------------------------------------------------------------------------------------
    def split_seeded_clusters(self, old_cpath):  # NOTE similarity to clusterpath.remove_unseeded_clusters()
        start = time.time()
        seeded_clusters, unseeded_clusters = utils.split_partition_with_criterion(old_cpath.partitions[old_cpath.i_best_minus_x], lambda cluster: self.args.seed_unique_id in cluster)
        self.unseeded_seqs = [uid for uclust in unseeded_clusters for uid in uclust]  # note that we no longer expect them to all be singletons, since we're merging queries with identical naive seqs before passing to glomerator.cc
        seeded_singleton_set = set([uid for sclust in seeded_clusters for uid in sclust])  # in case there's duplicates
        seeded_partition = utils.collapse_naive_seqs(self.synth_sw_info(seeded_singleton_set), split_by_cdr3=True)
        seeded_cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
        seeded_cpath.add_partition(seeded_partition, -1., 1)
        print('      removed %d sequences in unseeded clusters,' % len(self.unseeded_seqs), end=' ')
        print('split %d seeded clusters into %d singletons, and merged these into %d clusters with identical naive seqs (%.1f sec)' % (len(seeded_clusters), len(seeded_singleton_set), len(seeded_cpath.partitions[seeded_cpath.i_best_minus_x]), time.time() - start))

        return seeded_cpath

    # ----------------------------------------------------------------------------------------
    def remove_small_clusters(self, old_cpath):
        assert self.small_cluster_seqs is None  # at least for now, we want to call this only once (would otherwise need to modify it)
        big_clusters, small_clusters = utils.split_partition_with_criterion(old_cpath.partitions[old_cpath.i_best_minus_x], lambda cluster: len(cluster) not in self.args.small_clusters_to_ignore)
        if self.args.seed_unique_id is not None:  # should probably be implemented at some point
            print('  %s not specifically keeping --seed-unique-id sequence when removing small clusters' % utils.wrnstr())
        if self.args.queries_to_include is not None:
            kept_clusts = []
            for ism, sclust in enumerate(small_clusters):
                if any(q in sclust for q in self.args.queries_to_include):
                    small_clusters[ism] = None
                    big_clusters.append(sclust)
                    kept_clusts.append(sclust)
            small_clusters = [c for c in small_clusters if c is not None]
            if len(kept_clusts) > 0:
                print('    --queries-to-include: keeping %d small clusters that include specified queries with sizes: %s' % (len(kept_clusts), ' '.join(str(len(c)) for c in sorted(kept_clusts, reverse=True))))
        self.small_cluster_seqs = [sid for sclust in small_clusters for sid in sclust]
        new_cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
        new_cpath.add_partition(big_clusters, -1., 1)
        ntot = sum(len(c) for c in old_cpath.best())
        print('      --small-clusters-to-ignore: removing %d / %d (%.3f) sequences in %d / %d (%.3f) clusters (with sizes among %s)' % (len(self.small_cluster_seqs), ntot, len(self.small_cluster_seqs) / ntot, len(small_clusters), len(old_cpath.best()), len(small_clusters) / len(old_cpath.best()), ' '.join([str(sz) for sz in self.args.small_clusters_to_ignore])))
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
        print('        new n_procs %d (initial seqs/proc: %.2f   new seqs/proc: %.2f' % (new_n_procs, float(initial_nseqs) / initial_nprocs, float(new_n_clusters) / new_n_procs))
        return new_n_procs

    # ----------------------------------------------------------------------------------------
    def shall_we_reduce_n_procs(self, last_n_procs):
        if self.timing_info[-1]['total'] < self.args.min_hmm_step_time:  # mostly for when you're running on really small samples
            return True
        n_calcd_per_process = self.get_n_calculated_per_process()
        if n_calcd_per_process < self.args.n_max_to_calc_per_process and last_n_procs > 2:  # should be replaced by time requirement, since especially in later iterations, the larger clusters make this a crappy metric (2 is kind of a special case, becase, well, small integers and all)
            return True
        times_to_try_this_n_procs = max(4, last_n_procs)  # if we've already milked this number of procs for most of what it's worth (once you get down to 2 or 3, you don't want to go lower)
        if self.n_proc_list.count(last_n_procs) >= times_to_try_this_n_procs:
            return True

        return False

    # ----------------------------------------------------------------------------------------
    def prepare_next_iteration(self, cpath, initial_nseqs):
        last_n_procs = self.n_proc_list[-1]
        next_n_procs = last_n_procs

        factor = 1.3
        if self.shall_we_reduce_n_procs(last_n_procs):
            next_n_procs = int(next_n_procs / float(factor))

        def time_to_remove_some_seqs(n_proc_threshold):
            return len(self.n_proc_list) >= n_proc_threshold or next_n_procs == 1

        if self.args.small_clusters_to_ignore is not None and self.small_cluster_seqs is None and time_to_remove_some_seqs(self.args.n_steps_after_which_to_ignore_small_clusters):
            cpath = self.remove_small_clusters(cpath)
            next_n_procs = self.scale_n_procs_for_new_n_clusters(initial_nseqs, self.n_proc_list[0], cpath)
        if self.args.seed_unique_id is not None and self.unseeded_seqs is None and time_to_remove_some_seqs(3):  # if we didn't already remove the unseeded clusters in a previous partition step
            if (self.args.n_final_clusters is not None or self.args.min_largest_cluster_size is not None) and not self.set_force_args:  # need to add an additional iteration here with at least one of the force args set
                self.set_force_args = True
                next_n_procs = last_n_procs
            else:
                cpath = self.split_seeded_clusters(cpath)
                next_n_procs = self.scale_n_procs_for_new_n_clusters(initial_nseqs, self.n_proc_list[0], cpath)

        return next_n_procs, cpath

    # ----------------------------------------------------------------------------------------
    def get_n_calculated_per_process(self):
        assert self.bcrham_proc_info is not None

        total = 0.  # sum over each process
        for procinfo in self.bcrham_proc_info:
            if 'vtb' not in procinfo['calcd'] or 'fwd' not in procinfo['calcd']:
                print('%s couldn\'t find vtb/fwd in:\n%s' % (utils.color('red', 'warning'), procinfo['calcd']))  # may as well not fail, it probably just means we lost some stdout somewhere. Which, ok, is bad, but let's say it shouldn't be fatal.
                return 1.  # er, or something?
            if self.args.naive_hamming_cluster:  # make sure we didn't accidentally calculate some fwds
                assert procinfo['calcd']['fwd'] == 0.
            total += procinfo['calcd']['vtb'] + procinfo['calcd']['fwd']
        if self.args.debug:
            print('          vtb + fwd calcd: %d (%.1f per proc)' % (total, float(total) / len(self.bcrham_proc_info)))
        return float(total) / len(self.bcrham_proc_info)

    # ----------------------------------------------------------------------------------------
    def merge_shared_clusters(self, cpath, debug=False):  # replace the most likely partition with a new partition in which any clusters that share a sequence have been merged
        # cpath.partitions[cpath.i_best] = [['a', 'b', 'c', 'e'], ['d'], ['f', 'a'], ['g'], ['h'], ['i'], ['j', 'a'], ['x', 'y', 'z', 'd'], ['xx', 'x']]
        partition = cpath.partitions[cpath.i_best]

        if debug:
            print('merging shared clusters')
            cpath.print_partitions()

        # find every pair of clusters that has some overlap
        cluster_groups = []
        if debug:
            print(' making cluster_groups')
        for iclust in range(len(partition)):
            for jclust in range(iclust + 1, len(partition)):
                if len(set(partition[iclust]) & set(partition[jclust])) > 0:
                    if debug:
                        print('  %d %d' % (iclust, jclust))
                    cluster_groups.append(set([iclust, jclust]))

        # merge these pairs of clusters into groups
        while True:
            no_more_merges = True
            for cp1, cp2 in itertools.combinations(cluster_groups, 2):
                if len(cp1 & cp2) > 0:
                    if debug:
                        print('  merging %s and %s' % (cp1, cp2))
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
            print(' removing')
        for iclust in sorted([i for cgroup in cluster_groups for i in cgroup], reverse=True):
            if debug:
                print('    %d' % iclust)
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
            print('  stopping with %d (<= %d) clusters' % (len(cpath.partitions[cpath.i_best]), self.args.n_final_clusters))
            return True
        elif self.args.max_cluster_size is not None and max([len(c) for c in cpath.partitions[cpath.i_best]]) > self.args.max_cluster_size:  # NOTE I *think* I want the best, not best-minus-x here (hardish to be sure a.t.m., since I'm not really using the minus-x part right now)
            print('   --max-cluster-size (partitiondriver): stopping with a cluster of size %d (> %d)' % (max([len(c) for c in cpath.partitions[cpath.i_best]]), self.args.max_cluster_size))
            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------
    def synth_sw_info(self, queries):  # only used for passing info to utils.collapse_naive_seqs()
        # this uses the cached hmm naive seqs (since we have them and they're better) but then later we pass the hmm the sw annotations, so we have to make sure the sw cdr3 length is the same within each cluster (it's very rare that it isn't)
        synth_sw_info = {q : {'naive_seq' : s, 'cdr3_length' : self.sw_info[q]['cdr3_length']} for q, s in self.get_cached_hmm_naive_seqs(queries).items()}  # NOTE code duplication in cluster_with_bcrham()
        synth_sw_info['queries'] = list(synth_sw_info.keys())
        return synth_sw_info

    # ----------------------------------------------------------------------------------------
    def init_cpath(self, n_procs):
        initial_nseqs = len(self.sw_info['queries'])  # maybe this should be the number of clusters, now that we're doing some preclustering here?
        if self.input_partition is not None and self.args.continue_from_input_partition:
            print('      --continue-from-input-partition: using input partition for initial cpath')
            cpath = self.input_cpath
            print('             %d clusters (%d seqs)' % (len(cpath.best()), sum(len(c) for c in cpath.best())))
            # maybe i should split by cdr3?
            # nsets = utils.split_clusters_by_cdr3(nsets, self.sw_info, warn=True)
        else:
            initial_nsets = utils.collapse_naive_seqs(self.synth_sw_info(self.sw_info['queries']), split_by_cdr3=True, debug=True)
            cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
            cpath.add_partition(initial_nsets, logprob=0., n_procs=n_procs)  # NOTE sw info excludes failed sequences (and maybe also sequences with different cdr3 length)
        os.makedirs(self.cpath_progress_dir)
        if self.args.debug:
            print('    initial cpath:')
            cpath.print_partitions(abbreviate=self.args.abbreviate, reco_info=self.reco_info)
        return cpath, initial_nseqs

    # ----------------------------------------------------------------------------------------
    def merge_cpaths_from_previous_steps(self, final_cpath, debug=False):
        if debug:
            print('final (unmerged) cpath:')
            final_cpath.print_partitions(abbreviate=self.args.abbreviate)
            print('')

        n_before, n_after = self.args.n_partitions_to_write, self.args.n_partitions_to_write  # this takes more than we need, since --n-partitions-to-write is the *full* width, not half-width, but oh, well
        if self.args.debug or (self.args.calculate_alternative_annotations and self.args.subcluster_annotation_size is None) or self.args.get_selection_metrics:  # take all of 'em
            n_before, n_after = sys.maxsize, sys.maxsize
        elif self.args.write_additional_cluster_annotations is not None:
            n_before, n_after = [max(waca, n_) for waca, n_ in zip(self.args.write_additional_cluster_annotations, (n_before, n_after))]
        # NOTE we don't actually do anything with <n_after>, since we can't add any extra partitions here (well, we don't want to)

        cpfnames = self.get_all_cpath_progress_fnames()  # list of cpath files for each clustering step (last one corresponds to <final_cpath>)

        if final_cpath.i_best >= n_before or len(cpfnames) < 2:  # if we already have enough partitions, or if there was only one step, there's nothing to do
            if debug:
                print('   nothing to merge')
            return final_cpath

        icpfn = len(cpfnames) - 1
        merged_cp = ClusterPath(fname=cpfnames[icpfn], seed_unique_id=self.args.seed_unique_id)  # merged one is initially just the cp from the last step
        assert merged_cp.partitions[merged_cp.i_best] == final_cpath.partitions[final_cpath.i_best]  # shouldn't really be necessary, and is probably kind of slow
        while merged_cp.i_best < n_before and icpfn > 0:  # keep trying to add them until we have <n_before> of them
            icpfn -= 1
            previous_cp = ClusterPath(fname=cpfnames[icpfn], seed_unique_id=self.args.seed_unique_id)
            for ip in range(len(merged_cp.partitions)):
                if len(merged_cp.partitions[ip]) == len(previous_cp.partitions[-1]):  # skip identical partitions (for speed, first check if they have the same number of clusters, then whether the clusters are the same)
                    if set([tuple(c) for c in merged_cp.partitions[ip]]) == set([tuple(c) for c in previous_cp.partitions[-1]]):  # no, they're not always in the same order (I think because they get parcelled out to different processes, and then read back in random order)
                        if math.isinf(merged_cp.logprobs[ip]) and not math.isinf(previous_cp.logprobs[-1]):  # it should only be possible for the *later* partition to have non-infinite logprob, since we usually only calculate full logprobs in the last clustering step (which is why we're taking the later partition, from merged_cp), so print an error if the earlier one, that we're about to throw away, is the one that's non-infinite
                            print('%s earlier partition (that we\'re discarding) has non-infinite logprob %f, while later partition\'s is infinite %f' % (utils.color('red', 'error'), previous_cp.logprobs[-1], merged_cp.logprobs[ip]))
                        previous_cp.remove_partition(len(previous_cp.partitions) - 1)  # remove it from previous_cp, since we want the one that may have a logprob set (and which has smaller n_procs, although I don't think we care about that)
                previous_cp.add_partition(list(merged_cp.partitions[ip]), merged_cp.logprobs[ip], merged_cp.n_procs[ip])  # add each partition in the existing merged cp to the previous cp
            merged_cp = previous_cp
            assert merged_cp.partitions[merged_cp.i_best] == final_cpath.partitions[final_cpath.i_best]  # shouldn't really be necessary, and is probably kind of slow
            if debug:
                print('%s' % utils.color('red', str(icpfn)))
                merged_cp.print_partitions()

        if icpfn > 0:
            print('  %s not merging entire cpath history' % utils.color('yellow', 'note'))

        # this kind of sucks, and shouldn't be necessary, but someone reported that they're seeing cluster paths with all logprobs -inf (which is really bad since the actual ClusterPath code will know that the last partition should be the best, but someone reading the file by hand won't know). If I had a working example of this happening i could figure out a better place to put this check, but I don't
        if all(l == float('-inf') for l in merged_cp.logprobs):
            print('  %s all %d partitions in cluster path have log prob -inf, so setting last (best) to zero' % (utils.wrnstr(), len(merged_cp.partitions)))
            merged_cp.logprobs[-1] = 0.

        return merged_cp

    # ----------------------------------------------------------------------------------------
    def cluster_with_bcrham(self):
        tmpstart = time.time()
        self.set_force_args = False  # annoying shenanigans to make sure that if both --seed-unique-id and either of --n-final-clusters or --min-largest-cluster-size are set, that the "force" args are set in bcrham *before* we remove unseeded seqs
        n_procs = self.args.n_procs
        cpath, initial_nseqs = self.init_cpath(n_procs)
        self.n_proc_list = []
        self.istep = 0
        start = time.time()
        while n_procs > 0:
            if n_procs > len(cpath.bmx()):
                print('  reducing n procs to number of clusters: %d --> %d' % (n_procs, len(cpath.bmx())))
                n_procs = len(cpath.bmx())
            print('%s%d clusters with %d proc%s%s' % ('' if self.print_status else '  ', len(cpath.bmx()), n_procs, utils.plural(n_procs), '\n' if self.print_status else ''), end=' ')  # NOTE that a.t.m. i_best and i_best_minus_x are usually the same, since we're usually not calculating log probs of partitions (well, we're trying to avoid calculating any extra log probs, which means we usually don't know the log prob of the entire partition)
            cpath, _, _ = self.run_hmm('forward', self.sub_param_dir, n_procs=n_procs, partition=cpath.bmx(), shuffle_input=True)  # note that this annihilates the old <cpath>, which is a memory optimization (but we write all of them to the cpath progress dir)
            self.n_proc_list.append(n_procs)
            if self.are_we_finished_clustering(n_procs, cpath):
                break
            n_procs, cpath = self.prepare_next_iteration(cpath, initial_nseqs)
            self.istep += 1

        if self.args.max_cluster_size is not None:
            print('   --max-cluster-size (partitiondriver): merging shared clusters')
            self.merge_shared_clusters(cpath)

        cpath = self.merge_cpaths_from_previous_steps(cpath)

        print('      partition loop time: %.1f' % (time.time()-start))
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
            print('  %s %d queries missing from partition' % (utils.wrnstr(), len(missing_ids)))

    # # ----------------------------------------------------------------------------------------
    # def auto_nprocs(self, nseqs):
    #     if self.args.n_precache_procs is not None:  # command line override
    #         return self.args.n_precache_procs

    #     n_max_procs = 100

    #     if nseqs < 1000:
    #         seqs_per_proc = 250
    #     elif nseqs < 3000:
    #         seqs_per_proc = 500
    #     else:
    #         seqs_per_proc = 1000
    #     if self.args.batch_system is not None:  # if we're using a batch system, all the overhead (and priority issues) means it makes sense to have fewer processes
    #         seqs_per_proc *= 4
    #     n_precache_procs = int(math.ceil(float(nseqs) / seqs_per_proc))
    #     n_precache_procs = min(n_precache_procs, n_max_procs)  # I can't get more'n a few hundred slots at a time, so it isn't worth using too much more than that
    #     if self.args.batch_system is None:  # if we're not on a batch system, make sure it's less than the number of cpus
    #         n_precache_procs = min(n_precache_procs, multiprocessing.cpu_count())
    #     else:
    #         n_precache_procs = min(n_precache_procs, self.args.n_procs)  # aw, screw it, just limit it to --n-procs

    #     return n_precache_procs

    # ----------------------------------------------------------------------------------------
    # make new/fake annotations for <partition> using sw info (i.e. convert single-seq sw annotations to multi-seq annotations corresponding to clusters in <partition>)
    def convert_sw_annotations(self, partition=None):
        if partition is None:
            partition = self.get_nsets('viterbi', None)
        antn_dict, n_failed = OrderedDict(), 0
        for cluster in partition:
            antn = utils.synthesize_multi_seq_line_from_reco_info(cluster, self.sw_info, warn=False)  # can't really warn on different values from different single-seq annotations, there's way too many that do
            utils.remove_all_implicit_info(antn)  # gotta remove + re-add implicit info to get the naive seq the right length (I think just to remove padding)
            try:
                utils.add_implicit_info(self.sw_glfo, antn, reset_indel_genes=True)
            except:
                elines = traceback.format_exception(*sys.exc_info())
                print(utils.pad_lines(''.join(elines)))
                print('      couldn\'t convert sw annotation (implicit info adding failed, see previous lines): %s' % ':'.join(cluster))
                n_failed += 1
                continue
            antn_dict[':'.join(cluster)] = antn
        if n_failed > 0:
            print('    %s failed to convert %d/%d sw annotations' % (utils.wrnstr(), n_failed, len(partition)))
        self.glfo = self.sw_glfo  # this is ugly, but we do need to replace it (i'm just a bit worried about doing it so bluntly, but otoh the idea is that we don't really even use these annotations for anything except passing around the pair info)
        return antn_dict, set()  # maybe empty set for hmm failures makes sense since we're not actually running hmm?

    # ----------------------------------------------------------------------------------------
    def actually_get_annotations_for_clusters(self, clusters_to_annotate=None, n_procs=None, dont_print_annotations=False):  # yeah yeah this name sucks
        if self.args.use_sw_annotations: # or self.args.naive_vsearch:
            print('    using sw annotations (to make a multi-seq annotation for each cluster) instead of running hmm since either --use-sw-annotations or --naive-vsearch/--fast were set')
            all_annotations, hmm_failures = self.convert_sw_annotations(partition=clusters_to_annotate)
        else:
            if self.args.naive_vsearch:
                print('      calculating hmm annotations even though --fast/--naive-vsearch was set (you can set --use-sw-annotations to get faster (but substantially less accurate) annotations, or set --dont-calculate-annotations if you only need the partition)')
            _, all_annotations, hmm_failures = self.run_hmm('viterbi', self.sub_param_dir, count_parameters=self.args.count_parameters, parameter_out_dir=self.final_multi_paramdir, partition=clusters_to_annotate, n_procs=n_procs, dont_print_annotations=dont_print_annotations)
        return all_annotations, hmm_failures

    # ----------------------------------------------------------------------------------------
    def get_annotations_for_partitions(self, cpath):  # we used to try to have glomerator.cc try to guess what annoations to write (using ::WriteAnnotations()), but now we go back and rerun a separate bcrham process just to get the annotations we want, partly for that control, but also partly because we typically want all these annotations calculated without any uid translation.
        # ----------------------------------------------------------------------------------------
        def get_clusters_to_annotate():
            clusters_to_annotate = cpath.best()
            additional_clusters = set()  # note that we can in general add clusters to <additional_clusters> that also occur in the best partition (but then don't add duplicates to <clusters_to_annotate>)
            if self.args.write_additional_cluster_annotations is not None:
                istart = max(0, cpath.i_best - self.args.write_additional_cluster_annotations[0])
                istop = min(len(cpath.partitions), cpath.i_best + 1 + self.args.write_additional_cluster_annotations[1])
                for ip in range(istart, istop):
                    additional_clusters |= set([tuple(c) for c in cpath.partitions[ip]])
            if (self.args.calculate_alternative_annotations and self.args.subcluster_annotation_size is None):  # add every cluster from the entire clustering history NOTE stuff that's in self.hmm_cachefname (which we used to use for this), but not in the cpath progress files: translated clusters, clusters that we calculated but didn't merge (maybe? not sure)
                additional_clusters |= set([(uid,) for cluster in cpath.best() for uid in cluster])  # add the singletons separately, since we don't write a singleton partition before collapsing naive sequences before the first clustering step
                additional_clusters |= set([tuple(cluster) for partition in cpath.partitions for cluster in partition])  # kind of wasteful to re-add clusters from the best partition here, but oh well
            if self.args.n_final_clusters is not None or self.args.min_largest_cluster_size is not None:  # add the clusters from the last partition
                additional_clusters |= set([tuple(c) for c in cpath.last()])
            if len(additional_clusters) > 0 and any(list(c) not in clusters_to_annotate for c in additional_clusters):
                cluster_set = set([tuple(c) for c in clusters_to_annotate]) | additional_clusters
                clusters_to_annotate = [list(c) for c in cluster_set]
                actual_new_clusters = [c for c in clusters_to_annotate if c not in cpath.best()]
                self.added_extra_clusters_to_annotate = True
                print('    added %d cluster%s with size%s %s (in addition to the %d from the best partition) before running cluster annotations' % (len(actual_new_clusters), utils.plural(len(actual_new_clusters)), utils.plural(len(actual_new_clusters)),
                                                                                                                                                    ' '.join(str(len(c)) for c in actual_new_clusters), len(cpath.best())))
                if self.args.debug:
                    print('       %s these additional clusters will also be printed below, since --debug is greater than 0' % utils.color('yellow', 'note:'))
            return clusters_to_annotate

        # ----------------------------------------------------------------------------------------
        if self.args.dont_calculate_annotations:
            print('  not calculating annotations')
            if self.args.outfname is not None:
                self.write_output([], set(), cpath=cpath)
            if self.args.debug:
                self.print_results(cpath, [])  # used to just return, but now i want to at least see the cpath
            return
        self.added_extra_clusters_to_annotate = False
        clusters_to_annotate = get_clusters_to_annotate()
        if len(clusters_to_annotate) == 0:
            print('  no final clusters')
            return
        action_cache = self.current_action  # hackey, but probably not worth trying (more) to improve
        self.current_action = 'annotate'
        clusters_to_annotate = sorted(clusters_to_annotate, key=len, reverse=True)  # as opposed to in clusterpath, where we *don't* want to sort by size, it's nicer to have them sorted by size here, since then as you're scanning down a long list of cluster annotations you know once you get to the singletons you won't be missing something big
        n_procs = min(self.args.n_procs, len(clusters_to_annotate))  # we want as many procs as possible, since the large clusters can take a long time (depending on if we're translating...), but in general we treat <self.args.n_procs> as the maximum allowable number of processes
        print('getting annotations for final partition%s%s' % (' (including additional clusters)' if len(clusters_to_annotate) > len(cpath.best()) else '', ' (with star tree annotation since --subcluster-annotation-size is None)' if self.args.subcluster_annotation_size is None else ''))
        all_annotations, hmm_failures = self.actually_get_annotations_for_clusters(clusters_to_annotate=clusters_to_annotate, n_procs=n_procs, dont_print_annotations=True)  # have to print annotations below so we can also print the cpath
        if self.args.get_selection_metrics:
            self.calc_tree_metrics(all_annotations, cpath=cpath)  # adds tree metrics to <annotations>

        if self.args.calculate_alternative_annotations:
            for cluster in sorted(cpath.best(), key=len, reverse=True):
                if len(set(cluster) & hmm_failures) == len(cluster):
                    continue
                if len(cluster) < 5:
                    continue
                self.process_alternative_annotations(cluster, all_annotations, cpath=cpath, debug=self.args.debug)  # NOTE modifies the annotations (adds 'alternative-annotations' key)

        if self.args.outfname is not None:
            self.write_output(list(all_annotations.values()), hmm_failures, cpath=cpath)

        if self.args.count_parameters and not self.args.dont_write_parameters:  # not sure this is absolutely the most sensible place to put this, but I'm trying to kind of mimic where we write the hmms in self.cache_parameters()
            self.write_hmms(self.final_multi_paramdir)  # note that this modifies <self.glfo>

        self.current_action = action_cache

        if self.args.plotdir is not None and not self.args.no_partition_plots:
            ptn_to_use, annotation_list = self.restrict_ex_out_clusters(cpath, list(all_annotations.values()))
            partplotter = PartitionPlotter(self.args, glfo=self.glfo)
            partplotter.plot(self.args.plotdir + '/partitions', ptn_to_use, utils.get_annotation_dict(annotation_list), reco_info=self.reco_info, args=self.args)

        if self.args.seed_unique_id is not None:
            cpath.print_seed_cluster_size(queries_to_include=self.args.queries_to_include)

        if self.args.debug:
            print('final')
            self.print_results(cpath, list(all_annotations.values()))

    # ----------------------------------------------------------------------------------------
    def get_cached_hmm_naive_seqs(self, queries=None):
        if queries is not None:
            expected_queries = queries
        elif self.input_partition is not None:
            expected_queries = [':'.join(c) for c in self.input_partition]
        else:
            expected_queries = self.sw_info['queries']

        cached_naive_seqs = {}
        with open(self.hmm_cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                if ':' in line['unique_ids'] and self.input_partition is None:  # if it's a cache file left over from a previous partitioning, there'll be clusters in it, too
                    continue
                if self.args.persistent_cachefname is not None and line['unique_ids'] not in expected_queries:  # probably can only happen if self.args.persistent_cachefname is set, and it's slow on huge samples
                    continue
                cached_naive_seqs[line['unique_ids']] = line['naive_seq']
                if len(cached_naive_seqs) == len(expected_queries):  # already got everybody
                    break

        if set(cached_naive_seqs) != set(expected_queries):  # can happen if hmm can't find a path for a sequence for which sw *did* have an annotation (but in that case the annotation is almost certainly garbage)
            extra = set(cached_naive_seqs) - set(expected_queries)
            missing = set(expected_queries) - set(cached_naive_seqs)
            if len(extra) > 0:
                print('    %s read %d extra queries from hmm cache file %s' % (utils.color('yellow', 'warning:'), len(extra), ' '.join(extra)))
            if len(missing) > 0:
                print('    %s missing %d/%d queries from hmm cache file (using sw naive sequence instead): %s' % (utils.color('yellow', 'warning:'), len(missing), len(expected_queries), ' '.join(missing)))
                for ustr in missing:
                    cached_naive_seqs[ustr] = self.sw_info[ustr.split(':')[0]]['naive_seq']

        return cached_naive_seqs

    # ----------------------------------------------------------------------------------------
    def cluster_with_naive_vsearch_or_swarm(self, parameter_dir=None):
        start = time.time()

        naive_seq_list = []
        assert parameter_dir is not None
        threshold = self.get_hfrac_bounds(parameter_dir)[0]  # lo and hi are the same
        cached_naive_seqs = self.get_cached_hmm_naive_seqs()
        if self.input_partition is not None and self.args.continue_from_input_partition:
            tclusters = self.input_partition
            cdr3_info = {':'.join(c) : {'cdr3_length' : self.input_antn_dict[':'.join(c)]['cdr3_length']} for c in tclusters}  # have to make a 'fake' sw info to pass to the naive seq collapse fcn
            print('      --continue-from-input-partition: using input partition clusters to initialize vsearch')
        else:
            tclusters = [[u] for u in self.sw_info['queries']]
            cdr3_info = self.sw_info
        for tclust in tclusters:
            tkey = ':'.join(tclust)
            if tkey not in cached_naive_seqs:
                raise Exception('naive sequence for %s not found in %s' % (tkey, self.hmm_cachefname))
            tnseq = cached_naive_seqs[tkey]
            if tnseq == '':  # ugh, not sure why this is happening, but it seems to only be from a persistent cache file from a previous step from test.py, so maybe not important
                print('    %s empty naive seq in hmm cache for \'%s\', using sw info' % (utils.wrnstr(), tkey))
                tnseq = self.sw_info[tkey]['naive_seq']
            naive_seq_list.append((tkey, tnseq))

        nseq_map, nseq_hashes = utils.collapse_naive_seqs_with_hashes(naive_seq_list, cdr3_info)

        print('    using hfrac bound for vsearch %.3f' % threshold)

        partition = []
        print('    running vsearch %d times (once for each cdr3 length class):' % len(nseq_map), end=' ')
        for cdr3_length, sub_naive_seqs in nseq_map.items():
            if len(sub_naive_seqs) == 1:  # no need to run clustering for one sequence
                sub_hash_partition = [[list(sub_naive_seqs.keys())[0]]]
            else:
                sub_hash_partition = utils.run_vsearch('cluster', sub_naive_seqs, self.args.workdir + '/vsearch', threshold, vsearch_binary=self.args.vsearch_binary)
            sub_uid_partition = [[uid for hashstr in hashcluster for ustr in nseq_hashes[cdr3_length][hashstr] for uid in ustr.split(':')] for hashcluster in sub_hash_partition]
            partition += sub_uid_partition
            print('.', end=' ')
            sys.stdout.flush()
        print('')

        partition = utils.split_clusters_by_cdr3(partition, self.sw_info, warn=True)  # if we didn't use sw cdr3 lengths to start with here, there can be clusters with multiple sw cdr3 lengths, which then break when we go to get annotations later (it'd be too hard to avoid using sw info when getting annotations). And we don't really care which is right -- in general if we can get different cdr3 lengths for a sequence, it's absolute garbage to begin with
        # should really also re-pad after this i think (at least i'm getting naive seq len errors when merging paired partitions cause some look un/wrongly padded)
        # utils.re_pad_hmm_seqs(self.input_antn_list, self.input_glfo, self.sw_info)  # can't do this in the self init fcn since at that point we haven't yet read the sw cache file

        perf_metrics = None
        if not self.args.is_data:  # i'm no longer calculating ccfs here, so could remove/move some of this inside the pairwise clustering metric block
            queries_without_annotations = set(self.input_info) - set(self.sw_info['queries'])
            tmp_partition = copy.deepcopy(partition) + [[q, ] for q in queries_without_annotations]  # just add the missing ones as singletons
            self.check_partition(tmp_partition)
            if self.args.add_pairwise_clustering_metrics:
                true_partition = utils.get_partition_from_reco_info(self.reco_info)
                perf_metrics = {'pairwise' : utils.pairwise_cluster_metrics('pairwise', tmp_partition, true_partition)}
        cpath = ClusterPath(seed_unique_id=self.args.seed_unique_id)
        cpath.add_partition(partition, logprob=0.0, n_procs=1, perf_metrics=perf_metrics)
        if self.args.small_clusters_to_ignore is not None:  # maybe should do this a few lines higher? but this is where we're making the cpath atm
            cpath = self.remove_small_clusters(cpath)

        print('      vsearch time: %.1f' % (time.time()-start))
        sys.stdout.flush()
        return cpath

    # ----------------------------------------------------------------------------------------
    def get_hfrac_bounds(self, parameter_dir):  # parameterize the relationship between mutation frequency and naive sequence inaccuracy
        # ----------------------------------------------------------------------------------------
        def dbgstr():
            if self.args.naive_hamming_cluster:
                return utils.color('blue_bkg', '--naive-hamming-cluster')
            elif self.args.naive_vsearch:  # set lo and hi to the same thing, so we don't use log prob ratios, i.e. merge if less than this, don't merge if greater than this
                return utils.color('blue_bkg', '--naive-vsearch')
            else:
                return ''
        # ----------------------------------------------------------------------------------------
        if self.cached_naive_hamming_bounds is not None:  # only run the stuff below once
            if self.print_status:
                print('      %s naive hfrac bounds: %.3f %.3f' % (dbgstr(), self.cached_naive_hamming_bounds[0], self.cached_naive_hamming_bounds[1]))
            return self.cached_naive_hamming_bounds

        # this is a bit weird and messy because i just split out the guts into a fcn in utils (long after writing it), and i don't want to fiddle with this too much
        if self.args.naive_hamming_cluster:
            pmethod = 'naive-hamming'
        elif self.args.naive_vsearch:  # set lo and hi to the same thing, so we don't use log prob ratios, i.e. merge if less than this, don't merge if greater than this
            pmethod = 'naive-vsearch'
        else:  # these are a bit larger than the tight ones and should almost never merge non-clonal sequences, i.e. they're appropriate for naive hamming preclustering if you're going to run the full likelihood on nearby sequences
            pmethod = 'likelihood'
        self.cached_naive_hamming_bounds = utils.get_naive_hamming_bounds(pmethod, parameter_dir=parameter_dir)
        if self.print_status:
            print('      %s setting naive hfrac bounds: %.3f %.3f' % (dbgstr(), self.cached_naive_hamming_bounds[0], self.cached_naive_hamming_bounds[1]))
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
        cmd_str += ' --random-seed ' + str(self.args.random_seed)
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

                hfrac_bounds = self.get_hfrac_bounds(parameter_dir)
                if self.args.naive_hamming_cluster:  # shouldn't be able to happen, but...
                    assert hfrac_bounds[0] == hfrac_bounds[1]
                cmd_str += ' --hamming-fraction-bound-lo ' + str(hfrac_bounds[0])
                cmd_str += ' --hamming-fraction-bound-hi ' + str(hfrac_bounds[1])
                cmd_str += ' --logprob-ratio-threshold ' + str(self.args.logprob_ratio_threshold)
                cmd_str += ' --biggest-naive-seq-cluster-to-calculate ' + str(self.args.biggest_naive_seq_cluster_to_calculate)
                cmd_str += ' --biggest-logprob-cluster-to-calculate ' + str(self.args.biggest_logprob_cluster_to_calculate)
                cmd_str += ' --n-partitions-to-write ' + str(self.args.n_partitions_to_write)  # don't write too many, since calculating the extra logprobs is kind of expensive (note that in practice bcrham rarely has the default of 10 to work with any more, mostly because it's only looking at merges from one clustering step. Also note that we typically won't have the logprobs for all of them, for the same reason)
                if n_procs == 1:  # if this is the last time through, with one process, we want glomerator.cc to calculate the total logprob of each partition NOTE this is quite expensive, since we have to turn off translation entirely
                    cmd_str += '  --write-logprob-for-each-partition'

                if self.args.seed_unique_id is not None and self.unseeded_seqs is None:  # if we're in the last few cycles (i.e. we've removed unseeded clusters so self.unseeded_seqs is set) we want bcrham to *not* know about the seed (this gives more accurate clustering 'cause it means we're really doing hierarchical agglomeration)
                    cmd_str += ' --seed-unique-id ' + self.args.seed_unique_id

                if n_procs == 1 or self.set_force_args:
                    if self.args.n_final_clusters is not None:
                        cmd_str += ' --n-final-clusters ' + str(self.args.n_final_clusters)
                    if self.args.min_largest_cluster_size is not None:
                        cmd_str += ' --min-largest-cluster-size ' + str(self.args.min_largest_cluster_size)

                if self.args.max_cluster_size is not None:
                    cmd_str += ' --max-cluster-size ' + str(self.args.max_cluster_size)

        cmd_str += ' --ambig-base ' + utils.ambig_base

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
        if not self.print_status:
            return

        dbgstr = []
        pwidth = str(len(str(len(self.input_info))))  # close enough
        if 'read-cache' in utils.bcrham_dbgstrs[actionstr]:
            if sum(summaryfo['read-cache'].values()) == 0:
                dbgstr.append('                no/empty cache file')
            else:
                dbgstr.append(('          read from cache:  naive-seqs %' + pwidth + 'd   logprobs %' + pwidth + 'd') % (summaryfo['read-cache']['naive-seqs'], summaryfo['read-cache']['logprobs']))
        dbgstr.append(('                    calcd:         vtb %' + pwidth + 'd        fwd %' + pwidth + 'd') % (summaryfo['calcd']['vtb'], summaryfo['calcd']['fwd']))
        if 'merged' in utils.bcrham_dbgstrs[actionstr]:
            dbgstr.append(('                   merged:       hfrac %' + pwidth + 'd     lratio %' + pwidth + 'd') % (summaryfo['merged']['hfrac'], summaryfo['merged']['lratio']))

        if len(self.bcrham_proc_info) == 1:
            dbgstr.append('                     time:  %.1f sec' % summaryfo['time']['bcrham'][0])
        else:
            dbgstr.append('             min-max time:  %.1f - %.1f sec' % (summaryfo['time']['bcrham'][0], summaryfo['time']['bcrham'][1]))
        dbgstr = '\n'.join(dbgstr)
        print(dbgstr)

    # ----------------------------------------------------------------------------------------
    def check_wait_times(self, wait_time):
        max_bcrham_time = max([procinfo['time']['bcrham'] for procinfo in self.bcrham_proc_info])
        if max_bcrham_time > 0. and wait_time / float(max_bcrham_time) > 1.5 and wait_time > 30.:  # if we were waiting for a lot longer than the slowest process took, and if it took long enough for us to care
            print('    spent much longer waiting for bcrham (%.1fs) than bcrham reported taking (max per-proc time %.1fs)' % (wait_time, max_bcrham_time))

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

        if self.print_status:
            print('    running %d proc%s' % (n_procs, utils.plural(n_procs)))
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
    def subcl_split(self, n_seqs):  # return true if we want to split cluster of size <n_seqs> into subclusters for purposes of annotation accuracy
        return n_seqs > self.args.subcluster_annotation_size
        # return n_seqs >= 2 * self.args.subcluster_annotation_size  # used to do it this way, and there's pluses and minuses to both, but it turns out it's better to split smaller clusters

    # ----------------------------------------------------------------------------------------
    def run_subcluster_annotate(self, init_partition, parameter_in_dir, count_parameters=False, parameter_out_dir='', dont_print_annotations=False, debug=False):  # NOTE nothing to do with subcluster naive seqs above
        # ----------------------------------------------------------------------------------------
        def skey(c):
            return ':'.join(c)
        def skey_inverse(ustr):
            return ustr.split(':')
        # ----------------------------------------------------------------------------------------
        def hashstr(c):
            return 'subc_%s' % utils.uidhashstr(skey(c))
        # ----------------------------------------------------------------------------------------
        def getsubclusters(superclust, shuffle=False, sdbg=False):  # NOTE similar to code in bin/partis simulation fcn
            superclust = copy.deepcopy(superclust)  # should only need this if we shuffle, but you *really* don't want to modify it since if you change the clusters in <init_partition> you'll screw up the actual partition, and if you change the ones in <clusters_still_to_do> you'll end up losing them cause the uidstrs won't be right
            if shuffle:
                random.shuffle(superclust)
            # n_clusters = len(superclust) // self.args.subcluster_annotation_size  # old way: truncates the decimal (see note in subcl_split())
            n_clusters = int(math.ceil(len(superclust) / float(self.args.subcluster_annotation_size)))  # taking the ceiling keeps the max un-split cluster size equal to <self.args.subcluster_annotation_size> (rather than one less than twice that)
            if n_clusters < 2:  # initially we don't call this function unless superclust is big enough for two clusters of size self.args.subcluster_annotation_size, but on subsequent rounds the clusters fom naive_ancestor_hashes can be smaller than that
                return [superclust]

            if False: # self.args.kmeans_subclusters:  # this gives you clusters that are "tighter" -- i.e. clusters similar sequences together
                from . import mds  # this works fine, but it's not really different (on balance with kmeans is probably a bit worse) than the simple way. Which is weird, I would think it would help? but otoh it gives you non-equal-sized clusters, which sometimes i think is worse, although sometimes i also think is better
                seqfos = [{'name' : u, 'seq' : self.sw_info[u]['seqs'][0]} for u in superclust]
                return_clusts = mds.run_sklearn_mds(None, n_clusters, seqfos, self.args.random_seed, aligned=True)
            else:
                n_seq_list = [len(superclust) // n_clusters for _ in range(n_clusters)]  # start with the min possible number (without remainders) for each cluster
                if sdbg:
                    print('len %3d  sas %2d   r %2d  %s' % (len(superclust), self.args.subcluster_annotation_size, len(superclust) % n_clusters, n_seq_list), end=' ')
                for iextra in range(len(superclust) % n_clusters):  # spread out any extra ones
                    n_seq_list[iextra] += 1
                assert sum(n_seq_list) == len(superclust)
                ns_sums = [sum(n_seq_list[:i]) for i in range(len(n_seq_list))]
                return_clusts = [superclust[nsum : nsum + n] for n, nsum in zip(n_seq_list, ns_sums)]  # could do some cleverer tree-based clustering, but when partitioning they should be ordered by similarity (i.e. the order in which they got hierarchically agglomerated)

            if sdbg:
                print(' '.join([str(len(c)) for c in return_clusts]))
                if sdbg > 1:
                    # ptnprint(return_clusts)
                    for dbclust in return_clusts:
                        print('   len %d mean pw hfrac %.2f' % (len(dbclust), utils.mean_pairwise_hfrac([self.sw_info[u]['seqs'][0] for u in dbclust])))
                        naive_seq = self.reco_info[dbclust[0]]['naive_seq']  # it'd be nice to not have to use reco info, but at this point i think the only other naive seqs we have are the single-sequence sw ones, which won't be accurate enough (i guess i could use the consensus sequence or something)
                        for uid in dbclust:
                            utils.color_mutants(naive_seq, utils.remove_ambiguous_ends(self.sw_info[uid]['seqs'][0]), align_if_necessary=True, print_result=True, only_print_seq=True, extra_str='      ', post_str=' %22s'%uid)
            assert sum(len(c) for c in return_clusts) == len(superclust)
            return return_clusts
            # old simple way:
            # return [superclust[i : i + self.args.subcluster_annotation_size] for i in range(0, len(superclust), self.args.subcluster_annotation_size)]  # could do some cleverer tree-based clustering, but when partitioning they should be ordered by similarity (i.e. the order in which they got hierarchically agglomerated)

        # ----------------------------------------------------------------------------------------
        def add_hash_seq(line, hashid, naive_hash_ids):  # need to add it also to self.input_info and self.sw_info for self.combine_queries()
            subcluster_hash_seqs[hashid] = line['naive_seq']
            uid_to_copy = line['unique_ids'][0]  # we copy from just this one, but then OR (more or less) info from the other ones
            swfo_to_copy = self.sw_info[uid_to_copy]
            swhfo = copy.deepcopy(swfo_to_copy)  # <swhfo> is the new sw info for <hashid>
            swhfo.update({'unique_ids' : [hashid], 'seqs' : [line['naive_seq']], 'input_seqs' : [line['naive_seq']], 'cdr3_length' : line['cdr3_length'], 'naive_seq' : line['naive_seq']})
            assert len(swhfo['all_matches']) == 1
            for tid in line['unique_ids']:  # this duplicates code in combine_queries()
                for region in utils.regions:
                    swhfo['all_matches'][0][region].update(self.sw_info[tid]['all_matches'][0][region])  # this screws up/overwrites the match info (scores, etc.) but all we care about is adding any other genes
                    if region != 'j':  # also OR kbounds
                        kn = 'k_' + region
                        for mfcn in [min, max]:
                            mn = mfcn.__name__
                            swhfo[kn][mn] = mfcn(swhfo[kn][mn], self.sw_info[tid][kn][mn])
            if hashid in self.sw_info:
                # raise Exception('hashid %s already in sw info (i.e. the uids that made the hash were already read: %s)' % (hashid, ':'.join(line['unique_ids'])))
                if not self.added_extra_clusters_to_annotate:
                    print('  %s hashid %s already in sw info (i.e. the uids that made the hash were already read: %s)' % (utils.color('yellow', 'warning'), hashid, ':'.join(line['unique_ids'])))
                return
            self.sw_info[hashid] = swhfo
            if len(set(self.sw_info[u]['cdr3_length'] for u in naive_hash_ids)) > 1:  # the time this happened, it was because sw was still allowing conserved codon deletion, and now it (kind of) isn't so maybe it won't happen any more? ("kind of" because it does actually allow it, but it expands kbounds to let the hmm not delete the codon, and the hmm builder also by default doesn't allow them, so... it shouldn't happen)
                existing_cdr3_lengths = list(set([self.sw_info[u]['cdr3_length'] for u in naive_hash_ids[:-1]]))
                if len(existing_cdr3_lengths) > 1:
                    print('    %s multiple existing cdr3 lengths when adding hash seq in subcluster annotation: %s' % (utils.wrnstr(), existing_cdr3_lengths))
                    return
                # raise Exception('added hash seq %s with cdr3 len %d to list of len %d all with cdr3 of %d' % (hashid, self.sw_info[hashid]['cdr3_length'], len(naive_hash_ids), utils.get_single_entry(existing_cdr3_lengths)))
                print('  %s added hash seq %s with cdr3 len %d to list of len %d all with cdr3 of %d' % (utils.color('yellow', 'warning'), hashid, self.sw_info[hashid]['cdr3_length'], len(naive_hash_ids), utils.get_single_entry(existing_cdr3_lengths)))
                self.sw_info[hashid]['cdr3_length'] = utils.get_single_entry(existing_cdr3_lengths)  # I'm not totally sure that this fix makes sense

            self.input_info[hashid] = copy.deepcopy(self.input_info[uid_to_copy])  # i think i only need this for input meta info
            self.input_info[hashid].update({'unique_ids' : [hashid], 'seqs' : [line['naive_seq']]})

            if self.reco_info is not None: # and self.args_correct_multi_hmm_boundaries:  # atm this only gets used when correcting multi hmm boundaries, which we want to immediately stop doing as soon as this fcn is working, but oh well, it's might be nice to have the hash/naive seqs in reco info
                self.reco_info[hashid] = copy.deepcopy(self.reco_info[uid_to_copy])
                utils.replace_seqs_in_line(self.reco_info[hashid], [{'name' : hashid, 'seq' : line['naive_seq']}], self.simglfo, try_to_fix_padding=True)  # i wish i could set refuse_to_align=True since it's slow, but sometimes it needs to align (when j length to right of tryp differs)

        # ----------------------------------------------------------------------------------------
        subc_start = time.time()
        final_annotations = {}
        all_hmm_failures = set()
        subcluster_hash_seqs = {}  # all hash-named naive seqs, i.e. that we only made as intermediate steps, but don't care about afterwards (we keep track here just so we can remove them from input sw, and reco info afterwards)

        istep = 0
        print('  subcluster annotating %d cluster%s%s: %s%s' % (len(init_partition), utils.plural(len(init_partition)), '' if self.print_status else ' with steps', '' if not debug else ' '.join(utils.color('blue' if self.subcl_split(len(c)) else None, str(len(c))) for c in init_partition), '\n' if self.print_status else ''), end=' ')
        sys.stdout.flush()
        clusters_still_to_do = [copy.deepcopy(c) for c in init_partition]
        subd_clusters = {}  # keeps track of all the extra info for clusters that we actually had to subcluster: for each such cluster, stores a list where each entry is the subclusters for that round (i.e. the first entry has subclusters composed of the actual seqs in the cluster, and after that it's intermediate naives/hashid seqs)
        naive_ancestor_hashes = {}  # list of inferred naive "intermediate" (hashid) seqs, for each subclustered cluster, that we'll run on in the next step (if there is a next step)
        while len(clusters_still_to_do) > 0:
            if debug:
                print('   %s %d cluster%s left with size%s %s' % (utils.color('green', 'istep %d:'%istep), len(clusters_still_to_do), utils.plural(len(clusters_still_to_do)), utils.plural(len(clusters_still_to_do)), ' '.join(str(len(c)) for c in clusters_still_to_do)))
                print('        init     N prior  current    N     mean pw   subcluster')  # we call it "mean pw hfrac", but it's pairwise *within* each cluster, then the weighted average over clusters (i.e. the smaller it is the tighter/better the clusters are)
                print('        size     splits    size     subcl   hfrac     sizes')
            clusters_to_run = []
            for tclust in clusters_still_to_do:
                if not self.subcl_split(len(tclust)):  # if <tclust> is small enough we don't need to split it up
                    if tclust not in clusters_to_run:  # it should only be possible to have it already in there if we're adding <additional_clusters> from get_annotations_for_partitions() e.g. from --n-final-clusters or something
                        clusters_to_run.append(tclust)
                    if debug:
                        print('       %4d                   ' % len(tclust))
                else:
                    if skey(tclust) in subd_clusters:  # subsequent rounds: get subclusters from intermediate inferred naives (hashids) from the last round
                        subclusters = getsubclusters(naive_ancestor_hashes[skey(tclust)])
                        del naive_ancestor_hashes[skey(tclust)]
                    else:  # first time through: add <tclust> to subd_clusters and get initial subclusters
                        subd_clusters[skey(tclust)] = []
                        subclusters = getsubclusters(tclust, shuffle=self.reco_info is not None and istep==0)  # the order of uids in each cluster that comes out of simulation is the order of leaves in the tree (i.e. they're sorted by similarity), so it's *extremely* important to shuffle them to get a fair comparison
                    subd_clusters[skey(tclust)].append(subclusters)
                    clusters_to_run += [c for c in subclusters if c not in clusters_to_run]  # it should only be possible to have a <c> already in there if we're adding <additional_clusters> from get_annotations_for_partitions() e.g. from --n-final-clusters or something
                    if debug:
                        n_prev = len(subd_clusters[skey(tclust)]) - 1
                        mphfracs = [utils.mean_pairwise_hfrac([self.sw_info[u]['seqs'][0] for u in c]) for c in subclusters]
                        mean_weighted_pw_hfrac = numpy.average(mphfracs, weights=[len(c) / float(len(tclust)) for c in subclusters])
                        print('       %4d      %s      %s      %3d      %4.2f      %s' % (len(tclust), '   ' if n_prev==0 else '%3d'%n_prev, '   ' if n_prev==0 else '%3d'%sum(len(c) for c in subclusters), len(subclusters), mean_weighted_pw_hfrac, ' '.join(str(len(c)) for c in subclusters)))
            _, step_antns, step_failures = self.run_hmm('viterbi', parameter_in_dir, partition=clusters_to_run, is_subcluster_recursed=True)  # is_subcluster_recursed is really just a speed optimization so it doesn't have to check the length of every cluster
            if istep == 0 and self.args.calculate_alternative_annotations:  # could do this in one of the loops below, but it's nice to have it separate since it doesn't really have anything to do with subcluster annotation
                for tline in step_antns.values():
                    final_annotations[skey(tline['unique_ids'])] = tline
            # make sure all missing cluster are accounted for in <step_failures>
            missing_clusters = [c for c in clusters_to_run if skey(c) not in step_antns]
            if len(missing_clusters) > 0:
                print('    missing %d annotations' % len(missing_clusters))
            for mclust in missing_clusters:
                if len(set(mclust) - step_failures) > 0:
                    print('  %s cluster missing from output with uids not in hmm failures: %s\n    missing uids: %s' % (utils.wrnstr(), skey(mclust), ' '.join(set(mclust) - step_failures)))
                if mclust in clusters_still_to_do:  # i.e. if it's a simple/whole cluster (one way this happened was when the single-seq cluster had failed when caching hmm naive seqs)
                    print('      removing missing cluster from clusters_still_to_do: %s' % skey(mclust))
                    clusters_still_to_do.remove(mclust)
                    all_hmm_failures |= set(mclust)
            all_hmm_failures |= step_failures
            # process each annotation that we just got back from bcrham: store inferred naive/hash seq if it's a subcluster, or add to final annotations if it's a simple/whole cluster
            n_hashed, n_whole_finished = 0, 0
            failed_super_clusters = []
            for super_uidstr, subcluster_lists in subd_clusters.items():  # we used to loop over <step_antns>, which was cleaner in some ways, but then it was hard to make sure the hashid seqs stayed in the same order, which is important if we're not using kmeans
                if any(sclust in missing_clusters for sclust in subcluster_lists[-1]):
                    failed_super_clusters.append(super_uidstr)
                    continue
                for sclust in subcluster_lists[-1]:
                    sline = step_antns[skey(sclust)]
                    hashid = hashstr(sline['unique_ids'])
                    if super_uidstr not in naive_ancestor_hashes:
                        naive_ancestor_hashes[super_uidstr] = []
                    naive_ancestor_hashes[super_uidstr].append(hashid)  # add it to the list of inferred naives that we need to run next time through (if there is a next time)
                    add_hash_seq(sline, hashid, naive_ancestor_hashes[super_uidstr])  # add it to stuff so it can get run on
                    n_hashed += 1
            for fsclust in failed_super_clusters:
                print('    giving up on size %d cluster with failed seqs (was just split into %d subclusters): %s' % (len(skey_inverse(fsclust)), len(subd_clusters[fsclust][-1]), fsclust))
                all_hmm_failures |= set(skey_inverse(fsclust))
                clusters_still_to_do.remove(skey_inverse(fsclust))
                del subd_clusters[fsclust]
            for sclust in [c for c in clusters_still_to_do if skey(c) in step_antns]:  # non-subclustered/simple/whole clusters (only happens first time through) [there's simpler ways to get these, but we need to account for failures]
                final_annotations[skey(sclust)] = step_antns[skey(sclust)]
                clusters_still_to_do.remove(sclust)
                n_whole_finished += 1
            # finalize any subclustered clusters that're finished (i.e. that only had one subcluster this time through)
            n_sub_finished = 0
            for uidstr, subcluster_lists in list(subd_clusters.items()):  # handle the ones that're done (at this point they should be an annotation of about length self.args.subcluster_annotation_size consisting of just inferred subcluster naives)
                if len(subcluster_lists[-1]) > 1:  # not finished yet
                    continue
                sclust = skey_inverse(uidstr)
                if debug:
                    print('  %s cluster with original size %d and split history: %s' % (utils.color('blue', 'finishing'), len(sclust), '   '.join(' '.join(str(len(c)) for c in sclist) for sclist in subcluster_lists)))
                line = step_antns[skey(subcluster_lists[-1][0])]
                sfos_to_add = [{'name' : u, 'seq' : self.sw_info[u]['seqs'][0]} for u in sclust]  # original "leaf" seqs that we actually care about (i.e. not inferred hashid naive seqs)
                utils.replace_seqs_in_line(line, sfos_to_add, self.glfo, try_to_fix_padding=True, refuse_to_align=True)  # i think aligning should be unnecessary here, and it's really slow so we don't want to do it by accident (but trimming Ns off the ends seems pretty harmless and it might be enough)
                self.add_per_seq_sw_info(line)
                final_annotations[uidstr] = line
                clusters_still_to_do.remove(sclust)
                del subd_clusters[uidstr]
                n_sub_finished += 1
            if self.print_status:
                print('    read %d new subcluster annotation%s: added hashid %d   whole finished %d   subcl finished %d' % (len(step_antns), utils.plural(len(step_antns)), n_hashed, n_whole_finished, n_sub_finished))
            istep += 1

        for uid in subcluster_hash_seqs:
            if uid in self.sw_info[uid]:
                del self.sw_info[uid]
            if uid in self.input_info:
                del self.input_info[uid]
            if self.reco_info is not None and uid in self.reco_info:
                del self.reco_info[uid]

        annotation_list = [final_annotations[skey(c)] for c in init_partition if skey(c) in final_annotations]  # the missing ones should all be in all_hmm_failures
        if self.args.calculate_alternative_annotations:
            print('    --calculate-alternative-annotations: adding annotations for initial subcluster annotation step (each with size ~%d)' % self.args.subcluster_annotation_size)
            annotation_list += [l for l in final_annotations.values() if l not in annotation_list]
        self.process_annotation_output(annotation_list, all_hmm_failures, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir, print_annotations=self.args.debug and not dont_print_annotations)
        print('%s    subcluster annotation time %.1f' % ('' if self.print_status else '\n', time.time() - subc_start))

        return None, OrderedDict([(skey(l['unique_ids']), l) for l in annotation_list]), all_hmm_failures

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, algorithm, parameter_in_dir, parameter_out_dir='', count_parameters=False, n_procs=None, precache_all_naive_seqs=False, partition=None, shuffle_input=False, is_subcluster_recursed=False, dont_print_annotations=False):
        """ 
        Run bcrham, possibly with many processes, and parse and interpret the output.
        NOTE the local <n_procs>, which overrides the one from <self.args>
        """
        start = time.time()
        if len(self.sw_info['queries']) == 0:
            print('  %s no input queries for hmm' % utils.color('red', 'warning'))
            return

        nsets = self.get_nsets(algorithm, partition)
        if n_procs is None:
            n_procs = self.args.n_procs
            if len(nsets) < n_procs:
                # print '  note: reducing N procs to the number of nsets %d --> %d' % (n_procs, len(nsets))
                n_procs = len(nsets)

        if self.args.subcluster_annotation_size is not None and algorithm == 'viterbi' and not is_subcluster_recursed and not precache_all_naive_seqs and any(self.subcl_split(len(c)) for c in nsets):
            assert not shuffle_input
            return self.run_subcluster_annotate(nsets, parameter_in_dir, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir, dont_print_annotations=dont_print_annotations, debug=self.args.debug)

        self.write_to_single_input_file(self.hmm_infname, nsets, parameter_in_dir, shuffle_input=shuffle_input)  # single file gets split up later if we've got more than one process
        glutils.write_glfo(self.my_gldir, self.glfo)
        if self.print_status and time.time() - start > 0.1:
            print('        hmm prep time: %.1f' % (time.time() - start))

        cmd_str = self.get_hmm_cmd_str(algorithm, self.hmm_infname, self.hmm_outfname, parameter_dir=parameter_in_dir, precache_all_naive_seqs=precache_all_naive_seqs, n_procs=n_procs)

        if n_procs > 1:
            self.split_input(n_procs, self.hmm_infname)

        exec_start = time.time()
        self.execute(cmd_str, n_procs)
        exec_time = time.time() - exec_start

        glutils.remove_glfo_files(self.my_gldir, self.args.locus)

        cpath, annotations, hmm_failures = None, None, None
        if self.current_action == 'partition' or n_procs > 1:
            cpath = self.merge_all_hmm_outputs(n_procs, precache_all_naive_seqs)
            if cpath is not None:
                cpath.write(self.get_cpath_progress_fname(self.istep), self.args.is_data, reco_info=self.reco_info, true_partition=utils.get_partition_from_reco_info(self.reco_info) if not self.args.is_data else None)

        if algorithm == 'viterbi' and not precache_all_naive_seqs:
            annotations, hmm_failures = self.read_annotation_output(self.hmm_outfname, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir, print_annotations=self.args.debug and not dont_print_annotations, is_subcluster_recursed=is_subcluster_recursed)

        if os.path.exists(self.hmm_infname):
            os.remove(self.hmm_infname)

        step_time = time.time() - start
        if self.print_status or parameter_out_dir != '':
            if step_time - exec_time > 0.1:
                print('         infra time: %.1f' % (step_time - exec_time))  # i.e. time for non-executing, infrastructure time
            print('      hmm step time: %.1f' % step_time)
        else:
            print('(%.1fs)%s' % (step_time, '' if is_subcluster_recursed else '\n'), end=' ')
            sys.stdout.flush()
        self.timing_info.append({'exec' : exec_time, 'total' : step_time})  # NOTE in general, includes pre-cache step

        return cpath, annotations, hmm_failures

    # ----------------------------------------------------------------------------------------
    # compare annotations for all clusters with overlap with uids_of_interest (i.e., all the "alternative" annotations for uids_of_interest for which there's info in cluster_annotations)
    def process_alternative_annotations(self, uids_of_interest, cluster_annotations, cpath=None, debug=False):
        # NOTE that when seed partitioning, in the steps before we throw out non-seeded clusters, there's *tons* of very overlapping clusters, since I think we take the the biggest seeded cluster, and pass that to all the subprocs, so then if a different sequence gets added to it in each subproc you end up with lots of different versions (including potentially lots of different orderings of the exact same cluster)
        # ----------------------------------------------------------------------------------------
        def print_naive_line(other_genes, uid_str_list, naive_seq, n_independent_seqs, max_len_other_gene_str=20):
            gene_strs_of_interest = {r : utils.color_gene(genes_of_interest[r]) for r in utils.regions}  # have to do this up here so we know the width before we start looping
            # gswidth = str(utils.len_excluding_colors('  '.join(gene_strs_of_interest.values())))  # out of order, but doesn't matter since we only want the length
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

            pre_str, post_str = '', ''
            if uidstr_of_interest in uid_str_list:
                pre_str = utils.color('blue', '-->', width=5)
                post_str = utils.color('blue', ' <-- requested uids', width=5)
            true_str = ''
            if print_simu:
                true_naive_seqs = list(set([self.reco_info[u]['naive_seq'] for u in uids_of_interest]))
                if len(true_naive_seqs) > 1:  # the true_str will be the wrong width if there's ever actually more than one, but whatever
                    print('    %s multiple true naive seqs for uids in cluster of interest (they\'re probably not really clonal)' % utils.color('yellow', 'warning'))
                if any(len(ts) < len(naive_seq) for ts in true_naive_seqs):  # try to apply padding to left side of naive seqs (this is a pretty shitty way of doing this)
                    true_naive_seqs = [utils.ambig_base * max(0, (len(naive_seq) - len(ts))) + ts for ts in true_naive_seqs]
                true_hdists = [utils.hamming_distance(ts, naive_seq, align_if_necessary=True) for ts in true_naive_seqs]  # , align_if_necessary=True  # something's wrong if we need to align, but it probably just means they're not clonal *and* there's a clusterfuck of shm indels, so whatever
                def tcolor(d):
                    if d == 0: return 'blue_bkg'
                    elif d < 4: return None
                    elif d < 7: return 'yellow'
                    else: return 'red'
                true_str = ' '.join(utils.color(tcolor(d), str(d), width=2) for d in true_hdists)
                true_str = ' %s      ' % true_str
            if self.args.print_all_annotations:
                print('  printing annotations for all clusters supporting naive seq with fraction %.3f:' % final_info['naive-seqs'][naive_seq])
            csize_str = utils.cluster_size_str(uid_str_list, split_strs=True, clusters_to_emph=[uidstr_of_interest], short=True)
            print(('%5s %s  %s%4d  %4.2f     %s     %s    %s%s') % (pre_str, utils.color_mutants(line_of_interest['naive_seq'], naive_seq), true_str, n_independent_seqs, final_info['naive-seqs'][naive_seq], gene_str, other_gene_str, csize_str, post_str))
            if self.args.print_all_annotations:
                for uidstr in sorted(uid_str_list, key=lambda x: x.count(':'), reverse=True):
                    if uidstr in cluster_annotations:
                        utils.print_reco_event(cluster_annotations[uidstr], extra_str='      ', label='')
                    else:
                        print('    %s missing from cluster annotations' % uidstr)
                print('\n')
        # ----------------------------------------------------------------------------------------
        def print_header():
            print('  alternative annotations for cluster with %d seqs:\n    %s\n   (in %s below)' % (uidstr_of_interest.count(':') + 1, uidstr_of_interest, utils.color('blue', 'blue')))
            print('')
            utils.print_reco_event(utils.synthesize_single_seq_line(line_of_interest, iseq=0), extra_str='      ', label='annotation for a single (arbitrary) sequence from the cluster:')
            print('')
            print('')
            print('%s       %s   unique       inconsistent (%s: missing info)     cluster' % (' ' * len(line_of_interest['naive_seq']), 'hdist ' if print_simu else '', utils.color('blue', 'x')))  # 'missing info' means that we didn\'t have an annotation for at least one of the clusters in that line. This is probably because we're using an old hmm cache file, and at best passed in the sw cache file annotations, which are only single-sequence.
            print('%s       %s seqs  frac     regions         genes               sizes (+singletons)' % (' ' * len(line_of_interest['naive_seq']), 'to true ' if print_simu else ''))
        # ----------------------------------------------------------------------------------------
        def init_dicts(region, gene_call):  # ick
            unique_seqs_for_each_gene[region][gene_call] = set()
            cluster_sizes_for_each_gene[region][gene_call] = []
            uidstrs_for_each_gene[region][gene_call] = []
        # ----------------------------------------------------------------------------------------
        def print_gene_call_summary():
            print('')
            print('                        unique')
            print('                      seqs  frac          cluster sizes (+singletons)')
            for region in utils.regions:
                print('    %s' % utils.color('green', region))
                for gene_call, uids in sorted(list(unique_seqs_for_each_gene[region].items()), key=lambda s: len(s[1]), reverse=True):
                    csizes = cluster_sizes_for_each_gene[region][gene_call]
                    if self.args.print_all_annotations:
                        print('  printing annotations for all clusters supprting gene call with fraction %.3f:' % final_info['gene-calls'][region][gene_call])
                    print('      %s %4d  %4.2f             %s' % (utils.color_gene(gene_call, width=15), len(uids), final_info['gene-calls'][region][gene_call], utils.cluster_size_str(csizes, only_passing_lengths=True, short=True)))
                    if self.args.print_all_annotations:
                        for uidstr in sorted(uidstrs_for_each_gene[region][gene_call], key=lambda x: x.count(':'), reverse=True):
                            if uidstr in cluster_annotations:
                                utils.print_reco_event(cluster_annotations[uidstr], extra_str='      ', label='')
                            else:
                                print('    %s missing from cluster annotations' % uidstr)
                        print('\n')
        # ----------------------------------------------------------------------------------------
        print_simu = self.args.infname is not None and self.reco_info is not None
        if self.args.seed_unique_id is not None and debug and self.args.seed_unique_id not in uids_of_interest:
            print('  note: only printing alternative annotations for seed clusters (the rest are still written to disk, but if you want to print the others, don\'t set --seed-unique-id)')
            debug = False

        # <clusters_in_partitions> is by default None (it's only used for old-style, i.e. non-subcluster-annotation, alternative annotations) (cpath should actually always be set, but i want it a kw arg to make clear that it's basically optional)
        clusters_in_partitions = None if (self.args.subcluster_annotation_size is not None or cpath is None) else set([':'.join(c) for partition in cpath.partitions for c in partition])  # this is a list of the clusters that actually occur in one of the partitions in the clusterpath. It's so we can ignore annotations for clusters that were never actually formed, i.e. where we got the annotation for a potential cluster, but decided not to make the merge (at least, I think this is why some annotations don't correspond to clusters in the cluster path)
        uids_of_interest = set(uids_of_interest)
        uidstr_of_interest = None  # we don't yet know in what order the uids in <uids_of_interest> appear in the file
        sub_info = {}  # map from naive seq : sub uid strs
        final_info = {'naive-seqs' : OrderedDict(), 'gene-calls' : {r : OrderedDict() for r in utils.regions}}  # this is the only info that's persistent, it get's added to the cluster annotation corresponding to <uids_of_interest>

        # find all the clusters that have any overlap with <uids_of_interest>
        sum_of_used_cluster_sizes = 0
        if debug:
            print('  processing alternative annotations with %s' % ('partition path clusters' if self.args.subcluster_annotation_size is None else 'subcluster annotations'))  # NOTE this has no way of knowing if you're running 'view-alternative-annotations' with different options than what you used to write the file
            print('    looking among %d annotations for overlap with %d uids of interest' % (len(cluster_annotations), len(uids_of_interest)))
        for uidstr, info in cluster_annotations.items():
            if info['naive_seq'] == '':  # hmm cache file lines that only have logprobs UPDATE no longer support reading hmm cache files, but maybe this is also what failed queries look like sometimes?
                continue
            uid_set = set(info['unique_ids'])
            if len(uid_set & uids_of_interest) == 0:
                continue
            if clusters_in_partitions is not None and ':' in uidstr and uidstr not in clusters_in_partitions:  # first bit is because singletons are not in general all in a partition, since the first partition is from after initial naive seq merging (at least seems to be? I kind of thought I added the singleton partition, but I guess not)
                continue
            sum_of_used_cluster_sizes += len(uid_set)
            naive_seq = cluster_annotations[uidstr]['naive_seq']
            if naive_seq not in sub_info:
                sub_info[naive_seq] = []
            sub_info[naive_seq].append(uidstr)
            if uid_set == uids_of_interest:
                uidstr_of_interest = uidstr

        if uidstr_of_interest is None:
            raise Exception('alternative annotations: output file doesn\'t have a cluster with the exact requested uids (run without setting --queries to get a list of available clusters): %s' % uids_of_interest)
        if self.args.subcluster_annotation_size is not None and sum_of_used_cluster_sizes > 2.1 * len(uids_of_interest):  # if the annotations are from subcluster annotation, it should actually be exactly 2 times (since atm we just add each small cluster in the first round of subcluster annotation, plus the final cluster at the end)
            raise Exception('--subcluster-annotation-size is non-None (so we expect to see each uid only ~twice in alternative annotations) but sum of chosen cluster sizes is %d (vs %d uids of interest): you probably wrote the output file with --subcluster-annotation-size turned off (None), but didn\'t set that option when running \'view-alternative-annotations\'' % (sum_of_used_cluster_sizes, len(uids_of_interest)))

        line_of_interest = cluster_annotations[uidstr_of_interest]
        genes_of_interest = {r : line_of_interest[r + '_gene'] for r in utils.regions}
        if debug:
            print_header()

        # loop over each unique naive sequence that was inferred for at least one cluster, counting up various properties (and maybe printing dbg info)
        independent_seq_info = {naive_seq : set([uid for uidstr in uid_str_list for uid in uidstr.split(':')]) for naive_seq, uid_str_list in sub_info.items()}
        total_independent_seqs = sum(len(uids) for uids in independent_seq_info.values())
        unique_seqs_for_each_gene = {r : {} for r in utils.regions}  # for each gene, keeps track of the number of unique sequences that contributed to an annotation that used that gene
        cluster_sizes_for_each_gene = {r : {} for r in utils.regions}  # same, but number/sizes of different clusters (so we can tell if a gene is only supported by like one really large cluster, but all the smaller clusters point to another gene)
        uidstrs_for_each_gene = {r : {} for r in utils.regions}  # need the actual uid strs, but only use them for --print-all-annotations
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

        # divide by totals for gene call info (to get fractions)
        for region in utils.regions:
            total_unique_seqs_this_region = sum(len(uids) for uids in unique_seqs_for_each_gene[region].values())  # yes, it is in general different for each region
            for gene_call, uids in sorted(list(unique_seqs_for_each_gene[region].items()), key=lambda s: len(s[1]), reverse=True):
                final_info['gene-calls'][region][gene_call] = len(uids) / float(total_unique_seqs_this_region)
        if debug:
            print_gene_call_summary()

        # this is messy because we want to acces them as dicts in this fcn, for dbg printing, but want them in the output file as simple combinations of lists/tuples
        line_of_interest['alternative-annotations'] = {'naive-seqs' : list(final_info['naive-seqs'].items()), 'gene-calls' : {r : list(final_info['gene-calls'][r].items()) for r in utils.regions}}

    # ----------------------------------------------------------------------------------------
    def get_padded_true_naive_seq(self, qry):
        assert len(self.sw_info[qry]['padlefts']) == 1
        return self.sw_info[qry]['padlefts'][0] * utils.ambig_base + self.reco_info[qry]['naive_seq'] + self.sw_info[qry]['padrights'][0] * utils.ambig_base

    # ----------------------------------------------------------------------------------------
    def split_input(self, n_procs, infname):

        # should we pull out the seeded clusters, and carefully re-inject them into each process?
        separate_seeded_clusters = self.current_action == 'partition' and self.args.seed_unique_id is not None and self.unseeded_seqs is None  # I think it's no longer possible to have seed_unique_id set if we're not partitioning, but I'll leave it just to be safe (otherwise we get the seed seq sent to every process)

        # read single input file
        info = []
        seeded_clusters = {}
        with open(infname, 'r') as infile:
            reader = csv.DictReader(infile, delimiter=str(' '))
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
            return csv.DictWriter(sub_outfile, reader.fieldnames, delimiter=str(' '))
        def copy_cache_files(n_procs):
            cmdfos = [{'cmd_str' : 'cp ' + self.hmm_cachefname + ' ' + self.subworkdir(iproc, n_procs) + '/',
                       'workdir' : self.subworkdir(iproc, n_procs),
                       'outfname' : self.subworkdir(iproc, n_procs) + '/' + os.path.basename(self.hmm_cachefname)}
                      for iproc in range(n_procs)]
            utils.run_cmds(cmdfos)

        # initialize output/cache files
        for iproc in range(n_procs):
            utils.prep_dir(self.subworkdir(iproc, n_procs))
            sub_outfile = get_sub_outfile(iproc, utils.csv_wmode())
            get_writer(sub_outfile).writeheader()
            sub_outfile.close()  # can't leave 'em all open the whole time 'cause python has the thoroughly unreasonable idea that one oughtn't to have thousands of files open at once
        if self.current_action == 'partition' and os.path.exists(self.hmm_cachefname):  # copy cachefile to this subdir (first clause is just for so when we're getting cluster annotations we don't copy over the cache files)
            copy_cache_files(n_procs)

        seed_clusters_to_write = list(seeded_clusters.keys())  # the keys in <seeded_clusters> that we still need to write
        for iproc in range(n_procs):
            sub_outfile = get_sub_outfile(iproc, utils.csv_wmode('a'))
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
            outfile = open(outfname, utils.csv_wmode())
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
            # print '    nothing to merge into %s' % outfname
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
            pass
            # print '    nothing to merge into %s' % outfname
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
    def merge_all_hmm_outputs(self, n_procs, precache_all_naive_seqs):  # for search: read hmm read_hmm
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
        print('  writing hmms', end=' ')
        sys.stdout.flush()
        start = time.time()

        from .hmmwriter import HmmWriter
        hmm_dir = parameter_dir + '/hmms'
        utils.prep_dir(hmm_dir, '*.yaml')
        # hmglfo = copy.deepcopy(self.glfo)  # it might be better to not modify self.glfo here, but there's way too many potential downstream effects to change it at this point
        glutils.restrict_to_observed_genes(self.glfo, parameter_dir, debug=True)  # this is kind of a weird place to put this... it would make more sense to read the glfo from the parameter dir, but I don't want to mess around with changing that a.t.m.

        if self.args.debug:
            print('to %s' % parameter_dir + '/hmms', end=' ')

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

        print('(%.1f sec)' % (time.time()-start))
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
            uids_and_lengths = sorted(list(uids_and_lengths.items()), key=operator.itemgetter(1))
            uids, lengths = zip(*uids_and_lengths)
            raise Exception('cdr3 lengths not all the same (%s) for %s (probably need to add more criteria for call to utils.split_clusters_by_cdr3())' % (' '.join([str(c) for c in lengths]), ' '.join(uids)))
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
                regmatches = set(swfo['all_matches'][0][region])  # the best <n_max_per_region> matches for this query (note: currently not sorted by score, since now it's a dict keyed by gene that includes the score; *but* when reading old sw cache files it's a list, and it *is* sorted by score. Neither of which matters right here since we're just making a set from all of them, but still)
                genes_to_use |= regmatches & available_genes

            # OR this query's genes into the ones from previous queries
            combo['only_genes'] |= genes_to_use  # NOTE using the OR of all sets of genes (from all query seqs) like this *really* helps,

        combo['only_genes'] = list(combo['only_genes'])  # maybe I don't need to convert it to a list, but it used to be a list, and I don't want to make sure there wasn't a reason for that

        # if we don't have at least one gene for each region, add all available genes from that regions
        if set(utils.regions) > set(utils.get_region(g) for g in combo['only_genes']):  # this is *very* rare
            missing_regions = set(utils.regions) - set(utils.get_region(g) for g in combo['only_genes'])  # this is wasteful, but it hardly *ever* happens
            combo['only_genes'] += [g for g in available_genes if utils.get_region(g) in missing_regions]
            if self.args.debug:
                print('    %s no %s genes for query %s, so added in all the genes from these regions' % ('note:', ' or '.join(missing_regions), query_names))  # these other genes may not work (the kbounds are probably wrong), but they might work, and it's better than just giving up on the query

        for kb in ['k_v', 'k_d']:
            if combo[kb]['min'] <= 0 or combo[kb]['min'] >= combo[kb]['max']:
                raise Exception('nonsense k bounds for %s (v: %d %d  d: %d %d)' % (':'.join(query_names), combo['k_v']['min'], combo['k_v']['max'], combo['k_d']['min'], combo['k_d']['max']))

        return combo

    # ----------------------------------------------------------------------------------------
    def write_bcrham_cache_file(self, nsets, ctype='fake', antn_list=None):
        """
        If <ctype> is 'fake', write a cache file which, instead of the inferred naive sequences, has the *true* naive sequences. Used to generate 'synthetic' partitions (see papers).
        If <ctype> is 'input', use annotations from the input partition to write the cache file.
        """
        if ctype == 'fake':
            # NOTE this will be out of sync with sw info afterwards (e.g. indel info/cdr3 length), so things may break if you do extra things (e.g. subcluster annotate, run with debug set)
            if self.reco_info is None:
                raise Exception('can\'t write fake cache file for --synthetic-distance-based-partition unless --is-simu is specified (and there\'s sim info in the input file)')
            dbgstr = 'true'
            def naive_seq_fcn(query_name_list):
                return self.get_padded_true_naive_seq(query_name_list[0])  # NOTE just using the first one... but a.t.m. I think I'll only run this fcn the first time through when they're all singletons, anyway
        elif ctype == 'input':
            dbgstr = 'input annotation'
            # old version that used sw annotations (don't quite want to delete yet):
            # antn_dict, _ = self.convert_sw_annotations(partition=nsets)  # we *have* the annotations for the input partition, which would be better, but their seqs aren't correctly padded to the same length (since they were run in different subsets/procs), and propagating the new sw padding info (sw gets re-padded when we read the subset-merged sw info) to the input partition annotations would be hard, so we just use the sw annotations to get the naive seqs here
            # def naive_seq_fcn(query_name_list):
            #     return antn_dict.get(':'.join(query_name_list))['naive_seq']
            utils.re_pad_hmm_seqs(self.input_antn_list, self.input_glfo, self.sw_info)  # can't do this in the self init fcn since at that point we haven't yet read the sw cache file
            self.input_antn_dict = utils.get_annotation_dict(self.input_antn_list)  # fcn in previous line only replaces elements in list (rather than modifying those elements) in very rare cases, but for cases when it does, whe have to replace the associated dict
            def naive_seq_fcn(query_name_list):
                if ':'.join(query_name_list) not in self.input_antn_dict:
                    print('    %s annotation for %s not in input annotation dict, so looking for overlapping clusters (something is probably wrong, maybe to do with --queries-to-include duplicate removal)' % (utils.wrnstr(), ':'.join(query_name_list)))
                return self.input_antn_dict.get(':'.join(query_name_list))['naive_seq']
        else:
            assert False
        with open(self.hmm_cachefname, utils.csv_wmode()) as cachefile:
            writer = csv.DictWriter(cachefile, utils.partition_cachefile_headers)
            writer.writeheader()
            for query_name_list in nsets:
                writer.writerow({
                    'unique_ids' : ':'.join([qn for qn in query_name_list]),
                    'naive_seq' : naive_seq_fcn(query_name_list),
                })
        print('    %s %s naive seqs in bcrham cache file for %d clusters (%d seqs)' % (utils.color('blue_bkg', 'caching'), dbgstr, len(nsets), sum(len(c) for c in nsets)))

    # ----------------------------------------------------------------------------------------
    def write_to_single_input_file(self, fname, nsets, parameter_dir, shuffle_input=False):
        csvfile = open(fname, utils.csv_wmode())
        header = ['names', 'k_v_min', 'k_v_max', 'k_d_min', 'k_d_max', 'mut_freq', 'cdr3_length', 'only_genes', 'seqs']
        writer = csv.DictWriter(csvfile, header, delimiter=str(' '))
        writer.writeheader()

        if shuffle_input:  # shuffle nset order (this is absolutely critical when clustering with more than one process, in order to redistribute sequences among the several processes)
            random.shuffle(nsets)

        genes_with_hmm_files = self.get_existing_hmm_files(parameter_dir)
        genes_with_enough_counts = utils.get_genes_with_enough_counts(parameter_dir, self.args.min_allele_prevalence_fractions)  # it would be nice to do this at some earlier step, but then we have to rerun sw (note that there usually won't be any to remove for V, since the same threshold was already applied in alleleremover, but there we can't do d and j since we don't yet have annotations for them)
        glfo_genes = set([g for r in utils.regions for g in self.glfo['seqs'][r]])
        if self.args.only_genes is None and len(genes_with_hmm_files - glfo_genes) > 0:
            print('  %s hmm files for genes that aren\'t in glfo: %s' % (utils.color('red', 'warning'), utils.color_genes(genes_with_hmm_files - glfo_genes)))
        if len(glfo_genes - genes_with_hmm_files) > 0:
            print('    skipping matches from %d genes that don\'t have hmm files: %s' % (len(glfo_genes - genes_with_hmm_files), utils.color_genes(glfo_genes - genes_with_hmm_files)))
        if self.current_action == 'cache-parameters' and len(glfo_genes - genes_with_enough_counts) > 0:
            print('    skipping matches from %d genes without enough counts: %s' % (len(glfo_genes - genes_with_enough_counts), utils.color_genes(glfo_genes - genes_with_enough_counts)))
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
    # @utils.timeprinter
    def get_nsets(self, algorithm, partition):

        if partition is not None:
            if len(partition) < 10000 and any(partition.count(c) > 1 for c in partition):  # there's nothing really *wrong* with having duplicates, but a) it's wasteful and b) it typically means something is wrong/nonsensical in the code that decided to send you the same task twice (first clause is just cause this is slow on super large samples, and this whole check is only really likely to trigger if we're debugging new code)
                for tcount, tclust in set([(partition.count(c), ':'.join(c)) for c in partition if partition.count(c) > 1]):
                    print('  %s cluster occurs %d times in the <nsets> we\'re sending to bcrham: %s' % (utils.color('yellow', 'warning'), tcount, tclust))
            nsets = copy.deepcopy(partition)  # needs to be a deep copy so we can shuffle the order
            if self.input_partition is not None or self.args.simultaneous_true_clonal_seqs or self.args.annotation_clustering:  # this is hackey, but we absolutely cannot have different cdr3 lengths in the same cluster, and these are two cases where it can happen (in very rare cases, usually on really crappy sequences)
                nsets = utils.split_clusters_by_cdr3(nsets, self.sw_info, warn=True)
        else:
            qlist = self.sw_info['queries']  # shorthand

            if self.args.simultaneous_true_clonal_seqs:  # NOTE this arg can now also be set when partitioning, but it's dealt with elsewhere
                print('  --simultaneous-true-clonal-seqs: grouping seqs according to true partition')
                nsets = utils.get_partition_from_reco_info(self.reco_info, ids=qlist)
                nsets = utils.split_clusters_by_cdr3(nsets, self.sw_info, warn=True)  # arg, have to split some clusters apart by cdr3, for rare cases where we call an shm indel in j within the cdr3
            elif self.args.all_seqs_simultaneous:  # everybody together
                nsets = [qlist]
                nsets = utils.split_clusters_by_cdr3(nsets, self.sw_info, warn=True)  # ok, this shouldn't happen any more (with msa_vs_info)
            elif self.args.n_simultaneous_seqs is not None:  # set number of simultaneous seqs
                nlen = self.args.n_simultaneous_seqs  # shorthand
                # nsets = [qlist[iq : min(iq + nlen, len(qlist))] for iq in range(0, len(qlist), nlen)]  # this way works fine, but it's hard to get right 'cause it's hard to understand
                nsets = [list(group) for _, group in itertools.groupby(qlist, key=lambda q: qlist.index(q) // nlen)]  # integer division
            else:  # plain ol' singletons
                nsets = [[q] for q in qlist]

        return nsets

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
    def process_dummy_d_hack(self, line, debug=False):
        """
        a.t.m. we force bcrham to give us D of length one for loci with no D genes.
        Here, we delete the dummy D base, and give it to either V, J, or the insertion.
        """
        tmpline = copy.deepcopy(line)
        utils.add_implicit_info(self.glfo, tmpline, reset_indel_genes=True)
        if debug:
            print('')
            print('  dummy d hack for %s' % ' '.join(line['unique_ids']))
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
            print('    gl_j_base', gl_j_base)
            print('    gl_v_base', gl_v_base)

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

        sorted_votes = sorted(sorted(votes.items()), key=operator.itemgetter(1), reverse=True)  # first do default sort (will alphabetize keys), then reverse sort by votes (first step ensures repeatable order of ties)
        winner = sorted_votes[0][0]
        sorted_qr_base_votes = sorted(sorted(qr_base_votes.items()), key=operator.itemgetter(1), reverse=True)  # same comment ^ as previous sort
        qr_base_winner = sorted_qr_base_votes[0][0]
        if debug:
            print('   ', sorted_votes)
            print('   ', sorted_qr_base_votes)
            print('    winner', winner, qr_base_winner)

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
    def check_for_unexpectedly_missing_keys(self, annotation_list, hmm_failure_ids):
        missing_input_keys = set(self.input_info)
        missing_input_keys -= set([uid for line in annotation_list for uid in line['unique_ids']])  # set(self.sw_info['queries'])  # all the queries for which we had decent sw annotations (sw failures are accounted for below)
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
            print('  %s couldn\'t account for %d missing input uid%s%s' % (utils.color('red', 'warning'), len(missing_input_keys), utils.plural(len(missing_input_keys)), ': %s' % ' '.join(missing_input_keys) if len(missing_input_keys) < 15 else ''))

    # ----------------------------------------------------------------------------------------
    def process_annotation_output(self, annotation_list, hmm_failures, count_parameters=False, parameter_out_dir=None, print_annotations=False):
        # ----------------------------------------------------------------------------------------
        def print_bad_annotations(mname, bad_annotations):  # useful for parsing resulting log file: grep -v annotations log | grep -A3 'true\|inferred\|worst'
            worst_values, worst_annotations = zip(*sorted(bad_annotations[mname], key=lambda x: abs(x[0]), reverse=True))
            print('\n\n%s %d by %s: %s%s%s %s' % (utils.color('red', 'printing worst'), self.args.print_n_worst_annotations, mname, utils.color('blue', '['), ' '.join(str(v) for v in worst_values[:self.args.print_n_worst_annotations]), utils.color('blue', ']'), ' '.join(str(v) for v in worst_values[self.args.print_n_worst_annotations:])))
            def pstr(val, line):
                vstr = ('%d'%val) if 'hamming' in mname else '%d' % (len(l[mname]) if 'insertion' in mname else l[mname])
                return '  %s %s%s' % (mname, vstr, '' if 'hamming' in mname else ' (off by %d)'%val)
            perf_strs = [pstr(v, l) for v, l in zip(worst_values, worst_annotations)]
            self.print_results(None, worst_annotations[:self.args.print_n_worst_annotations], dont_sort=True, label_list=perf_strs[:self.args.print_n_worst_annotations], extra_str='  ')

        # ----------------------------------------------------------------------------------------
        pcounter = ParameterCounter(self.glfo, self.args, count_correlations=self.args.count_correlations) if count_parameters else None
        true_pcounter = ParameterCounter(self.simglfo, self.args, count_correlations=self.args.count_correlations) if (count_parameters and not self.args.is_data) else None
        perfplotter = PerformancePlotter('hmm') if self.args.plot_annotation_performance else None

        if self.args.print_n_worst_annotations:
            bad_annotations = OrderedDict([
                # NOTE that in the loop below we only count iseq=0, so (unless you change that) these have to be per-family (i.e. *not* per-sequence) variables
                ('v_hamming_to_true_naive', []), ('cdr3_hamming_to_true_naive', []), ('j_hamming_to_true_naive', []),
                ('v_3p_del', []), ('d_5p_del', []), ('d_3p_del', []), ('j_5p_del', []),
                ('vd_insertion', []), ('dj_insertion', []),
            ])

        for line_to_use in annotation_list:
            if pcounter is not None:
                pcounter.increment(line_to_use)

            if perfplotter is not None:
                for iseq in range(len(line_to_use['unique_ids'])):  # NOTE this counts rearrangement-level parameters once for every mature sequence, which is inconsistent with the pcounters... but I think might make more sense here?
                    pvals = perfplotter.evaluate(self.reco_info[line_to_use['unique_ids'][iseq]], utils.synthesize_single_seq_line(line_to_use, iseq), simglfo=self.simglfo)
                    if self.args.print_n_worst_annotations and iseq == 0 and pvals is not None:
                        for mname in bad_annotations:
                            bad_annotations[mname].append((pvals[mname], line_to_use))

        if true_pcounter is not None:
            for uids in utils.get_partition_from_reco_info(self.reco_info, ids=self.sw_info['queries']):  # NOTE this'll include queries that passed sw but failed the hmm... there aren't usually really any of those
                true_pcounter.increment(utils.synthesize_multi_seq_line_from_reco_info(uids, self.reco_info))

        # parameter and performance writing/plotting
        if pcounter is not None:
            path_str = 'hmm' if parameter_out_dir is None else os.path.basename(parameter_out_dir)  # at the moment, this is just 'hmm' for regular parameter caching, and 'multi-hmm' for parameter caching after partitioning
            if self.args.plotdir is not None:
                pcounter.plot('%s/%s' % (self.args.plotdir, path_str), only_csv=self.args.only_csv_plots, only_overall=not self.args.make_per_gene_plots, make_per_base_plots=self.args.make_per_gene_per_base_plots)
                if true_pcounter is not None:
                    true_pcounter.plot('%s/%s' % (self.args.plotdir, path_str.replace('hmm', 'true')), only_csv=self.args.only_csv_plots, only_overall=not self.args.make_per_gene_plots, make_per_base_plots=self.args.make_per_gene_per_base_plots)
            if not self.args.dont_write_parameters:
                pcounter.write(parameter_out_dir)
                if true_pcounter is not None:
                    true_pcounter.write('%s/%s' % (os.path.dirname(parameter_out_dir), path_str.replace('hmm', 'true')))

        if perfplotter is not None and self.args.plotdir is not None:
            perfplotter.plot(self.args.plotdir + '/hmm', only_csv=self.args.only_csv_plots)

        if self.args.align_constant_regions:
            utils.parse_constant_regions(self.args.species, self.args.locus, annotation_list, self.args.workdir, debug=self.args.debug)

        if print_annotations or self.args.print_n_worst_annotations is not None:
            if self.args.print_n_worst_annotations is None:
                self.print_results(None, annotation_list)
            else:
                for mname in bad_annotations:
                    print_bad_annotations(mname, bad_annotations)
        self.check_for_unexpectedly_missing_keys(annotation_list, hmm_failures)  # NOTE not sure if it's really correct to use <annotation_list>, [maybe since <hmm_failures> has ones that failed the conversion to eroded line (and maybe other reasons)]

    # ----------------------------------------------------------------------------------------
    def add_per_seq_sw_info(self, line):
        uids = line['unique_ids']
        line['indelfos'] = [self.sw_info['indels'].get(uid, indelutils.get_empty_indel()) for uid in uids]  # reminder: hmm was given a sequence that had any indels reversed (i.e. <self.sw_info['indels'][uid]['reverersed_seq']>)
        line['input_seqs'] = [self.sw_info[uid]['input_seqs'][0] for uid in uids]  # not in <line>, since the hmm doesn't know anything about the input (i.e. non-indel-reversed) sequences
        line['duplicates'] = [self.duplicates.get(uid, []) for uid in uids]
        def gv(lkey, uid): return self.sw_info[uid][lkey][0] if lkey in self.sw_info[uid] else utils.input_metafile_defaults(lkey)
        for lkey in list(utils.input_metafile_keys.values()) + ['leader_seqs', 'c_gene_seqs']:
            if any(lkey in self.sw_info[u] for u in uids):  # it used to be that if it was in one it had to be in all, but now no longer (see comments in input meta info reading in seqfileopener)
                line[lkey] = [gv(lkey, u) for u in uids]

    # ----------------------------------------------------------------------------------------
    def read_annotation_output(self, annotation_fname, count_parameters=False, parameter_out_dir=None, print_annotations=False, is_subcluster_recursed=False):  # for search: read hmm read_hmm
        """ Read bcrham annotation output """
        # ----------------------------------------------------------------------------------------
        def check_invalid(line, hmm_failures):
            if line['invalid']:
                counts['n_invalid_events'] += 1
                if self.args.debug:
                    print('      %s line invalid' % uidstr)
                    # utils.print_reco_event(line, extra_str='    ', label='invalid:')  # can't normally print invalid lines, so commenting this (although sometimes maybe you can?)
                hmm_failures |= set(line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                return True
            return False
        # ----------------------------------------------------------------------------------------
        def remove_insertions_and_deletions(line, hmm_failures):
            ends = ['v_3p', 'j_5p'] if not utils.has_d_gene(self.args.locus) else utils.real_erosions  # need 1-base d erosion for light chain
            del_lens = [line[e+'_del'] for e in ends]
            ins_lens = [len(line[b+'_insertion']) for b in utils.boundaries]
            if sum(ins_lens) != sum(del_lens):  # if the total insertion and deletion lengths aren't equal, there's no possible annotation without insertions or deletions
                print('  --no-insertions-or-deletions: total insertion len %d (%s) not equal to total del len %d (%s) so set annotation to failed for %s' % (sum(ins_lens), ins_lens, sum(del_lens), del_lens, line['unique_ids']))
                hmm_failures |= set(line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                return True
            for end in ends:
                line[end+'_del'] = 0
            for bound in utils.boundaries:
                line[bound+'_insertion'] = ''
            return False
        # ----------------------------------------------------------------------------------------
        if self.print_status or not is_subcluster_recursed:
            print('    reading output')
        sys.stdout.flush()

        counts = {n : 0 for n in ['n_lines_read', 'n_seqs_processed', 'n_events_processed', 'n_invalid_events']}
        eroded_annotations, padded_annotations = OrderedDict(), OrderedDict()
        hmm_failures = set()  # hm, does this duplicate info I'm already keeping track of in one of these other variables?
        errorfo = {}
        with open(annotation_fname, 'r') as hmm_csv_outfile:
            reader = csv.DictReader(hmm_csv_outfile)
            for padded_line in reader:  # line coming from hmm output is N-padded such that all the seqs are the same length
                utils.process_input_line(padded_line)
                counts['n_lines_read'] += 1

                failed = self.check_did_bcrham_fail(padded_line, errorfo)
                if failed:
                    hmm_failures |= set(padded_line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                    continue

                uids = padded_line['unique_ids']
                uidstr = ':'.join(uids)

                self.add_per_seq_sw_info(padded_line)
                if not utils.has_d_gene(self.args.locus):
                    self.process_dummy_d_hack(padded_line)
                if self.args.no_insertions_or_deletions:
                    failed = remove_insertions_and_deletions(padded_line, hmm_failures)
                    if failed:
                        continue

                try:
                    utils.add_implicit_info(self.glfo, padded_line, aligned_gl_seqs=self.aligned_gl_seqs, reset_indel_genes=True)
                except:  # I really don't like just swallowing it, but it's crashing deep in the new[ish] indel code on an extraordinarily rare and I think super screwed up sequence, and I can't replicate it without running on the entire stupid huge sample
                    lines = traceback.format_exception(*sys.exc_info())
                    print('      %s implicit info adding failed for %s when reading hmm output (so adding to failed queries):' % (utils.color('red', 'warning'), uidstr))
                    print(utils.pad_lines(''.join(lines)))
                    hmm_failures |= set(padded_line['unique_ids'])  # NOTE adds the ids individually (will have to be updated if we start accepting multi-seq input file)
                    continue

                utils.process_per_gene_support(padded_line)  # switch per-gene support from log space to normalized probabilities

                if check_invalid(padded_line, hmm_failures):
                    continue

                if uidstr in padded_annotations:  # this shouldn't happen, but it's more an indicator that something else has gone wrong than that in and of itself it's catastrophic
                    print('  %s uidstr %s already read from file %s' % (utils.color('yellow', 'warning'), uidstr, annotation_fname))
                padded_annotations[uidstr] = padded_line

                line_to_use = padded_line
                if self.args.mimic_data_read_length:  # this code is old and hasn't been tested recently so may be a bit dodgy
                    if len(uids) > 1:
                        print('  %s can\'t mimic data read length on multi-hmm annotations, since we need the padding to make lengths compatible (at least, I think it will crash just below here if you try)' % utils.color('red', 'error'))
                    # get a new dict in which we have edited the sequences to swap Ns on either end (after removing fv and jf insertions) for v_5p and j_3p deletions
                    eroded_line = utils.reset_effective_erosions_and_effective_insertions(self.glfo, padded_line, aligned_gl_seqs=self.aligned_gl_seqs)  #, padfo=self.sw_info)
                    if check_invalid(eroded_line, hmm_failures):
                        continue
                    line_to_use = eroded_line
                    eroded_annotations[uidstr] = eroded_line  # these only get used if there aren't any multi-seq lines, so it's ok that they don't all get added if there is a multi seq line

                counts['n_events_processed'] += 1
                counts['n_seqs_processed'] += len(uids)

        os.remove(annotation_fname)

        if self.print_status or not is_subcluster_recursed:
            print('        read %d hmm output lines with %d sequences in %d events  (%d failures)' % (counts['n_lines_read'], counts['n_seqs_processed'], counts['n_events_processed'], len(hmm_failures)))
        if counts['n_invalid_events'] > 0:
            print('            %s skipped %d invalid events' % (utils.color('red', 'warning'), counts['n_invalid_events']))
        for ecode in errorfo:
            if ecode == 'no_path':
                print('          %s no valid paths: %s' % (utils.color('red', 'warning'), ' '.join(sorted(errorfo[ecode]))))
            elif ecode == 'boundary':
                print('          %d boundary warnings' % len(errorfo[ecode]))
                if self.args.debug:
                    print('                %s' % ' '.join(errorfo[ecode]))
            else:
                print('          %s unknown ecode \'%s\': %s' % (utils.color('red', 'warning'), ecode, ' '.join(errorfo[ecode])))

        annotation_list = list(eroded_annotations.values()) if self.args.mimic_data_read_length else list(padded_annotations.values())
        seqfileopener.add_input_metafo(self.input_info, annotation_list, keys_not_to_overwrite=['multiplicities'])  # don't overwrite any info that's already in there (presumably multiplicities) since it will have been updated in waterer after collapsing duplicates NOTE/UPDATE if you screw something up though, this may end up not overwriting 'paired-uids' that you *do* want it to overwrite
        if not is_subcluster_recursed:
            self.process_annotation_output(annotation_list, hmm_failures, count_parameters=count_parameters, parameter_out_dir=parameter_out_dir, print_annotations=print_annotations)

        return utils.get_annotation_dict(annotation_list), hmm_failures

    # ----------------------------------------------------------------------------------------
    def write_output(self, annotation_list, hmm_failures, cpath=None, dont_write_failed_queries=False, write_sw=False, outfname=None, extra_headers=None):
        if outfname is None:
            outfname = self.args.outfname

        if write_sw:
            assert annotation_list is None
            annotation_list = [self.sw_info[q] for q in self.input_info if q in self.sw_info['queries']]

        failed_queries = None
        if not dont_write_failed_queries:  # write empty lines for seqs that failed either in sw or the hmm
            failed_queries = [{'unique_ids' : [uid], 'invalid' : True, 'input_seqs' : self.input_info[uid]['seqs']} for uid in self.sw_info['failed-queries'] | hmm_failures]  # <uid> *needs* to be single-sequence (but there shouldn't really be any way for it to not be)

        if self.args.presto_output:
            presto_annotation_fname = outfname
            if cpath is not None:  # note that if we're partitioning, we get here with <self.current_action> set to 'annotate'
                presto_annotation_fname = utils.getprefix(outfname) + '.tsv'
                cpath.write_presto_partitions(outfname, self.input_info)
            utils.write_presto_annotations(presto_annotation_fname, annotation_list, failed_queries=failed_queries)
            return
        elif self.args.airr_output:
            utils.write_airr_output(utils.replace_suffix(outfname, '.tsv'), annotation_list, cpath=cpath, failed_queries=failed_queries, extra_columns=self.args.extra_annotation_columns, args=self.args)  # suffix may already be .tsv, but that's fine
            if utils.getsuffix(outfname) == '.tsv':  # if it isn't .tsv, we also write the regular partis file
                return

        partition_lines = None
        if cpath is not None:
            true_partition = utils.get_partition_from_reco_info(self.reco_info) if not self.args.is_data else None
            partition_lines = cpath.get_partition_lines(reco_info=self.reco_info, true_partition=true_partition,
                                                        n_to_write=self.args.n_partitions_to_write, calc_missing_values=('all' if (len(annotation_list) < 500) else 'best'), fail_frac=self.args.max_ccf_fail_frac, add_pairwise_metrics=self.args.add_pairwise_clustering_metrics)

        if self.args.extra_annotation_columns is not None and 'linearham-info' in self.args.extra_annotation_columns:  # it would be nice to do this in utils.add_extra_column(), but it requires sw info, which would then have to be passed through all the output infrastructure
            utils.add_linearham_info(self.sw_info, annotation_list, self.glfo, min_cluster_size=5)  # NOTE not really worth trying to propagate through --cluster-indices or --seed/lineage-unique-ids here (i.e. propagate from linearham scons file) (setting hard coded 5 will fuck you if you want linearham on a tree with 4 seqs, but you probably don't really, and both not having a threshold here, and using --min-selection-metric-cluster-size suck [yes i tried both])

        headers = utils.sw_cache_headers if write_sw else utils.annotation_headers
        if extra_headers is not None:
            headers += extra_headers
        headers = utils.add_lists(headers, self.args.extra_annotation_columns)
        if utils.getsuffix(outfname) == '.csv':
            if cpath is not None:
                cpath.write(outfname, self.args.is_data, partition_lines=partition_lines)  # don't need to pass in reco_info/true_partition since we passed them when we got the partition lines
            annotation_fname = outfname if cpath is None else self.args.cluster_annotation_fname
            utils.write_annotations(annotation_fname, self.glfo, annotation_list, headers, failed_queries=failed_queries)
        elif utils.getsuffix(outfname) == '.yaml':
            utils.write_annotations(outfname, self.glfo, annotation_list, headers, failed_queries=failed_queries, partition_lines=partition_lines, use_pyyaml=self.args.write_full_yaml_output, dont_write_git_info=self.args.dont_write_git_info)
        else:
            raise Exception('unhandled annotation file suffix %s' % outfname)
