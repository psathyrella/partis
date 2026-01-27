#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import random
import numpy
import os
import csv
import glob
import math
import shutil
import time
from collections import OrderedDict
from subprocess import Popen, PIPE, check_call, check_output, CalledProcessError
import copy
import colored_traceback.always
import sys
from io import open
from pathlib import Path
import yaml

from partis.baseutils import get_extra_str
import partis.utils as utils
import partis.glutils as glutils
from partis.hist import Hist
from partis.clusterpath import ClusterPath
import partis.paircluster as paircluster

# ----------------------------------------------------------------------------------------
class Tester(object):
    # ----------------------------------------------------------------------------------------
    def dirs(self, tstr, force_paired=False):
        assert tstr in ['ref', 'new']
        bdir = 'test%s' % ('/paired' if (args.paired or force_paired) else '')
        if not os.path.exists(bdir):
            bdir = '%s/%s' % (utils.get_partis_dir(), bdir)
        return '%s/%s-results%s' % (bdir, tstr, '-slow' if args.slow else '')
    # ----------------------------------------------------------------------------------------
    def nqr(self, act):
        if act == 'quick':  # bears no relation to the others, so makes sense to handle it differently
            return 10
        nqdict = {'normal' : {'simu' :   50, 'data' :   100 if args.paired else 50},
                    'slow' : {'simu' : 1000, 'data' :   -1}}
        return nqdict['slow' if args.slow else 'normal'][act]
    # ----------------------------------------------------------------------------------------
    def get_stypes(self, ptest):
        namelist = ptest.split('-')
        if ptest == 'simulate':
            input_stype = 'new'
            input_dtype = None
        elif 'cache-parameters' in ptest:
            input_stype = 'new'
            input_dtype = namelist[-1]
        else:
            input_stype, input_dtype = namelist[-2:]

        assert input_stype in self.stypes + [None]
        assert input_dtype in self.dtypes + [None]
        return input_stype, input_dtype
    # ----------------------------------------------------------------------------------------
    def inpath(self, st, dt):
        if dt == 'data':
            return self.paired_datafname if args.paired else self.datafname
        else:
            spath = self.label + '/simu' + ('' if args.paired else '.yaml')
            if st is None:
                return spath
            return self.dirs(st) + '/' + spath
    # ----------------------------------------------------------------------------------------
    def paramdir(self, st, dt):
        pd = self.label + '/parameters/' + dt
        if st is None:  # if <st> isn't set, we want the subpath (without parent dir), e.g. when --dont-run/evaluating results
            return pd
        return self.dirs(st) + '/' + pd
    # ----------------------------------------------------------------------------------------
    def astr(self, inout, dt=None):
        if not args.paired or (args.paired and inout == 'in' and dt == 'data'):
            return '%sfname' % inout
        return 'paired-%sdir' % inout
    # ----------------------------------------------------------------------------------------
    # i'm adding this late (especially for the non-production tests) so there's probably some more places it could be used
    def opath(self, ptest, st=None, force_paired=False):  # don't set <st> if you just want the basename-type stuff (no base/parent dirs)
        if 'cache-parameters' in ptest:
            return self.paramdir(None, ptest.split('-')[2])
        elif ptest == 'simulate':
            return self.inpath(None, 'simu')
        else:
            op = '%s%s' % (ptest, '' if (args.paired or force_paired) else '.yaml')
            if (args.paired or force_paired) and 'get-selection-metrics' in ptest:
                op += '-chosen-abs.csv'
            if st is None:
                return op
            return '%s/%s' % (self.dirs(st, force_paired=force_paired), op)
    # ----------------------------------------------------------------------------------------
    def ptn_cachefn(self, st, for_cmd=False, lpair=None, locus=None):  # see note above for opath()
        assert st == 'new'  # i think?
        cfn = ''
        if for_cmd:
            if args.paired:
                return 'paired-outdir'
            else:
                cfn += self.dirs('new')
        if args.paired:
            assert lpair is not None and locus is not None
            cfn += '%s/persistent-cache-%s.csv' % ('+'.join(lpair), locus)  # duplicates code in bin/partis getofn()
        else:
            cfn = '%s%scache-%s-partition.csv' % (cfn, '' if cfn=='' else '/', st)
        return cfn
    # ----------------------------------------------------------------------------------------
    def all_ptn_cachefns(self):  # return all of them (ok atm it's juse the one, but we used to also have the 'ref' one, and maybe will want it in the future?)
        if args.paired:
            return [self.ptn_cachefn(s, lpair=lp, locus=l) for s in ['new'] for lp in utils.locus_pairs[args.ig_or_tr] for l in lp]
        else:
            return [self.ptn_cachefn(s) for s in ['new']]
    # ----------------------------------------------------------------------------------------
    def is_prod_test(self, ptest):
        return 'cache-parameters' in ptest or ptest == 'simulate'
    # ----------------------------------------------------------------------------------------
    def sclust_sizes(self):  # NOTE depends on self.n_simu_leaves
        return (15, 20) if args.slow else (5, 10)
    # ----------------------------------------------------------------------------------------
    def min_smetric_cluster_size(self):
        return 10 if args.slow else 3
    # ----------------------------------------------------------------------------------------
    def cluster_size_args(self):
        return ['--min-selection-metric-cluster-size', str(self.min_smetric_cluster_size()),
                '--min-paired-cluster-size-to-read', str(self.min_smetric_cluster_size())]
    # ----------------------------------------------------------------------------------------
    def __init__(self):
        self.partis_path = 'partis' if shutil.which('partis') else '%s/bin/partis' % utils.get_partis_dir()  # use version in PATH if it's there (pipx seems to leave two incompatible versions lying around)
        if args.prepend_coverage:
            self.partis_path = 'coverage3 run --append %s' % self.partis_path
        self.datafname = 'test/mishmash.fa'  # some data from adaptive, chaim, and vollmers
        # generate paired data dir with: UPDATE i cat'd the ig?.fa files into all-seqs.fa (in the same dir) so extract-pair-info and split-loci get run
        #  - ./bin/split-loci.py /fh/fast/matsen_e/data/10x-examples/data/sc5p_v2_hs_B_prevax_10k_5gex_B_vdj_b_filtered_contig.fasta --outdir test/paired-data --input-metafname /fh/fast/matsen_e/data/10x-examples/processed-data/v0/sc5p_v2_hs_B_prevax_10k_5gex_B_vdj_b_filtered_contig/meta.yaml --n-max-queries 100 >test/paired-data/split-loci.log
        #  - ./bin/extract-pairing-info.py /fh/fast/matsen_e/data/10x-examples/data/sc5p_v2_hs_B_prevax_10k_5gex_B_vdj_b_filtered_contig.fasta test/paired-data/meta.yaml --n-max-queries 100 >test/paired-data/extract.log
        self.paired_datafname = 'test/paired-data/all-seqs.fa'
        self.input_metafname = 'test/input-meta.yaml'

        self.stypes = ['ref', 'new']  # I don't know what the 's' stands for
        self.dtypes = ['data', 'simu']
        if not os.path.exists(self.dirs('new')):
            os.makedirs(self.dirs('new'))
        self.common_extras = ['--random-seed', '1', '--n-procs', '10']  # would be nice to set --n-procs based on the machine, but for some reason the order of things in the parameter files gets shuffled a bit if you change the number of procs
        self.label = 'test' #  i really don't think there's any reason to have this, but i don't feel like removing it atm since it's not really causing much complication

        self.tiny_eps = 1e-4
        self.run_times = {}
        self.eps_vals = {}  # fractional difference which we allow for each test type (these were generated with the code in get_typical_variances() above)
        self.eps_vals['v_call'] = 0.02  # hm, actually, I think I just made the annotation ones up
        self.eps_vals['d_call'] = 0.02
        self.eps_vals['j_call'] = 0.02
        self.eps_vals['mean_hamming']   = 0.1
        self.eps_vals['v_hamming'] = 0.1
        self.eps_vals['d_hamming'] = 0.1
        self.eps_vals['j_hamming'] = 0.1
        self.eps_vals['cdr3_hamming'] = 0.1
        self.eps_vals['purity']         = 0.08
        self.eps_vals['completeness']   = 0.08

        self.n_simu_leaves = 5

        self.selection_metrics = ['lbi', 'lbr', 'cons-dist-aa', 'aa-lbi', 'aa-lbr']  # NOTE kind of duplicates treeutils.selection_metrics, but I want to be able to change the latter
        self.pair_clean_metrics = ['correct', 'unpaired', 'mispaired'] if args.paired else []
        self.expected_trees = ['tree', 'aa-tree']

        self.logfname = self.dirs('new') + '/test.log'

        self.tests = OrderedDict()

        # ----------------------------------------------------------------------------------------
        def add_inference_tests(input_stype):  # if input_stype is 'ref', infer on old simulation and parameters, if it's 'new' use the new ones
            if not args.paired:
                self.tests['annotate-%s-simu'%input_stype] = {'extras' : ['--plot-annotation-performance', ]}
                self.tests['multi-annotate-%s-simu'%input_stype] = {'extras' : ['--plot-annotation-performance', '--simultaneous-true-clonal-seqs']}  # NOTE this is mostly different to the multi-seq annotations from the partition step because it uses the whole sample
            self.tests['partition-%s-simu'%input_stype] = {'extras' : ['--plot-annotation-performance', '--max-ccf-fail-frac', '0.10']} # '--biggest-logprob-cluster-to-calculate', '5', '--biggest-naive-seq-cluster-to-calculate', '5',
            if args.paired:
                self.tests['subset-partition-%s-simu'%input_stype] = {'extras' : ['--max-ccf-fail-frac', '0.15']} # '--biggest-logprob-cluster-to-calculate', '5', '--biggest-naive-seq-cluster-to-calculate', '5',
            # this runs ok, but i's need to modify some things so its output is actually checked self.tests['subset-annotate-%s-simu'%input_stype] = {'extras' : ['--max-ccf-fail-frac', '0.15']} # '--biggest-logprob-cluster-to-calculate', '5', '--biggest-naive-seq-cluster-to-calculate', '5',
            self.tests['seed-partition-%s-simu'%input_stype] = {'extras' : ['--max-ccf-fail-frac', '0.10']}
            if not args.paired:
                self.tests['vsearch-partition-%s-simu'%input_stype] = {'extras' : ['--naive-vsearch']}
            self.tests['get-selection-metrics-%s-simu'%input_stype] = {'extras' : ['--existing-output-run-cfg', 'paired'] + self.cluster_size_args()}  # NOTE this runs on simulation, but it's checking the inferred selection metrics

        # ----------------------------------------------------------------------------------------
        if args.quick:
            self.tests['cache-parameters-quick-new-simu'] =  {'extras' : ['--n-max-queries', str(self.nqr('quick'))]}
        else:
            pcache_data_args = {'extras' : ['--n-max-queries', str(self.nqr('data'))]}
            pcache_simu_args = {'extras' : []}
            n_events = int(self.nqr('simu') / float(self.n_simu_leaves))
            simulate_args = {'extras' : ['--n-sim-events', str(n_events), '--n-trees', str(n_events), '--n-leaf-distribution', 'geometric', '--n-leaves', str(self.n_simu_leaves)]}
            if args.paired:
                simulate_args['extras'] += ['--min-observations-per-gene', '5', '--mean-cells-per-droplet', '1.25', '--constant-cells-per-droplet', '--fraction-of-reads-to-remove', '0.15']  # it was crashing and this fixes it, i dunno if we should turn it on also for non-paired but whatever
            if args.bust_cache:  # if we're cache busting, we need to run these *first*, so that the inference tests run on a simulation file in the new dir that was just made (i.e. *not* whatever simulation file in the new dir happens to be there)
                self.tests['cache-parameters-data'] = pcache_data_args
                self.tests['simulate'] = simulate_args
            self.tests['cache-parameters-simu'] = pcache_simu_args
            add_inference_tests('new')
            if not args.bust_cache:  # normally (if we're not cache busting) we want to run these last, to make it super clear that the inference tests are running on the *reference* simulation file
                self.tests['cache-parameters-data'] = pcache_data_args
                if not args.no_simu:
                    self.tests['simulate'] = simulate_args

        self.quick_tests = ['cache-parameters-quick-new-simu']  # this is kind of dumb to keep track of what the quick tests are in two different ways, but I only just started not adding the non-quick tests if --quick is set, and I don't want to mess with all the other places that use <self.quick_tests>

        self.perfdirs = {}  # set in fiddle_with_arguments() NOTE these correspond only to annotation performance, whereas <self.perf_info> has also partition performance
        for ptest, argfo in self.tests.items():
            self.fiddle_with_arguments(ptest, argfo)

        self.perf_info = {version_stype : {} for version_stype in self.stypes}

    # ----------------------------------------------------------------------------------------
    def test(self, args):
        if not args.dont_run:
            self.run(args)
        if args.dry_run or args.bust_cache or args.quick:
            return
        self.compare_production_results(['cache-parameters-simu'])
        self.compare_stuff(input_stype='new')
        self.compare_production_results(['cache-parameters-data'])
        if not args.no_simu:
            self.compare_production_results(['simulate'])
        self.compare_run_times()

    # ----------------------------------------------------------------------------------------
    def run_coverage(self, args):
        # ----------------------------------------------------------------------------------------
        def run_cmd(cmd, shell=False, logfname=None, dont_prepend=False):
            if logfname is not None:
                print('        writing log to %s' % logfname)
            cov_str = 'coverage3 run --append'  # --data-file=%s/coverage/%d.cov (this doesn't seem to be supported in my version
            utils.simplerun('%s%s' % ('' if dont_prepend else cov_str+' ', cmd), shell=shell, logfname=logfname)
        # ----------------------------------------------------------------------------------------
        ivsn = 0
        while True:
            odir = '%s/vsn-%d' % (args.coverage_outdir, ivsn)
            if os.path.exists(odir):
                print('  coverage outdir %s exists, you may want to rm -r it by hand' % odir)
            else:
                break
            ivsn += 1
        cfn = '%s/.coverage' % os.getcwd()
        if os.path.exists(cfn):
            print('  removing existing coverage file %s' % cfn)
            os.remove(cfn)

        run_cmd('./test/test.py --prepend-coverage', dont_prepend=True)  # NOTE tests may fail because of the coverage stuff, which is fine (at the least they'll be way too slow)
        run_cmd('./test/test.py --prepend-coverage --paired', dont_prepend=True)  # also note that we have to put dont_prepend since recursive subprocs having coverage commands breaks things (at least before coverage 6.3)

        # cp output files so that working files (e.g. tree inference output files) don't get scattered around the normal test output dir
        if not os.path.exists(odir):
            os.makedirs(odir)
        ptnfn = '%s/%s' % (odir, utils.insert_before_suffix('-single', os.path.basename(self.opath('partition-new-simu', st='ref'))))
        utils.simplerun('cp %s %s' % (self.opath('partition-new-simu', st='ref'), ptnfn))
        pair_ptndir = '%s/%s-paired' % (odir, os.path.basename(self.opath('partition-new-simu', st='ref', force_paired=True)))
        utils.simplerun('cp -r %s %s' % (self.opath('partition-new-simu', st='ref', force_paired=True), pair_ptndir))

        for ft in ['csv', 'fa', 'yaml']:
            run_cmd('./bin/parse-output.py %s %s/parse-output.%s' % (ptnfn, odir, ft))
        run_cmd('./bin/parse-output.py %s %s/parse-output-paired --paired' % (pair_ptndir, odir))
        run_cmd('./bin/cf-alleles.py --bases all', logfname='%s/cf-alleles.log'%odir)
        run_cmd('./bin/cf-alleles.py --bases 8-51-1', logfname='%s/cf-alleles-8-51.log'%odir)

        run_cmd('./bin/partis view-output --outfname %s' % ptnfn, logfname='%s/view-output.log'%odir)
        run_cmd('./bin/partis view-output --paired-loci --paired-outdir %s' % pair_ptndir, logfname='%s/view-output-paired.log'%odir)

        run_cmd('./bin/cf-germlines.py %s/hmm/germline-sets %s/hmm/germline-sets' % (self.paramdir('ref', 'simu'), self.paramdir('ref', 'data')), logfname='%s/cf-germlines.log'%odir)
        run_cmd('./bin/compare-plotdirs.py --outdir %s/compare-plotdirs --plotdirs %s/hmm/mutation:%s/sw/mutation --names hmm:sw' % (odir, self.opath('annotate-new-simu-annotation-performance', st='ref').replace('.yaml', ''), self.opath('annotate-new-simu-annotation-performance', st='ref').replace('.yaml', '')))

        ptn_plot_cmd = './bin/partis plot-partitions --partition-plot-cfg mds:trees --tree-inference-method iqtree --cluster-indices 0:2 %s' % ' '.join(self.cluster_size_args())
        run_cmd('%s --outfname %s --plotdir %s/plot-partitions' % (ptn_plot_cmd, ptnfn, odir))
        run_cmd('%s --paired-loci --paired-outdir %s --plotdir %s/plot-partitions-paired' % (ptn_plot_cmd, pair_ptndir, odir))

        run_cmd('./bin/plot-hmms.py --outdir %s/plot-hmms --infiles %s' % (odir, ':'.join(glob.glob('%s/hmm/hmms/IGHD1*.yaml'%self.paramdir('ref', 'data')))))

        gct_sm_cmd = './bin/partis get-selection-metrics --tree-inference-method gctree %s --cluster-indices 0:2' % ' '.join(self.cluster_size_args())
        run_cmd('%s --outfname %s --plotdir %s/gctree-smetric-plots' % (gct_sm_cmd, ptnfn, odir))
        run_cmd('./bin/read-gctree-output.py --locus igh --species human --gctreedir %s/gctree/iclust-0 --outdir %s/read-gctree-output' % (utils.getprefix(ptnfn), odir))  # NOTE this uses output from the previous line
        run_cmd('%s --paired-loci --paired-outdir %s --plotdir %s/paired-gctree-smetric-plots' % (gct_sm_cmd, pair_ptndir, odir))

        run_cmd('./bin/bcr-phylo-run.py --base-outdir %s/bcr-phylo-run' % odir)  # don't use multiple gc rounds here, since we need the tree in the next line (and the tree isn't written for multiple gc rounds)
        run_cmd('./bin/smetric-run.py --infname %s/bcr-phylo-run/selection/simu/mutated-simu.yaml --base-plotdir %s/smetric-run --metric-method lbi' % (odir, odir))  # NOTE uses results from previous line
        run_cmd('./bin/bcr-phylo-run.py --base-outdir %s/bcr-phylo-run-paired --paired --n-gc-rounds 3 --obs-times 30:5,10:5' % odir)
        run_cmd('./test/cf-paired-loci.py --label coverage --version v0 --n-replicates 2 --n-sub-procs 10 --scratch-mute-freq-list 0.01:0.1 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs --mutate-stop-codons\" --final-plot-xvar scratch-mute-freq --n-leaves-list 3 --n-sim-events-list 100 --single-light-locus igk --base-outdir %s/cf-paired-loci --perf-metrics precision:sensitivity --actions simu:cache-parameters:partition:plot:combine-plots' % odir, shell=True)

        # it'd be nice to add germline inference to the normal (non-coverage) tests, but doing more than a trivial test like this requires lots of sequences, which is slower than I really want to add to testing
        run_cmd('./bin/test-germline-inference.py --prepend-coverage-command --n-sim-events 1000 --outdir %s/test-germline-inference --sim-v-genes=IGHV1-18*01 --inf-v-genes=IGHV1-18*01 --snp-positions 27,55,88 --mutation-multiplier 0.00001 --seed 1' % odir, dont_prepend=True)

        print('now run: coverage3 report --omit=python/__init__.py')  # could automate this, but whatever

    # ----------------------------------------------------------------------------------------
    def fiddle_with_arguments(self, ptest, argfo):
        input_stype, input_dtype = self.get_stypes(ptest)
        argfo['input_stype'] = input_stype
        argfo['bin'] = self.partis_path

        if ptest == 'simulate':
            argfo['parameter-dir'] = self.paramdir(input_stype, 'data')
            if args.no_tree_gen:
                argfo['extras'] += ['--input-simulation-treefname', '%s/test/trees.nwk'%utils.get_partis_dir()]
            if args.no_per_base_mutation:
                argfo['extras'] += ['--no-per-base-mutation']
        else:
            argfo['inpath'] = self.inpath('new' if args.bust_cache else 'ref', input_dtype)
            if ptest.find('subset-') != 0:
                argfo['parameter-dir'] = self.paramdir(input_stype, input_dtype)
            if not args.paired:
                argfo['sw-cachefname'] = self.paramdir(input_stype, input_dtype) + '/sw-cache.yaml'
        if ptest != 'simulate' and input_dtype == 'simu' and not args.quick:
            argfo['extras'] += ['--is-simu']

        if 'annotate' in ptest:
            argfo['action'] = 'annotate'
        elif 'partition' in ptest:
            argfo['action'] = 'partition'
            if not args.paired:  # eh, i don't think there's really a reason to do this for paired (although I partially implemented -- i got the files in single-chain/, but then getting the igh+igk/ etc. was going to be more work)
                argfo['extras'] += ['--persistent-cachefname', self.ptn_cachefn(input_stype, for_cmd=True)]
        elif 'get-selection-metrics' in ptest:
            argfo['action'] = 'get-selection-metrics'  # could really remove almost all of the arguments, mostly just need --outfname
        elif 'cache-parameters-' in ptest:
            argfo['action'] = 'cache-parameters'
        else:
            argfo['action'] = ptest
        if ptest.find('subset-') == 0:
            argfo['action'] = 'subset-%s' % argfo['action']
            argfo['extras'] += ['--n-subsets', '2']

        if '--plot-annotation-performance' in argfo['extras']:
            self.perfdirs[ptest] = ptest + '-annotation-performance'
            argfo['extras'] += ['--plotdir', self.dirs('new') + '/' + self.perfdirs[ptest]]

        if '--plotdir' in argfo['extras']:
            argfo['extras'] += ['--only-csv-plots']
            if 'partition' in ptest:
                argfo['extras'] += ['--no-partition-plots']

        if '-data' in ptest and not args.paired:  # would be cleaner to check that inpath is self.datafname
            argfo['extras'] += ['--input-metafnames', self.input_metafname]

    # ----------------------------------------------------------------------------------------
    def compare_stuff(self, input_stype):
        print('%s input' % input_stype)
        for version_stype in self.stypes:  # <version_stype> is the code version, i.e. 'ref' is the reference results, 'new' is the results we just made with the new code
            self.read_annotation_performance(version_stype, input_stype)
            self.read_partition_performance(version_stype, input_stype)  # NOTE also calls read_annotation_performance()
            self.read_selection_metric_performance(version_stype, input_stype)
        self.compare_performance(input_stype)
        if not args.paired:
            self.compare_partition_cachefiles(input_stype)

    # ----------------------------------------------------------------------------------------
    def prepare_to_run(self, args, name, info):
        """ Pre-run stuff that you don't want to do until *right* before you actually run. """
        # ----------------------------------------------------------------------------------------
        def clean_dir(sdir):  # rm whole seed dir to make sure the dir for the previous seed id doesn't hang around
            if args.dry_run:
                print('    would rm %s' % sdir)
            else:  # maybe i should just rm the output dir for every test before running? although it might get complicated since some of them i think share dirs
                print('    removing %s' % sdir)
                shutil.rmtree(sdir)
        # ----------------------------------------------------------------------------------------
        def rm_file(fn):  # for search: remove
            if args.dry_run:
                files_to_rm.append(fn)
            else:
                check_call(['rm', '-v', fn])
        # ----------------------------------------------------------------------------------------
        files_to_rm = []  # just for dbg
        # delete any old partition cache files
        if name == 'partition-' + info['input_stype'] + '-simu':
            cachefnames = ['%s/%s' % (self.dirs('new'), f) for f in self.all_ptn_cachefns()]
            for cfn in [f for f in cachefnames if os.path.exists(f)]:
                rm_file(cfn)
        # and any old tree inference files
        if name == 'get-selection-metrics-' + info['input_stype'] + '-simu':
            tfns = []
            for subd in ['', '/*+*/partition-*']:
                for tmeth in ['fasttree', 'iqtree']:
                    for ftp in ['fasttree.out', 'log*', 'input-seqs.fa']:
                        tfns += glob.glob('%s%s/%s/iclust-*/%s' % (self.opath('partition-new-simu', st='new'), subd, tmeth, ftp))
            for ffn in [f for f in tfns if os.path.exists(f)]:
                rm_file(ffn)
        if len(files_to_rm) > 0:
            print('    would rm %d tree inference working files' % len(files_to_rm))

        # choose a seed uid
        if name == 'seed-partition-' + info['input_stype'] + '-simu':
            ifn = info['inpath']
            seed_uid, _ = utils.choose_seed_unique_id(ifn, self.sclust_sizes()[0], self.sclust_sizes()[1], paired=args.paired)  # , n_max_queries=self.nqr('partition')
            if args.paired:
                seed_uid, seed_loci = seed_uid
                info['extras'] += ['--seed-unique-id', ':'.join(seed_uid), '--seed-loci', ':'.join(seed_loci)]
                sdir = '%s/seeds' % self.opath(name, st=info['input_stype'])
                if os.path.exists(sdir):
                    clean_dir(sdir)
            else:
                info['extras'] += ['--seed-unique-id', seed_uid]

        if name.find('subset-') == 0:
            if os.path.exists(self.opath(name, st='new')):
                clean_dir(self.opath(name, st='new'))

    # ----------------------------------------------------------------------------------------
    def run(self, args):
        if not args.dry_run:
            open(self.logfname, 'w').close()

        os.environ['PYTHONHASHSEED'] = '0'  # turn off hash seed randomization for repeatable results
        for name, info in self.tests.items():
            if args.quick and name not in self.quick_tests:
                continue

            self.prepare_to_run(args, name, info)

            ist, idt = self.get_stypes(name)
            action = info['action']
            cmd_str = info['bin'] + ' ' + action + ' --dont-write-git-info'
            if args.prepend_coverage:
                cmd_str += ' --prepend-coverage-command'
            if args.paired:
                cmd_str += ' --paired-loci'
            for tkey in ['inpath', 'parameter-dir', 'sw-cachefname']:
                if tkey in info:
                    cmd_str += ' --%s %s' % (tkey if tkey != 'inpath' else self.astr('in', dt=idt), info[tkey])
            cmd_str += ' %s' % ' '.join(info['extras'] + self.common_extras)

            if name == 'simulate':
                cmd_str += ' --%s %s' % (self.astr('out'), self.inpath('new', 'simu'))
                cmd_str += ' --indel-frequency 0.2'  # super high indel freq, partly because we want to make sure we have some even with the new smaller default number of seqs
            elif 'get-selection-metrics-' in name:
                cmd_str += ' --%s %s' % (self.astr('out'), self.opath(name.replace('get-selection-metrics-', 'partition-'), st='new'))
                cmd_str += ' --%s %s' % ('chosen-ab-fname' if args.paired else 'selection-metric-fname', self.opath(name, st='new'))
                if args.slow:
                    cmd_str += ' --ab-choice-cfg %s/test/ab-choice-slow.yaml' % utils.get_partis_dir()
                clist = cmd_str.split()
                utils.remove_from_arglist(clist, '--%s'%self.astr('in', dt=idt), has_arg=True)
                utils.remove_from_arglist(clist, '--parameter-dir', has_arg=True)
                utils.remove_from_arglist(clist, '--sw-cachefname', has_arg=True)
                utils.remove_from_arglist(clist, '--is-simu')
                cmd_str = ' '.join(clist)
            elif 'cache-parameters-' not in name:
                cmd_str += ' --%s %s' % (self.astr('out'), self.opath(name, st='new'))

            logstr = '%s   %s' % (utils.color('green', name, width=30, padside='right'), cmd_str)
            print(logstr if utils.len_excluding_colors(logstr) < args.print_width else logstr[:args.print_width] + '[...]')
            if args.dry_run:
                continue
            logfile = open(self.logfname, 'a')
            logfile.write(logstr + '\n')
            logfile.close()
            start = time.time()
            try:
                check_call(cmd_str + ' 1>>' + self.logfname + ' 2>>' + self.logfname, shell=True)
                if args.quick:
                    print('\n  %s' % 'ok')
            except CalledProcessError as err:
                # print err  # this just says it exited with code != 0
                print('  log tail: %s' % self.logfname)
                print(utils.pad_lines(check_output(['tail', self.logfname], universal_newlines=True)))
                sys.exit(1)  # raise Exception('exited with error')
            self.run_times[name] = time.time() - start  # seconds

        if not args.quick and not args.dry_run:
            self.write_run_times()

    # ----------------------------------------------------------------------------------------
    def remove_reference_results(self, expected_content):
        print('  removing ref files')
        dir_content = set([os.path.basename(f) for f in glob.glob(self.dirs('ref') + '/*')])
        if len(dir_content - expected_content) > 0 or len(expected_content - dir_content) > 0:
            if len(dir_content - expected_content) > 0:
                print('in ref dir but not expected\n    %s' % (utils.color('red', ' '.join(dir_content - expected_content))))
            if len(expected_content - dir_content) > 0:
                print('expected but not in ref dir\n    %s' % (utils.color('red', ' '.join(expected_content - dir_content))))
            raise Exception('unexpected or missing content in reference dir (see above)')
        for fname in [self.dirs('ref') + '/' + ec for ec in expected_content]:
            print('    %srm %s' % ('(would) ' if args.dry_run else '', fname))
            if args.dry_run:
                continue
            if os.path.isdir(fname):
                shutil.rmtree(fname)
            else:
                os.remove(fname)

    # ----------------------------------------------------------------------------------------
    def bust_cache(self):
        test_outputs = [self.opath(k) for k in self.tests if not self.is_prod_test(k)]
        expected_content = set(test_outputs + list(self.perfdirs.values()) + [os.path.basename(self.logfname), self.label])
        if not args.paired:
            expected_content |= set(self.all_ptn_cachefns())  # they're in the partition outdir if --paired is set, so don't need to be moved
        expected_content.add('run-times.csv')

        # remove (very, very gingerly) whole reference dir
        self.remove_reference_results(expected_content)

        # copy over parameters, simulation, and plots
        print('  mv new files to ref')
        for fname in expected_content:
            print('    mv %s   -->  %s/' % (fname, self.dirs('ref')))
            if args.dry_run:
                continue
            shutil.move(self.dirs('new') + '/' + fname, self.dirs('ref') + '/')

    # ----------------------------------------------------------------------------------------
    def read_annotation_performance(self, version_stype, input_stype, these_are_cluster_annotations=False):
        for sequence_multiplicity in ['single', 'multi']:
            self.read_each_annotation_performance(sequence_multiplicity, version_stype, input_stype, these_are_cluster_annotations=these_are_cluster_annotations)

    # ----------------------------------------------------------------------------------------
    def read_each_annotation_performance(self, sequence_multiplicity, version_stype, input_stype, these_are_cluster_annotations=False):  # <these_are_cluster_annotations> means this fcn is being called from within read_partition_performance()
        """ version_stype is the code version, while input_stype is the input data version, i.e. 'ref', 'new' is the reference code version (last commit) run on the then-new simulation and parameters"""
        if these_are_cluster_annotations:
            ptest = '-'.join(['partition', input_stype, 'simu'])
            methods = ['hmm']
        elif sequence_multiplicity == 'single':
            ptest = '-'.join(['annotate', input_stype, 'simu'])
            methods = ['sw', 'hmm']
        elif sequence_multiplicity == 'multi':
            ptest = '-'.join(['multi', 'annotate', input_stype, 'simu'])
            methods = ['hmm']
        else:
            assert False
        if (args.quick and ptest not in self.quick_tests) or ptest not in self.tests:  # ok i think now i'm adding the second clause i don't need the first, but not quite sure
            return
        if input_stype not in self.perf_info[version_stype]:
            self.perf_info[version_stype][input_stype] = OrderedDict()
        if ptest not in self.perf_info[version_stype][input_stype]:
            self.perf_info[version_stype][input_stype][ptest] = OrderedDict()
        perfdir = self.dirs(version_stype) + '/' + self.perfdirs[ptest]
        perffo = self.perf_info[version_stype][input_stype][ptest]
        for method in methods:
            perffo[method] = OrderedDict()  # arg, this is deeper than I'd like
            perffo[method]['mean_hamming'] = Hist(fname=perfdir + '/' + method + '/mutation/hamming_to_true_naive.csv').get_mean()
            for region in utils.regions + ['cdr3']:
                perffo[method][region + '_hamming'] = Hist(fname=perfdir + '/' + method + '/mutation/' + region + '_hamming_to_true_naive.csv').get_mean()
            for bound in utils.boundaries:
                perffo[method][bound + '_insertion'] = Hist(fname=perfdir + '/' + method + '/boundaries/' + bound + '_insertion.csv').get_mean(absval=True)
            for erosion in utils.real_erosions:
                perffo[method][erosion + '_del'] = Hist(fname=perfdir + '/' + method + '/boundaries/' + erosion + '_del.csv').get_mean(absval=True)

    # ----------------------------------------------------------------------------------------
    def do_this_test(self, tstr, input_stype, pt):
        if tstr not in pt:
            return False
        if input_stype not in pt:
            return False
        if args.quick and pt not in self.quick_tests:
            return False
        return True

    # ----------------------------------------------------------------------------------------
    def read_partition_performance(self, version_stype, input_stype, debug=False):
        """ Read new partitions from self.dirs('new'), and put the comparison numbers in self.perf_info (compare either to true, for simulation, or to the partition in reference dir, for data). """
        # ----------------------------------------------------------------------------------------
        def read_cpath(fname):
            _, _, cpath = utils.read_yaml_output(fname=fname, skip_annotations=True)
            ccfs = cpath.ccfs[cpath.i_best]
            if None in ccfs:
                raise Exception('none type ccf read from %s' % fname)
            if debug:
                print('    %5.2f          %5.2f      %-28s   to true partition' % (ccfs[0], ccfs[1], fname)) #os.path.basename(fname))
            return ccfs
        # ----------------------------------------------------------------------------------------
        ptest_list = [k for k in self.tests.keys() if self.do_this_test('partition', input_stype, k)]
        if len(ptest_list) == 0:
            return
        if input_stype not in self.perf_info[version_stype]:
            self.perf_info[version_stype][input_stype] = OrderedDict()
        pinfo = self.perf_info[version_stype][input_stype]
        if debug:
            print('  version %s input %s partitioning' % (version_stype, input_stype))
            print('  purity         completeness        test                    description')
        for ptest in ptest_list:
            if ptest not in pinfo:
                pinfo[ptest] = OrderedDict()
            if args.paired:
                l_ccfs = []
                for locus in utils.sub_loci(args.ig_or_tr):
                    ofn = '%s/partition-%s.yaml' % (self.opath(ptest, st=version_stype), locus)
                    if 'seed-' in ptest:
                        sfns = glob.glob(ofn.replace('/partition-', '/seeds/*/partition-'))
                        if len(sfns) == 0:  # non-seed light chain
                            continue
                        if len(sfns) > 1:
                            raise Exception('multiple seed subdirs in %s/seeds/, you probably changed the seed seq and need to delete the old dir (yes --bust-cache should do this, but doesn\'t atm)' % os.path.dirname(ofn) + '/seeds/')
                        ofn = sfns[0]
                    l_ccfs.append(read_cpath(ofn))
                pinfo[ptest]['purity'], pinfo[ptest]['completeness'] = [numpy.mean(lcfs) for lcfs in zip(*l_ccfs)]
                if 'seed-' not in ptest:
                    htmp = Hist(fname='%s/true-pair-clean-performance.csv' % self.opath(ptest, st=version_stype))
                    ttot = htmp.integral(False)
                    for pcat in self.pair_clean_metrics:
                        pinfo[ptest][pcat] = htmp.bin_contents[htmp.find_bin(None, label=pcat)] / float(ttot)
            else:
                ccfs = read_cpath(self.opath(ptest, st=version_stype))  # self.dirs(version_stype) + '/' + ptest + '.yaml')
                pinfo[ptest]['purity'], pinfo[ptest]['completeness'] = ccfs
                if ptest in self.perfdirs:
                    self.read_each_annotation_performance('single', version_stype, input_stype, these_are_cluster_annotations=True)

    # ----------------------------------------------------------------------------------------
    def read_selection_metric_performance(self, version_stype, input_stype, debug=False):
        # ----------------------------------------------------------------------------------------
        def read_smfile(fname, smfo):
            if not os.path.exists(fname):  # probably e.g. igh+igl for a sample with only igh+igk
                print('  %s selection metric output file doesn\'t exist: %s' % (utils.wrnstr(), fname))
                return
            with open(fname) as yfile:
                lbfos = yaml.load(yfile, Loader=yaml.CLoader)
            for metric in self.selection_metrics:
                for lbfo in lbfos:  # one lbfo for each cluster
                    smfo[metric] += list(lbfo['lb'][metric].values())
            if debug:
                print('      read lbfos for %d cluster%s from %s' % (len(lbfos), utils.plural(len(lbfos)), fname))
        # ----------------------------------------------------------------------------------------
        def read_chosen_abs(fname):
            with open(fname) as chfile:
                chlines = list(csv.DictReader(chfile))
            return set((l['h_id'], l['l_id']) for l in chlines)  # not really proper to key only by h_id, but the pairing info shouldn't be able to change, right?
        # ----------------------------------------------------------------------------------------
        pinfo = self.perf_info[version_stype][input_stype]
        if debug:
            print('  version %s input %s selection metrics' % (version_stype, input_stype))
        ptest_list = [k for k in self.tests.keys() if self.do_this_test('get-selection', input_stype, k)]
        for ptest in ptest_list:
            if ptest not in pinfo:  # perf_info should already have all the parent keys cause we run read_partition_performance() first
                pinfo[ptest] = OrderedDict([(m, []) for m in self.selection_metrics])
            if args.paired:
                for lpair in utils.locus_pairs[args.ig_or_tr]:
                    for locus in lpair:
                        smfname = '%s/%s/partition-%s-selection-metrics.yaml' % (self.opath(ptest.replace('get-selection-metrics', 'partition'), st=version_stype), '+'.join(lpair), locus)
                        read_smfile(smfname, pinfo[ptest])
                if debug:
                    print('  total values read: %s' % '  '.join('%s %d'%(m, len(pinfo[ptest][m])) for m in self.selection_metrics))
                pinfo[ptest]['chosen-abs'] = read_chosen_abs(self.opath(ptest, st=version_stype))
            else:
                read_smfile(self.opath(ptest, st=version_stype), pinfo[ptest])

    # ----------------------------------------------------------------------------------------
    def compare_performance(self, input_stype):
        # ----------------------------------------------------------------------------------------
        def print_comparison_str(ref_val, new_val, epsval, fw=7, dp=3, pm=False):
            fractional_change = 0. if ref_val == 0. else (new_val - ref_val) / float(ref_val)  # NOTE not the abs value yet
            if abs(fractional_change) > epsval:
                color = 'red'
            elif abs(fractional_change) > self.tiny_eps:
                color = 'yellow'
            else:
                color = None
            def floatstr(v):
                fmstr = '%%-%d.%df' % (fw, dp)
                if pm:
                    fmstr = fmstr.replace('%', '%+')
                return fmstr % v
            print('  %s%s  ' % (floatstr(ref_val), (fw+4)*' ' if color is None else utils.color(color, '--> %s'%floatstr(new_val))), end=' ')

        # ----------------------------------------------------------------------------------------
        print('  performance with %s simulation and parameters (smaller is better for all annotation metrics)' % input_stype)
        all_annotation_ptests = ['annotate-' + input_stype + '-simu', 'multi-annotate-' + input_stype + '-simu', 'partition-' + input_stype + '-simu']  # hard code for order
        all_partition_ptests = [flavor + 'partition-' + input_stype + '-simu' for flavor in ['', 'vsearch-', 'seed-', 'subset-']]
        annotation_ptests = [pt for pt in all_annotation_ptests if pt in self.perf_info['ref'][input_stype]]
        partition_ptests = [pt for pt in all_partition_ptests if pt in self.perf_info['ref'][input_stype]]
        selection_metric_tests = ['get-selection-metrics-'+input_stype+'-simu']
        metricstrs = {
            'mean_hamming' : 'hamming',
            'v_hamming' : 'v  ',
            'd_hamming' : 'd  ',
            'j_hamming' : 'j  ',
            'cdr3_hamming' : 'cdr3  ',
            'vd_insertion' : 'vd insert',
            'dj_insertion' : 'dj insert',
            'd_call' : 'd  ',
            'j_call' : 'j  ',
            'completeness' : 'compl.',
            'cons-dist-aa' : 'aa-cdist',
            'correct' : 'pair clean correct',
            'mispaired' : '         mispaired',
            'unpaired' : '          unpaired',
        }
        refpfo, newpfo = [self.perf_info[st][input_stype] for st in ['ref', 'new']]


        # print annotation header
        print('%8s %9s' % ('', ''), end=' ')
        for ptest in annotation_ptests:
            for method in [m for m in refpfo[ptest] if m in ['sw', 'hmm']]:  # 'if' is just to skip purity and completeness
                printstr = method
                if 'multi-annotate' in ptest:
                    printstr = 'multi %s' % method
                if 'partition' in ptest:
                    printstr = 'partition %s' % method
                print('    %-15s' % printstr, end=' ')
        print('')

        # print values
        if 'hmm' in refpfo[annotation_ptests[0]]:  # it's not in there for paired partition test, since (at least atm) we don't do annotation tests for it
            allmetrics = [m for m in refpfo[annotation_ptests[0]]['hmm']]
            for metric in allmetrics:
                alignstr = '' if len(metricstrs.get(metric, metric).strip()) < 5 else '-'
                print(('%8s %' + alignstr + '9s') % ('', metricstrs.get(metric, metric)), end=' ')
                for ptest in annotation_ptests:
                    for method in [m for m in refpfo[ptest] if m in ['sw', 'hmm']]:  # 'if' is just to skip purity and completeness
                        if set(refpfo[ptest]) != set(newpfo[ptest]):
                            raise Exception('different metrics in ref vs new:\n  %s\n  %s' % (sorted(refpfo[ptest]), sorted(newpfo[ptest])))
                        print_comparison_str(refpfo[ptest][method][metric], newpfo[ptest][method][metric], self.eps_vals.get(metric, 0.1))
                print('')

        # print partition header
        print('%8s %5s' % ('', ''), end=' ')
        for ptest in partition_ptests:
            print('    %-18s' % ptest.split('-')[0], end=' ')
        print('')
        for metric in ['purity', 'completeness'] + self.pair_clean_metrics:
            alignstr = '' if len(metricstrs.get(metric, metric).strip()) < 5 else '-'
            print(('%8s %' + alignstr + '9s') % ('', metricstrs.get(metric, metric)), end=' ')
            for ptest in partition_ptests:
                if 'seed-' in ptest and metric in self.pair_clean_metrics:  # ick
                    continue
                if set(refpfo[ptest]) != set(newpfo[ptest]):
                    raise Exception('different metrics in ref vs new:\n  %s\n  %s' % (sorted(refpfo[ptest]), sorted(newpfo[ptest])))
                method = ptest.split('-')[0]
                if metric != 'purity':
                    method = ''
                print_comparison_str(refpfo[ptest][metric], newpfo[ptest][metric], self.eps_vals.get(metric, 0.1))
            print('')

        # selection metrics
        print('                      %s' % ''.join(['%-23s'%metricstrs.get(m, m) for m in self.selection_metrics]))
        for mfname, mfcn in [('mean', numpy.mean), ('min', min), ('max', max), ('len', len)]:
            print('             %5s' % mfname, end=' ')
            for metric in self.selection_metrics:
                for ptest in selection_metric_tests:  # this'll break if there's more than one selection metric ptest
                    if set(refpfo[ptest]) != set(newpfo[ptest]):
                        raise Exception('different metrics in ref vs new:\n  %s\n  %s' % (sorted(refpfo[ptest]), sorted(newpfo[ptest])))
                    ref_list, new_list = [self.perf_info[rn][input_stype][ptest][metric] for rn in ['ref', 'new']]
                    dp = 1 if metric=='cons-dist-aa' else 3
                    if mfname=='len': dp = 0
                    print_comparison_str(mfcn(ref_list), mfcn(new_list), self.eps_vals.get(metric, 0.1), dp=dp, pm=metric=='cons-dist-aa')
            print('')
        if args.paired:
            ptest = utils.get_single_entry(selection_metric_tests)
            ref_abs, new_abs = refpfo[ptest]['chosen-abs'], newpfo[ptest]['chosen-abs']
            n_ref, n_new = len(ref_abs), len(new_abs)
            n_common = len(ref_abs & new_abs)
            n_only_ref, n_only_new = len(ref_abs - new_abs), len(new_abs - ref_abs)
            diffstr = '    ok'
            if n_ref != n_new or n_only_ref > 0 or n_only_new > 0:
                diffstr = '       %s in common, %s only in ref, %s only in new' % (utils.color(None if n_common==n_ref else 'red', str(n_common)), utils.color(None if n_only_ref==0 else 'red',  str(n_only_ref)), utils.color(None if n_only_new==0 else 'red',  str(n_only_new)))
            print('    chose %d abs %s%s' % (n_ref, '' if n_new==n_ref else utils.color('red', '--> %d'%n_new), diffstr))

    # ----------------------------------------------------------------------------------------
    def compare_production_results(self, ptests):
        print('diffing production results')
        for ptest in ptests:
            if args.quick and ptest not in self.quick_tests:
                continue
            fname = self.opath(ptest)  # sometimes a dir rather than a file
            print('    %-30s' % fname, end=' ')
            cmd = 'diff -qbr ' + ' '.join(self.dirs(st) + '/' + fname for st in self.stypes)
            proc = Popen(cmd.split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
            out, err = proc.communicate()
            if proc.returncode == 0:
                print('       ok')
            else:
                differlines = [ l for l in out.split('\n') if 'differ' in l]
                onlylines = [ l for l in out.split('\n') if 'Only' in l]
                print('')
                if len(differlines) > 0:
                    n_total_files = int(check_output('find ' + self.dirs('ref') + '/' + fname + ' -type f | wc -l', shell=True, universal_newlines=True))
                    if n_total_files == 1:
                        assert len(differlines) == 1
                        print(utils.color('red', '      file differs'), end=' ')
                    else:
                        print(utils.color('red', '      %d / %d files differ' % (len(differlines), n_total_files)), end=' ')
                if len(onlylines) > 0:
                    for st in self.stypes:
                        theseonlylines = [l for l in onlylines if self.dirs(st) + '/' + fname in l]
                        if len(theseonlylines) > 0:
                            print(utils.color('red', '      %d files only in %s' % (len(theseonlylines), st)), end=' ')
                if differlines == 0 and onlylines == 0:
                    print(utils.color('red', '      not sure why, but diff returned %d' % proc.returncode), end=' ')
                print('  (%s)' % cmd)
                if err != '':
                    print(err)

    # ----------------------------------------------------------------------------------------
    def write_run_times(self):
        with open(self.dirs('new') + '/run-times.csv', utils.csv_wmode()) as newfile:
            writer = csv.DictWriter(newfile, ('name', 'seconds'))
            writer.writeheader()
            for name, seconds in self.run_times.items():
                writer.writerow({'name' : name, 'seconds' : '%.1f'%seconds})

    # ----------------------------------------------------------------------------------------
    def compare_run_times(self):
        print('checking run times')

        def read_run_times(stype):
            times[stype] = {}
            with open(self.dirs(stype) + '/run-times.csv') as timefile:
                reader = csv.DictReader(timefile)
                for line in reader:
                    times[stype][line['name']] = float(line['seconds'])
        times = {}
        for stype in self.stypes:
            read_run_times(stype)

        for name in self.tests:
            if args.quick and name not in self.quick_tests:
                continue
            print('  %30s   %7.1f' % (name, times['ref'][name]), end=' ')
            if name not in times['new']:
                print('  no new time for %s' % utils.color('red', name))
                continue
            fractional_change = (times['new'][name] - times['ref'][name]) / float(times['ref'][name])
            if abs(fractional_change) > 0.2:
                print('--> %-5.1f %s' % (times['new'][name], utils.color('red', '(%+.3f)' % fractional_change)), end=' ')
            elif abs(fractional_change) > 0.1:
                print('--> %-5.1f %s' % (times['new'][name], utils.color('yellow', '(%+.3f)' % fractional_change)), end=' ')
            else:
                print('    ok   ', end=' ')
            print('')

    # ----------------------------------------------------------------------------------------
    def compare_partition_cachefiles(self, input_stype, debug=False):
        # ----------------------------------------------------------------------------------------
        def print_key_differences(vtype, refkeys, newkeys):
            print('    %s keys' % vtype)
            if len(refkeys - newkeys) > 0 or len(newkeys - refkeys) > 0:
                if len(refkeys - newkeys) > 0:
                    print(utils.color('red', '      %d only in ref version' % len(refkeys - newkeys)))
                if len(newkeys - refkeys) > 0:
                    print(utils.color('red', '      %d only in new version' % len(newkeys - refkeys)))
                print('      %d in common' % len(refkeys & newkeys))
            else:
                print('        %d identical keys in new and ref cache' % len(refkeys))
        # ----------------------------------------------------------------------------------------
        def readcache(fname):
            if debug: print('      reading partition cache from %s' % fname)
            cache = {'naive_seqs' : {}, 'logprobs' : {}}
            with open(fname) as cachefile:
                reader = csv.DictReader(cachefile)
                for line in reader:
                    if line['naive_seq'] != '':
                        cache['naive_seqs'][line['unique_ids']] = line['naive_seq']
                    if line['logprob'] != '':
                        cache['logprobs'][line['unique_ids']] = float(line['logprob'])
            return cache
        # ----------------------------------------------------------------------------------------
        def compare_files(fname):
            print('  %s input partition cache file' % input_stype)
            refcache = readcache(self.dirs('ref') + '/' + fname)
            newcache = readcache(self.dirs('new') + '/' + fname)

            # work out intersection and complement
            refkeys, newkeys = {}, {}
            for vtype in ['naive_seqs', 'logprobs']:
                refkeys[vtype] = set(refcache[vtype].keys())
                newkeys[vtype] = set(newcache[vtype].keys())
                print_key_differences(vtype, refkeys[vtype], newkeys[vtype])

            hammings = []
            n_hammings = 0
            n_different_length, n_big_hammings = 0, 0
            hamming_eps = 0.
            vtype = 'naive_seqs'
            for uids in refkeys[vtype] & newkeys[vtype]:
                refseq = refcache[vtype][uids]
                newseq = newcache[vtype][uids]
                n_hammings += 1
                if len(refseq) == len(newseq):
                    hamming_fraction = utils.hamming_fraction(refseq, newseq)
                    if hamming_fraction > hamming_eps:
                        n_big_hammings += 1
                        hammings.append(hamming_fraction)
                else:
                    n_different_length += 1

            diff_hfracs_str = '%3d / %4d' % (n_big_hammings, n_hammings)
            mean_hfrac_str = '%.3f' % (numpy.average(hammings) if len(hammings) > 0 else 0.)
            if n_big_hammings > 0:
                diff_hfracs_str = utils.color('red', diff_hfracs_str)
                mean_hfrac_str = utils.color('red', mean_hfrac_str)

            abs_delta_logprobs = []
            n_delta_logprobs = 0
            n_big_delta_logprobs = 0
            logprob_eps = 1e-5
            vtype = 'logprobs'
            for uids in refkeys[vtype] & newkeys[vtype]:
                refval = refcache[vtype][uids]
                newval = newcache[vtype][uids]
                n_delta_logprobs += 1
                abs_delta_logprob = abs(refval - newval)
                if abs_delta_logprob > logprob_eps:
                    # print '%s  %s  ref  %f  new %f' % (vtype, uids, refval, newval)
                    n_big_delta_logprobs += 1
                    abs_delta_logprobs.append(abs_delta_logprob)

            diff_logprob_str = '%3d / %4d' % (n_big_delta_logprobs, n_delta_logprobs)
            mean_logprob_str = '%.3f' % (numpy.average(abs_delta_logprobs) if len(abs_delta_logprobs) > 0 else 0.)
            if n_big_delta_logprobs > 0:
                diff_logprob_str = utils.color('red', diff_logprob_str)
                mean_logprob_str = utils.color('red', mean_logprob_str)
            print('                  fraction different     mean abs difference among differents')
            print('      naive seqs     %s                      %s      (hamming fraction)' % (diff_hfracs_str, mean_hfrac_str))
            print('      log probs      %s                      %s' % (diff_logprob_str, mean_logprob_str))
            if n_different_length > 0:
                print(utils.color('red', '        %d different length' % n_different_length))

        # ----------------------------------------------------------------------------------------
        ptest = 'partition-' + input_stype + '-simu'
        if args.quick and ptest not in self.quick_tests:
            return

        if args.paired:
            assert False  # eh, probably not really any point
            # for locus in XXX:
            # compare_files(self.ptn_cachefn(input_stype))
        else:
            compare_files(self.ptn_cachefn(input_stype))

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--dont-run', action='store_true', help='don\'t actually run anything, just check the results')
parser.add_argument('--dry-run', action='store_true', help='do all preparations to run, but don\'t actually run the commands, and don\'t check results')
parser.add_argument('--quick', action='store_true', help='only run one command: cache-parameters on a small numbrer of simulation events')
parser.add_argument('--slow', action='store_true', help='by default, we run tests on a fairly small number of sequences, which is sufficient for checking that *nothing* has changed. But --slow is for cases where you\'ve made changes that you know will affect results, and you want to look at the details of how they\'re affected, for which you need to run on more sequences. Note that whether --slow is set or not (runs all tests with more or less sequences) is separate from --quick (which only runs one test).')
parser.add_argument('--run-coverage', action='store_true', help='instead of running the normal suite of tests (which compare results to make sure they haven\'t changed), instead run a series of commands that\'s designed to execute as many of the lines as possible (without comparing results).')
parser.add_argument('--prepend-coverage', action='store_true', help='run normal tests, but prepending coverage append commands')
parser.add_argument('--coverage-outdir', default='%s/partis/tmp/coverage' % os.getenv('fs', default=os.getenv('HOME')))
parser.add_argument('--bust-cache', action='store_true', help='overwrite current ref info, i.e. run this when things have changed, but you\'ve decided they\'re fine')
parser.add_argument('--only-bust-current', action='store_true', help='only bust cache for current command line args (as opposed to the default of busting caches for both slow and non-slow, paired and non-paired)')
parser.add_argument('--paired', action='store_true', help='run paired tests, i.e. with --paired-loci. Note that this doesn\'t test all the things (e.g. seed partitioning) that non-paired does.')
parser.add_argument('--run-all', action='store_true', help='run all four combinations of tests: paired/non-paired and slow/non-slow (by default only runs one). *Not* for use with --bust-cache, which runs all of them by default.')
parser.add_argument('--no-simu', action='store_true', help='don\'t run simulation, e.g. if using a minimal install')
parser.add_argument('--no-tree-gen', action='store_true', help='if set, use the pre-generated newick file instead of generating trees with R/TreeSim (to avoid R installation)')
parser.add_argument('--no-per-base-mutation', action='store_true', help='use simpler, non-per-base mutation model (to avoid bpp-newlik compilation)')
parser.add_argument('--ig-or-tr', default='ig')
parser.add_argument('--print-width', type=int, default=300, help='set to 0 for infinite')

parser.add_argument('--glfo-dir', default='data/germlines/human')
parser.add_argument('--locus', default='igh')
args = parser.parse_args()
assert not (args.quick and args.slow)  # it just doesn't make sense
assert not (args.quick and args.paired)  # --quick ignores --paired, which is confusing

random.seed(0)
numpy.random.seed(0)

if args.print_width == 0:
    args.print_width = 99999

if args.run_all or (args.bust_cache and not args.only_bust_current):  # run all four combos
    for slowval in [False, True]:
        for pairedval in [False, True]:
            clist = copy.deepcopy(sys.argv)
            utils.remove_from_arglist(clist, '--slow')
            utils.remove_from_arglist(clist, '--paired')
            if args.bust_cache:
                assert not args.run_all
                clist += ['--only-bust-current']
            else:
                utils.remove_from_arglist(clist, '--run-all')
            cmd_str = ' '.join(clist)
            if slowval:
                cmd_str += ' --slow'
            if pairedval:
                cmd_str += ' --paired'
            utils.simplerun(cmd_str, dryrun=args.dry_run)
    sys.exit(0)

tester = Tester()
if args.run_coverage:
    tester.run_coverage(args)
    sys.exit(0)
tester.test(args)
if args.bust_cache:
    tester.bust_cache()

# ----------------------------------------------------------------------------------------
def get_typical_variances():
    # NOTE don't delete this, since it was used (and might be needed again) to get the expected variances hardcoded above
    raise Exception('needs updating to work as a function')
    # cp = ClusterPath(fname='tmp.csv')
    # cp.print_partitions()
    # sys.exit()
    # cps = []
    # adj_mis, ccf_unders, ccf_overs = [], [], []
    # for iseed in range(6):
    #     # print 'seed %d' % iseed
    #     cp = ClusterPath(fname='%d.csv' % iseed)
    #     cp.print_partitions()  #(cp.i_best)  #, abbreviate=False)
    #     adj_mis.append(cp.adj_mis[cp.i_best])
    #     ccf_unders.append(cp.ccfs[cp.i_best][0])
    #     ccf_overs.append(cp.ccfs[cp.i_best][1])
    #     cps.append(cp)
    # def print_mean_variance(vals):
    #     mean = numpy.average(vals)
    #     variance = numpy.average((vals - mean)**2)  #, weights=wgts)
    #     print 'mean %.2f   std dev %.3f   (%.1f%%)' % (mean, math.sqrt(variance), 100. * math.sqrt(variance) / mean)

    # # mean/var for six random seeds
    # print_mean_variance(adj_mis)     # mean 0.61   std dev 0.053   (8.7%)
    # print_mean_variance(ccf_unders)  # mean 0.74   std dev 0.026   (3.5%)
    # print_mean_variance(ccf_overs)   # mean 0.90   std dev 0.015   (1.7%)
    # # for iseed in range(len(cps)):
    # #     icp = cps[iseed]
    # #     for jseed in range(iseed, len(cps)):
    # #         jcp = cps[jseed]
    # #         print '  %d %d   %.3f' % (iseed, jseed, utils.adjusted_mutual_information(icp.partitions[icp.i_best], jcp.partitions[jcp.i_best]))

