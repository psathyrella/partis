#!/usr/bin/env python
import argparse
import numpy
import os
import csv
import glob
import math
import shutil
import time
from collections import OrderedDict
from subprocess import Popen, PIPE, check_call, check_output, CalledProcessError
import colored_traceback.always
import sys
sys.path.insert(1, './python')

from baseutils import get_extra_str
import utils
import glutils
from hist import Hist
from clusterpath import ClusterPath

# ----------------------------------------------------------------------------------------
class Tester(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self):
        self.partis = './bin/partis'
        self.datafname = 'test/mishmash.fa'  # some data from adaptive, chaim, and vollmers
        self.label = 'test'

        self.stypes = ['ref', 'new']  # I don't know what the 's' stands for
        self.dtypes = ['data', 'simu']
        self.dirs = {'ref' : 'test/reference-results', 'new' : 'test/new-results'}
        if not os.path.exists(self.dirs['new']):
            os.makedirs(self.dirs['new'])
        self.infnames = {st : {dt : self.datafname if dt == 'data' else self.dirs[st] + '/' + self.label + '/simu.yaml' for dt in self.dtypes} for st in self.stypes}
        self.param_dirs = { st : { dt : self.dirs[st] + '/' + self.label + '/parameters/' + dt for dt in self.dtypes} for st in self.stypes}  # muddafuggincomprehensiongansta
        self.common_extras = ['--seed', '1', '--n-procs', '10']  # would be nice to set --n-procs based on the machine, but for some reason the order of things in the parameter files gets shuffled a bit if you changed the number of procs

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

        self.n_partition_queries = '500'
        self.n_sim_events = '500'
        self.logfname = self.dirs['new'] + '/test.log'
        self.sw_cache_paths = {st : {dt : self.param_dirs[st][dt] + '/sw-cache' for dt in self.dtypes} for st in self.stypes}  # don't yet know the 'new' ones (they'll be the same only if the simulation is the same) #self.stypes}
        self.cachefnames = { st : 'cache-' + st + '-partition.csv' for st in ['new']}  # self.stypes

        self.tests = OrderedDict()

        def add_inference_tests(input_stype):  # if input_stype is 'ref', infer on old simulation and parameters, if it's 'new' use the new ones
            self.tests['annotate-' + input_stype + '-simu']          = {'extras' : ['--plot-annotation-performance', ]}
            self.tests['multi-annotate-' + input_stype + '-simu']    = {'extras' : ['--plot-annotation-performance', '--simultaneous-true-clonal-seqs']}  # NOTE this is mostly different to the multi-seq annotations from the partition step because it uses the whole sample
            self.tests['partition-' + input_stype + '-simu']         = {'extras' : [
                '--n-max-queries', self.n_partition_queries,
                '--n-precache-procs', '10',
                '--plot-annotation-performance',
                # '--biggest-logprob-cluster-to-calculate', '5', '--biggest-naive-seq-cluster-to-calculate', '5',
            ]}
            self.tests['seed-partition-' + input_stype + '-simu']    = {'extras' : ['--n-max-queries', self.n_partition_queries]}
            self.tests['vsearch-partition-' + input_stype + '-simu'] = {'extras' : ['--naive-vsearch', '--n-max-queries', self.n_partition_queries]}

        self.tests['cache-parameters-simu']  = {'extras' : [], 'output-path' : 'test/parameters/simu'}
        add_inference_tests('new')
        self.tests['cache-parameters-data']  = {'extras' : [], 'output-path' : 'test/parameters/data'}
        self.tests['simulate']  = {'extras' : ['--n-sim-events', self.n_sim_events, '--n-trees', self.n_sim_events, '--n-leaves', '5'], 'output-path' : 'test/simu.yaml'}

        self.quick_tests = ['cache-parameters-simu', 'annotate-new-simu']

        self.perfdirs = {}  # set in fiddle_with_arguments() NOTE these correspond only to annotation performance, whereas <self.perf_info> has also partition performance
        for ptest, argfo in self.tests.items():
            self.fiddle_with_arguments(ptest, argfo)

        self.perf_info = {version_stype : {} for version_stype in self.stypes}

    # ----------------------------------------------------------------------------------------
    def test(self, args):
        if args.comparison_plots:
            assert False  # needs updating
            self.make_comparison_plots()
            return
        if not args.dont_run:
            self.run(args)
        if args.dryrun or args.bust_cache:
            return
        self.compare_production_results(['cache-parameters-simu'])
        self.compare_stuff(input_stype='new')
        self.compare_production_results(['cache-parameters-data', 'simulate'])
        self.compare_run_times()

    # ----------------------------------------------------------------------------------------
    def fiddle_with_arguments(self, ptest, argfo):
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
        argfo['input_stype'] = input_stype
        argfo['bin'] = self.partis

        if ptest == 'simulate':
            argfo['parameter-dir'] = self.param_dirs[input_stype]['data']
        else:
            argfo['infname'] = self.infnames['new' if args.bust_cache else 'ref'][input_dtype]
            argfo['parameter-dir'] = self.param_dirs[input_stype][input_dtype]
            if input_dtype == 'simu':
                argfo['extras'] += ['--is-simu']
            argfo['sw-cachefname'] = self.sw_cache_paths[input_stype][input_dtype] + '.yaml'

        if 'annotate' in ptest:
            argfo['action'] = 'annotate'
        elif 'partition' in ptest:
            argfo['action'] = 'partition'
            argfo['extras'] += ['--persistent-cachefname', self.dirs['new'] + '/' + self.cachefnames[input_stype]]
        elif 'cache-parameters-' in ptest:
            argfo['action'] = 'cache-parameters'
            # if True:  #args.make_plots:
            #     argfo['extras'] += ['--plotdir', self.dirs['new'] + '/' + self.label + '/plots/' + input_dtype]
        else:
            argfo['action'] = ptest

        if '--plot-annotation-performance' in argfo['extras']:
            self.perfdirs[ptest] = ptest + '-annotation-performance'
            argfo['extras'] += ['--plotdir', self.dirs['new'] + '/' + self.perfdirs[ptest]]

        if '--plotdir' in argfo['extras']:
            argfo['extras'] += ['--only-csv-plots', '--only-overall-plots']

    # ----------------------------------------------------------------------------------------
    def compare_stuff(self, input_stype):
        print '%s input' % input_stype
        for version_stype in self.stypes:  # <version_stype> is the code version, i.e. 'ref' is the reference results, 'new' is the results we just made with the new code
            self.read_annotation_performance(version_stype, input_stype)
            self.read_partition_performance(version_stype, input_stype)  # NOTE also calls read_annotation_performance()
        self.compare_performance(input_stype)
        self.compare_partition_cachefiles(input_stype)

    # ----------------------------------------------------------------------------------------
    def prepare_to_run(self, args, name, info):
        """ Pre-run stuff that you don't want to do until *right* before you actually run. """

        # delete old partition cache file
        if name == 'partition-' + info['input_stype'] + '-simu':
            this_cachefname = self.dirs['new'] + '/' + self.cachefnames[info['input_stype']]
            if os.path.exists(this_cachefname):
                if args.dryrun:
                    print '   would remove %s' % this_cachefname
                else:
                    check_call(['rm', '-v', this_cachefname])

        # choose a seed uid
        if name == 'seed-partition-' + info['input_stype'] + '-simu':
            seed_uid, _ = utils.choose_seed_unique_id(info['infname'], 5, 8, n_max_queries=int(self.n_partition_queries), debug=False)
            info['extras'] += ['--seed-unique-id', seed_uid]

    # ----------------------------------------------------------------------------------------
    def run(self, args):
        if not args.dryrun:
            open(self.logfname, 'w').close()

        for name, info in self.tests.items():
            if args.quick and name not in self.quick_tests:
                continue

            self.prepare_to_run(args, name, info)

            action = info['action']
            cmd_str = info['bin'] + ' ' + action
            if 'infname' in info:
                cmd_str += ' --infname %s' % info['infname']
            cmd_str += ' --parameter-dir %s' % info['parameter-dir']
            if 'sw-cachefname' in info:
                cmd_str += ' --sw-cachefname %s' % info['sw-cachefname']
            cmd_str += ' %s' % ' '.join(info['extras'] + self.common_extras)

            if name == 'simulate':
                cmd_str += ' --outfname ' + self.infnames['new']['simu']
                cmd_str += ' --indel-frequency 0.01 --indel-location v'
            elif 'cache-parameters-' not in name:
                cmd_str += ' --outfname ' + self.dirs['new'] + '/' + name + '.yaml'

            logstr = '%s   %s' % (utils.color('green', name, width=30, padside='right'), cmd_str)
            print logstr if utils.len_excluding_colors(logstr) < args.print_width else logstr[:args.print_width] + '[...]'
            if args.dryrun:
                continue
            logfile = open(self.logfname, 'a')
            logfile.write(logstr + '\n')
            logfile.close()
            start = time.time()
            try:
                check_call(cmd_str + ' 1>>' + self.logfname + ' 2>>' + self.logfname, shell=True)
            except CalledProcessError, err:
                # print err  # this just says it exited with code != 0
                print '  log tail:'
                print utils.pad_lines(check_output(['tail', self.logfname]))
                sys.exit(1)  # raise Exception('exited with error')
            self.run_times[name] = time.time() - start  # seconds

        self.write_run_times()

    # ----------------------------------------------------------------------------------------
    def remove_reference_results(self, expected_content):
        print '  removing ref files'
        dir_content = set([os.path.basename(f) for f in glob.glob(self.dirs['ref'] + '/*')])
        if len(dir_content - expected_content) > 0 or len(expected_content - dir_content) > 0:
            if len(dir_content - expected_content) > 0:
                print 'in ref dir but not expected\n    %s' % (utils.color('red', ' '.join(dir_content - expected_content)))
            if len(expected_content - dir_content) > 0:
                print 'expected but not in ref dir\n    %s' % (utils.color('red', ' '.join(expected_content - dir_content)))
            raise Exception('unexpected or missing content in reference dir (see above)')
        for fname in [self.dirs['ref'] + '/' + ec for ec in expected_content]:
            print '    rm %s' % fname
            if args.dryrun:
                continue
            if os.path.isdir(fname):
                shutil.rmtree(fname)
            else:
                os.remove(fname)

    # ----------------------------------------------------------------------------------------
    def bust_cache(self):
        test_outputs = [k + '.yaml' for k, tfo in self.tests.items() if 'output-path' not in tfo]
        expected_content = set(test_outputs + self.perfdirs.values() + self.cachefnames.values() + [os.path.basename(self.logfname), self.label])
        expected_content.add('run-times.csv')

        # remove (very, very gingerly) whole reference dir
        self.remove_reference_results(expected_content)

        # copy over parameters, simulation, and plots
        print '  copy new files to ref'
        for fname in expected_content:
            print '    mv %s   -->  %s/' % (fname, self.dirs['ref'])
            if args.dryrun:
                continue
            shutil.move(self.dirs['new'] + '/' + fname, self.dirs['ref'] + '/')

    # ----------------------------------------------------------------------------------------
    def read_annotation_performance(self, version_stype, input_stype, these_are_cluster_annotations=False):  # <these_are_cluster_annotations> means this fcn is being called from within read_partition_performance()
        for sequence_multiplicity in ['single', 'multi']:
            self.read_each_annotation_performance(sequence_multiplicity, version_stype, input_stype, these_are_cluster_annotations=these_are_cluster_annotations)

    # ----------------------------------------------------------------------------------------
    def read_each_annotation_performance(self, sequence_multiplicity, version_stype, input_stype, these_are_cluster_annotations=False):  # <these_are_cluster_annotations> means this fcn is being called from within read_partition_performance()
        """ version_stype is the code version, while input_stype is the input data version, i.e. 'ref', 'new' is the reference code version (last commit) run on the then-new simulation and parameters"""
        def read_performance_file(fname, column, only_ibin=None):
            values = []
            with open(fname) as csvfile:
                reader = csv.DictReader(csvfile)
                ibin = 0
                for line in reader:
                    if only_ibin is not None and ibin != only_ibin:
                        ibin += 1
                        continue
                    values.append(float(line[column]))
                    ibin += 1
            if len(values) == 1:
                return values[0]
            else:
                return values

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
        if args.quick and ptest not in self.quick_tests:
            return
        if input_stype not in self.perf_info[version_stype]:
            self.perf_info[version_stype][input_stype] = OrderedDict()
        if ptest not in self.perf_info[version_stype][input_stype]:
            self.perf_info[version_stype][input_stype][ptest] = OrderedDict()
        perfdir = self.dirs[version_stype] + '/' + self.perfdirs[ptest]
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
            # for region in utils.regions:
            #     fraction_correct = read_performance_file(perfdir + '/' + method + '/gene-call/' + region + '_gene.csv', 'contents', only_ibin=1)
            #     perffo[method][region + '_call'] = fraction_correct

    # ----------------------------------------------------------------------------------------
    def read_partition_performance(self, version_stype, input_stype, debug=False):
        """ Read new partitions from self.dirs['new'], and put the comparison numbers in self.perf_info (compare either to true, for simulation, or to the partition in reference dir, for data). """
        def do_this_test(pt):
            if 'partition' not in pt:
                return False
            if input_stype not in pt:
                return False
            if args.quick and pt not in self.quick_tests:
                return False
            return True

        ptest_list = [k for k in self.tests.keys() if do_this_test(k)]
        if len(ptest_list) == 0:
            return
        if input_stype not in self.perf_info[version_stype]:
            self.perf_info[version_stype][input_stype] = OrderedDict()
        if debug:
            print '  version %s input %s partitioning' % (version_stype, input_stype)
            print '  purity         completeness        test                    description'
        for ptest in ptest_list:
            if ptest not in self.perf_info[version_stype][input_stype]:
                self.perf_info[version_stype][input_stype][ptest] = OrderedDict()
            _, _, cpath = utils.read_yaml_output(fname=self.dirs[version_stype] + '/' + ptest + '.yaml', skip_annotations=True)
            ccfs = cpath.ccfs[cpath.i_best]
            if None in ccfs:
                raise Exception('none type ccf read from %s' % self.dirs[version_stype] + '/' + ptest + '.yaml')
            self.perf_info[version_stype][input_stype][ptest]['purity'], self.perf_info[version_stype][input_stype][ptest]['completeness'] = ccfs
            if debug:
                print '    %5.2f          %5.2f      %-28s   to true partition' % (self.perf_info[version_stype][input_stype][ptest]['purity'], self.perf_info[version_stype][input_stype][ptest]['completeness'], ptest)

            if ptest in self.perfdirs:
                self.read_each_annotation_performance('single', version_stype, input_stype, these_are_cluster_annotations=True)

    # ----------------------------------------------------------------------------------------
    def compare_performance(self, input_stype):

        def print_comparison_str(ref_val, new_val, epsval):
            printstr = '%-5.3f' % ref_val
            fractional_change = 0. if ref_val == 0. else (new_val - ref_val) / ref_val  # NOTE not the abs value yet
            color = None
            if abs(fractional_change) > epsval:
                color = 'red'
            elif abs(fractional_change) > self.tiny_eps:
                color = 'yellow'
            if color is None:
                printstr += '          '
            else:
                printstr += utils.color(color, ' --> %-5.3f' % new_val)
            print '    %-s' % printstr,

        print '  performance with %s simulation and parameters (smaller is better for all annotation metrics)' % input_stype
        all_annotation_ptests = ['annotate-' + input_stype + '-simu', 'multi-annotate-' + input_stype + '-simu', 'partition-' + input_stype + '-simu']  # hard code for order
        all_partition_ptests = [flavor + 'partition-' + input_stype + '-simu' for flavor in ['', 'vsearch-', 'seed-']]
        annotation_ptests = [pt for pt in all_annotation_ptests if pt in self.perf_info['ref'][input_stype]]
        partition_ptests = [pt for pt in all_partition_ptests if pt in self.perf_info['ref'][input_stype]]
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
        }


        # print annotation header
        print '%8s %9s' % ('', ''),
        for ptest in annotation_ptests:
            for method in [m for m in self.perf_info['ref'][input_stype][ptest] if m in ['sw', 'hmm']]:  # 'if' is just to skip purity and completeness
                printstr = method
                if 'multi-annotate' in ptest:
                    printstr = 'multi %s' % method
                if 'partition' in ptest:
                    printstr = 'partition %s' % method
                print '    %-15s' % printstr,
        print ''

        # print values
        allmetrics = [m for m in self.perf_info['ref'][input_stype][annotation_ptests[0]]['hmm']]
        for metric in allmetrics:
            alignstr = '' if len(metricstrs.get(metric, metric).strip()) < 5 else '-'
            print ('%8s %' + alignstr + '9s') % ('', metricstrs.get(metric, metric)),
            for ptest in annotation_ptests:
                for method in [m for m in self.perf_info['ref'][input_stype][ptest] if m in ['sw', 'hmm']]:  # 'if' is just to skip purity and completeness
                    if set(self.perf_info['ref'][input_stype][ptest]) != set(self.perf_info['new'][input_stype][ptest]):
                        raise Exception('different metrics in ref vs new:\n  %s\n  %s' % (sorted(self.perf_info['ref'][input_stype][ptest]), sorted(self.perf_info['new'][input_stype][ptest])))
                    print_comparison_str(self.perf_info['ref'][input_stype][ptest][method][metric], self.perf_info['new'][input_stype][ptest][method][metric], self.eps_vals.get(metric, 0.1))
            print ''

        # print partition header
        print '%8s %7s' % ('', ''),
        for ptest in partition_ptests:
            print '    %-15s' % ptest.split('-')[0],
        print ''
        for metric in ['purity', 'completeness']:
            alignstr = '' if len(metricstrs.get(metric, metric).strip()) < 5 else '-'
            print ('%8s %' + alignstr + '9s') % ('', metricstrs.get(metric, metric)),
            for ptest in partition_ptests:
                if set(self.perf_info['ref'][input_stype][ptest]) != set(self.perf_info['new'][input_stype][ptest]):
                    raise Exception('different metrics in ref vs new:\n  %s\n  %s' % (sorted(self.perf_info['ref'][input_stype][ptest]), sorted(self.perf_info['new'][input_stype][ptest])))
                method = ptest.split('-')[0]
                if metric != 'purity':
                    method = ''
                print_comparison_str(self.perf_info['ref'][input_stype][ptest][metric], self.perf_info['new'][input_stype][ptest][metric], self.eps_vals.get(metric, 0.1))
            print ''

    # ----------------------------------------------------------------------------------------
    def compare_production_results(self, ptests):
        print 'diffing production results'
        for ptest in ptests:
            if args.quick and ptest not in self.quick_tests:
                continue
            fname = self.tests[ptest]['output-path']  # this is kind of messy
            print '    %-30s' % fname,
            cmd = 'diff -qbr ' + ' '.join(self.dirs[st] + '/' + fname for st in self.stypes)
            proc = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate()
            if proc.returncode == 0:
                print '       ok'
            else:
                differlines = [ l for l in out.split('\n') if 'differ' in l]
                onlylines = [ l for l in out.split('\n') if 'Only' in l]
                print ''
                if len(differlines) > 0:
                    n_total_files = int(check_output('find ' + self.dirs['ref'] + '/' + fname + ' -type f | wc -l', shell=True))
                    if n_total_files == 1:
                        assert len(differlines) == 1
                        print utils.color('red', '      file differs'),
                    else:
                        print utils.color('red', '      %d / %d files differ' % (len(differlines), n_total_files)),
                if len(onlylines) > 0:
                    for st in self.stypes:
                        theseonlylines = [l for l in onlylines if self.dirs[st] + '/' + fname in l]
                        if len(theseonlylines) > 0:
                            print utils.color('red', '      %d files only in %s' % (len(theseonlylines), st)),
                if differlines == 0 and onlylines == 0:
                    print utils.color('red', '      not sure why, but diff returned %d' % proc.returncode),
                print '  (%s)' % cmd
                if err != '':
                    print err

    # ----------------------------------------------------------------------------------------
    def write_run_times(self):
        with open(self.dirs['new'] + '/run-times.csv', 'w') as newfile:
            writer = csv.DictWriter(newfile, ('name', 'seconds'))
            writer.writeheader()
            for name, seconds in self.run_times.items():
                writer.writerow({'name' : name, 'seconds' : seconds})

    # ----------------------------------------------------------------------------------------
    def compare_run_times(self):
        print 'checking run times'

        def read_run_times(stype):
            times[stype] = {}
            with open(self.dirs[stype] + '/run-times.csv') as timefile:
                reader = csv.DictReader(timefile)
                for line in reader:
                    times[stype][line['name']] = float(line['seconds'])
        times = {}
        for stype in self.stypes:
            read_run_times(stype)

        for name in self.tests:
            if args.quick and name not in self.quick_tests:
                continue
            print '  %30s   %7.1f' % (name, times['ref'][name]),
            if name not in times['new']:
                print '  no new time for %s' % utils.color('red', name)
                continue
            fractional_change = (times['new'][name] - times['ref'][name]) / times['ref'][name]
            if abs(fractional_change) > 0.2:
                print '--> %-5.1f %s' % (times['new'][name], utils.color('red', '(%+.3f)' % fractional_change)),
            elif abs(fractional_change) > 0.1:
                print '--> %-5.1f %s' % (times['new'][name], utils.color('yellow', '(%+.3f)' % fractional_change)),
            else:
                print '    ok   ',
            print ''

    # ----------------------------------------------------------------------------------------
    def make_comparison_plots(self):
        assert False  # self.perfdirs treatment (among probably lots of other things) needs work
        plotdirs = [
            # self.perfdirs['ref'] + '/sw',  # ref sw performance
            # self.perfdirs['ref'] + '/hmm', # ref hmm performance
            'test/plots/data/sw/overall',          # sw data parameters
            'test/plots/data/hmm/overall',         # hmm data parameters
            'test/plots/data/hmm/mute-freqs/overall',
            'test/plots/data/sw/mute-freqs/overall'
            # 'test/plots/data/hmm/mute-freqs'
            # 'test/plots/simu/hmm-true',    # true simulation parameters
            # 'test/plots/simu/hmm-true/mute-freqs'
        ]
        base_check_cmd = './bin/compare.py --linewidths 7:2 --alphas 0.55:1 --str-colors #006600:#990012 --markersizes 3:1 --names reference:new'  # --colors 595:807:834 --scale-errors 1.414 "
        if os.getenv('www') is None:
            www_dir = '_test-plots'
        else:
            www_dir = os.getenv('www') + '/partis/test'

        # # if you want to do *all* the subdirs use this:
        # recursive_subdirs = []
        # for plotdir in plotdirs:
        #     find_plotdirs_cmd = 'find ' + self.dirs['ref'] + '/' + plotdir + ' -name "*.csv" -exec dirname {} \;|sed \'s@/plots$@@\' | sort | uniq'
        #     recursive_subdirs += check_output(find_plotdirs_cmd, shell=True).split()

        for plotdir in plotdirs:
            print plotdir
            check_cmd = base_check_cmd + ' --plotdirs '  + self.dirs['ref'] + '/' + plotdir + ':' + self.dirs['new'] + '/' + plotdir
            check_cmd += ' --outdir ' + www_dir + '/' + plotdir
            check_call(check_cmd.split())
        
#            # env.Command('test/_results/%s.passed' % name, out,
#            #             './bin/diff-parameters.py --arg1 test/regression/parameters/' + actions[name]['target'] + ' --arg2 ' + stashdir + '/test/' + actions[name]['target'] + ' && touch $TARGET')
        
            # ----------------------------------------------------------------------------------------
    def compare_partition_cachefiles(self, input_stype):
        """ NOTE only writing this for the ref input_stype a.t.m. """

        # ----------------------------------------------------------------------------------------
        def print_key_differences(vtype, refkeys, newkeys):
            print '    %s keys' % vtype
            if len(refkeys - newkeys) > 0 or len(newkeys - refkeys) > 0:
                if len(refkeys - newkeys) > 0:
                    print utils.color('red', '      %d only in ref version' % len(refkeys - newkeys))
                if len(newkeys - refkeys) > 0:
                    print utils.color('red', '      %d only in new version' % len(newkeys - refkeys))
                print '      %d in common' % len(refkeys & newkeys)
            else:
                print '        %d identical keys in new and ref cache' % len(refkeys)

        ptest = 'partition-' + input_stype + '-simu'
        if args.quick and ptest not in self.quick_tests:
            return

        # ----------------------------------------------------------------------------------------
        print '  %s input partition cache file' % input_stype
        def readcache(fname):
            cache = {'naive_seqs' : {}, 'logprobs' : {}}
            with open(fname) as cachefile:
                reader = csv.DictReader(cachefile)
                for line in reader:
                    if line['naive_seq'] != '':
                        cache['naive_seqs'][line['unique_ids']] = line['naive_seq']
                    if line['logprob'] != '':
                        cache['logprobs'][line['unique_ids']] = float(line['logprob'])
            return cache

        refcache = readcache(self.dirs['ref'] + '/' + self.cachefnames[input_stype])
        newcache = readcache(self.dirs['new'] + '/' + self.cachefnames[input_stype])

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
        print '                  fraction different     mean abs difference among differents'
        print '      naive seqs     %s                      %s      (hamming fraction)' % (diff_hfracs_str, mean_hfrac_str)
        print '      log probs      %s                      %s' % (diff_logprob_str, mean_logprob_str)
        if n_different_length > 0:
            print utils.color('red', '        %d different length' % n_different_length)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--dont-run', action='store_true', help='don\'t actually run anything, just check the results')
parser.add_argument('--dryrun', action='store_true', help='do all preparations to run, but don\'t actually run the commands, and don\'t check results')
parser.add_argument('--quick', action='store_true')
parser.add_argument('--bust-cache', action='store_true', help='copy info from new dir to reference dir, i.e. overwrite old test info')
parser.add_argument('--comparison-plots', action='store_true')
parser.add_argument('--print-width', type=int, default=300)
# parser.add_argument('--make-plots', action='store_true')  # needs updating
# example to make comparison plots:
#   ./bin/compare-plotdirs.py --plotdirs test/reference-results/simu-new-performance/sw:test/new-results/simu-new-performance/sw --names ref:new --outdir $www/partis/tmp/test-plots

parser.add_argument('--glfo-dir', default='data/germlines/human')
parser.add_argument('--locus', default='igh')
args = parser.parse_args()

tester = Tester()
if args.bust_cache:
    tester.test(args)
    tester.bust_cache()
else:
    tester.test(args)

def get_typical_variances():
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

