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
from subprocess import Popen, PIPE, check_call, check_output
import sys
sys.path.insert(1, './python')
from baseutils import get_extra_str
import utils
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
        self.dirs = {'ref' : 'test/reference-results', 'new' : 'test/_new-results'}
        self.perfdirs = {st : 'simu-' + st + '-performance' for st in self.stypes}
        if not os.path.exists(self.dirs['new']):
            os.makedirs(self.dirs['new'])
        self.simfnames = {st : self.dirs[st] + '/' + self.label + '/simu.csv' for st in self.stypes}
        self.param_dirs = { st : { dt : self.dirs[st] + '/' + self.label + '/parameters/' + dt for dt in ['simu', 'data']} for st in self.stypes}  # muddafuggincomprehensiongansta
        self.common_extras = ['--seed', '1', '--n-procs', '10', '--only-genes', 'TEST', '--only-csv-plots']

        self.perf_info = { version_stype : OrderedDict() for version_stype in self.stypes }

        # check against reference csv file
        self.tiny_eps = 1e-4
        self.run_times = {}
        self.eps_vals = {}  # fractional difference which we allow for each test type (these were generated with the code in get_typical_variances() above)
        self.eps_vals['v_gene_correct'] = 0.02  # hm, actually, I think I just made the annotation ones up
        self.eps_vals['d_gene_correct'] = 0.02
        self.eps_vals['j_gene_correct'] = 0.02
        self.eps_vals['mean_hamming']   = 0.1
        self.eps_vals['precision']      = 0.08
        self.eps_vals['sensitivity']    = 0.08

        self.n_partition_queries = '250'
        n_data_inference_queries = '50'
        self.logfname = self.dirs['new'] + '/test.log'
        self.cachefnames = { st : 'cache-' + st + '-partition.csv' for st in self.stypes }

        self.quick_tests = ['annotate-ref-simu']
        self.production_tests = ['cache-parameters-data', 'simulate', 'cache-parameters-simu']  # vs "inference" tests. Kind of crappy names, but we need to distinguish these three from all the other ones

        self.tests = OrderedDict()

        def add_inference_tests(input_stype):  # if input_stype is 'ref', infer on old simulation and parameters, if it's 'new' use the new ones
            self.tests['annotate-' + input_stype + '-simu']          = {'extras' : ['--plotdir', self.dirs['new'] + '/' + self.perfdirs[input_stype], '--plot-performance']}
            self.tests['annotate-' + input_stype + '-data']          = {'extras' : ['--n-max-queries', n_data_inference_queries]}
            self.tests['partition-' + input_stype + '-simu']         = {'extras' : ['--n-max-queries', self.n_partition_queries, '--persistent-cachefname', self.dirs['new'] + '/' + self.cachefnames[input_stype], '--n-precache-procs', '10', '--biggest-logprob-cluster-to-calculate', '2', '--biggest-naive-seq-cluster-to-calculate', '2']}
            self.tests['seed-partition-' + input_stype + '-simu']    = {'extras' : ['--n-max-queries', '-1', '--n-precache-procs', '10']}
            self.tests['vsearch-partition-' + input_stype + '-simu'] = {'extras' : ['--naive-vsearch', '--n-max-queries', self.n_partition_queries, '--n-precache-procs', '10']}

        if not args.skip_ref:
            add_inference_tests('ref')
        if not args.only_ref:
            self.tests['cache-parameters-data']  = {'extras' : []}
            self.tests['simulate']  = {'extras' : ['--n-sim-events', '500', '--n-trees', '500', '--n-leaves', '2', '--mimic-data-read-length']}
            self.tests['cache-parameters-simu']  = {'extras' : []}
            add_inference_tests('new')

        # add some arguments
        for ptest, argfo in self.tests.items():
            namelist = ptest.split('-')
            argfo['bin'] = self.partis
            if 'annotate' in ptest:
                argfo['action'] = 'run-viterbi'
            elif 'partition' in ptest:
                argfo['action'] = 'partition'
            elif 'cache-parameters-' in ptest:
                argfo['action'] = 'cache-parameters'
                dtype = namelist[-1]
                assert dtype in self.dtypes
                argfo['extras'] += ['--plotdir', self.dirs['new'] + '/' + self.label + '/plots/' + dtype]
            else:
                argfo['action'] = ptest

            if ptest in self.production_tests:
                input_stype = 'new'  # sort of...
            else:
                input_stype = namelist[-2]
                assert input_stype in self.stypes
            argfo['input_stype'] = input_stype
            if ptest == 'simulate':
                argfo['extras'] += ['--parameter-dir', self.param_dirs[input_stype]['data']]
            elif namelist[-1] == 'simu':
                argfo['extras'] += ['--is-simu', ]
                argfo['extras'] += ['--infname', self.simfnames[input_stype]]
                argfo['extras'] += ['--parameter-dir', self.param_dirs[input_stype]['simu']]
            elif namelist[-1] == 'data':
                argfo['extras'] += ['--infname', self.datafname]
                argfo['extras'] += ['--parameter-dir', self.param_dirs[input_stype]['data']]
            else:
                raise Exception('-'.join(namelist))

    # ----------------------------------------------------------------------------------------
    def test(self, args):
        if args.make_plots:
            self.make_comparison_plots()
            return
        if not args.dont_run:
            self.run(args)
        if not args.skip_ref:
            self.compare_stuff(input_stype='ref')
        if not args.only_ref and not args.quick:
            self.compare_production_results()
            self.compare_stuff(input_stype='new')

        self.compare_run_times()

    # ----------------------------------------------------------------------------------------
    def compare_stuff(self, input_stype):
        print '%s input' % input_stype
        for version_stype in self.stypes:
            self.read_annotation_performance(version_stype, input_stype)
            self.read_partition_performance(version_stype, input_stype)
        self.compare_performance(input_stype)
        self.compare_partition_cachefiles(input_stype)
        self.compare_data_annotation(input_stype)

    # ----------------------------------------------------------------------------------------
    def prepare_to_run(self, args, name, info):
        """ Pre-run stuff that you don't want to do until *right* before you actually run. """

        # delete old partition cache file
        if name == 'partition-' + info['input_stype'] + '-simu':
            this_cachefname = self.dirs['new'] + '/' + self.cachefnames[info['input_stype']]
            if os.path.exists(this_cachefname):
                check_call(['rm', '-v', this_cachefname])

        ref_globfnames = [fn for dtype in self.dtypes for fn in glob.glob(self.param_dirs['ref'][dtype] + '/sw-cache-*')]
        if len(ref_globfnames) > 0:
            raise Exception('found reference sw cache files %s -- but you really want ref sw to run from scratch' % ' '.join(ref_globfnames))

        # delete any old sw cache files
        for dtype in self.dtypes:
            globfnames = glob.glob(self.param_dirs['new'][dtype] + '/sw-cache-*.csv')
            if len(globfnames) == 0:  # not there
                continue
            elif len(globfnames) != 1:
                raise Exception('unexpected sw cache files: %s' % ' '.join(globfnames))
            check_call(['rm', '-v', globfnames[0]])
            sw_cache_gldir = globfnames[0].replace('.csv', '-glfo')
            glutils.remove_glfo_files(sw_cache_gldir, args.chain)
            os.rmdir(sw_cache_gldir)

        # choose a seed uid
        if name == 'seed-partition-' + info['input_stype'] + '-simu':
            seed_uid, _ = utils.choose_seed_unique_id(args.glfo_dir, args.chain, self.simfnames[info['input_stype']], 5, 8, n_max_queries=int(self.n_partition_queries), debug=False)
            info['extras'] += ['--seed-unique-id', seed_uid]

    # ----------------------------------------------------------------------------------------
    def run(self, args):
        open(self.logfname, 'w').close()

        for name, info in self.tests.items():
            if args.quick and name not in self.quick_tests:
                continue

            self.prepare_to_run(args, name, info)

            action = info['action']
            cmd_str = info['bin'] + ' ' + action
            cmd_str += ' ' + ' '.join(info['extras'] + self.common_extras)
            if name == 'simulate':
                cmd_str += ' --outfname ' + self.simfnames['new']
            elif 'cache-parameters-' not in name:
                cmd_str += ' --outfname ' + self.dirs['new'] + '/' + name + '.csv'

            logstr = 'TEST %30s   %s' % (name, cmd_str)
            print logstr
            logfile = open(self.logfname, 'a')
            logfile.write(logstr + '\n')
            logfile.close()
            start = time.time()
            check_call(cmd_str + ' 1>>' + self.logfname + ' 2>>' + self.logfname, shell=True)
            self.run_times[name] = time.time() - start  # seconds

        self.write_run_times()

    # ----------------------------------------------------------------------------------------
    def remove_reference_results(self, expected_content):
        print '  remove ref files'
        dir_content = set([os.path.basename(f) for f in glob.glob(self.dirs['ref'] + '/*')])
        if len(dir_content - expected_content) > 0 or len(expected_content - dir_content) > 0:
            if len(dir_content - expected_content) > 0:
                print 'in ref dir but not expected\n    %s' % (utils.color('red', ' '.join(dir_content - expected_content)))
            if len(expected_content - dir_content) > 0:
                print 'expected but not in ref dir\n    %s' % (utils.color('red', ' '.join(expected_content - dir_content)))
            raise Exception('unexpected or missing content in reference dir')
        for fname in [self.dirs['ref'] + '/' + ec for ec in expected_content]:
            print '    rm %s' % fname
            if os.path.isdir(fname):
                shutil.rmtree(fname)
            else:
                os.remove(fname)

    # ----------------------------------------------------------------------------------------
    def bust_cache(self):
        test_outputs = [k + '.csv' for k in self.tests.keys() if k not in self.production_tests]
        expected_content = set(test_outputs + self.perfdirs.values() + self.cachefnames.values() + [os.path.basename(self.logfname), self.label])
        expected_content.add('run-times.csv')

        # remove (very, very gingerly) whole reference dir
        self.remove_reference_results(expected_content)

        # copy over parameters, simulation, and plots
        # NOTE in the ref dir, the ref in new input_stype files are always identical, since we'd always make sure we actually pass the tests before committing
        print '  copy new files to ref'
        for fname in expected_content:
            source_fname = self.dirs['new'] + '/' + fname
            if '-ref-' in fname:  # these correspond to stuff run on the old reference simulation and parameters, so we no longer need it
                if os.path.isdir(source_fname):
                    shutil.rmtree(source_fname)
                else:
                    os.remove(source_fname)
                continue

            if '-new-' in fname:  # whereas the stuff that's on the new simulation and parameters replaces both the ref and new stuff
                print '    cp %s   -->  %s/ (new --> ref)' % (fname, self.dirs['ref'])
                if os.path.isdir(source_fname):
                    shutil.copytree(source_fname, self.dirs['ref'] + '/' + fname.replace('-new-', '-ref-'))
                else:
                    shutil.copy(source_fname, self.dirs['ref'] + '/' + fname.replace('-new-', '-ref-'))
            print '    mv %s   -->  %s/' % (fname, self.dirs['ref'])
            shutil.move(source_fname, self.dirs['ref'] + '/')

    # ----------------------------------------------------------------------------------------
    def read_annotation_performance(self, version_stype, input_stype, debug=False):
        """ version_stype is the code version, while input_stype is the input data version, i.e. 'ref', 'new' is the reference code version (last commit) run on the then-new simulation and parameters"""
        ptest = 'annotate-' + input_stype + '-simu'
        if args.quick and ptest not in self.quick_tests:
            return
        if debug:
            print '  version %s input %s annotation' % (version_stype, input_stype)

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

        perfdir = self.dirs[version_stype] + '/' + self.perfdirs[input_stype]
        for method in ['sw', 'hmm']:
            if debug:
                print '   ', method

            # fraction of genes correct
            for region in utils.regions:
                fraction_correct = read_performance_file(perfdir + '/' + method + '/' + region + '_gene.csv', 'contents', only_ibin=1)
                if debug:
                    print '      %s %.3f' % (region, fraction_correct)
                self.perf_info[version_stype][input_stype + '-' + method + '-' + region + '_gene_correct'] = fraction_correct

            # hamming fraction
            hamming_hist = Hist(fname=perfdir + '/' + method + '/hamming_to_true_naive.csv')
            if debug:
                print '      mean hamming %.2f' % hamming_hist.get_mean()
            self.perf_info[version_stype][input_stype + '-' + method + '-mean_hamming'] = hamming_hist.get_mean()

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
        if debug:
            print '  version %s input %s partitioning' % (version_stype, input_stype)
            print '  precision      sensitivity        test                    description'
        for ptest in ptest_list:
            cp = ClusterPath(-1)
            cp.readfile(self.dirs[version_stype] + '/' + ptest + '.csv')
            ccfs = cp.ccfs[cp.i_best]
            if None in ccfs:
                raise Exception('none type ccf read from %s' % self.dirs[version_stype] + '/' + ptest + '.csv')
            self.perf_info[version_stype][ptest + '-precision'], self.perf_info[version_stype][ptest + '-sensitivity'] = ccfs
            if debug:
                print '    %5.2f          %5.2f      %-28s   to true partition' % (self.perf_info[version_stype][ptest + '-precision'], self.perf_info[version_stype][ptest + '-sensitivity'], ptest)

    # ----------------------------------------------------------------------------------------
    def compare_data_annotation(self, input_stype):
        # NOTE don't really need to do this for simulation, since for simulation we already compare the performance info
        ptest = 'annotate-' + input_stype + '-data'
        if args.quick and ptest not in self.quick_tests:
            return
        print '  %s data annotation' % input_stype
        infnames = [self.dirs[version_stype] + '/' + ptest + '.csv' for version_stype in self.stypes]
        cmd = 'diff -u ' + ' '.join(infnames) + ' | grep "^+[^+]" | wc -l'
        n_diff_lines = int(check_output(cmd, shell=True))
        if n_diff_lines == 0:
            print '      ok'
        else:
            n_total_lines = int(check_output(['wc', '-l', infnames[0]]).split()[0])
            print utils.color('red', '      %d / %d lines differ' % (n_diff_lines, n_total_lines)),
            print '   (%s)' % cmd

    # ----------------------------------------------------------------------------------------
    def compare_performance(self, input_stype):
        performance_metric_list = [n for n in self.perf_info['ref'] if input_stype in n]
        if len(performance_metric_list) == 0:
            return

        print '  performance with %s simulation and parameters' % input_stype

        # make sure there's a new performance value for each reference one, and vice versa
        refkeys = set(self.perf_info['ref'].keys())
        newkeys = set(self.perf_info['new'].keys())
        if len(refkeys - newkeys) > 0 or len(newkeys - refkeys) > 0:
            print '  %d keys only in ref' % len(refkeys - newkeys)
            print '  %d keys only in new' % len(newkeys - refkeys)
            print '  %d in common' % len(refkeys & newkeys)
            raise Exception('')

        for name in performance_metric_list:  # don't use the sets above so we get the nice ordering
            ref_val = self.perf_info['ref'][name]
            new_val = self.perf_info['new'][name]
            val_type = name.split('-')[-1]
            print '    %-28s %-15s       %-5.3f' % (name.replace('-' + val_type, ''), val_type, ref_val),
            fractional_change = (new_val - ref_val) / ref_val  # NOTE not the abs value yet
            if abs(fractional_change) > self.eps_vals[val_type]:
                print '--> %-5.3f %s' % (new_val, utils.color('red', '(%+.3f)' % fractional_change)),
            elif abs(fractional_change) > self.tiny_eps:
                print '--> %-5.3f %s' % (new_val, utils.color('yellow', '(%+.3f)' % fractional_change)),
            else:
                print '    ok   ',
            print ''

    # ----------------------------------------------------------------------------------------
    def compare_production_results(self):
        if args.quick:
            return
        print 'diffing production results'
        for fname in ['test/parameters/data', 'test/simu.csv', 'test/parameters/simu']:
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

        for name in times['ref']:
            if args.quick and name not in self.quick_tests:
                continue
            if args.only_ref and '-ref-' not in name:
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
        # check_call(['chmod', '664', htmlfname])
        check_call(['./bin/permissify-www', www_dir])
        
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
parser.add_argument('--dont-run', action='store_true', help='don\'t actually run the tests, presumably so you can just check the results ')
parser.add_argument('--quick', action='store_true')
parser.add_argument('--only-ref', action='store_true', help='only run with input_stype of \'ref\'')
parser.add_argument('--skip-ref', action='store_true', help='skip stuff that\'s run by --only-ref')
parser.add_argument('--bust-cache', action='store_true', help='copy info from new dir to reference dir, i.e. overwrite old test info')
parser.add_argument('--make-plots', action='store_true')
parser.add_argument('--glfo-dir', default='data/germlines/human')
parser.add_argument('--chain', default='h')
args = parser.parse_args()

tester = Tester()
if args.bust_cache:
    tester.bust_cache()
else:
    tester.test(args)

def get_typical_variances():
    raise Exception('needs updating to work as a function')
    # cp = ClusterPath()
    # cp.readfile('tmp.csv')
    # cp.print_partitions()
    # sys.exit()
    # cps = []
    # adj_mis, ccf_unders, ccf_overs = [], [], []
    # for iseed in range(6):
    #     # print 'seed %d' % iseed
    #     cp = ClusterPath()
    #     cp.readfile('%d.csv' % iseed)
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

