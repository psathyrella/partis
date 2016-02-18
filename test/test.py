#!/usr/bin/env python
import argparse
import numpy
import os
import csv
import glob
import math
import shutil
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
        self.partis = './bin/partis.py'
        self.datafname = 'test/mishmash.fa'  # some data from adaptive, chaim, and vollmers
        self.label = 'test'

        self.stypes = ['ref', 'new']  # I don't know what the 's' stands for
        self.dirs = {'ref' : 'test/reference-results', 'new' : 'test/_new-results'}
        self.perfdirs = {st : 'simu-' + st + '-performance' for st in self.stypes}
        if not os.path.exists(self.dirs['new']):
            os.makedirs(self.dirs['new'])
        simfnames = {st : self.dirs[st] + '/' + self.label + '/simu.csv' for st in self.stypes}
        param_dirs = { st : { dt : self.dirs[st] + '/' + self.label + '/parameters/' + dt + '/hmm' for dt in ['simu', 'data']} for st in self.stypes}  # muddafuggincomprehensiongansta
        run_driver = './bin/run-driver.py --label ' + self.label + ' --stashdir ' + self.dirs['new']
        self.common_extras = ['--seed', '1', '--n-procs', '10', '--only-genes', utils.test_only_genes, '--only-csv-plots']

        # check against reference csv file
        self.tiny_eps = 1e-4
        self.eps_vals = {}  # fractional difference which we allow for each test type (these were generated with the code in get_typical_variances() above)
        self.eps_vals['v_gene_correct'] = 0.02  # hm, actually, I think I just made the annotation ones up
        self.eps_vals['d_gene_correct'] = 0.02
        self.eps_vals['j_gene_correct'] = 0.02
        self.eps_vals['mean_hamming']   = 0.1
        self.eps_vals['precision']      = 0.08
        self.eps_vals['sensitivity']    = 0.08

        n_partition_queries = '250'
        self.logfname = self.dirs['new'] + '/test.log'
        self.cachefnames = { st : 'cache-' + st + '-partition.csv' for st in self.stypes }

        self.quick_tests = ['annotate-ref-simu']
        self.production_tests = ['cache-data-parameters', 'simulate', 'cache-simu-parameters']  # vs "inference" tests. Kind of crappy names, but it's to distinguish these three from all the other ones

        self.tests = OrderedDict()

        def add_inference_tests(input_stype):  # if input_stype is 'ref', infer on old simulation and parameters, if it's 'new' use the new ones
            self.tests['annotate-' + input_stype + '-simu']          = {'bin' : self.partis, 'action' : 'run-viterbi', 'extras' : ['--seqfile', simfnames[input_stype], '--parameter-dir', param_dirs[input_stype]['simu'], '--plotdir', self.dirs['new'] + '/' + self.perfdirs[input_stype], '--plot-performance']}
            self.tests['partition-' + input_stype + '-simu']         = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--seqfile', simfnames[input_stype], '--parameter-dir', param_dirs[input_stype]['simu'], '--n-max-queries', n_partition_queries, '--persistent-cachefname', self.dirs['new'] + '/' + self.cachefnames[input_stype], '--n-precache-procs', '10']}
            # self.tests['partition-' + input_stype + '-data']         = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--seqfile', self.datafname, '--parameter-dir', param_dirs[input_stype]['data'], '--is-data', '--skip-unproductive', '--n-max-queries', n_partition_queries, '--n-precache-procs', '10']}
            self.tests['point-partition-' + input_stype + '-simu']   = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--naive-hamming', '--seqfile', simfnames[input_stype], '--parameter-dir', param_dirs[input_stype]['simu'], '--n-max-queries', n_partition_queries, '--n-precache-procs', '10']}
            self.tests['vsearch-partition-' + input_stype + '-simu'] = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--naive-vsearch', '--seqfile', simfnames[input_stype], '--parameter-dir', param_dirs[input_stype]['simu'], '--n-max-queries', n_partition_queries, '--n-precache-procs', '10']}

        add_inference_tests('ref')
        self.tests['cache-data-parameters']  = {'bin' : run_driver, 'extras' : []}  # ['--skip-unproductive']}
        self.tests['simulate']  = {'bin' : run_driver, 'extras' : ['--n-sim-events', 500, '--n-leaves', 2, '--mimic-data-read-length']}
        self.tests['cache-simu-parameters']  = {'bin' : run_driver, 'extras' : []}
        add_inference_tests('new')

        self.perf_info = { version_stype : OrderedDict() for version_stype in self.stypes }

    # ----------------------------------------------------------------------------------------
    def test(self, args):
        if args.make_plots:
            self.make_comparison_plots()
            sys.exit(0)
        if not args.dont_run:
            self.run(args)
        for version_stype in self.stypes:
            self.read_performance_info(version_stype)
        self.compare_performance(input_stype='ref')
        self.compare_partition_cachefiles(input_stype='ref')
        self.compare_production_results()
        self.compare_performance(input_stype='new')
        self.compare_partition_cachefiles(input_stype='new')

    # ----------------------------------------------------------------------------------------
    def run(self, args):
        open(self.logfname, 'w').close()

        for name, info in self.tests.items():
            if args.quick and name not in self.quick_tests:
                continue

            input_stype = None
            if name not in self.production_tests:
                input_stype = 'ref' if '-ref-' in name else 'new'
                assert '-' + input_stype + '-' in name
                if name == 'partition-' + input_stype + '-simu' and os.path.exists(self.dirs['new'] + '/' + self.cachefnames[input_stype]):
                    check_call(['rm', '-v', self.dirs['new'] + '/' + self.cachefnames[input_stype]])

            action = info['action'] if 'action' in info else name
            cmd_str = info['bin'] + ' --action ' + action
            if info['bin'] == self.partis:
                cmd_str += ' ' + ' '.join(info['extras'] + self.common_extras)
                cmd_str += ' --outfname ' + self.dirs['new'] + '/' + name + '.csv'
            else:
                cmd_str += get_extra_str(info['extras'] + self.common_extras)
                if action == 'cache-data-parameters':
                    cmd_str += ' --datafname ' + self.datafname
            logstr = 'TEST %30s   %s' % (name, cmd_str)
            print logstr
            logfile = open(self.logfname, 'a')
            logfile.write(logstr + '\n')
            logfile.close()
            check_call(cmd_str + ' 1>>' + self.logfname + ' 2>>' + self.logfname, shell=True)

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
    # collect summary performance info from a few places in stashdir
    def read_performance_info(self, version_stype):
        for input_stype in self.stypes:
            self.read_annotation_performance(version_stype, input_stype)
            self.read_partition_performance(version_stype, input_stype)

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
        ptest = 'partition-' + input_stype + '-simu'
        if args.quick and ptest not in self.quick_tests:
            return
        if debug:
            print '  version %s input %s partitioning' % (version_stype, input_stype)
            print '  precision      sensitivity        test                    description'
        for ptest in [k for k in self.tests.keys() if 'partition' in k and input_stype in k]:
            if args.quick and ptest not in self.quick_tests:
                continue
            cp = ClusterPath(-1)
            cp.readfile(self.dirs[version_stype] + '/' + ptest + '.csv')
            if 'data' in ptest:
                raise Exception('needs fixing')
                # ref_cp = ClusterPath(-1)
                # ref_cp.readfile(self.dirs['xxxref'] + '/' + ptest + '.csv')
                # self.perf_info['xxx'][ptest] = utils.adjusted_mutual_information(cp.partitions[cp.i_best], ref_cp.partitions[ref_cp.i_best])  # adj mi between the reference and the new data partitions
                # if debug:
                #     print '    %5.2f   %-28s   to reference partition' % (self.perf_info['xxx'][ptest], ptest)
            else:
                self.perf_info[version_stype][ptest + '-precision'], self.perf_info[version_stype][ptest + '-sensitivity'] = cp.ccfs[cp.i_best]
                if debug:
                    print '    %5.2f          %5.2f      %-28s   to true partition' % (self.perf_info[version_stype][ptest + '-precision'], self.perf_info[version_stype][ptest + '-sensitivity'], ptest)

    # ----------------------------------------------------------------------------------------
    def compare_performance(self, input_stype):
        print 'performance with %s simulation and parameters' % input_stype

        # make sure there's a new performance value for each reference one, and vice versa
        refkeys = set(self.perf_info['ref'].keys())
        newkeys = set(self.perf_info['new'].keys())
        if len(refkeys - newkeys) > 0 or len(newkeys - refkeys) > 0:
            print '  %d keys only in ref' % len(refkeys - newkeys)
            print '  %d keys only in new' % len(newkeys - refkeys)
            print '  %d in common' % len(refkeys & newkeys)
            raise Exception('')

        for name in self.perf_info['ref']:  # don't use the sets above so we get the nice ordering
            if input_stype not in name:
                continue
            ref_val = self.perf_info['ref'][name]
            new_val = self.perf_info['new'][name]
            val_type = name.split('-')[-1]
            print '  %-28s %-15s       %-5.3f' % (name.replace('-' + val_type, ''), val_type, ref_val),
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
        for fname in ['test/parameters/data', 'test/simu.csv', 'test/parameters/simu/hmm-true', 'test/parameters/simu/sw', 'test/parameters/simu/hmm']:
            print '    %s' % fname
            cmd = 'diff -qbr ' + ' '.join(self.dirs[st] + '/' + fname for st in self.stypes)
            proc = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                outlines = [ l for l in out.split('\n') if 'differ' in l ]
                n_total_files = int(check_output('find ' + self.dirs['ref'] + '/' + fname + ' -type f | wc -l', shell=True))
                print utils.color('red', '      %d / %d files differ' % (len(outlines), n_total_files)),
                print '  (%s)' % cmd
                if err != '':
                    print err

    # ----------------------------------------------------------------------------------------
    def make_comparison_plots(self):
        plotdirs = [
            self.perfdirs['ref'] + '/sw',  # ref sw performance
            self.perfdirs['ref'] + '/hmm', # ref hmm performance
            # 'test/plots/data/sw',          # sw data parameters
            # 'test/plots/data/hmm',         # hmm data parameters
            # 'test/plots/data/hmm/mute-freqs/v',
            # 'test/plots/data/sw/mute-freqs',
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
        ptest = 'partition-' + input_stype + '-simu'
        if args.quick and ptest not in self.quick_tests:
            return

        print '%s input partition cache file' % input_stype

        def readcache(fname):
            cache = {}
            with open(fname) as cachefile:
                reader = csv.DictReader(cachefile)
                for line in reader:
                    cache[line['unique_ids']] = {'naive_seq' : line['naive_seq'], 'logprob' : float(line['logprob'])}
            return cache

        refcache = readcache(self.dirs['ref'] + '/' + self.cachefnames[input_stype])
        newcache = readcache(self.dirs['new'] + '/' + self.cachefnames[input_stype])

        # work out intersection and complement
        refkeys = set(refcache.keys())
        newkeys = set(newcache.keys())
        if len(refkeys - newkeys) > 0 or len(newkeys - refkeys) > 0:
            if len(refkeys - newkeys) > 0:
                print utils.color('red', '  %d only in ref version' % len(refkeys - newkeys))
            if len(newkeys - refkeys) > 0:
                print utils.color('red', '  %d only in new version' % len(newkeys - refkeys))
            print '  %d in common' % len(refkeys & newkeys)
        else:
            print '    %d identical keys in new and ref cache' % len(refkeys)

        hammings, delta_logprobs = [], []
        n_hammings, n_delta_logprobs = 0, 0
        n_different_length, n_big_hammings, n_big_delta_logprobs = 0, 0, 0
        hamming_eps = 0.
        logprob_eps = 1e-5
        for uids in refkeys & newkeys:
            refline = refcache[uids]
            newline = newcache[uids]
            if refline['naive_seq'] != '':
                n_hammings += 1
                if len(refline['naive_seq']) == len(newline['naive_seq']):
                    hamming_fraction = utils.hamming_fraction(refline['naive_seq'], newline['naive_seq'])
                    if hamming_fraction > hamming_eps:
                        n_big_hammings += 1
                        hammings.append(hamming_fraction)
                else:
                    n_different_length += 1
            if refline['logprob'] != '':
                n_delta_logprobs += 1
                delta_logprob = abs(float(refline['logprob']) - float(newline['logprob']))
                if delta_logprob > logprob_eps:
                    n_big_delta_logprobs += 1
                    delta_logprobs.append(delta_logprob)

        diff_hfracs_str = '%3d / %4d' % (n_big_hammings, n_hammings)
        mean_hfrac_str = '%.3f' % (numpy.average(hammings) if len(hammings) > 0 else 0.)
        if n_big_hammings > 0:
            diff_hfracs_str = utils.color('red', diff_hfracs_str)
            mean_hfrac_str = utils.color('red', mean_hfrac_str)

        diff_logprob_str = '%3d / %4d' % (n_big_delta_logprobs, n_delta_logprobs)
        mean_logprob_str = '%.6f' % (numpy.average(delta_logprobs) if len(delta_logprobs) > 0 else 0.)
        if n_big_delta_logprobs > 0:
            diff_logprob_str = utils.color('red', diff_logprob_str)
            mean_logprob_str = utils.color('red', mean_logprob_str)
        print '                fraction different     mean difference among differents'
        print '    naive seqs     %s                      %s      (hamming fraction)' % (diff_hfracs_str, mean_hfrac_str)
        print '    log probs      %s                      %s' % (diff_logprob_str, mean_logprob_str)
        if n_different_length > 0:
            print utils.color('red', '      %d different length' % n_different_length)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--dont-run', action='store_true', help='don\'t actually run the tests, presumably so you can just check the results ')
parser.add_argument('--dont-plot', action='store_true', help='don\'t make all the comparison plots')
parser.add_argument('--quick', action='store_true')
parser.add_argument('--bust-cache', action='store_true', help='copy info from new dir to reference dir, i.e. overwrite old test info')
parser.add_argument('--make-plots', action='store_true')
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

