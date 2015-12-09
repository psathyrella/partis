#!/usr/bin/env python
import argparse
import numpy
import os
import csv
import glob
import math
from collections import OrderedDict
from subprocess import check_call, check_output
import sys
sys.path.insert(1, './python')
from baseutils import get_extra_str
import utils
from hist import Hist
from clusterpath import ClusterPath

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

# ----------------------------------------------------------------------------------------
class Tester(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self):
        self.partis = './bin/partis.py'
        self.datafname = 'test/mishmash.csv'  # some data from adaptive, chaim, and vollmers
        self.label = 'test'
        self.common_extras = ['--seed', '1', '--n-procs', '10', '--only-genes', utils.test_only_genes]

        stypes = ['ref', 'new']  # I don't know what the 's' stands for
        self.dirs = {'ref' : 'test/reference-results', 'new' : 'test/_new-results'}
        if not os.path.exists(self.dirs['new']):
            os.makedirs(self.dirs['new'])
        simfnames = {st : self.dirs[st] + '/' + self.label + '/simu.csv' for st in stypes}
        param_dirs = { st : { dt : self.dirs[st] + '/' + self.label + '/parameters/' + dt + '/hmm' for dt in ['simu', 'data']} for st in stypes}  # muddafuggincomprehensiongansta
        run_driver = './bin/run-driver.py --label ' + self.label + ' --stashdir ' + self.dirs['new']

        n_partition_queries = '250'
        self.cachefname = 'hmm_cached_info.csv'
        self.logfname = self.dirs['new'] + '/test.log'
        open(self.logfname, 'w').close()

        self.quick_tests = ['annotate-ref-simu']

        self.tests = OrderedDict()

        def add_inference_tests(stype):  # if stype is 'ref', infer on old simulation and parameters, if it's 'new' use the new ones
            self.tests['annotate-' + stype + '-simu']          = {'bin' : self.partis, 'action' : 'run-viterbi', 'extras' : ['--seqfile', simfnames[stype], '--parameter-dir', param_dirs[stype]['simu'], '--plotdir', self.dirs['new'] + '/' + self.label + '/plots/' + stype + '-simu-performance', '--plot-performance']}
            self.tests['partition-' + stype + '-simu']         = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--seqfile', simfnames[stype], '--parameter-dir', param_dirs[stype]['simu'], '--n-max-queries', n_partition_queries, '--persistent-cachefname', self.dirs['new'] + '/' + self.cachefname]}
            # self.tests['partition-' + stype + '-data']         = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--seqfile', self.datafname, '--parameter-dir', param_dirs[stype]['data'], '--is-data', '--skip-unproductive', '--n-max-queries', n_partition_queries]}
            self.tests['point-partition-' + stype + '-simu']   = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--naive-hamming', '--seqfile', simfnames[stype], '--parameter-dir', param_dirs[stype]['simu'], '--n-max-queries', n_partition_queries]}
            self.tests['vsearch-partition-' + stype + '-simu'] = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--naive-vsearch', '--seqfile', simfnames[stype], '--parameter-dir', param_dirs[stype]['simu'], '--n-max-queries', n_partition_queries]}

        add_inference_tests('ref')

        # # then infer new parameters, and make new simulation
        # self.tests['cache-data-parameters']  = {'bin' : run_driver, 'extras' : ['--skip-unproductive']}
        # self.tests['simulate']  = {'bin' : run_driver, 'extras' : ['--n-sim-events', 500, '--n-leaves', 2, '--mimic-data-read-length']}
        # self.tests['cache-simu-parameters']  = {'bin' : run_driver, 'extras' : []}

        # # then test performance on the new simulation, with the new parameter values
        # self.tests['annotate-new-simu']          = {'bin' : self.partis, 'action' : 'run-viterbi', 'extras' : ['--seqfile', new_simfname, '--parameter-dir', param_dirs['new']['simu'], '--plotdir', self.dirs['new'] + '/' + self.label + '/plots/new-simu-performance', '--plot-performance']}
        # # self.tests['annotate-new-simu']  = {'bin' : run_driver, 'action' : 'plot-performance', 'extras' : []}
        # # self.tests['single-point-estimate']      = {'bin' : self.partis, 'action' : 'run-viterbi', 'extras' : ['--seqfile', new_simfname, '--parameter-dir', param_dirs['new']['simu'], '--n-max-queries', '10']}
        # self.tests['partition-data']             = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--seqfile', self.datafname, '--parameter-dir', param_dirs['new']['data'], '--is-data', '--skip-unproductive', '--n-max-queries', n_partition_queries]}
        # self.tests['partition-new-simu']         = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--seqfile', new_simfname, '--parameter-dir', param_dirs['new']['simu'], '--n-max-queries', n_partition_queries]}
        # self.tests['point-partition-new-simu']   = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--seqfile', new_simfname, '--parameter-dir', param_dirs['new']['simu'], '--naive-hamming', '--n-max-queries', n_partition_queries]}
        # self.tests['vsearch-partition-new-simu'] = {'bin' : self.partis, 'action' : 'partition',   'extras' : ['--seqfile', new_simfname, '--parameter-dir', param_dirs['new']['simu'], '--naive-vsearch', '--n-max-queries', n_partition_queries]}

        self.perf_info = OrderedDict()

    # ----------------------------------------------------------------------------------------
    def test(self, args):
        if not args.dont_run:
            self.run(args)
        self.read_new_performance_info()
        self.compare_to_reference_performance_file()
        self.compare_partition_cachefiles()
        self.write_new_performance_file()

    # ----------------------------------------------------------------------------------------
    def run(self, args):
        for name, info in self.tests.items():
            if args.quick and name not in self.quick_tests:
                continue
            action = info['action'] if 'action' in info else name
            cmd_str = info['bin'] + ' --action ' + action
            if info['bin'] == self.partis:
                cmd_str += ' ' + ' '.join(info['extras'] + self.common_extras)
                cmd_str += ' --outfname ' + self.dirs['new'] + '/' + name + '.csv'
            else:
                cmd_str += get_extra_str(info['extras'] + self.common_extras)
                if action == 'cache-data-parameters':
                    cmd_str += ' --datafname ' + self.datafname
            if self.cachefname in cmd_str and os.path.exists(self.dirs['new'] + '/' + self.cachefname):
                check_call(['rm', '-v', self.dirs['new'] + '/' + self.cachefname])
            logstr = 'TEST %30s   %s' % (name, cmd_str)
            print logstr
            logfile = open(self.logfname, 'a')
            logfile.write(logstr + '\n')
            logfile.close()
            check_call(cmd_str + ' 1>>' + self.logfname + ' 2>>' + self.logfname, shell=True)

    # ----------------------------------------------------------------------------------------
    # collect summary performance info from a few places in stashdir
    def read_new_performance_info(self):
        for simtype in ['ref']:  #, 'new']
            self.read_annotation_performance(simtype)
        self.read_partition_performance()

    # ----------------------------------------------------------------------------------------
    def read_annotation_performance(self, simtype):
        assert simtype in ['ref', 'new']
        print simtype, 'annotation'

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

        perfdir = self.dirs['new'] + '/' + self.label + '/plots/' + simtype + '-simu-performance'
        for method in ['sw', 'hmm']:
            print '   ', method

            # fraction of genes correct
            for region in utils.regions:
                fraction_correct = read_performance_file(perfdir + '/' + method + '/plots/' + region + '_gene.csv', 'contents', only_ibin=1)
                print '      %s %.3f' % (region, fraction_correct)
                self.perf_info[simtype + '-' + method + '-' + region + '_gene_correct'] = fraction_correct

            # hamming fraction
            hamming_hist = Hist(fname=perfdir + '/' + method + '/plots/hamming_to_true_naive.csv')
            print '      mean hamming %.2f' % hamming_hist.get_mean()
            self.perf_info[simtype + '-' + method + '-mean_hamming'] = hamming_hist.get_mean()

    # ----------------------------------------------------------------------------------------
    def read_partition_performance(self):
        """ Read new partitions from self.dirs['new'], and put the comparison numbers in self.perf_info (compare either to true, for simulation, or to the partition in reference dir, for data). """
        print 'partitioning'
        print '    adj mi   ccf under/over        test                    description'
        for ptest in [k for k in self.tests.keys() if 'partition' in k]:
            if args.quick and ptest not in self.quick_tests:
                continue
            cp = ClusterPath(-1)
            cp.readfile(self.dirs['new'] + '/' + ptest + '.csv')
            if 'data' in ptest:
                ref_cp = ClusterPath(-1)
                ref_cp.readfile(self.dirs['ref'] + '/' + ptest + '.csv')
                self.perf_info[ptest] = utils.adjusted_mutual_information(cp.partitions[cp.i_best], ref_cp.partitions[ref_cp.i_best])  # adj mi between the reference and the new data partitions
                print '    %5.2f   %-28s   to reference partition' % (self.perf_info[ptest], ptest)
            else:
                self.perf_info[ptest + '-adj_mi'] = cp.adj_mis[cp.i_best]  # adj mi to true partition
                self.perf_info[ptest + '-ccf_under'], self.perf_info[ptest + '-ccf_over'] = cp.ccfs[cp.i_best]
                print '    %5.2f    %5.2f %5.2f      %-28s   to true partition' % (self.perf_info[ptest + '-adj_mi'], self.perf_info[ptest + '-ccf_under'], self.perf_info[ptest + '-ccf_over'], ptest)

    # ----------------------------------------------------------------------------------------
    def compare_to_reference_performance_file(self):
        # NOTE does *not* regenerate the reference performance file based on the reference outputs
        print 'comparing to reference file'

        # check against reference csv file
        tiny_eps = 1e-4
        eps_vals = {}  # fractional difference which we allow for each test type (these were generated with the code in get_typical_variances() above)
        eps_vals['v_gene_correct'] = 0.001  # hm, actually, I think I just made the annotation ones up
        eps_vals['d_gene_correct'] = 0.001
        eps_vals['j_gene_correct'] = 0.001
        eps_vals['mean_hamming']   = 0.001
        eps_vals['adj_mi']         = 0.2  # the three partitioning ones are roughly two sigma
        eps_vals['ccf_under']      = 0.08
        eps_vals['ccf_over']       = 0.08

        with open(self.dirs['ref'] + '/performance-info.csv') as perf_file:
            reader = csv.DictReader(perf_file)
            line = reader.next()
            for name, new_val in self.perf_info.items():
                ref_val = float(line[name])
                val_type = name.split('-')[-1]
                print '  %-28s %-15s       %-5.3f' % (name.replace('-' + val_type, ''), val_type, ref_val),
                fractional_change = (new_val - ref_val) / ref_val  # NOTE not the abs value yet
                if abs(fractional_change) > eps_vals[val_type]:
                    print '--> %-5.3f %s' % (new_val, utils.color('red', '(%+.3f)' % fractional_change)),
                elif abs(fractional_change) > tiny_eps:
                    print '--> %-5.3f %s' % (new_val, utils.color('yellow', '(%+.3f)' % fractional_change)),
                else:
                    print '    ok   ',
                print ''

    # ----------------------------------------------------------------------------------------
    def compare_partition_cachefiles(self):
        print '\npartition cache file'

        def readcache(fname):
            cache = {}
            with open(fname) as cachefile:
                reader = csv.DictReader(cachefile)
                for line in reader:
                    cache[line['unique_ids']] = {'naive_seq' : line['naive_seq'], 'logprob' : float(line['logprob'])}
            return cache

        refcache = readcache(self.dirs['ref'] + '/' + self.cachefname)
        newcache = readcache(self.dirs['new'] + '/' + self.cachefname)

        # work out intersection and complement
        refkeys = set(refcache.keys())
        newkeys = set(newcache.keys())
        if len(refkeys - newkeys) > 0 or len(newkeys - refkeys) > 0:
            print '  %d in ref but not in new cache' % len(refkeys - newkeys)
            print '  %d in new but not in ref cache' % len(newkeys - refkeys)
            print '  %d cached uids in common' % len(refkeys & newkeys)
        else:
            print '  all keys in new and ref cache are the same'

        hammings, delta_logprobs = [], []
        n_different_length, n_big_hammings, n_big_delta_logprobs = 0, 0, 0
        hamming_eps = 0.
        logprob_eps = 1e-5
        for uids in refkeys & newkeys:
            refline = refcache[uids]
            newline = newcache[uids]
            if refline['naive_seq'] != '':
                if len(refline['naive_seq']) == len(newline['naive_seq']):
                    hamming_fraction = utils.hamming_fraction(refline['naive_seq'], newline['naive_seq'])
                    if hamming_fraction > hamming_eps:
                        n_big_hammings += 1
                        hammings.append(hamming_fraction)
                else:
                    n_different_length += 1
            if refline['logprob'] != '':
                delta_logprob = abs(float(refline['logprob']) - float(newline['logprob']))
                if delta_logprob > logprob_eps:
                    n_big_delta_logprobs += 1
                    delta_logprobs.append(delta_logprob)

        print '              fraction different     mean difference among differents'
        print '  naive seqs     %d / %d                      %.3f' % (n_big_hammings, len(hammings) + n_different_length, numpy.average(hammings) if len(hammings) > 0 else 0.)
        print '  log probs      %d / %d                      %.3f' % (n_big_delta_logprobs, len(delta_logprobs), numpy.average(delta_logprobs) if len(delta_logprobs) > 0 else 0.)
        if n_different_length > 0:
            print '    %d different length' % n_different_length

    # ----------------------------------------------------------------------------------------
    def write_new_performance_file(self):
        with open(self.dirs['new'] + '/performance-info.csv', 'w') as perf_file:
            writer = csv.DictWriter(perf_file, self.perf_info.keys())
            writer.writeheader()
            writer.writerow(self.perf_info)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--dont-run', action='store_true', help='don\'t actually run the tests, presumably so you can just check the results ')
parser.add_argument('--dont-plot', action='store_true', help='don\'t make all the comparison plots')
parser.add_argument('--quick', action='store_true')
args = parser.parse_args()
print 'TODO it\'s confusing that the otuput annotation and partition files are in the reference dir, but we only look in the performance csv file'
print 'TODO add in cache-bust option'

tester = Tester()
tester.test(args)

# ----------------------------------------------------------------------------------------
# make a bunch of comparison plots
print 'skipping plots for the moment'
sys.exit()
plotdirs = ['test/plots/data/sw', 'test/plots/data/hmm', 'test/plots/simu/hmm-true', 'test/plots/ref-simu-performance/sw', 'test/plots/ref-simu-performance/hmm']
base_check_cmd = './bin/compare.py --dont-calculate-mean-info --graphify --linewidth 1 --markersizes 2:1 --names reference:new'  # --colors 595:807:834 --scale-errors 1.414
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
    if args.dont_plot:
        continue
    plotdirstr = plotdir.replace(self.dirs['ref'] + '/', '')
    check_cmd = base_check_cmd + ' --plotdirs '  + self.dirs['ref'] + '/' + plotdirstr + ':' + self.dirs['new'] + '/' + plotdirstr
    check_cmd += ' --outdir ' + www_dir + '/' + plotdirstr
    check_call(check_cmd.split())
check_call(['./bin/permissify-www', www_dir])

#     env.Command('test/_results/%s.passed' % name, out,
#                 './bin/diff-parameters.py --arg1 test/regression/parameters/' + actions[name]['target'] + ' --arg2 ' + stashdir + '/test/' + actions[name]['target'] + ' && touch $TARGET')
