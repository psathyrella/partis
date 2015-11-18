#!/usr/bin/env python
import argparse
import os
import csv
import glob
from collections import OrderedDict
from subprocess import check_call, check_output
import sys
sys.path.insert(1, './python')
from baseutils import get_extra_str
import utils
from hist import Hist
from clusterpath import ClusterPath

parser = argparse.ArgumentParser()
parser.add_argument('--dont-run', action='store_true', help='don\'t actually run the tests, presumably so you can just check the results ')
parser.add_argument('--dont-plot', action='store_true', help='don\'t make all the comparison plots')
args = parser.parse_args()

stashdir = 'test/_new-results'
referencedir = 'test/reference-results'
datafname = 'test/mishmash.csv'
label = 'test'
# these are the top 10 v and d genes, and top six js, from mishmash.csv. Restricting to these should make testing much more stable and much faster.
only_genes = 'IGHV4-61*08:IGHV3-48*01:IGHV5-51*02:IGHV3-69-1*02:IGHV1/OR15-1*04:IGHV3-66*03:IGHV3-23D*01:IGHV3-71*03:IGHV1-2*04:IGHV1-2*02:IGHD3-16*02:IGHD2-2*03:IGHD2-8*01:IGHD3-22*01:IGHD6-13*01:IGHD4-17*01:IGHD6-19*01:IGHD3-10*01:IGHD2-15*01:IGHD2-21*02:IGHJ5*02:IGHJ3*02:IGHJ2*01:IGHJ1*01:IGHJ6*03:IGHJ4*02'
common_extras = ['--seed', '1', '--n-procs', '10', '--only-genes', only_genes]
param_dir = stashdir + '/' + label  # change this to <referencedir> if you want to see if these results change without checking if parameter stashing has changed
ref_simfname = referencedir + '/' + label + '/simu.csv'
new_simfname = param_dir + '/simu.csv'
simu_param_dir = param_dir + '/parameters/simu/hmm'
data_param_dir = param_dir + '/parameters/data/hmm'
run_driver = './bin/run-driver.py --label ' + label + ' --stashdir ' + stashdir
partis = './bin/partis.py'

tests = OrderedDict()
n_partition_queries = '100'
# tests['partition-ref-simu']         = {'bin' : partis, 'action' : 'partition',   'extras' : ['--seqfile', ref_simfname, '--parameter-dir', simu_param_dir, '--n-max-queries', n_partition_queries]}
# tests['point-partition-ref-simu']   = {'bin' : partis, 'action' : 'partition',   'extras' : ['--naive-hamming', '--seqfile', ref_simfname, '--parameter-dir', simu_param_dir, '--n-max-queries', n_partition_queries]}
# tests['vsearch-partition-ref-simu'] = {'bin' : partis, 'action' : 'partition',   'extras' : ['--naive-vsearch', '--seqfile', ref_simfname, '--parameter-dir', simu_param_dir, '--n-max-queries', n_partition_queries]}

tests['cache-data-parameters']  = {'bin' : run_driver, 'extras' : ['--skip-unproductive']}
tests['simulate']  = {'bin' : run_driver, 'extras' : ['--n-sim-events', 100, '--mimic-data-read-length']}
tests['cache-simu-parameters']  = {'bin' : run_driver, 'extras' : []}

tests['annotate-simu']  = {'bin' : run_driver, 'action' : 'plot-performance', 'extras' : []}
tests['single-point-estimate']      = {'bin' : partis, 'action' : 'run-viterbi', 'extras' : ['--seqfile', new_simfname, '--parameter-dir', simu_param_dir, '--n-max-queries', '10']}
tests['partition-data']             = {'bin' : partis, 'action' : 'partition',   'extras' : ['--seqfile', datafname, '--parameter-dir', data_param_dir, '--is-data', '--skip-unproductive', '--n-max-queries', n_partition_queries]}
tests['partition-new-simu']         = {'bin' : partis, 'action' : 'partition',   'extras' : ['--seqfile', new_simfname, '--parameter-dir', simu_param_dir, '--n-max-queries', n_partition_queries]}
tests['point-partition-new-simu']   = {'bin' : partis, 'action' : 'partition',   'extras' : ['--seqfile', new_simfname, '--parameter-dir', simu_param_dir, '--naive-hamming', '--n-max-queries', n_partition_queries]}
tests['vsearch-partition-new-simu'] = {'bin' : partis, 'action' : 'partition',   'extras' : ['--seqfile', new_simfname, '--parameter-dir', simu_param_dir, '--naive-vsearch', '--n-max-queries', n_partition_queries]}
for name, info in tests.items():
    if args.dont_run:
        continue
    action = info['action'] if 'action' in info else name
    cmd_str = info['bin'] + ' --action ' + action
    if info['bin'] == partis:
        cmd_str += ' ' + ' '.join(info['extras'] + common_extras)
        cmd_str += ' --outfname ' + stashdir + '/' + name + '.csv'
    else:
        cmd_str += get_extra_str(info['extras'] + common_extras)
        if action == 'cache-data-parameters':
            cmd_str += ' --datafname ' + datafname
    print 'TEST %30s   %s' % (name, cmd_str)
    check_call(cmd_str.split())

# ----------------------------------------------------------------------------------------
# collect summary performance info from a few places
print 'reading performance info'
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

perf_info = OrderedDict()
perfdir = param_dir + '/plots/simu-performance'
for method in ['sw', 'hmm']:
    print '  ', method
    for region in utils.regions:
        fraction_correct = read_performance_file(perfdir + '/' + method + '/plots/' + region + '_gene.csv', 'contents', only_ibin=1)
        print '    %s %.3f' % (region, fraction_correct)
        perf_info[method + '-' + region + '_gene_correct'] = fraction_correct
    hamming_hist = Hist(fname=perfdir + '/' + method + '/plots/hamming_to_true_naive.csv')
    print '    mean hamming %.2f' % hamming_hist.get_mean()
    perf_info[method + '-mean_hamming'] = hamming_hist.get_mean()

print '    adj mi      test                    description'
for ptest in [k for k in tests.keys() if 'partition' in k]:
    cp = ClusterPath(-1)
    cp.readfile(stashdir + '/' + ptest + '.csv')
    if 'data' in ptest:
        ref_cp = ClusterPath(-1)
        ref_cp.readfile(referencedir + '/' + ptest + '.csv')
        # adj mi between the reference and the new data partitions
        perf_info[ptest] = utils.adjusted_mutual_information(cp.partitions[cp.i_best], ref_cp.partitions[ref_cp.i_best])
        print '    %5.2f   %-28s   to reference partition' % (perf_info[ptest], ptest)
    else:
        perf_info[ptest] = cp.adj_mis[cp.i_best]  # adj mi to true partition
        print '    %5.2f   %-28s   to true partition' % (perf_info[ptest], ptest)

# ----------------------------------------------------------------------------------------
# check against reference file
eps = 0.05
print 'comparing to reference file'
with open(referencedir + '/performance-info.csv') as perf_file:
    reader = csv.DictReader(perf_file)
    line = reader.next()
    for name, new_val in perf_info.items():
        ref_val = float(line[name])
        if abs(ref_val - new_val) / ref_val > eps:
            print '  %s %-28s       ref %-5.2f  new %-5.2f' % (utils.color('red', 'different:'), name, ref_val, new_val)
        else:
            print '   ok   %-28s' % name

# ----------------------------------------------------------------------------------------
# write performance info file
with open(stashdir + '/performance-info.csv', 'w') as perf_file:
    writer = csv.DictWriter(perf_file, perf_info.keys())
    writer.writeheader()
    writer.writerow(perf_info)

# ----------------------------------------------------------------------------------------
# make a bunch of comparison plots
plotdirs = ['test/plots/data/sw', 'test/plots/data/hmm', 'test/plots/simu/hmm-true', 'test/plots/simu-performance/sw', 'test/plots/simu-performance/hmm']
base_check_cmd = './bin/compare.py --dont-calculate-mean-info --graphify --linewidth 1 --markersizes 2:1 --names reference:new'  # --colors 595:807:834 --scale-errors 1.414
if os.getenv('www') is None:
    www_dir = '_test-plots'
else:
    www_dir = os.getenv('www') + '/partis/test'
# # if you want to do *all* the subdirs use this:
# recursive_subdirs = []
# for plotdir in plotdirs:
#     find_plotdirs_cmd = 'find ' + referencedir + '/' + plotdir + ' -name "*.csv" -exec dirname {} \;|sed \'s@/plots$@@\' | sort | uniq'
#     recursive_subdirs += check_output(find_plotdirs_cmd, shell=True).split()
for plotdir in plotdirs:
    if args.dont_plot:
        continue
    plotdirstr = plotdir.replace(referencedir + '/', '')
    check_cmd = base_check_cmd + ' --plotdirs '  + referencedir + '/' + plotdirstr + ':' + stashdir + '/' + plotdirstr
    check_cmd += ' --outdir ' + www_dir + '/' + plotdirstr
    check_call(check_cmd.split())
check_call(['./bin/permissify-www', www_dir])

#     env.Command('test/_results/%s.passed' % name, out,
#                 './bin/diff-parameters.py --arg1 test/regression/parameters/' + actions[name]['target'] + ' --arg2 ' + stashdir + '/test/' + actions[name]['target'] + ' && touch $TARGET')
