#!/usr/bin/env python
import argparse
import os
import glob
from collections import OrderedDict
from subprocess import check_call, check_output
import sys
sys.path.insert(1, './python')
from baseutils import get_extra_str

# parser = argparse.ArgumentParser()
# parser.add_argument('--plotdir')
# args = parser.parse_args()

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
sys.exit()

# ----------------------------------------------------------------------------------------
# make some comparison plots
checks = OrderedDict()
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
    plotdirstr = plotdir.replace(referencedir + '/', '')
    check_cmd = base_check_cmd + ' --plotdirs '  + referencedir + '/' + plotdirstr + ':' + stashdir + '/' + plotdirstr
    check_cmd += ' --outdir ' + www_dir + '/' + plotdirstr
    check_call(check_cmd.split())
check_call(['./bin/permissify-www', www_dir])

print 'TODO check that fraction correct genes, naive hamming distnace, and adj mi (& co)  don\'t change too much'
print 'TODO work out what happened to adj mi (and make sure true partition knows which queries we skipped in sw)'


# # ----------------------------------------------------------------------------------------
# for name, cmd_str in tests.items():
#     print '%30s   %s' % (name, cmd_str)
#     check_call(cmd_str.split())
    # sys.exit()

    # out = 'test/_results/%s.out' % name
    # Depends(out, glob.glob('python/*.py') + ['packages/ham/bcrham',])  # this is a woefully inadequate description of dependencies, but it's probably not worth trying to improve it
    # if name in actions:
    #     # print cmd_str
    #     # continue
    #     env.Command(out, script, cmd_str + ' && touch $TARGET')  # it's kind of silly to put partis.py as the SOURCE, but you've got to put *something*, and we've already got the deps covered...
    #     env.Command('test/_results/%s.passed' % name, out,
    #                 './bin/diff-parameters.py --arg1 test/regression/parameters/' + actions[name]['target'] + ' --arg2 ' + stashdir + '/test/' + actions[name]['target'] + ' && touch $TARGET')
    # else:
    #     # print cmd_str
    #     # continue
    #     env.Command(out, script, cmd_str + ' --outfname $TARGET')
    #     # touch a sentinel `passed` file if we get what we expect
    #     env.Command('test/_results/%s.passed' % name,
    #             ['test/regression/%s.out' % name, out],
    #             'diff -ub ${SOURCES[0]} ${SOURCES[1]} && touch $TARGET')  # ignore the lines telling you how long things took

# # Set up sentinel dependency of all passed on the individual_passed sentinels.
# Command(all_passed,
#         individual_passed,
#         'cat $SOURCES > $TARGET')
