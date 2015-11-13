#!/usr/bin/env python
import os
import glob
from collections import OrderedDict
import sys
sys.path.insert(1, './python')
from baseutils import get_extra_str

stashdir = '_output'
datafname = 'test/mishmash.csv'  #test/adaptive-A-250.tsv.bz2

common_extras = ['--seed', '1', '--n-procs', '5']

actions = OrderedDict()
# key is name, value is target (note that the target corresponds to a directory or file in <stashdir>
actions['cache-data-parameters'] = {'target' : 'data', 'extras' : []}
actions['simulate'] = {'target' : 'simu.csv', 'extras' : ['--n-sim-events', '150']}
actions['cache-simu-parameters'] = {'target' : 'simu', 'extras' : []}
actions['plot-performance'] = {'target' : 'simu-performance', 'extras' : []}

tests = OrderedDict()
# first add the tests that run over the framework (using run-driver.py)
for action, config in actions.items():
    tests[action] = './bin/run-driver.py --label test --datafname ' + datafname + ' --stashdir _output' + ' --action ' + action + get_extra_str(config['extras'] + common_extras)
    if action == 'plot-performance':
        tests[action] += ' --plotdir ' + stashdir

simu_parameter_dir = 'test/regression/parameters/simu/hmm'
data_parameter_dir = 'test/regression/parameters/data/hmm'
# then add the simple, few-sequence tests (using partis.py)
script = './bin/partis.py'
tests['single-point-estimate'] = script + ' --action run-viterbi --seqfile test/regression/parameters/simu.csv --parameter-dir ' + simu_parameter_dir + ' --n-max-queries 3 --debug 1 ' + ' '.join(common_extras)
tests['viterbi-pair'] = script + ' --action run-viterbi --n-sets 2 --all-combinations --seqfile test/regression/parameters/simu.csv --parameter-dir ' + simu_parameter_dir + ' --debug 1 --n-max-queries 3 ' + ' '.join(common_extras)
tests['partition-data'] = script + ' --action partition --seqfile ' + datafname + ' --is-data --random-divvy --parameter-dir ' + data_parameter_dir + ' --n-max-queries 30 --n-procs 5 --debug 1 ' + ' '.join(common_extras)
tests['partition-simu'] = script + ' --action partition --seqfile test/regression/parameters/simu.csv --random-divvy --parameter-dir ' + simu_parameter_dir + ' --n-max-queries 30 --n-procs 5 --debug 1 ' + ' '.join(common_extras)
tests['naive-hamming-partition-simu'] = script + ' --action partition --naive-hamming --seqfile test/regression/parameters/simu.csv --random-divvy --parameter-dir ' + simu_parameter_dir + ' --n-max-queries 30 --n-procs 5 --debug 1 ' + ' '.join(common_extras)
tests['vsearch-hamming-partition-simu'] = script + ' --action partition --naive-vsearch --seqfile test/regression/parameters/simu.csv --random-divvy --parameter-dir ' + simu_parameter_dir + ' --n-max-queries 30 --n-procs 5 --debug 1 ' + ' '.join(common_extras)

# ----------------------------------------------------------------------------------------
all_passed = 'test/_results/ALL.passed'
individual_passed = ['test/_results/{}.passed'.format(name) for name in tests.keys()]

for path in individual_passed + [all_passed]:
    if os.path.exists(path):
        print 'removing', path
        os.remove(path)

for name, cmd_str in tests.items():
    print '%30s   %s' % (name, cmd_str)
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
