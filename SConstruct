import os
import glob
from collections import OrderedDict
import sys
from SCons.Script import Command, Depends

env = Environment(ENV=os.environ, SHELL='/bin/bash')
sys.path.append(os.getenv('HOME') + '/bin')

Alias('validate', '_output/validation/valid.out')
env.Command('_output/validation/valid.out', './bin/run-driver.py', './bin/run-driver.py --label validation --plotdir _output/validation/plots --datafname test/adaptive-A-250.tsv.bz2 && touch $TARGET')

# ----------------------------------------------------------------------------------------
# scons test

Alias('test', 'test/_results/ALL.passed')


def get_extra_str(extra_list):
    modified_list = [ex.replace(':', '.').replace('--', '__') for ex in extra_list]
    return ' --extra-args ' + ':'.join(modified_list)

stashdir = '_output'
datafname = 'test/mishmash.csv'  #test/adaptive-A-250.tsv.bz2

common_extras = ['--seed', '1', '--no-plot']
semi_common_extras = ['--match-mismatch', '5:1']
base_cmd = './bin/run-driver.py --label test --datafname ' + datafname + ' --stashdir _output --n-procs 5 --plotdir ' + stashdir + ' --no-skip-unproductive'

actions = OrderedDict()
# key is name, value is target (note that the target corresponds to a directory or file in <stashdir>
actions['cache-data-parameters'] = {'target' : 'data', 'extras' : semi_common_extras}
actions['simulate'] = {'target' : 'simu.csv', 'extras' : ['--random-number-of-leaves', '--n-sim-events', '150']}
actions['cache-simu-parameters'] = {'target' : 'simu', 'extras' : semi_common_extras}
actions['plot-performance'] = {'target' : 'performance', 'extras' : semi_common_extras}

tests = OrderedDict()
# first add the tests that run over the framework (using run-driver.py)
for action, config in actions.items():
    tests[action] = base_cmd + ' --action ' + action + get_extra_str(config['extras'] + common_extras)

simu_parameter_dir = 'test/regression/parameters/simu/hmm'
data_parameter_dir = 'test/regression/parameters/data/hmm'
# then add the simple, few-sequence tests (using partis.py)
cmd = './bin/partis.py'
tests['single-point-estimate'] = cmd + ' --action run-viterbi --seqfile test/regression/parameters/simu.csv --parameter-dir ' + simu_parameter_dir + ' --n-max-queries 3 --debug 1 ' + ' '.join(common_extras + semi_common_extras)
tests['partition-data'] = cmd + ' --action partition --seqfile ' + datafname + ' --is-data --parameter-dir ' + data_parameter_dir + ' --n-max-queries 30 --n-procs 5 --debug 1 ' + ' '.join(common_extras + semi_common_extras)
tests['partition-simu'] = cmd + ' --action partition --seqfile test/regression/parameters/simu.csv --parameter-dir ' + simu_parameter_dir + ' --n-max-queries 30 --n-procs 5 --debug 1 ' + ' '.join(common_extras + semi_common_extras)
tests['viterbi-pair'] = cmd + ' --action run-viterbi --n-sets 2 --all-combinations --seqfile test/regression/parameters/simu.csv --parameter-dir ' + simu_parameter_dir + ' --debug 1 --n-max-queries 3 ' + ' '.join(common_extras + semi_common_extras)

# ----------------------------------------------------------------------------------------
all_passed = 'test/_results/ALL.passed'
individual_passed = ['test/_results/{}.passed'.format(name) for name in tests.keys()]

for path in individual_passed + [all_passed]:
    if os.path.exists(path):
        print 'removing', path
        os.remove(path)

for name, test_cmd in tests.items():
    out = 'test/_results/%s.out' % name
    Depends(out, glob.glob('python/*.py') + ['packages/ham/bcrham',])
    if name in actions:
        env.Command(out, cmd, test_cmd + ' && touch $TARGET')  # it's kind of silly to put partis.py as the SOURCE, but you've got to put *something*, and we've already got the deps covered...
        env.Command('test/_results/%s.passed' % name, out,
                    './bin/diff-parameters.py --arg1 test/regression/parameters/' + actions[name]['target'] + ' --arg2 ' + stashdir + '/test/' + actions[name]['target'] + ' && touch $TARGET')
    else:
        env.Command(out, cmd, test_cmd + ' --outfname $TARGET')
        # touch a sentinel `passed` file if we get what we expect
        # NOTE [vdj]: regex is a hack. I can't figure out a.t.m. why the missing genes come up in a different order each time
        env.Command('test/_results/%s.passed' % name,
                ['test/regression/%s.out' % name, out],
                'diff -ub ${SOURCES[0]} ${SOURCES[1]} && touch $TARGET')  # ignore the lines telling you how long things took

# Set up sentinel dependency of all passed on the individual_passed sentinels.
Command(all_passed,
        individual_passed,
        'cat $SOURCES > $TARGET')
