#!/usr/bin/env python
# run sequence of commands to build models and plot performance
from subprocess import check_call
import os
import sys
import argparse
import multiprocessing
sys.path.insert(1, './python')
import utils

# ----------------------------------------------------------------------------------------
def run_command(cmd_str):
    cmd_str += ' --print-git-commit'
    print 'RUN', cmd + cmd_str
    sys.stdout.flush()
    check_call([cmd,] + cmd_str.split())

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
if os.getenv('USER') is None:  # if we're in docker
    fsdir = '/partis'  # shouldn't need it in this case
else:
    fsdir = '/fh/fast/matsen_e/' + os.getenv('USER') + '/work/partis-dev'
parser.add_argument('--label', required=True)  # label for this test run. e.g. results are written to dirs with this name
parser.add_argument('--stashdir', required=True)  #default=fsdir + '/_output')
parser.add_argument('--extra-args')  # args to pass on to commands (colon-separated) NOTE have to add space and quote like so: --extra-args __option (NOTE replaces __ with --, and . with :)
parser.add_argument('--datafname')
parser.add_argument('--is-simu', action='store_true')
parser.add_argument('--old-style-dir-structure', action='store_true')
parser.add_argument('--simfname')
parser.add_argument('--outfname')
parser.add_argument('--plotdir')

all_actions = ('cache-data-parameters', 'simulate', 'cache-simu-parameters', 'plot-performance', 'partition', 'run-viterbi')
default_actions = all_actions[ : -2]
parser.add_argument('--actions', default=':'.join(default_actions), choices=all_actions, help='Colon-separated list of actions to perform')
# parser.add_argument('--n-procs', default=str(max(1, multiprocessing.cpu_count() / 2)))

args = parser.parse_args()
args.extra_args = utils.get_arg_list(args.extra_args)
args.actions = utils.get_arg_list(args.actions)
if args.plotdir is None:
    args.plotdir = args.stashdir + '/' + args.label + '/plots'
# assert '--n-procs' not in args.extra_args

cmd = './bin/partis.py'
# common_args = ' --n-procs ' + str(args.n_procs)
common_args = ''

if args.extra_args is not None:
    common_args += ' ' + ' '.join(args.extra_args).replace('__', '--').replace(',', ':').replace('+', ' ')
if args.simfname is None:
    args.simfname = args.stashdir + '/' + args.label + '/simu.csv'
if args.old_style_dir_structure:  # oh, backwards compatibility, you're such a craven old bitch
    param_dir = args.stashdir + '/' + args.label
else:
    param_dir = args.stashdir + '/' + args.label + '/parameters'

assert '--outfname' not in common_args

if 'cache-data-parameters' in args.actions:
    if args.datafname is None or not os.path.exists(args.datafname):
        raise Exception('ERROR datafname d.n.e.: ' + str(args.datafname))
    # cache parameters from data
    cmd_str = ' cache-parameters --seqfile ' + args.datafname + common_args
    # cmd_str += ' --skip-unproductive'
    cmd_str += ' --parameter-dir ' + param_dir + '/data'
    if args.plotdir is not None:
        cmd_str += ' --plotdir ' + args.plotdir + '/data'
    run_command(cmd_str)

if 'simulate' in args.actions:
    # simulate based on data parameters
    cmd_str = ' simulate --outfname ' + args.simfname + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/data/hmm'
    run_command(cmd_str)

sim_name = os.path.basename(args.simfname).replace('.csv', '')  # 'sim', if simfname is just 'simu.csv'
if 'cache-simu-parameters' in args.actions:
    # cache parameters from simulation
    cmd_str = ' cache-parameters --seqfile ' + args.simfname + ' --is-simu' + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/' + sim_name
    if args.plotdir is not None:
        cmd_str += ' --plotdir ' + args.plotdir + '/' + sim_name
    run_command(cmd_str)

if 'plot-performance' in args.actions:  # run point estimation on simulation
    cmd_str = ' run-viterbi --plot-performance --seqfile ' + args.simfname + ' --is-simu' + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/' + sim_name + '/hmm'
    cmd_str += ' --plotdir ' + args.plotdir + '/' + sim_name + '-performance'
    run_command(cmd_str)

if 'partition' in args.actions or 'run-viterbi' in args.actions:
    assert len(args.actions) == 1
    action = args.actions[0]

    cmd_str = ' ' + action + common_args
    if not args.is_simu:
        seqfile = args.datafname
        pdir = param_dir + '/data/hmm'
        if args.outfname is None:
            args.outfname =  param_dir + '/data-' + action + '.csv'
    else:
        cmd_str += ' --is-simu'
        seqfile = args.simfname
        pdir = param_dir + '/' + sim_name + '/hmm'
        if args.outfname is None:
            args.outfname = args.simfname.replace('.csv', '-' + action + '.csv')
    cmd_str += ' --seqfile ' + seqfile
    cmd_str += ' --parameter-dir ' + pdir
    cmd_str += ' --outfname ' + args.outfname
    run_command(cmd_str)
