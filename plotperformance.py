#!/usr/bin/env python
# run sequence of commands to build models and plot performance
from subprocess import check_call
import os
import sys
import argparse

sys.path.insert(1, './python')
import utils

# ----------------------------------------------------------------------------------------
def run_that_shit(cmd_str):
    print 'RUN', cmd + cmd_str
    check_call([cmd,] + cmd_str.split())

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
# parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
parser.add_argument('--label', required=True)  # label for this test run. e.g. results are written to dirs with this name
parser.add_argument('--n-queries', default='-1')  # label for this test run. e.g. results are written to dirs with this name
parser.add_argument('--n-sim-seqs', default='10000')
parser.add_argument('--datadir', required=True)  #default='data/imgt')
parser.add_argument('--extra-args')  # args to pass on to commands (colon-separated) NOTE have to add space and quote like so: --extra-args ' --option'
parser.add_argument('--datafname', default='test/every-hundredth-data.tsv.bz2')
parser.add_argument('--simfname')
parser.add_argument('--plotdir')
parser.add_argument('--actions', default='cache-data-parameters:simulate:cache-simu-parameters:plot-performance', help='Colon-separated list of actions to perform')

args = parser.parse_args()
args.extra_args = utils.get_arg_list(args.extra_args)
args.actions = utils.get_arg_list(args.actions)

cmd = './runpart.py'
common_args = ' --n-procs 10 --datadir ' + args.datadir  #+ ' --only-genes \'IGHV1-18*01:IGHD3-10*01:IGHJ6*03\''
if args.extra_args != None:
    common_args += ' ' + ' '.join(args.extra_args)
if args.simfname == None:
    args.simfname = 'caches/recombinator/performance/' + args.label + '/simu.csv'
param_dir = 'caches/performance/' + args.label
if args.plotdir == None:
    args.plotdir = os.getenv('www') + '/partis/performance/' + args.label

if 'cache-data-parameters' in args.actions:
    # cache parameters from data
    cmd_str = ' --cache-parameters --seqfile ' + args.datafname + ' --is-data --skip-unproductive' + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/data'
    cmd_str += ' --plotdir ' + args.plotdir + '/params/data'
    cmd_str += ' --n-max-queries ' + args.n_queries
    run_that_shit(cmd_str)

if 'simulate' in args.actions:
    n_reco_events = float(args.n_sim_seqs) / 5  # a.t.m. I'm just hard coding five seqs per reco event
    assert n_reco_events > 0
    # simulate based on data parameters
    cmd_str = ' --simulate --outfname ' + args.simfname + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/data/hmm_parameters'
    cmd_str += ' --n-max-queries ' + str(int(n_reco_events))
    run_that_shit(cmd_str)

if 'cache-simu-parameters' in args.actions:
    # cache parameters from simulation
    cmd_str = ' --cache-parameters --seqfile ' + args.simfname + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/simu'
    cmd_str += ' --plotdir ' + args.plotdir + '/params/simu'
    cmd_str += ' --n-max-queries ' + args.n_queries
    run_that_shit(cmd_str)

if 'plot-performance' in args.actions:  # run point estimation on simulation
    cmd_str = ' --point-estimate --plot-performance --seqfile ' + args.simfname + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/simu/hmm_parameters'
    cmd_str += ' --plotdir ' + args.plotdir
    cmd_str += ' --n-max-queries ' + args.n_queries
    run_that_shit(cmd_str)
