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
parser.add_argument('--n-queries', default='10000')  # label for this test run. e.g. results are written to dirs with this name
parser.add_argument('--datadir', required=True)  #default='data/imgt')
parser.add_argument('--extra-args')  # args to pass on to commands (colon-separated) NOTE have to add space and quote like so: --extra-args ' --option'
parser.add_argument('--skip-simulation', action='store_true')  # skip param caching on data and subsequent simulation
parser.add_argument('--datafname', default='test/every-hundredth-data.tsv.bz2')
parser.add_argument('--simfname')
parser.add_argument('--only-plot-performance', action='store_true')  # also skip parameter caching on simulation

args = parser.parse_args()
args.extra_args = utils.get_arg_list(args.extra_args)
if args.only_plot_performance:
    args.skip_simulation = True

cmd = './runpart.py'
common_args = ' --n-procs 10 --datadir ' + args.datadir  #+ ' --only-genes \'IGHV1-18*01:IGHD3-10*01:IGHJ6*03\''
if args.extra_args != None:
    common_args += ' ' + ' '.join(args.extra_args)
if args.simfname == None:
    args.simfname = 'caches/recombinator/performance/' + args.label + '/simu.csv'
param_dir = 'caches/performance/' + args.label
plot_dir = os.getenv('www') + '/partis/performance/' + args.label

if not args.skip_simulation:
    # cache parameters from data
    cmd_str = ' --cache-parameters --seqfile ' + args.datafname + ' --is-data --skip-unproductive' + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/data'
    cmd_str += ' --plotdir ' + plot_dir + '/params/data'
    cmd_str += ' --n-max-queries ' + args.n_queries
    run_that_shit(cmd_str)
    
    # simulate based on data parameters
    cmd_str = ' --simulate --outfname ' + args.simfname + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/data/hmm_parameters'
    cmd_str += ' --n-max-queries ' + str(int(float(args.n_queries) / 5))
    run_that_shit(cmd_str)

if not args.only_plot_performance:
    # cache parameters from simulation
    cmd_str = ' --cache-parameters --seqfile ' + args.simfname + common_args
    cmd_str += ' --parameter-dir ' + param_dir + '/simu'
    cmd_str += ' --plotdir ' + plot_dir + '/params/simu'
    cmd_str += ' --n-max-queries ' + args.n_queries
    run_that_shit(cmd_str)

# run point estimation on simulation
cmd_str = ' --point-estimate --plot-performance --seqfile ' + args.simfname + common_args
cmd_str += ' --parameter-dir ' + param_dir + '/simu/hmm_parameters'
cmd_str += ' --plotdir ' + plot_dir
cmd_str += ' --n-max-queries ' + args.n_queries
run_that_shit(cmd_str)
