#!/usr/bin/env python
# run sequence of commands to build models and plot performance
from subprocess import check_call
import os
import argparse

# ----------------------------------------------------------------------------------------
def run_that_shit(cmd_str):
    print 'RUN', cmd + cmd_str
    check_call([cmd,] + cmd_str.split())

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
# parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
parser.add_argument('--label', required=True)  # label for this test run. e.g. results are written to dirs with this name
parser.add_argument('--n-queries', default='50000')  # label for this test run. e.g. results are written to dirs with this name
# parser.add_argument('--skip-simulation', action='store_true')  # skip param caching on data and subsequent simulation?
parser.add_argument('--datadir', default='data/many')
args = parser.parse_args()

cmd = './runpart.py'
common_args = ' --n_procs 10 --datadir ' + args.datadir
simu_file = 'caches/recombinator/performance/' + args.label + '/simu.csv'
param_dir = 'caches/performance/' + args.label
plot_dir = os.getenv('www') + '/partis/performance/' + args.label

# cache parameters from data
cmd_str = ' --cache_parameters --seqfile test/every-hundredth-data.tsv.bz2 --is_data --skip_unproductive' + common_args
cmd_str += ' --parameter_dir ' + param_dir + '/data'
cmd_str += ' --plotdir ' + plot_dir + '/params/data'
cmd_str += ' --n_max_queries ' + args.n_queries
run_that_shit(cmd_str)

# simulate based on data parameters
cmd_str = ' --simulate --outfname ' + simu_file + common_args
cmd_str += ' --parameter_dir ' + param_dir + '/data/hmm_parameters'
cmd_str += ' --n_max_queries ' + str(int(float(args.n_queries) / 5))
run_that_shit(cmd_str)

# cache parameters from this simulation
cmd_str = ' --cache_parameters --seqfile ' + simu_file + common_args
cmd_str += ' --parameter_dir ' + param_dir + '/simu'
cmd_str += ' --plotdir ' + plot_dir + '/params/simu'
run_that_shit(cmd_str)

# run point estimation on this simulation (should split the simulation into training and testing sets)
cmd_str = ' --point_estimate --plot_performance --seqfile ' + simu_file + common_args
cmd_str += ' --parameter_dir ' + param_dir + '/simu/hmm_parameters'
cmd_str += ' --plotdir ' + plot_dir
run_that_shit(cmd_str)
