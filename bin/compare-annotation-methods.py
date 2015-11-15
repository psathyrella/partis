#!/usr/bin/env python

import sys
import os
import argparse
from subprocess import check_call
import multiprocessing
sys.path.insert(1, './python')
import utils

all_codes = ['partis', 'multi-partis', 'ihhhmmm', 'imgt', 'igblast']

parser = argparse.ArgumentParser()
parser.add_argument('--run-codes', help='Which codes to check? (if none, we just plot preexisting results)')
parser.add_argument('--n-queries', type=int)
args = parser.parse_args()
if args.run_codes == 'all':
    args.run_codes = all_codes
elif args.run_codes == 'none':
    args.run_codes = []
else:
    args.run_codes = utils.get_arg_list(args.run_codes)

outdir = '/fh/fast/matsen_e/dralph/work/partis-dev'
simfname = outdir + '/data/performance/simu.csv'
base_plotdir = outdir + '/_compare'

for code in args.run_codes:
    if code not in all_codes:
        raise Exception('ERROR bad code: ' + code)
    if 'partis' in code:
        n_procs = max(1, multiprocessing.cpu_count() / 2)
        for action in ('cache-simu-parameters', 'plot-performance'):
            cmd = './bin/run-driver.py --label comparisons-' + code + ' --action ' + action + ' --plotdir ' + base_plotdir + '/' + code + ' --simfname ' + simfname + ' --n-procs 50:10'  # + str(n_procs)
            if 'multi-' in code:
                cmd += ' --extra-args __n-sets:5:--slurm:--workdir:_tmp/foop'
            else:
                cmd += ' --extra-args __slurm:--workdir:_tmp/foop'
            print cmd
            check_call(cmd.split())
    else:
        cmd = './python/' + code + 'parser.py --plotdir ' + base_plotdir + '/' + code + ' --simfname ' + simfname
        if args.n_queries is not None:
            cmd += ' --n-queries ' + str(args.n_queries)
        if code == 'imgt':
            cmd += ' --indir ' + os.path.dirname(simfname) + '/imgt/smaller_repertoire_higher_mutation/IMGT_HighV-QUEST_individual_files_folder'
        elif code == 'igblast':
            cmd += ' --infname ' + os.path.dirname(simfname) + '/igblast/simu.html'
        elif code == 'ihhhmmm':
            cmd += ' --indir ' + os.path.dirname(simfname) + '/ihhhmmm'
        print cmd
        check_call(cmd.split())

# ----------------------------------------------------------------------------------------
plotdirs = [ 'multi-partis/hmm/performance',
             'partis/hmm/performance',
             'partis/sw/performance',
             'ihhhmmm/',
             'igblast/',
             'imgt/' ]
plotdirs = [ base_plotdir + '/' + pd for pd in plotdirs ]
names = 'partis@(k=5):partis@(k=1):ighutil:iHMMunealign:igblast:imgt'

plot_cmd = './bin/compare.py --plot-performance --no-errors --normalize --markersizes 0 --linestyles 2:1:1:1:1:1 --outdir ' + base_plotdir + '/all-vs-all' \
           + ' --graphify --plotdirs ' + ':'.join(plotdirs) + ' --names ' + names + ' --colors 810:634:596:418:798:869'
print plot_cmd
check_call(plot_cmd.split())
