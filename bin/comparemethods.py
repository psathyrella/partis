#!/usr/bin/env python

import sys
import argparse
from subprocess import check_call
import multiprocessing
sys.path.insert(1, './python')
import utils

all_codes = ['partis', 'multi-partis', 'ihhhmmm', 'imgt', 'igblast']

parser = argparse.ArgumentParser()
parser.add_argument('--run-codes', default='all', help='Which codes to check? (if none, we just plot preexisting results)')
args = parser.parse_args()
if args.run_codes == 'all':
    args.run_codes = all_codes
elif args.run_codes == 'none':
    args.run_codes = []
else:
    args.run_codes = utils.get_arg_list(args.run_codes)

base_plotdir = '_compare'

for code in args.run_codes:
    if code not in all_codes:
        raise Exception('ERROR bad code: ' + code)
    if 'partis' in code:
        n_procs = max(1, multiprocessing.cpu_count() / 2)
        for action in ('cache-simu-parameters', 'plot-performance'):
            cmd = './plotperformance.py --label comparisons-' + code + ' --action ' + action + ' --plotdir ' + base_plotdir + '/' + code + ' --simfname data/performance/simu.csv --n-procs ' + str(n_procs)
            if 'multi-' in code and action == 'plot-performance':
                cmd += ' --extra-args __n-sets:5'
            print cmd
            check_call(cmd.split())
    else:
        cmd = './python/' + code + 'parser.py --plotdir ' + base_plotdir + '/' + code
        print cmd
        check_call(cmd.split())

plotdirs = [ 'multi-partis/hmm/performance',
             'partis/hmm/performance',
             'partis/sw/performance',
             'ihhhmmm/',
             'igblast/',
             'imgt/' ]
plotdirs = [ base_plotdir + '/' + pd for pd in plotdirs ]
names = 'partis@\(k=5\):partis@\(k=1\):ighutil:iHMMunealign:igblast:imgt'

plot_cmd = './python/compare.py --plot-performance --no-errors --normalize --markersizes 0 --linestyles 2:1:1:1:1:1 --outdir ' + base_plotdir + '/all-vs-all' \
           + ' --graphify --plotdirs ' + ':'.join(plotdirs) + ' --names ' + names + ' --colors 810:634:596:418:798:869'
check_call(plot_cmd.split())
