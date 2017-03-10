#!/usr/bin/env python
import sys
from subprocess import check_call
sys.path.insert(1, './python')

import utils
import glutils

base_cmd = './bin/test-allele-finding.py'

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print '%s %s' % (utils.color('red', 'run'), cmd_str)
    sys.stdout.flush()
    check_call(cmd_str.split())

# ----------------------------------------------------------------------------------------
def cf_nsnps():
    for nsnp in args.nsnplist:
        cmd = base_cmd + ' --n-procs 3 --n-sim-events 300'
        cmd += ' --mut-mult ' + str(mut_mult)
        cmd += ' --outdir /fh/fast/matsen_e/dralph/partis/allele-finder/' + 'mut-mult-' + str(mut_mult)
        run(cmd)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['nsnp'])
parser.add_argument('--nsnplist', default='1:2')
args = parser.parse_args()

args.nsnplist = utils.get_arg_list(args.nsnplist)

# ----------------------------------------------------------------------------------------
if args.action == 'nsnp':
    cf_nsnps(args)
