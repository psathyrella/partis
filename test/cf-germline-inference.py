#!/usr/bin/env python
import random
import argparse
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
def cf_nsnps(args, glfo):
    v_gene = args.v_genes[0]
    for nsnp in args.nsnp_list:
        cmd = base_cmd + ' --n-procs 1 --n-tests 3 --n-sim-events 10'
        cmd += ' --sim-v-genes ' + v_gene
        cmd += ' --inf-v-genes ' + v_gene
        cmd += ' --nsnp-list ' + str(nsnp)
        cmd += ' --outdir /fh/fast/matsen_e/dralph/partis/allele-finder/' + 'nsnp-' + str(nsnp)
        run(cmd)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['nsnp'])
parser.add_argument('--nsnp-list', default='1')
parser.add_argument('--v-genes', default='IGHV4-39*01')
args = parser.parse_args()

args.nsnp_list = utils.get_arg_list(args.nsnp_list, intify=True)
args.v_genes = utils.get_arg_list(args.v_genes)

locus = 'igh'
glfo = glutils.read_glfo('data/germlines/human', locus=locus)

# ----------------------------------------------------------------------------------------
if args.action == 'nsnp':
    cf_nsnps(args, glfo)
