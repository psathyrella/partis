#!/usr/bin/env python
import sys
import os
from subprocess import check_call
sys.path.insert(1, './python')

import utils

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print 'RUN', cmd_str
    sys.stdout.flush()
    check_call(cmd_str.split())

outdir = '_tmp/my-tigger'
param_dir = os.getcwd() + '/test/reference-results/test/parameters/simu/hmm'
base_cmd = './bin/partis.py --parameter-dir ' + param_dir

existing_genes = 'IGHV3-69-1*02:IGHD4-17*01:IGHJ4*02'
# existing_genes = utils.test_only_genes

# simulate
cmd_str = base_cmd + ' --action simulate'
cmd_str += ' --n-sim-events 100'
cmd_str += ' --n-procs 5'
cmd_str += ' --only-genes ' + existing_genes
cmd_str += ' --outfname ' + outdir + '/simu.csv'
# run(cmd_str)

snps_to_add = {'IGHV3-69-1*02' : 1}
utils.rewrite_germline_fasta('data/imgt', outdir + '/germlines', only_genes=existing_genes.split(':'), snps_to_add=snps_to_add)
sys.exit()
# cache-parameters
cmd_str = 'python -m cProfile -s tottime -o prof.out ' + base_cmd + ' --action cache-parameters'
cmd_str += ' --seqfile ' + outdir + '/simu.csv'
cmd_str += ' --n-procs 10'
cmd_str += ' --datadir ' + outdir + '/germlines'
cmd_str += ' --only-genes ' + existing_genes
cmd_str += ' --parameter-dir ' + outdir + '/simu'
cmd_str += ' --plotdir ' + os.getenv('www') + '/tmp/test/plots'
run(cmd_str)
sys.exit()
# annotate
cmd_str = base_cmd + ' --action run-viterbi'
cmd_str += ' --seqfile ' + outdir + '/simu.csv'
cmd_str += ' --n-procs 10'
# cmd_str += ' --presto-output'
# cmd_str += ' --only-genes ' + existing_genes
cmd_str += ' --datadir ' + outdir + '/germlines'
cmd_str += ' --outfname ' + outdir + '/run-viterbi.csv'
# cmd_str += ' --debug 1 --n-max-queries 1'  # --is-simu
# run(cmd_str)

sys.exit()
rfname = '_tmp/tigger/run.R'
write_tigger_cmd(rfname)
check_call(['R', '--slave', '-f', rfname])
os.remove(rfname)
