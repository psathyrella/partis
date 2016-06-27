#!/usr/bin/env python
import sys
import os
import random
from subprocess import check_call
sys.path.insert(1, './python')

import utils
import glutils

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print 'RUN', cmd_str
    sys.stdout.flush()
    check_call(cmd_str.split())

outdir = '_tmp/light-chain'
base_cmd = './bin/partis'

dj_genes = 'IGHD6-19*01:IGHJ4*02'
v_genes = 'IGHV3-71*01' #:IGHV1-18*01'
all_genes= v_genes + ':' + dj_genes

# glutils.write_glfo(outdir + '/germline-set', input_dir=outdir + '/imgt', only_genes=all_genes.split(':'), debug=True)

# glutils.write_glfo('imgt', chain='l', input_dir='data/imgt', debug=True, generate_new_alignment=True)
assert False  # needs updating
glutils.write_glfo('imgt', chain='k', input_dir='imgt', debug=True, generate_new_alignment=True)
# glutils.write_glfo(outdir + '/germline-set', input_dir=outdir + '/imgt', debug=True)
sys.exit()

# simulate
cmd_str = base_cmd + ' simulate --n-sim-events 10 --simulate-partially-from-scratch --mutation-multiplier 0.5 --debug 1 --n-trees 10'
cmd_str += ' --initial-datadir ' + outdir + '/germline-set'
cmd_str += ' --outfname ' + outdir + '/simu.csv'
# run(cmd_str)

# cache parameters
cmd_str = base_cmd + ' cache-parameters --chain-weight light --light-chain-locus kappa --debug 1'
cmd_str += ' --infname ' + outdir + '/simu.csv'
cmd_str += ' --initial-datadir ' + outdir + '/germline-set'
cmd_str += ' --outfname ' + 'tmp.csv'
run(cmd_str)
