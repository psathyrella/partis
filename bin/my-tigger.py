#!/usr/bin/env python
import sys
import os
from subprocess import check_call
sys.path.insert(1, './python')

import utils
# original = 'GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGTTTCACCTTCAGTGACTACTACATGAGCTGGGTCCGCCAGGCTCCCGGGAAGGGGCTGGAGTGGGTAGGTTTCATTAGAAACAAAGCTAATGGTGGGACAACAGAATAGACCACGTCTGTGAAAGGCAGATTCACAATCTCAAGAGATGATTCCAAAAGCATCACCTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGA'
# snpd = 'GAGGTGCAGCTGGTGGAGTCGGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGTTTCACCTTCAGTGACTACTACATGAGCTCGGTCCGCCAGGCTCCCGGGAAGGGGCTGGAGTGGGTAGGTTTCATTAGAAACAAAGCTAATGGTGAGACAACAGAATAGACCACGTCTGTGAAAGGCATATTCACAATCTCAAGAGATGATTCCAAAAGCATCACCTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGA'
# oh_one = 'GAGGTGCAGCTGGTGGAGTCCGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTGACTACTACATGAGCTGGGTCCGCCAGGCTCCCGGGAAGGGGCTGGAGTGGGTAGGTTTCATTAGAAACAAAGCTAATGGTGGGACAACAGAATAGACCACGTCTGTGAAAGGCAGATTCACAATCTCAAGAGATGATTCCAAAAGCATCACCTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTGTATTACTGTGCGAGAGA'
# utils.color_mutants(original, snpd, print_result=True)
# print utils.color_mutants(original, oh_one)
# sys.exit()
# # python -m cProfile -s tottime -o prof.out ' + 

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print 'RUN', cmd_str
    sys.stdout.flush()
    check_call(cmd_str.split())

outdir = '_tmp/allele-finder'
# param_dir = os.getcwd() + '/test/reference-results/test/parameters/simu/hmm'
original_param_dir = '/fh/fast/matsen_e/dralph/work/partis-dev/_output/021-018/simu-3-leaves-1.0-mutate/hmm'
base_cmd = './bin/partis'

dj_genes = 'IGHD6-19*01:IGHJ4*02'
existing_genes = 'IGHV3-71*01' + ':' + dj_genes
new_allele = 'IGHV3-71*03'

# simulate
cmd_str = base_cmd + ' simulate --parameter-dir ' + original_param_dir + ' --n-sim-events 100 --n-procs 5'
cmd_str += ' --only-genes ' + existing_genes + ':' + new_allele
cmd_str += ' --uniform-vj-choice-probs'
cmd_str += ' --outfname ' + outdir + '/simu.csv'
run(cmd_str)

snps_to_add = None #{'IGHV3-71*03' : 4} #{'IGHV3-69-1*02' : 3}
utils.rewrite_germline_fasta('data/imgt', outdir + '/germlines', only_genes=existing_genes.split(':'), snps_to_add=snps_to_add)

# cache-parameters
cmd_str = base_cmd + ' cache-parameters --infname ' + outdir + '/simu.csv' + ' --n-procs 10'
cmd_str += ' --datadir ' + outdir + '/germlines'
cmd_str += ' --only-genes ' + existing_genes
cmd_str += ' --parameter-dir ' + outdir + '/simu'
cmd_str += ' --plotdir ' + os.getenv('www') + '/partis/tmp/tigger'
run(cmd_str)

# annotate
cmd_str = base_cmd + ' run-viterbi --infname ' + outdir + '/simu.csv --n-procs 10'
# cmd_str += ' --presto-output'
# cmd_str += ' --only-genes ' + existing_genes
cmd_str += ' --parameter-dir ' + outdir + '/simu'
cmd_str += ' --datadir ' + outdir + '/germlines'
cmd_str += ' --outfname ' + outdir + '/run-viterbi.csv'
# cmd_str += ' --debug 1 --n-max-queries 1'  # --is-simu
# run(cmd_str)
