#!/usr/bin/env python
import sys
import os
from subprocess import check_call
current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
if not os.path.exists(current_script_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % current_script_dir
sys.path.insert(1, current_script_dir)

import utils

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print 'RUN', cmd_str
    sys.stdout.flush()
    check_call(cmd_str.split())

# ----------------------------------------------------------------------------------------
def write_tigger_cmd(fname):
    cmdfile = open(fname, 'w')

    cmdfile.write("library(tigger)\n")

    # # ----------------------------------------------------------------------------------------
    # # example in docs:
    # cmdfile.write("data(sample_db, germline_ighv); annotations = sample_db; germlines = germline_ighv\n")
    # # write sample_db to a csv file:
    # # cmdfile.write("write.csv(sample_db, file='_tmp/tigger/sample_db.csv')\n")
    
    # # ----------------------------------------------------------------------------------------
    # # example in docs, but with my annotations:
    # # ./bin/partis.py --action cache-parameters --seqfile _tmp/tigger/sample_db.csv --parameter-dir _tmp/tigger/sample_db --n-procs 15 --name-column SEQUENCE_ID --seq-column SEQUENCE_INPUT
    # # ./bin/partis.py --action run-viterbi --seqfile _tmp/tigger/sample_db.csv --parameter-dir _tmp/tigger/sample_db/hmm --outfname _tmp/tigger/sample_db/run-viterbi.csv --n-procs 30 --slurm --workdir $fs/_tmp/$RANDOM --name-column SEQUENCE_ID --seq-column SEQUENCE_INPUT --presto-output
    # cmdfile.write("data(germline_ighv); germlines = germline_ighv\n")
    # cmdfile.write("annotations = read.csv('_tmp/tigger/sample_db/run-viterbi.csv')\n")

    # ----------------------------------------------------------------------------------------
    # my simulation with one gene:
    cmdfile.write("annotations = read.csv('_tmp/tigger/run-viterbi.csv')\n")
    cmdfile.write("germlines = readIgFasta('_tmp/tigger/germlines/ighv-aligned.fasta')\n")
    # cmdfile.write("germlines = readIgFasta('data/imgt/ighv-aligned.fasta')\n")

    cmdfile.write("novel_df = findNovelAlleles(annotations, germlines)\n")   #, nproc=1, germline_min=2)
    cmdfile.write("novel = selectNovel(novel_df)\n")
    cmdfile.write("glimpse(novel)\n")
    cmdfile.close()
    

outdir = '_tmp/tigger'
param_dir = os.getcwd() + '/test/reference-results/test/parameters/simu/hmm'
base_cmd = './bin/partis.py --parameter-dir ' + param_dir

# # existing_genes = 'IGHV1-2*02:IGHD6-19*01:IGHJ5*02'
existing_genes = utils.test_only_genes
# # new_gene = 'IGHV1-2*04'
# new_gene = 'IGHV4-61*08'
# all_genes = existing_genes + ':' + new_gene

# simulate
cmd_str = base_cmd + ' --action simulate'
cmd_str += ' --n-sim-events 500'
cmd_str += ' --n-procs 10'
cmd_str += ' --only-genes ' + existing_genes
cmd_str += ' --outfname ' + outdir + '/simu.csv'
# run(cmd_str)

# utils.write_germline_fasta(outdir + '/germlines', input_dir='data/imgt', only_genes=existing_genes.split(':'))

# cache-parameters
cmd_str = base_cmd + ' --action cache-parameters'
cmd_str += ' --seqfile ' + outdir + '/simu.csv'
cmd_str += ' --n-procs 10'
cmd_str += ' --datadir ' + outdir + '/germlines'
cmd_str += ' --parameter-dir ' + outdir + '/simu'
cmd_str += ' --plotdir ' + os.getenv('www') + '/tmp/test/plots'
run(cmd_str)

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
