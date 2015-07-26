#!/usr/bin/env python
from subprocess import check_call

fsdir = '/fh/fast/matsen_e/dralph/work/changeo'
label = '3_leaves_1_mutate'

def run(cmdstr):
    print 'RUN %s' % cmdstr
    check_call(cmdstr.split())

cmd = './MakeDb.py imgt -i ' + fsdir + '/' + label + ' -s ' + fsdir + '/simulation/head-simu-' + label.replace('_', '-') + '.fasta'
# apparently deprecated: ./SplitDb.py group -d IB_db-pass.tab -f FUNCTIONAL, instead use:
cmd = './ParseDb.py select -d ' + fsdir + '/' + label + '_db-pass.tab -f FUNCTIONAL -u T'
# Rscript tigger.R
cmd = './DefineClones.py bygroup -d ' + fsdir + '/' + label + '_db-pass_parse-select.tab --act first --model m1n --dist 7'
# ./CreateGermlines.py -d IB_db-pass_tigger-pass_clone-pass.tab --vfield V_CALL_GENOTYPED --cloned
run(cmd)
