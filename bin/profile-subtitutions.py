#!/usr/bin/env python3
import sys
import os
import argparse
import colored_traceback.always

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils

parser = argparse.ArgumentParser()
parser.add_argument('--fname', required=True, help='fasta file with input sequences')
parser.add_argument('--outdir', required=True, help='output directory')
parser.add_argument('--only-genes')
# parser.add_argument('--locus', default='igh')
# parser.add_argument('--plotdir', help='if set, plot annotation parameters from --fname to --plotdir and exit')
args = parser.parse_args()

args.only_genes = utils.get_arg_list(args.only_genes)

annotationfname = '%s/annotations.yaml' % args.outdir

cmd = './bin/partis cache-parameters --infname %s --outfname %s' % (args.fname, annotationfname)
if args.only_genes is not None:
    cmd += ' --only-genes %s' % ':'.join(args.only_genes)
utils.simplerun(cmd)

annotations = utils.read_output(annotationfname)
for line in annotations

# this is all for qa255-synth v18
# except 157-Vk no-indels

# seqfos = utils.read_fastx(args.fname)
