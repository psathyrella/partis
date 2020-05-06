#!/usr/bin/env python
import csv
import os
import sys
import argparse
import colored_traceback.always

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
# from clusterpath import ClusterPath

# parser = argparse.ArgumentParser()
# parser.add_argument('--fname', default=partis_dir + '/test/reference-results/partition-new-simu.yaml')
# parser.add_argument('--glfo-dir', default=partis_dir + '/data/germlines/human')
# parser.add_argument('--locus', default='igh')
# parser.add_argument('--plotdir', help='if set, plot annotation parameters from --fname to --plotdir and exit')
# args = parser.parse_args()

outdir = '_output/per-base-test'

# simcmd = './bin/partis simulate --n-sim-events 100 --n-procs 8 --only-genes IGHV1-2*02:IGHJ4*02 --n-leaves 100 --constant-number-of-leaves --mutation-multiplier 3'
# cfglist = [
#     ('old', ['--simulate-from-scratch', '--force-dont-generate-germline-set']),
#     ('old-flat', [' --simulate-from-scratch', '--force-dont-generate-germline-set', '--flat-mute-freq']),
#     ('old-non-scratch', ['--parameter-dir', 'test/reference-results/test/parameters/simu']),
# ]
# simcmd = './bin/partis simulate --n-sim-events 100 --n-procs 8 --only-genes IGHV1-2*02:IGHJ4*02 --n-leaves 100 --constant-number-of-leaves --seed 1'
simcmd = './bin/partis simulate --n-sim-events 1 --only-genes IGHV1-2*02:IGHJ4*02 --n-leaves 5 --debug 1 --constant-number-of-leaves --seed 4'
cfglist = [
    ('newlik-old', ['--parameter-dir', 'test/reference-results/test/parameters/simu']),
    ('newlik-v0', ['--per-base-mutation', '--parameter-dir', 'test/reference-results/test/parameters/simu']),
]

for label, xtra_args in cfglist:
    utils.simplerun('%s %s --outfname %s/%s/simu.yaml' % (simcmd, ' '.join(xtra_args), outdir, label)) #, dryrun=True)
    # utils.simplerun('./bin/example-parse-output.py --fname %s/%s/simu.yaml --plotdir %s/%s/plots' % (outdir, label, outdir, label))
    pass
sys.exit()


subd = 'mute-freqs/per-gene-per-position/v'
# subd = 'mute-freqs/overall'
pdirs = ':'.join('%s/%s/plots/%s' % (outdir, l, subd) for l, _ in cfglist)
cmd = './bin/compare-plotdirs.py --outdir _output/tmp --plotdirs %s --names %s --normalize' % (pdirs, ':'.join(l for l, _ in cfglist))
utils.simplerun(cmd)
