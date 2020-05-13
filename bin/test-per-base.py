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

cfglist = [
    # simcmd = './bin/partis simulate --n-sim-events 100 --n-procs 8 --only-genes IGHV1-2*02:IGHJ4*02 --n-leaves 100 --constant-number-of-leaves --mutation-multiplier 3'
    # ('old', ['--simulate-from-scratch', '--force-dont-generate-germline-set']),
    # ('old-flat', [' --simulate-from-scratch', '--force-dont-generate-germline-set', '--flat-mute-freq']),
    # ('old-non-scratch', ['--parameter-dir', '%s/test/reference-results/test/parameters/simu' % os.getenv('PWD')]),

    # # simcmd = './bin/partis simulate --n-sim-events 100 --n-procs 8 --only-genes IGHV1-2*02:IGHJ4*02 --n-leaves 100 --constant-number-of-leaves --seed 1 --check-tree-depths'
    # # ('newlik-old', ['--parameter-dir', '%s/test/reference-results/test/parameters/simu' % os.getenv('PWD')]),  # check out d1b5cb7 to recapitulate this
    # # ('newlik-v0', ['--per-base-mutation', '--parameter-dir', '%s/test/reference-results/test/parameters/simu' % os.getenv('PWD')]),  # this was wrong -- it had the bug where I had the germline base as 0.9 in the equilibrium freqs, where it should've been near 0
    # ('newlik-rewritten-old', ['--parameter-dir', '%s/test/reference-results/test/parameters/simu' % os.getenv('PWD')]),  # check this is almost the same as 'newlik-old' (this is after switching to one tree/bpp run for the whole sequence, so it should be similar to changing the random seed (and it is, yay)
    # ('newlik-rewritten-v1', ['--per-base-mutation', '--parameter-dir', '%s/test/reference-results/test/parameters/simu' % os.getenv('PWD')]),  # both adding 'rewritten' and incrementing version is kind of confusing, but it beats the alternatives

    # simcmd = './bin/partis simulate --n-sim-events 100 --n-procs 8 --only-genes IGHV1-2*02:IGHJ4*02 --n-leaves 10 --constant-number-of-leaves --seed 1 --mutation-multiplier 5' #  --check-tree-depths
    # ('tmp-test-old', ['--parameter-dir', '%s/test/reference-results/test/parameters/simu' % os.getenv('PWD')]),
    # ('tmp-test-v0', ['--per-base-mutation', '--parameter-dir', '%s/test/reference-results/test/parameters/simu' % os.getenv('PWD')]),
# added these lines to recombinator at l444 and l473, then cached parameters, then simulated with those parameters, then plotted. It looked as expected, per-base is definitely per-basing
# pbcounts = {n : 0.97 if n=='C' else 0.01 for n in utils.nukes}
# freq = 1. if rseq[inuke] != 'A' else 0.001
# ./bin/partis cache-parameters --infname _output/per-base-test/tmp-test-v0/simu.yaml --parameter-dir _output/per-base-test/cache-parameters-from-tmp-test-v0 --only-genes IGHV1-2*02:IGHJ4*02 --n-procs 8 --plotdir _output/per-base-test/cache-parameters-from-tmp-test-v0/plots --make-per-gene-per-base-plots
# ./bin/partis simulate --n-sim-events 100 --n-procs 8 --only-genes IGHV1-2*02:IGHJ4*02 --n-leaves 10 --constant-number-of-leaves --seed 1 --mutation-multiplier 1 --per-base-mutation --parameter-dir _output/per-base-test/cache-parameters-from-tmp-test-v0 --outfname _output/per-base-test/simu-from-cache-parameters-from-tmp-test-v0/simu.yaml
# ./bin/example-parse-output.py --fname _output/per-base-test/simu-from-cache-parameters-from-tmp-test-v0/simu.yaml --plotdir _output/per-base-test/simu-from-cache-parameters-from-tmp-test-v0/plots

    ('qa013-mar-8-old', ['--parameter-dir', '_output/per-base-test/qa013-mar-8/QA013-g-merged']),  # all the stuff using test/reference-results has a big super-high-mutation tail, which isn't a part of parameter space i care a lot about
    ('qa013-mar-8-v0', ['--per-base-mutation', '--parameter-dir', '_output/per-base-test/qa013-mar-8/QA013-g-merged']),
    ('qa013-mar-8-v1', ['--per-base-mutation', '--parameter-dir', '_output/per-base-test/qa013-mar-8/QA013-g-merged']),  # add self.per_base_mutation_multiplier
]
simcmd = './bin/partis simulate --n-sim-events 250 --n-procs 8 --only-genes IGHV1-2*02:IGHJ4*02 --n-leaves 10 --constant-number-of-leaves --seed 1'  #  --check-tree-depths


# for label, xtra_args in cfglist:
#     utils.simplerun('%s %s --outfname %s/%s/simu.yaml' % (simcmd, ' '.join(xtra_args), outdir, label)) #, dryrun=True)
#     utils.simplerun('./bin/example-parse-output.py --fname %s/%s/simu.yaml --plotdir %s/%s/plots' % (outdir, label, outdir, label))
#     pass
# sys.exit()


# subd = 'mute-freqs/per-gene-per-position/v'
subd = 'mute-freqs/overall'

# cfglist.insert(0, ('original', [None]))  # this is from just running example-parse-output plotting on the reference simulation file (most everybody's trying to look like this)
cfglist.insert(0, ('original-qa013-mar-8-sw-cache', [None]))

pdirs = ['%s/%s/plots/%s' % (outdir, l, subd) for l, _ in cfglist]
names = [l for l, _ in cfglist]
cmd = './bin/compare-plotdirs.py --outdir _output/tmp --plotdirs %s --names %s --normalize --translegend=-0.4:0 --extra-stats mean' % (':'.join(pdirs), ':'.join(names))  # 
utils.simplerun(cmd)
