#!/usr/bin/env python
import sys
import os
partis_dir = os.getcwd()
sys.path.insert(1, partis_dir + '/python')

import utils

label = 'subcluster-annotation-v0'
bdir = '%s/partis/fix-super-long-insertions/%s' % ('/fh/local/dralph', label) # os.getenv('fs')
n_trees = 300 #10 #50
n_procs = 20 #10
ntseq = 3000 #10000
mmstr = '--mutation-multiplier 3' #1' # 10

dryrun = False #True

for npg in [3, 25, 150]: #[3, 7]: #[25]: #3 7 25 50; do
    outdir = '%s/npg-%d' % (bdir, npg)
    cmd = './bin/bcr-phylo-run.py --base-outdir %s/bcr-phylo --actions simu --n-sim-seqs-per-generation %d --carry-cap 500 --obs-times 100 --n-sim-events %d --n-procs %d --only-csv-plots' % (outdir, npg, n_trees, n_procs)
    utils.simplerun(cmd, logfname='%s/bcr-phylo.log'%outdir, dryrun=dryrun)
    cmd = 'cat %s/bcr-phylo/selection/simu/event-*/simu.nwk |sed \'s/;/;\\n/g\' >%s/all-trees.nwk' % (outdir, outdir)
    utils.simplerun(cmd, shell=True, dryrun=dryrun) #, logfname='%s/tree-cat.log'%outdir)
    cmd = './bin/partis simulate --n-sim-events %d --outfname %s/simu.yaml --simulate-from-scratch --mutate-conserved-codons %s --input-simulation-treefname %s/all-trees.nwk --n-procs %d' % (ntseq / npg, outdir, mmstr, outdir, n_procs)
    utils.simplerun(cmd, logfname='%s/simulate.log'%outdir, dryrun=dryrun)
    common = ' --is-simu --infname %s/simu.yaml --parameter-dir %s/parameters --only-csv-plots --n-procs %d' % (outdir, outdir, n_procs)
    cmd = './bin/partis cache-parameters %s --plotdir %s/parameter-plots' % (common, outdir)
    utils.simplerun(cmd, logfname='%s/cache-parameters.log'%(outdir), dryrun=dryrun)
    for bcorr, acsize in [(True, 999999), (False, 999999), (False, 5), (False, 10), (False, 15)]:
        bstr = 'bcorr-%s-subc-%d' % (bcorr, acsize)
        cmd = './bin/partis annotate --simultaneous-true-clonal-seqs --plot-annotation-performance %s --outfname %s/%s/annotations.yaml --plotdir %s/%s/plots' % (common, outdir, bstr, outdir, bstr)
        if not bcorr:
            cmd += ' --dont-correct-multi-hmm-boundaries'
        cmd += ' --subcluster-annotation-size %d' % acsize
        utils.simplerun(cmd, logfname='%s/annotate-%s.log'%(outdir, bstr), dryrun=dryrun)
    break

sys.exit()

npg = 10
p_strs = ['erode-codons', 'preserve-codons']
b_strs = ['old-boundaries', 'correct-boundaries-%.1f'%1]
# p_strs = ['preserve-codons']
# b_strs = ['correct-boundaries-%.1f'%bcf for bcf in [0.5, 1, 2]]
# pt=mute-freqs/overall; subd=parameter-plots/true
subd = 'plots/hmm'
for pt in ['gene-call', 'boundaries', 'mutation']:
    cmd = './bin/compare-plotdirs.py --outdir %s/cf-plots/npg-%d/%s --normalize --extra-stats mean --translegend=-0.6:-0.1 ' % (bdir, npg, pt)
    names = ['%s_%s'%(p.split('-')[0], b.split('-')[-1]) for p in p_strs for b in b_strs]
    plotdirs = ['%s/npg-%d/%s/%s/%s/%s' % (bdir, npg, p, b, subd, pt) for p in p_strs for b in b_strs]
    cmd += ' --names %s --plotdirs %s' % (':'.join(names), ':'.join(plotdirs))
    utils.simplerun(cmd)

# # pt=mute-freqs/overall; subd=parameter-plots/true
# subd=corrected/plots/hmm; pt=boundaries #gene-call # #boundaries mutation
# ./bin/compare-plotdirs.py --outdir $bdir/cf-plots/$pt --names 3:7:25:50 --normalize --extra-stats mean --translegend=-0.6:0 \
# 			   --plotdirs $bdir/npg-3/$subd/$pt:$bdir/npg-7/$subd/$pt:$bdir/npg-25/$subd/$pt:$bdir/npg-50/$subd/$pt
