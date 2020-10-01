#!/usr/bin/env python
import sys
import os
partis_dir = os.getcwd()
sys.path.insert(1, partis_dir + '/python')

import utils

label = 'cleanup-check-old'
bdir = '%s/partis/fix-super-long-insertions/%s' % (os.getenv('fs'), label)
n_trees = 10 #300 #50
n_procs = 10 #20
ntseq = 1000 #10000
mmstr = '--mutation-multiplier 3' #1' # 10

# for npg in [10]: #[3, 7]: #[25]: #3 7 25 50; do
#     outdir = '%s/npg-%d' % (bdir, npg)
#     cmd = './bin/bcr-phylo-run.py --base-outdir %s/bcr-phylo --actions simu --n-sim-seqs-per-generation %d --carry-cap 500 --obs-times 100 --n-sim-events %d --n-procs %d --only-csv-plots' % (outdir, npg, n_trees, n_procs)
#     # utils.simplerun(cmd, logfname='%s/bcr-phylo.log'%outdir)
#     cmd = 'cat %s/bcr-phylo/selection/simu/event-*/simu.nwk |sed \'s/;/;\\n/g\' >%s/all-trees.nwk' % (outdir, outdir)
#     # utils.simplerun(cmd, shell=True) #, logfname='%s/tree-cat.log'%outdir)
#     cmd = './bin/partis simulate --n-sim-events %d --outfname %s/simu.yaml --simulate-from-scratch --mutate-conserved-codons %s --input-simulation-treefname %s/all-trees.nwk --n-procs %d' % (ntseq / npg, outdir, mmstr, outdir, n_procs)
#     # utils.simplerun(cmd, logfname='%s/simulate.log'%outdir)
#     for erode_codons in [True, False]:
#         pstr = 'erode-codons' if erode_codons else 'preserve-codons'
# 	common = ' --is-simu --infname %s/simu.yaml --parameter-dir %s/%s/parameters --only-csv-plots --n-procs %d' % (outdir, outdir, pstr, n_procs)
#         cmd = './bin/partis cache-parameters %s --plotdir %s/%s/parameter-plots' % (common, outdir, pstr)
#         # if not erode_codons:  # old version
#         #     cmd += ' --dont-erode-conserved-codons'
#         if erode_codons:
#             cmd += ' --allow-conserved-codon-deletion'
# 	utils.simplerun(cmd, logfname='%s/%s/cache-parameters.log'%(outdir, pstr)) #, dryrun=True)
#         for bcf in [None, 1]: #[None, 2]:
#             bstr = 'old-boundaries' if bcf is None else 'correct-boundaries-%.1f' % bcf
#             cmd = './bin/partis annotate --simultaneous-true-clonal-seqs --plot-annotation-performance %s --outfname %s/%s/%s/annotations.yaml --plotdir %s/%s/%s/plots' % (common, outdir, pstr, bstr, outdir, pstr, bstr)
#             if bcf is not None:
#                 # cmd += ' --boundary-correction-factor %.1f' % bcf  # old version
#                 assert bcf == 1  # would have to go back and rejigger bin/partis and partitiondriver.py
#             else:
#                 cmd += ' --dont-correct-multi-hmm-boundaries'
# 	    utils.simplerun(cmd, logfname='%s/%s/annotate-%s.log'%(outdir, pstr, bstr)) #, dryrun=True)

# sys.exit()

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
