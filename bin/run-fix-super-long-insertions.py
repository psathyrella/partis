#!/usr/bin/env python
import sys
import os
partis_dir = os.getcwd()
sys.path.insert(1, partis_dir + '/python')

import utils

label = 'subcluster-annotation-v1'  # NOTE already made -v2
bdir = '%s/partis/fix-super-long-insertions/%s' % ('/fh/local/dralph', label) # os.getenv('fs')
n_trees = 500 #300 #10 #50
n_procs = 25 #10
ntseq = 10000 #3000 #
mmstr = '--mutation-multiplier 2' #7' #1' # 10

dryrun = False #True

# for npg, npgvar in [(50, 98), ]: #[(50, 98), (300, 598)]:
#     outdir = '%s/npg-%d-width-%d' % (bdir, npg, npgvar)
#     common = '--is-simu --infname %s/simu.yaml --parameter-dir %s/parameters --only-csv-plots --n-procs %d' % (outdir, outdir, n_procs)
#     cmd = './bin/bcr-phylo-run.py --base-outdir %s/bcr-phylo --actions simu --n-sim-seqs-per-generation %d --carry-cap 1000 --obs-times 100 --n-sim-events %d --n-procs %d --only-csv-plots' % (outdir, npg, n_trees, n_procs)
#     cmd += ' --parameter-variances n-sim-seqs-per-generation,%d' % npgvar
#     # utils.simplerun(cmd, logfname='%s/bcr-phylo.log'%outdir, dryrun=dryrun)
#     # cmd = 'cat %s/bcr-phylo/selection/simu/event-*/simu.nwk |sed \'s/;/;\\n/g\' >%s/all-trees.nwk' % (outdir, outdir)
#     # utils.simplerun(cmd, shell=True, dryrun=dryrun) #, logfname='%s/tree-cat.log'%outdir)
#     # cmd = './bin/partis simulate --n-sim-events %d --outfname %s/simu.yaml --simulate-from-scratch --mutate-conserved-codons %s --input-simulation-treefname %s/all-trees.nwk --n-procs %d' % (ntseq / npg, outdir, mmstr, outdir, n_procs)
#     # utils.simplerun(cmd, logfname='%s/simulate.log'%outdir, dryrun=dryrun)
#     # cmd = './bin/partis cache-parameters %s --plotdir %s/parameter-plots' % (common, outdir)
#     # utils.simplerun(cmd, logfname='%s/cache-parameters.log'%(outdir), dryrun=dryrun)
#     for use_kmean in [False, True]:
#         for acsize in [None, 3, 5, 10]: #[None, 2, 3, 5, 10, 15]:
#             bstr = 'subc-%s_kmean-%s' % (acsize, use_kmean)
#             cmd = './bin/partis annotate --simultaneous-true-clonal-seqs --plot-annotation-performance %s --outfname %s/%s/annotations.yaml --plotdir %s/%s/plots' % (common, outdir, bstr, outdir, bstr)
#             if acsize is not None:
#                 cmd += ' --subcluster-annotation-size %d' % acsize
#             if use_kmean:
#                 cmd += ' --kmeans-subclusters'
#             utils.simplerun(cmd, logfname='%s/annotate-%s.log'%(outdir, bstr), dryrun=dryrun)
# sys.exit()

npg = 50; npgvar = 98
# npg = 300; npgvar = 598
p_strs = ['subc-%s'%sc for sc in ['None', 10, 5, 3]]
b_strs = ['kmean-False'] #, 'kmean-True']
xtrastr = '-kmean-False'
# pt=mute-freqs/overall; subd=parameter-plots/true
subd = 'plots/hmm'
for pt in ['gene-call', 'boundaries', 'mutation']:
    cmd = './bin/compare-plotdirs.py --outdir %s/cf-plots/npg-%d-width-%d%s/%s --normalize --performance-plots --extra-stats absmean --translegend=-0.6:-0.1 ' % (bdir.replace('/fh/local', '/fh/fast/matsen_e'), npg, npgvar, xtrastr, pt)
    names = ['%s%s%s'%(p, '' if b=='' else '_', b) for p in p_strs for b in b_strs]
    plotdirs = ['%s/npg-%d-width-%d/%s%s%s/%s/%s' % (bdir, npg, npgvar, p, '' if b=='' else '_', b, subd, pt) for p in p_strs for b in b_strs]
    cmd += ' --names %s --plotdirs %s' % (':'.join(names), ':'.join(plotdirs))
    utils.simplerun(cmd)

# # pt=mute-freqs/overall; subd=parameter-plots/true
# subd=corrected/plots/hmm; pt=boundaries #gene-call # #boundaries mutation
# ./bin/compare-plotdirs.py --outdir $bdir/cf-plots/$pt --names 3:7:25:50 --normalize --extra-stats mean --translegend=-0.6:0 \
# 			   --plotdirs $bdir/npg-3/$subd/$pt:$bdir/npg-7/$subd/$pt:$bdir/npg-25/$subd/$pt:$bdir/npg-50/$subd/$pt
