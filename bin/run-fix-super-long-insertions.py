#!/usr/bin/env python
import sys
import os
partis_dir = os.getcwd()
sys.path.insert(1, partis_dir + '/python')

import utils

label = 'subcluster-annotation-v5'
bdir = '%s/partis/fix-super-long-insertions/%s' % ('/fh/local/dralph', label) # os.getenv('fs')
n_trees = 300 #500 #10 #50
n_procs = 25 #10
ntseq = 3000 #10000 #
mmstr = '--mutation-multiplier 7' #7' #1' # 10

dryrun = False #True

# # for npg, npgvar in [(5, 2), ]: #[(50, 98), (300, 598)]:
# npgvar = None
# for use_bcr_phylo in [False]: #[True, False]:
#     for npg in [50]:
#         # outdir = '%s/npg-%d-width-%d' % (bdir, npg, npgvar)
#         outdir = '%s/npg-%d-bcrp-%s' % (bdir, npg, use_bcr_phylo)
#         common = '--is-simu --infname %s/simu.yaml --parameter-dir %s/parameters --only-csv-plots --n-procs %d' % (outdir, outdir, n_procs)
#         if use_bcr_phylo:
#             cmd = './bin/bcr-phylo-run.py --base-outdir %s/bcr-phylo --actions simu --n-sim-seqs-per-generation %d --carry-cap 1000 --obs-times 100 --n-sim-events %d --n-procs %d --only-csv-plots' % (outdir, npg, n_trees, n_procs)
#             # cmd += ' --parameter-variances n-sim-seqs-per-generation,%d' % npgvar
#             utils.simplerun(cmd, logfname='%s/bcr-phylo.log'%outdir, dryrun=dryrun)
#             cmd = 'cat %s/bcr-phylo/selection/simu/event-*/simu.nwk |sed \'s/;/;\\n/g\' >%s/all-trees.nwk' % (outdir, outdir)
#             utils.simplerun(cmd, shell=True, dryrun=dryrun) #, logfname='%s/tree-cat.log'%outdir)
#         cmd = './bin/partis simulate --n-sim-events %d --outfname %s/simu.yaml --simulate-from-scratch --mutate-conserved-codons %s --n-procs %d' % (ntseq / npg, outdir, mmstr, n_procs)
#         if use_bcr_phylo:
#             cmd += ' --input-simulation-treefname %s/all-trees.nwk' % outdir
#         else:
#             assert npgvar is None
#             cmd += ' --constant-number-of-leaves'
#         utils.simplerun(cmd, logfname='%s/simulate.log'%outdir, dryrun=dryrun)
#         cmd = './bin/partis cache-parameters %s --plotdir %s/parameter-plots' % (common, outdir)
#         utils.simplerun(cmd, logfname='%s/cache-parameters.log'%(outdir), dryrun=dryrun)
#         for acsize in [3, 999]: #[None, 3, 5, 10]:
#             bstr = 'subc-%s' % (acsize)
#             cmd = './bin/partis annotate --simultaneous-true-clonal-seqs --plot-annotation-performance %s --outfname %s/%s/annotations.yaml --plotdir %s/%s/plots' % (common, outdir, bstr, outdir, bstr)
#             if acsize is not None:
#                 cmd += ' --subcluster-annotation-size %d' % acsize
#             utils.simplerun(cmd, logfname='%s/annotate-%s.log'%(outdir, bstr), dryrun=dryrun)
# sys.exit()

# npg = 5; npgvar = 2
npg = 50
p_strs = ['subc-%s'%sc for sc in [999, 3]] #['None', 10, 5, 3]]
b_strs = ['bcrp-False'] #, 'bcrp-True']
xtrastr = '' #'-ceil' #'-' + str(b_strs[0])
# pt=mute-freqs/overall; subd=parameter-plots/true
subd = 'plots/hmm'
for pt in ['gene-call', 'boundaries', 'mutation']:
    npgstr = 'npg-%d%s' % (npg, xtrastr)
    cmd = './bin/compare-plotdirs.py --outdir %s/cf-plots/%s/%s --normalize --performance-plots --extra-stats absmean --translegend=-0.6:-0.1 ' % (bdir.replace('/fh/local', '/fh/fast/matsen_e'), npgstr, pt)
    names = ['%s%s%s'%(p, '' if b=='' else '-', b) for p in p_strs for b in b_strs]
    # plotdirs = ['%s/%s/%s%s%s/%s/%s' % (bdir, npgstr, b, p, '' if b=='' else '-', b, subd, pt) for p in p_strs for b in b_strs]
    plotdirs = ['%s/%s-%s/%s/%s/%s' % (bdir, npgstr, b, p, subd, pt) for p in p_strs for b in b_strs]
    cmd += ' --names %s --plotdirs %s' % (':'.join(names), ':'.join(plotdirs))
    utils.simplerun(cmd)

# # pt=mute-freqs/overall; subd=parameter-plots/true
# subd=corrected/plots/hmm; pt=boundaries #gene-call # #boundaries mutation
# ./bin/compare-plotdirs.py --outdir $bdir/cf-plots/$pt --names 3:7:25:50 --normalize --extra-stats mean --translegend=-0.6:0 \
# 			   --plotdirs $bdir/npg-3/$subd/$pt:$bdir/npg-7/$subd/$pt:$bdir/npg-25/$subd/$pt:$bdir/npg-50/$subd/$pt
