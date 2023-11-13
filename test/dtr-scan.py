#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import sys
import json
import colored_traceback.always
import os
import yaml
import argparse
import numpy
import operator
import math
from io import open

sys.path.insert(1, './python')
import utils
import treeutils
import lbplotting

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('action', choices=['train', 'test', 'plot'])
parser.add_argument('--label', default='test-dtr-v0')
parser.add_argument('--training-label', help='if set, use the dtr that was trained on this label (wheras --label is used for the sample we test on)')
parser.add_argument('--base-outdir', default='/fh/fast/matsen_e/%s/partis/tree-metrics'%os.getenv('USER'))
parser.add_argument('--training-seed', type=int, default=0)
parser.add_argument('--plot-seed', type=int)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--n-max-queries', type=int)
parser.add_argument('--n-max-procs', type=int, default=utils.auto_n_procs())
parser.add_argument('--cgroup')
parser.add_argument('--tvar')
args = parser.parse_args()

if args.training_label is not None:  # only makes sense for --training-label to be set if we're testing or plotting
    assert args.action in ['test', 'plot']
    if args.training_label == args.label:
        # print '  setting --training-label to None since it\'s the same as --label'
        args.training_label = None
cglist = treeutils.cgroups if args.cgroup is None else [args.cgroup]
def tvlist(cg):
    return treeutils.dtr_targets[cg] if args.tvar is None else [args.tvar]

# vsets = {  # NOTE if you change the allowed vars in treeutils, these will of course not reflect that
#     'all' : {
#         'within-families' : ['lbi', 'cons-dist-nuc', 'cons-dist-aa', 'edge-dist', 'lbr', 'shm', 'shm-aa'],
#         'among-families' : ['lbi', 'cons-dist-nuc', 'cons-dist-aa', 'edge-dist', 'lbr', 'shm', 'shm-aa', 'fay-wu-h', 'cons-seq-shm-nuc', 'cons-seq-shm-aa', 'mean-shm', 'max-lbi', 'max-lbr']},
#     # 'med' : {
#     #     'within-families' : ['lbi', 'cons-dist-nuc', 'cons-dist-aa', 'shm'],
#     #     'among-families' : ['lbi', 'cons-dist-nuc', 'cons-dist-aa', 'shm', 'fay-wu-h', 'max-lbi', 'max-lbr']},
#     # 'min' : {
#     #     'within-families' : ['lbi', 'cons-dist-nuc', 'cons-dist-aa'],
#     #     'among-families' : ['lbi', 'cons-dist-nuc', 'cons-dist-aa', 'mean-shm', 'max-lbi']}
# }

baseworkdir = '%s/_tmp' % os.getcwd()

basecmds = {
    'test-dtr-v0' : './test/cf-tree-metrics.py --label test-dtr-v0 --n-replicates 2 --n-sim-events-per-proc 5 --carry-cap-list 750 --obs-times-list 75 --n-sim-seqs-per-gen-list 80 --lb-tau-list 0.0025',
    'choose-among-families-v3' : './test/cf-tree-metrics.py --label choose-among-families-v3 --n-replicates 10 --n-sim-events-per-proc 30 --slurm --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --lb-tau-list 0.0025 --dont-observe-common-ancestors',
    'dtr-train-v0' : './test/cf-tree-metrics.py --label dtr-train-v0 --n-replicates 5 --n-sim-events-per-proc 1000   --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --selection-strength 0.75 --lb-tau-list 0.0025 --parameter-variances carry-cap,2000:obs-times,150:n-sim-seqs-per-generation,200:selection-strength,0.5',
    'dtr-train-v1' : './test/cf-tree-metrics.py --label dtr-train-v1 --n-replicates 4 --n-sim-events-per-proc 50000  --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 30  --selection-strength 0.75 --lb-tau-list 0.0025 --parameter-variances carry-cap,2000:obs-times,150:n-sim-seqs-per-generation,15:selection-strength,0.5',
    'dtr-train-v2' : './test/cf-tree-metrics.py --label dtr-train-v2 --n-replicates 2 --n-sim-events-per-proc 300000 --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 20  --selection-strength 0.75 --lb-tau-list 0.0025 --parameter-variances carry-cap,2000:obs-times,150:n-sim-seqs-per-generation,15:selection-strength,0.5',
    # NOTE not actually two replicates, but want to get the dir structure the same as the others UPDATE made a second replicate for testing, with only 1000 families
    'dtr-train-v3' : './test/cf-tree-metrics.py --label dtr-train-v3 --n-replicates 2 --n-sim-events-per-proc 50000  --carry-cap-list=-1   --obs-times-list=-1  --n-sim-seqs-per-gen-list=-1  --selection-strength=-1.  --lb-tau-list 0.0025 --parameter-variances carry-cap,250..500..900..1000..1100..1500..5000:obs-times,75..100..150..200..1000:n-sim-seqs-per-generation,15..30..75..150..500:selection-strength,0.5..0.9..0.95..1.0',

    #  --slurm
    # --n-sub-procs 30
}
basecmd = basecmds[args.label]

basecmd += ' --actions get-tree-metrics --metric-method dtr'
log_seed = 0  # for the latter, put our log file here, and also (only) check for existing output in this seed NOTE wait I'm not sure this gets used for checking for existing output, maybe it's only used to see if the job finished successfully

cmdfos, plotfo = [], []
for n_estimators in [30, 100]:
    for max_depth in [5, 10]:  # tried max depth 50, but gave up waiting for training to finish (waited 18hr or so on 2.4m seqs/50000 families)
        cfgfo = {'n_estimators' : n_estimators, 'max_depth' : max_depth}
        cmd = basecmd
        paramstr = 'n-estimators_%d_max-depth_%d' % (n_estimators, max_depth)
        xtrastrs = {s : '%s_%s' % (s, paramstr) for s in ['train', 'test']}
        trainstr = '_trained-on-%s' % args.training_label if args.training_label is not None else ''
        modeldir = '%s/%s/seed-%d/dtr/%s-dtr-models' % (args.base_outdir, args.label if args.training_label is None else args.training_label, args.training_seed, xtrastrs['train'])
        workdir = '%s/%s' % (baseworkdir, paramstr)
        cfgfname = '%s/cfg.yaml' % workdir
        cmd += ' --dtr-cfg %s' % cfgfname
        cmd += ' --base-outdir %s' % args.base_outdir
        if args.overwrite:
            cmd += ' --overwrite'
        if args.action == 'train':
            assert args.training_label is None
            if not os.path.exists(workdir):
                os.makedirs(workdir)
            with open(cfgfname, 'w') as cfile:
                yaml.dump(cfgfo, cfile, width=200)
            cmd += ' --train-dtr --iseed %d --extra-plotstr %s' % (args.training_seed, xtrastrs['train'])
            modelfnames = [treeutils.dtrfname(modeldir, cg, tv) for cg in cglist for tv in tvlist(cg)]
            outfname = modelfnames[-1]
            # logdir = os.path.dirname(outfname)  # NOTE has to be the '%s-dtr-models' dir (or at least can't be the '%s-plots' dir, since cf-tree-metrics uses the latter, and the /{out,err} files overwrite each other
            logdir = '%s/%s/seed-%d/dtr/%s-plots/dtr-scan-logs' % (args.base_outdir, args.label if args.training_label is None else args.training_label, args.training_seed, xtrastrs['train'])
        elif args.action == 'test':
            cmd += ' --dtr-path %s --extra-plotstr %s%s' % (modeldir, xtrastrs['test'], trainstr)
            plotdir = '%s/%s/seed-%d/dtr/%s%s-plots' % (args.base_outdir, args.label, log_seed, xtrastrs['test'], trainstr)
            allfns = [treeutils.tmfname(plotdir, 'dtr', lbplotting.getptvar(tv), cg=cg, tv=tv) for cg in cglist for tv in tvlist(cg)]  # i think this should really include all seeds, but the only problem it'd cause is output would be missing, which I'd immediately notice, so whatever
            outfname = allfns[-1]
            logdir = '%s/%s/seed-%d/dtr/%s%s-plots/dtr-scan-logs' % (args.base_outdir, args.label, log_seed, xtrastrs['test'], trainstr)
        elif args.action == 'plot':
            def pdfcn(tmpseed): return '%s/%s/seed-%d/dtr/%s%s-plots' % (args.base_outdir, args.label, tmpseed, xtrastrs['test'], trainstr)
            seed_iter = [args.plot_seed] if args.plot_seed is not None else list(range(int(utils.get_val_from_arglist(cmd.split(), '--n-replicates'))))
            plotfo.append({'cfg' : cfgfo, 'plotdirs' : [(s, pdfcn(s)) for s in seed_iter if os.path.exists(pdfcn(s))]})
            continue
        else:
            assert False

        if args.n_max_queries is not None:
            cmd += ' --n-max-queries %d' % args.n_max_queries

        print('    %s' % logdir)
        print(cmd)
        if not os.path.exists(logdir):
            os.makedirs(logdir)

        # utils.simplerun(cmd, debug=True) #, dryrun=True)
        cmdfo = {'cmd_str' : cmd,
                 'outfname' : outfname,
                 'logdir' : logdir,
                 'workdir' : workdir,
        }  # I don't think I actually use the work dir
        cmdfos.append(cmdfo)

if args.action == 'plot':
    print_all_seeds = False
    # print '%s --training-seed has no effect when we\'re plotting, i.e. we aren\'t checking what seed was used for training here' % utils.color('yellow', 'note')
    # ----------------------------------------------------------------------------------------
    def getptvals(yfname, choice_grouping, target_var, per_x='per-seq', iclust=None, min_ptile_to_plot=75.):  # NOTE duplicates code in cf-tree-metrics.py
        with open(yfname) as yfile:
            yamlfo = json.load(yfile)  # too slow with yaml
        ytmpfo = yamlfo['percentiles'][per_x]
        if per_x == 'per-seq':
            assert iclust is None
            ytmpfo = ytmpfo['all-clusters' if choice_grouping == 'among-families' else 'within-families-mean']
        else:
            ytmpfo = ytmpfo[choice_grouping]
        return [abs(pafp - afp) for lbp, afp, pafp in zip(ytmpfo['lb_ptiles'], ytmpfo['mean_%s_ptiles' % lbplotting.getptvar(target_var)], ytmpfo['perfect_vals']) if lbp > min_ptile_to_plot]

    # ----------------------------------------------------------------------------------------
    for cg in cglist:
        for tv in tvlist(cg):
            print('  %s %s   train: %s   test: %s' % (cg, tv, args.training_label if args.training_label is not None else args.label, args.label))
            print('    %s    diff to perfect   seeds' % ' '.join(tuple('%20s'%k for k in sorted(plotfo[0]['cfg']))))
            for pfo in plotfo:
                paramstr = '    '.join(tuple('%15d'%v for k, v in sorted(list(pfo['cfg'].items()), key=operator.itemgetter(0))))
                fstr = '1' if tv == 'affinity' else '2'
                estr = str(int(fstr)+1)
                diff_val_list, seed_list = [], []  # for each seed
                for ipd, (pseed, pdir) in enumerate(pfo['plotdirs']):
                    is_train_seed = args.training_label is None and pseed == args.training_seed  # even if it was set on command line, we set it to None if it's the same as --label
                    yfn = treeutils.tmfname(pdir, 'dtr', lbplotting.getptvar(tv), cg=cg, tv=tv)
                    if not os.path.exists(yfn):
                        print('  %s missing %s' % (utils.color('yellow', 'warning'), yfn))
                        continue
                    diff_vals = getptvals(yfn, cg, tv)  # for each percentile
                    diff_to_perfect = numpy.mean(diff_vals)
                    if not is_train_seed:
                        diff_val_list.append(diff_to_perfect)
                        seed_list.append(pseed)
                    if print_all_seeds:
                        print(('     %s        %3d   %6.'+fstr+'f  %s') % (paramstr if ipd == 0 else len(paramstr)*' ', pseed, diff_to_perfect, utils.color('red', 'training') if is_train_seed else ''))
                if not print_all_seeds:
                    print(('     %s        %6.'+fstr+'f +/-%-6.'+estr+'f     %s') % (paramstr, numpy.mean(diff_val_list), 0. if len(diff_val_list) < 2 else numpy.std(diff_val_list, ddof=1) / math.sqrt(len(diff_val_list)), ' '.join(str(s) for s in seed_list)))
else:
    print('  starting %d jobs' % len(cmdfos))
    utils.run_cmds(cmdfos, n_max_procs=args.n_max_procs, debug='write:cf-tree-metrics.log')
