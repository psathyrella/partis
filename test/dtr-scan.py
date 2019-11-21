#!/usr/bin/env python
import sys
import colored_traceback.always
import os
import yaml
import argparse

sys.path.insert(1, './python')
import utils
import treeutils

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('action', choices=['train', 'test'])
parser.add_argument('--label', default='test-dtr-v0')
args = parser.parse_args()

vsets = {
    'all' : {
        'within-families' : ['lbi', 'cons-dist', 'edge-dist', 'lbr', 'shm'],
        'among-families' : ['lbi', 'cons-dist', 'edge-dist', 'lbr', 'shm', 'fay-wu-h', 'cons-seq-shm', 'mean-shm', 'max-lbi', 'max-lbr']},
    # 'med' : {
    #     'within-families' : ['lbi', 'cons-dist', 'shm'],
    #     'among-families' : ['lbi', 'cons-dist', 'shm', 'fay-wu-h', 'max-lbi', 'max-lbr']},
    # 'min' : {
    #     'within-families' : ['lbi', 'cons-dist'],
    #     'among-families' : ['lbi', 'cons-dist', 'mean-shm', 'max-lbi']}
}

baseworkdir = '%s/_tmp' % os.getcwd()

basecmds = {
    'test-dtr-v0' : './test/cf-tree-metrics.py --label test-dtr-v0 --n-replicates 2 --n-sim-events-per-proc 5 --carry-cap-list 750 --obs-times-list 75 --n-sim-seqs-per-gen-list 80 --lb-tau-list 0.0025',
    'choose-among-families-v3' : './test/cf-tree-metrics.py --label choose-among-families-v3 --n-replicates 10 --n-sim-events-per-proc 30 --slurm --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --lb-tau-list 0.0025 --dont-observe-common-ancestors',
    'dtr-train-v0' : './test/cf-tree-metrics.py --label dtr-train-v0 --n-replicates 5 --n-sim-events-per-proc 1000 --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --selection-strength 0.75 --lb-tau-list 0.0025 --parameter-variances carry-cap,2000:obs-times,150:n-sim-seqs-per-generation,200:selection-strength,0.5',
    #  --slurm
}
basecmd = basecmds[args.label]

basecmd += ' --actions get-tree-metrics --metric-method dtr'
basepath = '/fh/fast/matsen_e/dralph/partis/tree-metrics/%s' % args.label
training_seed = 0  # just for output file names, I don't really want to keep track of here, but utils.run_cmds() requires it

cmdfos = []
for ensemble in ['grad-boost']: #, 'ada-boost', 'forest']: #, 'bag']:
    for n_estimators in [2]: #[30, 100]:
        for max_depth in [3]: #, 6]:
            for n_train_per_family in [1]: #, 3]:
                for vsname, varset in vsets.items():
                    cmd = basecmd
                    paramstr = 'ensemble_%s_n-estimators_%d_vars_%s_max-depth_%d_n-train-per-family_%d' % (ensemble, n_estimators, vsname, max_depth, n_train_per_family)
                    xtrastrs = {s : '%s_%s' % (s, paramstr) for s in ['train', 'test']}
                    modeldir = '%s/seed-%d/dtr/%s-dtr-models' % (basepath, training_seed, xtrastrs['train'])
                    workdir = '%s/%s' % (baseworkdir, paramstr)
                    cfgfname = '%s/cfg.yaml' % workdir
                    cmd += ' --dtr-cfg %s' % cfgfname
                    if args.action == 'train':
                        if not os.path.exists(workdir):
                            os.makedirs(workdir)
                        with open(cfgfname, 'w') as cfile:
                            yaml.dump({'ensemble' : ensemble, 'n_estimators' : n_estimators, 'vars' : varset, 'max_depth' : max_depth, 'n_train_per_family' : n_train_per_family}, cfile, width=200)
                        cmd += ' --iseed %d --extra-plotstr %s' % (training_seed, xtrastrs['train'])
                        modelfnames = [treeutils.dtrfname(modeldir, cg) for cg in treeutils.cgroups]
                        outfname = modelfnames[0]
                        # logdir = os.path.dirname(outfname)  # NOTE has to be the '%s-dtr-models' dir (or at least can't be the '%s-plots' dir, since cf-tree-metrics uses the latter, and the /{out,err} files overwrite each other
                        logdir = '%s/seed-%d/dtr/%s-plots/dtr-scan-logs' % (basepath, training_seed, xtrastrs['train'])
                    elif args.action == 'test':
                        cmd += ' --dtr-path %s --extra-plotstr %s' % (modeldir, xtrastrs['test'])
                        log_seed = 0  # put our log file here, and also (only) check for existing output in this seed NOTE wait I'm not sure this gets used for checking for existing output, maybe it's only used to see if the job finished successfully
                        outfname = ['%s/seed-%d/dtr/%s-plots/true-tree-metrics/%s-dtr/%s-dtr-vs-affinity-ptiles/%s-dtr-vs-affinity-true-tree-ptiles-all-clusters.yaml' % (basepath, log_seed, xtrastrs['test'], cg, cg, cg) for cg in treeutils.cgroups][0]
                        logdir = '%s/seed-%d/dtr/%s-plots/dtr-scan-logs' % (basepath, log_seed, xtrastrs['test'])
                    else:
                        assert False

                    print '    %s' % logdir
                    if not os.path.exists(logdir):
                        os.makedirs(logdir)
                    print cmd

                    # utils.simplerun(cmd, debug=True) #, dryrun=True)
                    cmdfo = {'cmd_str' : cmd,
                             'outfname' : outfname,
                             'logdir' : logdir,
                             'workdir' : workdir,
                    }  # I don't think I actually use the work dir
                    cmdfos.append(cmdfo)

print '  starting %d jobs' % len(cmdfos)
n_max_procs = 25  # utils.auto_n_procs()
utils.run_cmds(cmdfos, n_max_procs=n_max_procs, proc_limit_str='test/cf-tree-metrics', debug='write:cf-tree-metrics.log')
