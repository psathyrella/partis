#!/usr/bin/env python
import sys
import os
import yaml

sys.path.insert(1, './python')
import utils
import treeutils

vsets = {
    'all' : {
        'within-families' : ['lbi', 'cons-dist', 'edge-dist', 'lbr', 'shm'],
        'among-families' : ['lbi', 'cons-dist', 'edge-dist', 'lbr', 'shm', 'fay-wu-h', 'cons-seq-shm', 'mean-shm', 'max-lbi', 'max-lbr']},
    'med' : {
        'within-families' : ['lbi', 'cons-dist', 'shm'],
        'among-families' : ['lbi', 'cons-dist', 'shm', 'fay-wu-h', 'max-lbi', 'max-lbr']},
    'min' : {
        'within-families' : ['lbi', 'cons-dist'],
        'among-families' : ['lbi', 'cons-dist', 'mean-shm', 'max-lbi']}
}

train = True
predict = False

label = 'test-dtr-v0'
baseworkdir = '%s/_tmp' % os.getcwd()
basecmd = './test/cf-tree-metrics.py --label %s --n-replicates 2 --n-sim-events-per-proc 5 --carry-cap-list 750 --obs-times-list 75 --n-sim-seqs-per-gen-list 80 --lb-tau-list 0.0025 --actions get-tree-metrics --metric-method dtr' % label
basepath = '/fh/fast/matsen_e/dralph/partis/tree-metrics/%s' % label

cmdfos = []
for ensemble in treeutils.dtr_cfg_options['ensemble']:
    for n_estimators in [10, 30, 100, 500]:
        for vsname, varset in vsets.items():
            cmd = basecmd
            xtra_label = 'ensemble_%s_n-estimators_%d_vars_%s' % (ensemble, n_estimators, vsname)

            if train:
                workdir = '%s/%s' % (baseworkdir, xtra_label)
                if not os.path.exists(workdir):
                    os.makedirs(workdir)
                cfgfname = '%s/cfg.yaml' % workdir
                with open(cfgfname, 'w') as cfile:
                    yaml.dump({'ensemble' : ensemble, 'n_estimators' : n_estimators, 'vars' : varset}, cfile, width=200)
                cmd += ' --dtr-cfg %s --extra-plotstr train_%s' % (cfgfname, xtra_label)
                outfname = 'x'

            # if predict:
            #     x = '--dtr-path %s/seed-%d/dtr/dtr-models' % (basepath, iseed)

            cmdfo = {'cmd_str' : cmd,
                     'outfname' : outfname,
                     'workdir' : workdir}  # I don't think I actually need the work dir
            cmdfos.append(cmdfo)

print '  staring %d jobs' % len(cmdfos)
utils.run_cmds(cmdfos, n_max_procs=utils.auto_n_procs(), proc_limit_str='test/cf-tree-metrics')
