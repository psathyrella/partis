#!/usr/bin/env python
from subprocess import check_call
import sys
sys.path.insert(1, './python')
import os
import csv

from clusterpath import ClusterPath
import utils

def run_cluster(cl):
    cmd = './bin/run-driver.py --label ' + label + ' --action run-viterbi --is-data --datafname VRC01_heavy_chains-dealigned.fasta'  # --simfname ' + os.path.dirname(infname) + '/simu-foo-bar.csv'
    extras = ['--n-sets', len(cl), '--queries', ':'.join(cl), '--debug', 1, '--sw-debug', 0, '--n-procs', 1, '--n-best-events', 1]  # the space is needed in case query name starts with a '-'
    cmd += utils.get_extra_str(extras)
    check_call(cmd.split())

label = 'chaim-test'
infname = '/fh/fast/matsen_e/dralph/work/partis-dev/_output/' + label + '/partitions.csv'

cp = ClusterPath(-1)
cp.readfile(infname)

# print '---> annotations for clusters in best partition:'
# for cluster in cp.partitions[cp.i_best]:
#     run_cluster(cluster)

for ipart in range(cp.i_best, cp.i_best + 10):
    if ipart >= len(cp.partitions):
        break
    print '---> annotation for intersection of partitions %d and %d:' % (ipart - 1, ipart)
    cp.print_partition(ipart - 1)
    cp.print_partition(ipart)

    parents = cp.get_parent_clusters(ipart)
    if parents is None:
        print '   skipping synthetic rewind step'
        continue
    run_cluster(parents[0] + parents[1])
