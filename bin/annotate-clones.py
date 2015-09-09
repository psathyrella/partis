#!/usr/bin/env python
from subprocess import check_call, Popen
import sys
sys.path.insert(1, './python')
import argparse
import os
import csv

from clusterpath import ClusterPath
import baseutils
parser = argparse.ArgumentParser()
parser.add_argument('--infname')
parser.add_argument('--label', default='A')  #021-019')  # chaim-test
parser.add_argument('--istartstop', default='0-10000')
parser.add_argument('--method', default='-naive-hamming')
args = parser.parse_args()

datafname = '/fh/fast/matsen_e/dralph/work/partis-dev/_output/' + args.label + '/istartstop-' + args.istartstop + '/data.csv'
args.infname = datafname.replace('.csv', args.method + '-partition.csv')

def run_cluster(cl, iclust=None):
    cmd = './bin/run-driver.py --label ' + args.label + ' --action run-viterbi --is-data --datafname ' + datafname #VRC01_heavy_chains-dealigned.fasta'  # --simfname ' + os.path.dirname(infname) + '/simu-foo-bar.csv'
    extras = ['--n-sets', len(cl), '--queries', ':'.join(cl), '--debug', 1, '--sw-debug', 0, '--n-procs', 1, '--n-best-events', 1]
    cmd += baseutils.get_extra_str(extras)
    print cmd
    check_call(cmd.split())
    # Popen(cmd + ' >_tmp/' + str(iclust) + '.csv', shell=True)

cp = ClusterPath(-1)
cp.readfile(args.infname)

print '---> annotations for clusters in best partition:'
iclust = 0
for cluster in sorted(cp.partitions[cp.i_best], key=len, reverse=True):
    print len(cluster)
    run_cluster(cluster, iclust)
    # sys.exit()
    if iclust > 3:
        break
    iclust += 1

sys.exit()

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
