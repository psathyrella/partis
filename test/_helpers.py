import os
import subprocess
import sys

import yaml

REPO_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
TEST_DIR = os.path.join(REPO_DIR, 'test')

# import partis from this checkout, not from whatever is installed
sys.path.insert(0, REPO_DIR)

# tracked paired simulation and the parameters inferred on it
PAIRED_SIMU_DIR = os.path.join(TEST_DIR, 'paired', 'ref-results', 'test', 'simu')
PAIRED_PARAM_DIR = os.path.join(TEST_DIR, 'paired', 'ref-results', 'test', 'parameters', 'simu')
LOCI = ['igh', 'igk', 'igl']

BACKEND = os.environ.get('PARTIS_TEST_BACKEND', 'cpp')
if BACKEND not in ('cpp', 'zig'):
    raise Exception('PARTIS_TEST_BACKEND must be cpp or zig, got %s' % BACKEND)

# small enough that at least one cdr3 group per locus in the paired simulation splits into sub-groups
HFRAC_ARGS = ['--hfrac', '--hfrac-min-seqs', '5', '--hfrac-max-bin-size', '3']


def run_partis(args, logfname):
    cmd = [sys.executable, os.path.join(REPO_DIR, 'bin', 'partis')] + args + ['--n-procs', '2', '--random-seed', '1', '--dont-write-git-info']
    if BACKEND == 'zig':
        cmd.append('--zig')
    env = dict(os.environ, PATH='%s/bin:%s' % (REPO_DIR, os.environ['PATH']))  # subset jobs call partis by name
    with open(logfname, 'w') as logfile:
        retcode = subprocess.call(cmd, stdout=logfile, stderr=subprocess.STDOUT, cwd=REPO_DIR, env=env)
    if retcode != 0:
        with open(logfname) as logfile:
            tail = ''.join(logfile.readlines()[-20:])
        raise Exception('partis exited %d: %s\nlast lines of %s:\n%s' % (retcode, ' '.join(cmd), logfname, tail))


def read_yaml(fname):
    with open(fname) as yfile:
        return yaml.safe_load(yfile)


def input_uids(locus):
    return set(u for u, info in read_yaml(os.path.join(PAIRED_SIMU_DIR, 'meta.yaml')).items() if info['locus'] == locus)


def partition_uids(fname):
    from partis.clusterpath import ClusterPath
    return set(u for cluster in ClusterPath(fname=str(fname)).best() for u in cluster)
