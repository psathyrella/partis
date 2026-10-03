import os
import re
import shutil
import subprocess
import sys

import yaml

REPO_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
TEST_DIR = os.path.join(REPO_DIR, 'test')

# import partis from this checkout, not from whatever is installed
sys.path.insert(0, REPO_DIR)
from partis import utils
from partis import disjointgrouper as dg
from partis.clusterpath import ClusterPath
from partis.paircluster import paired_fn

# tracked paired simulation and the parameters inferred on it
PAIRED_SIMU_DIR = os.path.join(TEST_DIR, 'paired', 'ref-results', 'test', 'simu')
PAIRED_PARAM_DIR = os.path.join(TEST_DIR, 'paired', 'ref-results', 'test', 'parameters', 'simu')
# tracked 10x data, split by locus
PAIRED_DATA_DIR = os.path.join(TEST_DIR, 'paired-data')
UNPAIRED_IGH_FNAME = os.path.join(PAIRED_DATA_DIR, 'igh.fa')
LOCI = ['igh', 'igk', 'igl']
N_SUBSETS = 2
N_PROCS = 2  # passed to every partis run
MIN_ANNOTATED_FRAC = 0.9  # real data loses only the seqs that fail annotation
PARTIS_TIMEOUT = 900  # seconds per partis run

BACKEND = os.environ.get('PARTIS_TEST_BACKEND', 'cpp')
if BACKEND not in ('cpp', 'zig'):
    raise Exception('PARTIS_TEST_BACKEND must be cpp or zig, got %s' % BACKEND)

# small enough that at least one cdr3 group per locus in the paired simulation splits into sub-groups
HFRAC_ARGS = ['--hfrac', '--hfrac-min-seqs', '5', '--hfrac-max-bin-size', '3']


class PartisFailed(Exception):
    pass


def run_partis(args, logfname):
    cmd = [sys.executable, os.path.join(REPO_DIR, 'bin', 'partis')] + args + ['--n-procs', str(N_PROCS), '--random-seed', '1', '--dont-write-git-info']
    if BACKEND == 'zig':
        cmd.append('--zig')
    env = dict(os.environ, PATH='%s/bin:%s' % (REPO_DIR, os.environ['PATH']))  # subset jobs call partis by name
    with open(logfname, 'w') as logfile:
        retcode = subprocess.run(cmd, stdout=logfile, stderr=subprocess.STDOUT, cwd=REPO_DIR, env=env, timeout=PARTIS_TIMEOUT).returncode
    if retcode != 0:
        with open(logfname) as logfile:
            tail = ''.join(logfile.readlines()[-20:])
        raise PartisFailed('partis exited %d: %s\nlast lines of %s:\n%s' % (retcode, ' '.join(cmd), logfname, tail))


def read_yaml(fname):
    with open(fname) as yfile:
        return yaml.safe_load(yfile)


def simu_meta():
    return read_yaml(os.path.join(PAIRED_SIMU_DIR, 'meta.yaml'))


def input_uids(locus):
    return set(u for u, info in simu_meta().items() if info['locus'] == locus)


def well_paired_uids(locus):
    # paired to exactly one other-chain uid, which is paired back to only it
    meta = simu_meta()
    return set(u for u, info in meta.items() if info['locus'] == locus and len(info['paired-uids']) == 1 and meta[info['paired-uids'][0]]['paired-uids'] == [u])


def partition_uids(fname):
    return set(u for cluster in ClusterPath(fname=str(fname)).best() for u in cluster)


def output_partition_uids(fname):
    # also reads the multifile dir that replaces <fname> for a large locus
    _, _, cpath = utils.read_output(str(fname), skip_annotations=True)
    return set(u for cluster in cpath.best() for u in cluster)


def run_partis_fails(args, logfname):
    # run partis expecting a nonzero exit, return its log text
    try:
        run_partis(args, logfname)
    except PartisFailed:
        return read_log(logfname)
    raise Exception('partis succeeded but was expected to fail: %s' % ' '.join(args))


def read_log(logfname):
    # log text without the terminal color codes
    with open(logfname) as logfile:
        return re.sub(r'\x1b\[[0-9;]*m', '', logfile.read())


def disjoint_dir(outdir, locus):
    return os.path.join(str(outdir), 'single-chain', 'disjoint-groups', locus)


def copy_dir(srcdir, dstdir):
    # copy of a session fixture, so a test can change it without touching the fixture
    shutil.copytree(str(srcdir), str(dstdir), symlinks=True)
    return dstdir


def write_yaml(fname, data):
    with open(fname, 'w') as yfile:
        yaml.dump(data, yfile, width=400, default_flow_style=False)


def fixture_log(outdir):
    # log of the partis run that made a session fixture's <outdir>
    return outdir.parent / 'partis.log'


def logged_cmd(logfname):
    # partis subprocess logs start with their command line
    with open(logfname) as logfile:
        return logfile.readline().split()


def locus_manifest(outdir, locus):
    ddir = disjoint_dir(outdir, locus)
    return ddir, dg.read_manifest('%s/%s' % (ddir, dg.MANIFEST_FNAME))


def partition_fname(outdir, locus, single_chain=False):
    return paired_fn(str(outdir), locus, single_chain=single_chain, actstr='partition', suffix='.yaml')


def fasta_uids(fname):
    return set(sfo['name'] for sfo in utils.read_fastx(str(fname)))


def sw_cache_uids(pdir):
    return set(u for line in read_yaml('%s/%s' % (pdir, dg.SW_CACHE_FNAME))['events'] for u in line['unique_ids'])


def check_subsets_complete(basedir):
    assert utils.read_subset_index(str(basedir))['n_subsets'] == N_SUBSETS
    for isub in range(N_SUBSETS):
        assert utils.subset_is_marked_complete(utils.parameter_subset_dir(str(basedir), isub))
