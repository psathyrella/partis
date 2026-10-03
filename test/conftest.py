import pytest

from _helpers import HFRAC_ARGS, N_SUBSETS, PAIRED_PARAM_DIR, PAIRED_SIMU_DIR, run_partis


def run_chunked_cache_parameters(tmp_path_factory, name, extra_args):
    basedir = tmp_path_factory.mktemp(name)
    run_partis(['cache-parameters', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--n-subsets', str(N_SUBSETS), '--paired-outdir', str(basedir / 'out')] + extra_args, str(basedir / 'partis.log'))
    return basedir / 'out'


@pytest.fixture(scope='session')
def chunked_param_dir(tmp_path_factory):
    """paired cache-parameters on the tracked simulation, split into subsets, merged sw caches"""
    return run_chunked_cache_parameters(tmp_path_factory, 'chunked-merged', [])


@pytest.fixture(scope='session')
def chunked_param_index_dir(tmp_path_factory):
    """same as chunked_param_dir, but per-subset sw caches plus sw-cache-index.yaml"""
    return run_chunked_cache_parameters(tmp_path_factory, 'chunked-index', ['--no-merged-sw-cache'])


@pytest.fixture(scope='session')
def hfrac_partition_dir(tmp_path_factory):
    """paired partition --disjoint-groups --hfrac on the tracked simulation and its parameters"""
    basedir = tmp_path_factory.mktemp('hfrac-partition')
    run_partis(['partition', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--parameter-dir', PAIRED_PARAM_DIR, '--disjoint-groups'] + HFRAC_ARGS + ['--paired-outdir', str(basedir / 'out')], str(basedir / 'partis.log'))
    return basedir / 'out'
