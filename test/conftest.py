import os

import pytest

from _helpers import HFRAC_ARGS, N_SUBSETS, PAIRED_DATA_DIR, PAIRED_PARAM_DIR, PAIRED_SIMU_DIR, UNPAIRED_IGH_FNAME, paired_chunked_args, paired_partition_args, run_partis


def run_chunked_cache_parameters(tmp_path_factory, name, extra_args):
    basedir = tmp_path_factory.mktemp(name)
    run_partis(paired_chunked_args(basedir / 'out', extra_args=extra_args), str(basedir / 'partis.log'))
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


def run_fixture(tmp_path_factory, name, args, outdir_arg='--paired-outdir'):
    # one partis run into <basedir>/out, log in <basedir>/partis.log
    basedir = tmp_path_factory.mktemp(name)
    run_partis(args + [outdir_arg, str(basedir / 'out')], str(basedir / 'partis.log'))
    return basedir / 'out'


@pytest.fixture(scope='session')
def plain_multifile_partition_dir(tmp_path_factory):
    """paired partition --disjoint-groups without hfrac, multifile output forced by low caps"""
    return run_fixture(tmp_path_factory, 'plain-multifile', paired_partition_args(['--is-simu', '--multifile-min-seqs', '1', '--multifile-max-seqs-per-file', '10']))


@pytest.fixture(scope='session')
def unpaired_auto_partition_dir(tmp_path_factory):
    """partition --disjoint-groups --hfrac on one locus file, with no parameter dir so parameters are cached first"""
    basedir = tmp_path_factory.mktemp('unpaired-auto')
    run_partis(['partition', '--disjoint-groups'] + HFRAC_ARGS + ['--infname', UNPAIRED_IGH_FNAME, '--outfname', str(basedir / 'partition.yaml')], str(basedir / 'partis.log'))
    return basedir / 'partition'  # --outfname becomes a --paired-outdir named for its prefix


@pytest.fixture(scope='session')
def subset_hfrac_dir(tmp_path_factory):
    """paired subset-partition --disjoint-groups --hfrac"""
    return run_fixture(tmp_path_factory, 'subset-hfrac', ['subset-partition', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--disjoint-groups'] + HFRAC_ARGS + ['--n-subsets', str(N_SUBSETS)])


@pytest.fixture(scope='session')
def subset_unpaired_dir(tmp_path_factory):
    """subset-partition --disjoint-groups on real 10x data with its pairing info ignored"""
    return run_fixture(tmp_path_factory, 'subset-unpaired', ['subset-partition', '--paired-loci', '--paired-indir', PAIRED_DATA_DIR, '--no-pairing-info', '--keep-all-unpaired-seqs', '--disjoint-groups', '--n-subsets', str(N_SUBSETS)])


@pytest.fixture(scope='session')
def single_locus_chunked_dir(tmp_path_factory):
    """single locus cache-parameters split into subsets, merged into --parameter-dir"""
    return run_fixture(tmp_path_factory, 'single-locus-chunked', ['cache-parameters', '--infname', UNPAIRED_IGH_FNAME, '--n-subsets', str(N_SUBSETS)], outdir_arg='--parameter-dir')


@pytest.fixture(scope='session')
def unpaired_chunked_dir(tmp_path_factory):
    """chunked cache-parameters on one unpaired multi-locus file, through the paired infrastructure"""
    return run_fixture(tmp_path_factory, 'unpaired-chunked', ['cache-parameters', '--paired-loci', '--no-pairing-info', '--infname', os.path.join(PAIRED_DATA_DIR, 'all-seqs.fa'), '--n-subsets', str(N_SUBSETS)])


@pytest.fixture(scope='session')
def ha_refine_dir(tmp_path_factory):
    """paired partition --disjoint-groups --hfrac --ha-repartition --partition-refine"""
    return run_fixture(tmp_path_factory, 'ha-refine', paired_partition_args(HFRAC_ARGS + ['--ha-repartition', '--partition-refine']))


@pytest.fixture(scope='session')
def ha_only_dir(tmp_path_factory):
    """paired partition --disjoint-groups --ha-repartition"""
    return run_fixture(tmp_path_factory, 'ha-only', paired_partition_args(['--ha-repartition']))


@pytest.fixture(scope='session')
def refine_only_dir(tmp_path_factory):
    """paired partition --disjoint-groups --partition-refine"""
    return run_fixture(tmp_path_factory, 'refine-only', paired_partition_args(['--partition-refine']))
