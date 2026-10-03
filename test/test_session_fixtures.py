import collections
import os

import pytest

from _helpers import BACKEND, LOCI, N_SUBSETS, input_uids, partition_uids, read_yaml, well_paired_uids
from partis import utils
from partis.disjointgrouper import MANIFEST_FNAME, SW_CACHE_INDEX_FNAME, read_manifest
from partis.paircluster import paired_fn


@pytest.mark.parametrize('locus', LOCI)
def test_chunked_params_merged_form(chunked_param_dir, locus):
    pdir = chunked_param_dir / 'parameters' / locus
    assert not (pdir / SW_CACHE_INDEX_FNAME).exists()
    events = read_yaml(pdir / 'sw-cache.yaml')['events']
    assert set(u for line in events for u in line['unique_ids']) == input_uids(locus)


@pytest.mark.parametrize('locus', LOCI)
def test_chunked_params_index_form(chunked_param_index_dir, locus):
    pdir = chunked_param_index_dir / 'parameters' / locus
    assert not (pdir / 'sw-cache.yaml').exists()
    index = read_yaml(pdir / SW_CACHE_INDEX_FNAME)
    assert index['n_subsets'] == N_SUBSETS
    assert len(index['sw_caches']) == N_SUBSETS
    assert all((pdir / entry['path']).is_file() for entry in index['sw_caches'])
    assert sum(entry['n_sequences'] for entry in index['sw_caches']) == len(input_uids(locus))


@pytest.mark.parametrize('fixture_name', ['chunked_param_dir', 'chunked_param_index_dir'])
def test_chunked_params_subsets_complete(request, fixture_name):
    basedir = str(request.getfixturevalue(fixture_name))
    assert len(read_yaml(utils.subset_index_fname(basedir))['subsets']) == N_SUBSETS
    for isub in range(N_SUBSETS):
        assert os.path.isfile('%s/%s' % (utils.parameter_subset_dir(basedir, isub), utils.SUBSET_COMPLETE_FNAME))


def test_backend_forwarded_to_subset_jobs(chunked_param_dir):
    for isub in range(N_SUBSETS):
        with open('%s/log' % utils.parameter_subset_dir(str(chunked_param_dir), isub)) as logfile:
            assert ('--zig' in logfile.read().split()) == (BACKEND == 'zig')


@pytest.mark.parametrize('locus', LOCI)
def test_hfrac_partition_splits_a_group(hfrac_partition_dir, locus):
    manifest = read_manifest(str(hfrac_partition_dir / 'single-chain' / 'disjoint-groups' / locus / MANIFEST_FNAME))
    assert manifest['grouping-info']['method'] == 'cdr3-length+hfrac'
    assert all('/sub-groups/sub-' in group['fasta_path'] for group in manifest['groups'])
    n_sub_groups = collections.Counter(group['cdr3_length'] for group in manifest['groups'])
    assert max(n_sub_groups.values()) > 1


@pytest.mark.parametrize('locus', LOCI)
def test_hfrac_partition_keeps_uids(hfrac_partition_dir, locus):
    single_chain_uids = partition_uids(paired_fn(str(hfrac_partition_dir), locus, single_chain=True, actstr='partition', suffix='.yaml'))
    assert single_chain_uids == input_uids(locus)
    merged_uids = partition_uids(paired_fn(str(hfrac_partition_dir), locus, actstr='partition', suffix='.yaml'))
    assert well_paired_uids(locus) <= merged_uids <= single_chain_uids  # pair cleaning may drop only badly paired or unpaired seqs
