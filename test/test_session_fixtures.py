import collections

import pytest

from _helpers import BACKEND, LOCI, input_uids, partition_uids, read_yaml


@pytest.mark.parametrize('locus', LOCI)
def test_chunked_params_merged_form(chunked_param_dir, locus):
    pdir = chunked_param_dir / 'parameters' / locus
    assert (pdir / 'sw-cache.yaml').is_file()
    assert not (pdir / 'sw-cache-index.yaml').exists()


@pytest.mark.parametrize('locus', LOCI)
def test_chunked_params_index_form(chunked_param_index_dir, locus):
    pdir = chunked_param_index_dir / 'parameters' / locus
    assert (pdir / 'sw-cache-index.yaml').is_file()
    assert not (pdir / 'sw-cache.yaml').exists()


def test_chunked_params_subsets_complete(chunked_param_dir):
    sdir = chunked_param_dir / 'parameter-subsets'
    assert (sdir / 'subset-index.yaml').is_file()
    for isub in range(2):
        assert (sdir / ('subset-%d' % isub) / 'subset-complete.yaml').is_file()


def test_backend_forwarded_to_subset_jobs(chunked_param_dir):
    for isub in range(2):
        log_words = (chunked_param_dir / 'parameter-subsets' / ('subset-%d' % isub) / 'log').read_text().split()
        assert ('--zig' in log_words) == (BACKEND == 'zig')


@pytest.mark.parametrize('locus', LOCI)
def test_hfrac_partition_splits_a_group(hfrac_partition_dir, locus):
    manifest = read_yaml(hfrac_partition_dir / 'single-chain' / 'disjoint-groups' / locus / 'manifest.yaml')
    assert manifest['grouping-info']['method'] == 'cdr3-length+hfrac'
    assert all('/sub-groups/sub-' in group['fasta_path'] for group in manifest['groups'])
    n_sub_groups = collections.Counter(group['cdr3_length'] for group in manifest['groups'])
    assert max(n_sub_groups.values()) > 1


@pytest.mark.parametrize('locus', LOCI)
def test_hfrac_partition_keeps_uids(hfrac_partition_dir, locus):
    single_chain_uids = partition_uids(hfrac_partition_dir / 'single-chain' / ('partition-%s.yaml' % locus))
    assert single_chain_uids == input_uids(locus)
    assert partition_uids(hfrac_partition_dir / ('partition-%s.yaml' % locus)) <= single_chain_uids  # pair cleaning may drop unpaired seqs
