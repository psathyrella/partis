import glob
import os

import pytest

from _helpers import BACKEND, LOCI, UNPAIRED_IGH_FNAME, fasta_uids, fixture_log, input_uids, locus_manifest, logged_cmd, output_partition_uids, partition_fname, partition_uids, read_log, well_paired_uids
from partis import disjointgrouper as dg

# ----------------------------------------------------------------------------------------
# partition --disjoint-groups, integrated route


def group_cmds(outdir, locus):
    ddir, _ = locus_manifest(outdir, locus)
    return [logged_cmd(f) for f in glob.glob('%s/groups/**/log' % ddir, recursive=True)]


@pytest.mark.parametrize('locus', LOCI)
def test_plain_grouping_one_group_per_cdr3(plain_multifile_partition_dir, locus):
    _, manifest = locus_manifest(plain_multifile_partition_dir, locus)
    assert manifest['grouping-info']['method'] == 'cdr3-length'
    assert len(manifest['groups']) > 1
    assert len(set(g['cdr3_length'] for g in manifest['groups'])) == len(manifest['groups'])
    assert all('/sub-groups/' not in g['fasta_path'] for g in manifest['groups'])


@pytest.mark.parametrize('locus', LOCI)
def test_multifile_output_keeps_uids(plain_multifile_partition_dir, locus):
    sc_fname = partition_fname(plain_multifile_partition_dir, locus, single_chain=True)
    assert not os.path.exists(sc_fname)
    index = dg.read_multifile_index('%s/%s' % (dg.multifile_dir_path(sc_fname), dg.MULTIFILE_INDEX_FNAME))
    assert index['assembly']['n_cdr3_groups'] == len(locus_manifest(plain_multifile_partition_dir, locus)[1]['groups'])
    assert output_partition_uids(sc_fname) == input_uids(locus)


@pytest.mark.parametrize('locus', LOCI)
def test_multifile_feeds_paired_merge(plain_multifile_partition_dir, locus):
    merged_uids = partition_uids(partition_fname(plain_multifile_partition_dir, locus))
    assert well_paired_uids(locus) <= merged_uids <= input_uids(locus)


@pytest.mark.parametrize('locus', LOCI)
def test_group_command_construction(plain_multifile_partition_dir, locus):
    cmds = group_cmds(plain_multifile_partition_dir, locus)
    assert len(cmds) == len(locus_manifest(plain_multifile_partition_dir, locus)[1]['groups'])
    for cmd in cmds:
        assert cmd[cmd.index('partis') + 1] == 'partition'
        assert cmd[cmd.index('--locus') + 1] == locus
        assert cmd.count('--n-procs') == 1
        for arg in ['--refuse-to-cache-parameters', '--crash-on-duplicate-uids', '--naive-vsearch', '--sw-cachefname']:
            assert arg in cmd
        for arg in ['--disjoint-groups', '--paired-loci', '--paired-indir', '--paired-outdir', '--is-simu']:
            assert arg not in cmd
        assert os.path.basename(cmd[cmd.index('--parameter-dir') + 1]) == locus  # locus-level parameter dir
        assert ('--zig' in cmd) == (BACKEND == 'zig')


@pytest.mark.parametrize('locus', LOCI)
def test_grouping_from_tracked_sw_cache(hfrac_partition_dir, locus):
    _, manifest = locus_manifest(hfrac_partition_dir, locus)
    assert manifest['grouping-info']['total_input_sequences'] == len(input_uids(locus))
    assert manifest['grouping-info']['failed_sequences'] == 0


# ----------------------------------------------------------------------------------------
# one locus file, no parameter dir: auto-enabled paired infrastructure and parameter caching


def test_unpaired_auto_caches_parameters(unpaired_auto_partition_dir):
    pdir = unpaired_auto_partition_dir / 'parameters' / 'igh'
    assert (pdir / dg.SW_CACHE_FNAME).is_file()
    assert (pdir / 'hmm' / 'all-mean-mute-freqs.csv').is_file()
    assert 'turning on --paired-loci' in read_log(fixture_log(unpaired_auto_partition_dir))


def test_unpaired_auto_hfrac_partition(unpaired_auto_partition_dir):
    _, manifest = locus_manifest(unpaired_auto_partition_dir, 'igh')
    assert manifest['grouping-info']['method'] == 'cdr3-length+hfrac'
    assert all('/sub-groups/sub-' in g['fasta_path'] for g in manifest['groups'])
    sc_uids = partition_uids(partition_fname(unpaired_auto_partition_dir, 'igh', single_chain=True))
    assert len(sc_uids) == manifest['grouping-info']['total_grouped_sequences']
    assert sc_uids <= fasta_uids(UNPAIRED_IGH_FNAME)
    assert partition_uids(partition_fname(unpaired_auto_partition_dir, 'igh')) == sc_uids  # no pairing, so linked through
