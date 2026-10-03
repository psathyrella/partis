import os

import pytest

from _helpers import LOCI, MIN_ANNOTATED_FRAC, N_SUBSETS, PAIRED_DATA_DIR, PAIRED_SIMU_DIR, fasta_uids, fixture_log, input_uids, locus_manifest, logged_cmd, partition_fname, partition_uids, read_log, run_partis, run_partis_fails, well_paired_uids
from partis import utils

# ----------------------------------------------------------------------------------------
# subset-partition --disjoint-groups


def subset_dir(outdir, isub):
    return outdir / ('isub-%d' % isub)


def merged_run_cmd(outdir):
    # the final merged partition command, as subset-partition printed it
    lines = [l for l in read_log(fixture_log(outdir)).split('\n') if '--input-partition-fname' in l]
    assert len(lines) == 1
    return lines[0].split()


@pytest.mark.parametrize('isub', range(N_SUBSETS))
def test_subsets_run_disjoint_hfrac(subset_hfrac_dir, isub):
    sdir = subset_dir(subset_hfrac_dir, isub)
    assert utils.subset_is_marked_complete(str(sdir))
    cmd = logged_cmd(sdir / 'log')
    assert '--disjoint-groups' in cmd and '--hfrac' in cmd
    for locus in LOCI:
        assert locus_manifest(sdir, locus)[1]['grouping-info']['method'] == 'cdr3-length+hfrac'


def test_merged_run_drops_disjoint_args(subset_hfrac_dir):
    cmd = merged_run_cmd(subset_hfrac_dir)
    for arg in ['--disjoint-groups', '--hfrac', '--hfrac-min-seqs', '--hfrac-max-bin-size', '--hfrac-merge-factor', '--ha-repartition', '--partition-refine', '--infname']:
        assert arg not in cmd
    for arg in ['--continue-from-input-partition', '--ignore-sw-pair-info', '--refuse-to-cache-parameters', '--ignore-default-input-metafile']:
        assert arg in cmd
    assert cmd[cmd.index('--parameter-dir') + 1] == '%s/merged-subsets/parameters' % subset_hfrac_dir
    assert cmd[cmd.index('--input-partition-fname') + 1] == '%s/merged-subsets' % subset_hfrac_dir


@pytest.mark.parametrize('locus', LOCI)
def test_subset_hfrac_keeps_uids(subset_hfrac_dir, locus):
    # the final single-chain output holds only what survived pair cleaning, so check each subset's
    subset_uids = [partition_uids(partition_fname(subset_dir(subset_hfrac_dir, i), locus, single_chain=True)) for i in range(N_SUBSETS)]
    assert set().union(*subset_uids) == input_uids(locus)
    assert sum(len(u) for u in subset_uids) == len(input_uids(locus))
    merged_uids = partition_uids(partition_fname(subset_hfrac_dir, locus))
    assert well_paired_uids(locus) <= merged_uids <= input_uids(locus)


@pytest.mark.parametrize('locus', LOCI)
def test_subset_unpaired_keeps_uids(subset_unpaired_dir, locus):
    in_uids = fasta_uids(os.path.join(PAIRED_DATA_DIR, '%s.fa' % locus))
    out_uids = partition_uids(partition_fname(subset_unpaired_dir, locus))
    assert out_uids <= in_uids
    assert len(out_uids) >= MIN_ANNOTATED_FRAC * len(in_uids)
    locus_manifest(subset_dir(subset_unpaired_dir, 0), locus)  # raises if grouping did not run


def test_write_subsets_only(tmp_path):
    outdir = tmp_path / 'out'
    run_partis(['subset-partition', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--disjoint-groups', '--n-subsets', str(N_SUBSETS), '--write-subsets-only', '--paired-outdir', str(outdir)], str(tmp_path / 'partis.log'))
    n_written = 0
    for isub in range(N_SUBSETS):
        sdir = subset_dir(outdir, isub)
        n_written += len(fasta_uids(sdir / 'input-seqs.fa'))
        assert not utils.subset_is_marked_complete(str(sdir))
        assert not os.path.exists(partition_fname(sdir, 'igh'))
    assert n_written == sum(len(input_uids(l)) for l in LOCI)
    assert not (outdir / 'merged-subsets').exists()


@pytest.mark.parametrize('flag', ['--ha-repartition', '--partition-refine'])
def test_refinement_flags_rejected(tmp_path, flag):
    log = run_partis_fails(['subset-partition', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--disjoint-groups', flag, '--n-subsets', str(N_SUBSETS), '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert 'do not work with subset-partition' in log
