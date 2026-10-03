import json
import os
import shutil

import pytest

from _helpers import LOCI, PAIRED_PARAM_DIR, PAIRED_SIMU_DIR, copy_dir, disjoint_dir, fixture_log, input_uids, locus_manifest, partition_fname, partition_uids, read_log, read_yaml, run_partis, run_partis_fails, write_yaml
from partis import disjointgrouper as dg
from partis import partition_refinement as prf

# ----------------------------------------------------------------------------------------
# ha-repartition and partition-refine, integrated and standalone (single chain only)


def group_stage_exists(outdir, locus, stage):
    ddir, manifest = locus_manifest(outdir, locus)
    return [os.path.exists('%s/%s/%s' % (ddir, os.path.dirname(g['fasta_path']), dg.stage_fname(stage, locus))) for g in manifest['groups']]


def single_chain_uids(outdir, locus):
    return partition_uids(partition_fname(outdir, locus, single_chain=True))


@pytest.mark.parametrize('locus', LOCI)
def test_ha_alone(ha_only_dir, locus):
    assert all(group_stage_exists(ha_only_dir, locus, dg.STAGE_HAREP))
    assert not any(group_stage_exists(ha_only_dir, locus, dg.STAGE_REFINE))
    assert single_chain_uids(ha_only_dir, locus) == input_uids(locus)


@pytest.mark.parametrize('locus', LOCI)
def test_refine_alone(refine_only_dir, locus):
    assert all(group_stage_exists(refine_only_dir, locus, dg.STAGE_REFINE))
    assert not any(group_stage_exists(refine_only_dir, locus, dg.STAGE_HAREP))
    assert single_chain_uids(refine_only_dir, locus) == input_uids(locus)


def test_refine_alone_warns_on_vsearch_input(refine_only_dir):
    assert 'refine input: 0 groups from ha-repartition' in read_log(fixture_log(refine_only_dir))


@pytest.mark.parametrize('locus', LOCI)
def test_hfrac_ha_refine(ha_refine_dir, locus):
    assert all(group_stage_exists(ha_refine_dir, locus, dg.STAGE_HAREP))
    assert all(group_stage_exists(ha_refine_dir, locus, dg.STAGE_REFINE))
    _, manifest = locus_manifest(ha_refine_dir, locus)
    assert all(dg.stage_from_path(g['partition_path']) == dg.STAGE_REFINE for g in manifest['groups'])
    assert single_chain_uids(ha_refine_dir, locus) == input_uids(locus)


@pytest.mark.parametrize('locus', LOCI)
def test_refined_output_has_no_paired_merge(ha_refine_dir, locus):
    assert not os.path.exists(partition_fname(ha_refine_dir, locus))
    assert 'not running paired clustering' in read_log(fixture_log(ha_refine_dir))


def test_heavy_light_fork(ha_refine_dir):
    # only the heavy split reads the locus-wide naive threshold
    for locus in LOCI:
        assert os.path.exists(prf.locuswide_threshold_fname(disjoint_dir(ha_refine_dir, locus), locus)) == (locus == 'igh')


def test_merge_paired_refuses_refined_input(tmp_path, ha_refine_dir):
    outdir = copy_dir(ha_refine_dir, tmp_path / 'out')
    log = run_partis_fails(['merge-paired-partitions', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--paired-outdir', str(outdir)], str(tmp_path / 'merge.log'))
    assert 'came from partition-refine or ha-repartition, whose output is single-chain only' in log


@pytest.mark.parametrize('flag', ['--ha-repartition', '--partition-refine'])
def test_refinement_flags_need_disjoint_groups(tmp_path, flag):
    log = run_partis_fails(['partition', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--parameter-dir', PAIRED_PARAM_DIR, flag, '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert '--ha-repartition/--partition-refine require --disjoint-groups' in log


def test_refine_needs_length_tables(tmp_path):
    pdir = str(copy_dir(PAIRED_PARAM_DIR, tmp_path / 'params'))
    os.remove('%s/igh/hmm/%s' % (pdir, prf.LENGTH_TABLE_SPECS[0][1]))
    log = run_partis_fails(['partition', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--parameter-dir', pdir, '--disjoint-groups', '--partition-refine', '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert 'couldn\'t read length table' in log


def refine_job_args(outdir, locus, pdir=None, job_start=0, job_count=None):
    pdir = os.path.join(PAIRED_PARAM_DIR, locus) if pdir is None else pdir
    args = ['run-partition-refine-jobs', '--locus', locus, '--parameter-dir', pdir, '--paired-outdir', str(outdir), '--job-start', str(job_start)]
    return args + ([] if job_count is None else ['--job-count', str(job_count)])


def test_refine_needs_mute_freq_tables(tmp_path, plain_multifile_partition_dir):
    outdir = copy_dir(plain_multifile_partition_dir, tmp_path / 'out')
    pdir = str(copy_dir(os.path.join(PAIRED_PARAM_DIR, 'igh'), tmp_path / 'igh-params'))
    shutil.rmtree('%s/hmm/mute-freqs' % pdir)
    log = run_partis_fails(refine_job_args(outdir, 'igh', pdir=pdir), str(tmp_path / 'refine.log'))
    assert 'no per-gene mutation frequency tables' in log


@pytest.mark.parametrize('locus', ['igh', 'igk'])
def test_standalone_refine_jobs_in_slices(tmp_path, plain_multifile_partition_dir, locus):
    outdir = copy_dir(plain_multifile_partition_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, locus)
    n_groups = len(manifest['groups'])
    assert n_groups > 1
    run_partis(refine_job_args(outdir, locus, job_count=1), str(tmp_path / 'refine-0.log'))
    run_partis(refine_job_args(outdir, locus, job_start=1), str(tmp_path / 'refine-1.log'))
    assert all(group_stage_exists(outdir, locus, dg.STAGE_REFINE))
    assert all(os.path.exists(prf.bundle_marker_fname(ddir, s)) for s in [0, 1])
    outfname = tmp_path / ('assembled-%s.yaml' % locus)
    run_partis(['assemble-groups', '--locus', locus, '--paired-outdir', str(outdir), '--outfname', str(outfname)], str(tmp_path / 'assemble.log'))
    assert partition_uids(outfname) == input_uids(locus)
    assert 'all %d groups from %s' % (n_groups, dg.STAGE_REFINE) in read_log(tmp_path / 'assemble.log')


# ----------------------------------------------------------------------------------------
# stage files


def test_mixed_stages_refused_at_assemble(tmp_path, ha_refine_dir):
    outdir = copy_dir(ha_refine_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    ginfo = manifest['groups'][0]
    fasta_dir = os.path.dirname(ginfo['fasta_path'])
    os.remove('%s/%s/%s' % (ddir, fasta_dir, dg.stage_fname(dg.STAGE_REFINE, 'igh')))
    ginfo['partition_path'] = '%s/%s' % (fasta_dir, dg.stage_fname(dg.STAGE_HAREP, 'igh'))  # as if this group's refine never ran
    write_yaml('%s/%s' % (ddir, dg.MANIFEST_FNAME), manifest)
    with pytest.raises(Exception, match='cannot assemble a mix of stages'):
        dg.assemble_groups('igh', ddir, str(tmp_path / 'assembled-igh.yaml'))


def test_truncated_partition_refused_as_refine_input(tmp_path, plain_multifile_partition_dir):
    outdir = copy_dir(plain_multifile_partition_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    spec = max(prf.group_specs(ddir, manifest['groups'], 'igh'), key=lambda s: s['group']['sequence_count'])
    yinfo = read_yaml(spec['input'])
    best = yinfo['partitions'][-1]['partition']  # a block yaml file cut at a line boundary loses the tail of this list
    best[-1] = best[-1][:-1] if len(best[-1]) > 1 else None
    yinfo['partitions'][-1]['partition'] = [c for c in best if c is not None]
    write_yaml(spec['input'], yinfo)
    with pytest.raises(Exception, match='incomplete partition file'):
        prf.read_refine_inputs(spec['input'], spec['sw_cache'])


def write_stage_file(ddir, fasta_dir, stage):
    fname = '%s/%s/%s' % (ddir, fasta_dir, dg.stage_fname(stage, 'igh'))
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    open(fname, 'w').close()


def test_stage_precedence(tmp_path):
    ddir = str(tmp_path)
    ginfo = {'fasta_path': 'groups/cdr3-42/igh.fa', 'locus': 'igh'}
    assert dg.resolve_partition_path(ginfo, ddir) == (None, None, False)
    write_stage_file(ddir, 'groups/cdr3-42', dg.STAGE_VSEARCH)
    assert dg.resolve_partition_path(ginfo, ddir) == ('groups/cdr3-42/partition-igh.yaml', dg.STAGE_VSEARCH, False)
    write_stage_file(ddir, 'groups/cdr3-42', dg.STAGE_HAREP)
    assert dg.discover_partition_path(ginfo, ddir) == ('groups/cdr3-42/ha-repartition-igh.yaml', dg.STAGE_HAREP)
    # standalone actions leave the manifest's name stale, so the more refined file on disk wins
    ginfo['partition_path'] = 'groups/cdr3-42/partition-igh.yaml'
    assert dg.resolve_partition_path(ginfo, ddir) == ('groups/cdr3-42/ha-repartition-igh.yaml', dg.STAGE_HAREP, True)
    write_stage_file(ddir, 'groups/cdr3-42', dg.STAGE_REFINE)
    ginfo['partition_path'] = 'groups/cdr3-42/partition-refine-igh.yaml'
    assert dg.resolve_partition_path(ginfo, ddir) == ('groups/cdr3-42/partition-refine-igh.yaml', dg.STAGE_REFINE, False)


# ----------------------------------------------------------------------------------------
# locus-wide naive threshold cache


@pytest.fixture
def counted_threshold(monkeypatch):
    calls = []
    def fake_estimate(specs):
        calls.append(len(specs))
        return 0.125 * len(calls)
    monkeypatch.setattr(prf, 'estimate_locuswide_threshold', fake_estimate)
    return calls


def threshold_specs(tmp_path, n_specs):
    specs = []
    for ispec in range(n_specs):
        fname = tmp_path / ('input-%d.yaml' % ispec)
        if not fname.exists():
            fname.write_text('{}\n')
        specs.append({'input': str(fname)})
    return specs


def test_threshold_cache_reused_until_inputs_change(tmp_path, counted_threshold):
    specs = threshold_specs(tmp_path, 2)
    first = prf.locuswide_threshold(str(tmp_path), specs, 'igh')
    assert prf.locuswide_threshold(str(tmp_path), specs, 'igh') == first
    assert len(counted_threshold) == 1
    with open(prf.locuswide_threshold_fname(str(tmp_path), 'igh')) as tfile:
        cached = json.load(tfile)
    assert float(cached['threshold']) == first and cached['n_specs'] == 2
    os.utime(specs[0]['input'], (0, 12345))
    second = prf.locuswide_threshold(str(tmp_path), specs, 'igh')
    assert second != first  # an input changed
    assert prf.locuswide_threshold(str(tmp_path), threshold_specs(tmp_path, 3), 'igh') != second  # spec count changed
    assert len(counted_threshold) == 3
    prf.locuswide_threshold(str(tmp_path), threshold_specs(tmp_path, 3), 'igh', overwrite=True)
    assert len(counted_threshold) == 4


def test_threshold_none_without_d_gene_or_specs(tmp_path, counted_threshold):
    assert prf.locuswide_threshold(str(tmp_path), threshold_specs(tmp_path, 2), 'igk') is None
    assert prf.locuswide_threshold(str(tmp_path), [], 'igh') is None
    assert len(counted_threshold) == 0
    assert not os.path.exists(prf.locuswide_threshold_fname(str(tmp_path), 'igk'))
