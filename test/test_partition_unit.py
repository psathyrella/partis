import os
import re

import pytest

from _helpers import HFRAC_ARGS, LOCI, PAIRED_PARAM_DIR, copy_dir, disjoint_dir, group_cmds, group_stage_path, locus_manifest, paired_partition_args, partition_fname, read_log, read_yaml, run_partis, run_partis_fails
from partis import disjointgrouper as dg
from partis import ha_repartition
from partis import partition_refinement as prf
from partis import utils

# ----------------------------------------------------------------------------------------
# the partition unit (grouping, vsearch, HA, refine, assemble) re-run on an existing output dir


def stage_mtimes(outdir, locus, stage):
    ddir, manifest = locus_manifest(outdir, locus)
    return [os.path.getmtime(group_stage_path(ddir, g, stage, locus)) for g in manifest['groups']]


def create_groups_args(outdir, extra_args=None):
    return ['create-disjoint-groups', '--locus', 'igh', '--parameter-dir', PAIRED_PARAM_DIR, '--paired-outdir', str(outdir)] + (extra_args or [])


@pytest.mark.parametrize('fixture, extra_args', [('hfrac_partition_dir', HFRAC_ARGS), ('ha_only_dir', ['--ha-repartition']), ('refine_only_dir', ['--partition-refine'])], ids=['hfrac', 'ha', 'refine'])
def test_unchanged_rerun_resumes(request, tmp_path, fixture, extra_args):
    outdir = copy_dir(request.getfixturevalue(fixture), tmp_path / 'out')
    before = {l: locus_manifest(outdir, l)[1]['groups'] for l in LOCI}
    run_partis(paired_partition_args(extra_args, outdir=outdir), str(tmp_path / 'rerun.log'))
    logstr = read_log(str(tmp_path / 'rerun.log'))
    assert logstr.count('manifest with the same inputs exists, skipping grouping') == len(LOCI)
    assert not re.search(r'running \d+ group partition jobs', logstr)
    assert all(locus_manifest(outdir, l)[1]['groups'] == before[l] for l in LOCI)


@pytest.mark.parametrize('extra_args, errstr', [([], '--hfrac False but manifest has True'),
                                                (['--hfrac', '--hfrac-min-seqs', '5', '--hfrac-max-bin-size', '4'], '--hfrac-max-bin-size 4 but manifest has 3')])
def test_changed_grouping_args_raise(hfrac_partition_dir, tmp_path, extra_args, errstr):
    outdir = copy_dir(hfrac_partition_dir, tmp_path / 'out')
    logstr = run_partis_fails(paired_partition_args(extra_args, outdir=outdir), str(tmp_path / 'rerun.log'))
    assert 'disjoint dir %s was grouped from different inputs' % disjoint_dir(outdir, 'igh') in logstr
    assert errstr in logstr


@pytest.mark.parametrize('edited, recorded', [(dg.SW_CACHE_FNAME, dg.SW_CACHE_FNAME), ('sw/%s' % dg.MEAN_MFREQ_FNAME, 'sw/%s' % dg.MEAN_MFREQ_FNAME), ('hmm/%s' % dg.MEAN_MFREQ_FNAME, 'hmm')], ids=['sw-cache', 'hfrac-mfreq', 'hmm-dir'])
def test_changed_input_file_raises(hfrac_partition_dir, tmp_path, edited, recorded):
    outdir = copy_dir(hfrac_partition_dir, tmp_path / 'out')
    pdir = copy_dir(PAIRED_PARAM_DIR, tmp_path / 'parameters')
    with open('%s/igh/%s' % (pdir, edited), 'a') as efile:
        efile.write('\n')
    logstr = run_partis_fails(paired_partition_args(HFRAC_ARGS, param_dir=pdir, outdir=outdir), str(tmp_path / 'rerun.log'))
    assert re.search(r'grouped from different inputs .*: input files .*%s/igh/%s .* but manifest has' % (pdir, recorded), logstr)


@pytest.mark.parametrize('fixture, stage', [('ha_only_dir', dg.STAGE_HAREP), ('refine_only_dir', dg.STAGE_REFINE)], ids=['ha', 'refine'])
def test_dropped_stage_raises(request, tmp_path, fixture, stage):
    outdir = copy_dir(request.getfixturevalue(fixture), tmp_path / 'out')
    logstr = run_partis_fails(paired_partition_args([], outdir=outdir), str(tmp_path / 'rerun.log'))
    n_groups = len(locus_manifest(outdir, 'igh')[1]['groups'])
    assert 'disjoint dir %s has output from stages not requested in this run' % disjoint_dir(outdir, 'igh') in logstr
    assert '%d groups with %s' % (n_groups, stage) in logstr


def test_dir_hash(tmp_path):
    for fname, text in [('a.csv', 'x'), ('sub/b.csv', 'y')]:
        utils.mkdir(str(tmp_path / fname), isfile=True)
        (tmp_path / fname).write_text(text)
    first = utils.xxh3_dir_hash(str(tmp_path))
    assert utils.xxh3_dir_hash(str(tmp_path)) == first
    (tmp_path / 'sub/b.csv').write_text('z')
    assert utils.xxh3_dir_hash(str(tmp_path)) != first
    (tmp_path / 'sub/b.csv').write_text('y')
    (tmp_path / 'sub/b.csv').rename(tmp_path / 'sub/c.csv')
    assert utils.xxh3_dir_hash(str(tmp_path)) != first


def test_overwrite_redoes_unit(hfrac_partition_dir, tmp_path):
    outdir = copy_dir(hfrac_partition_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    stale = [group_stage_path(ddir, manifest['groups'][0], dg.STAGE_REFINE, 'igh'), ha_repartition.bundle_marker_fname(ddir, 'oldplan', 0), prf.bundle_marker_fname(ddir, 0)]
    for fname in stale:
        open(fname, 'w').close()
    run_partis(paired_partition_args(['--overwrite'], outdir=outdir), str(tmp_path / 'rerun.log'))
    assert not any(os.path.exists(f) for f in stale)
    cmds = group_cmds(outdir, 'igh')
    assert len(cmds) > 0 and all('--overwrite' not in cmd for cmd in cmds)
    assert 'removing %s' % partition_fname(outdir, 'igh', single_chain=True) in read_log(str(tmp_path / 'rerun.log'))
    ginfo = locus_manifest(outdir, 'igh')[1]['grouping-info']
    assert ginfo['method'] == 'cdr3-length'
    assert ginfo['inputs']['args'] == {'hfrac': False}


def test_overwrite_removes_multifile_output(plain_multifile_partition_dir, tmp_path):
    outdir = copy_dir(plain_multifile_partition_dir, tmp_path / 'out')
    run_partis(paired_partition_args(['--is-simu', '--multifile-min-seqs', '1', '--multifile-max-seqs-per-file', '10', '--overwrite'], outdir=outdir), str(tmp_path / 'rerun.log'))
    mfdir = dg.multifile_dir_path(partition_fname(outdir, 'igh', single_chain=True))
    assert 'removing %s' % mfdir in read_log(str(tmp_path / 'rerun.log'))
    assert os.path.exists('%s/%s' % (mfdir, dg.MULTIFILE_INDEX_FNAME))


def test_standalone_rerun(tmp_path):
    ddir = disjoint_dir(tmp_path, 'igh')
    stale = '%s/groups/cdr3-999/%s' % (ddir, dg.stage_fname(dg.STAGE_VSEARCH, 'igh'))  # left by a grouping that never wrote its manifest
    utils.mkdir(stale, isfile=True)
    open(stale, 'w').close()
    run_partis(create_groups_args(tmp_path), str(tmp_path / 'first.log'))
    assert not os.path.exists(stale)
    run_partis(create_groups_args(tmp_path), str(tmp_path / 'same.log'))
    assert 'manifest with the same inputs exists, skipping grouping' in read_log(str(tmp_path / 'same.log'))
    logstr = run_partis_fails(create_groups_args(tmp_path, HFRAC_ARGS), str(tmp_path / 'changed.log'))
    assert 'disjoint dir %s was grouped from different inputs' % ddir in logstr
    run_partis(create_groups_args(tmp_path, HFRAC_ARGS + ['--overwrite']), str(tmp_path / 'overwrite.log'))
    assert locus_manifest(tmp_path, 'igh')[1]['grouping-info']['method'] == 'cdr3-length+hfrac'


def test_hfrac_knobs_ignored_without_hfrac():
    swcache = '%s/igh/%s' % (PAIRED_PARAM_DIR, dg.SW_CACHE_FNAME)
    assert dg.grouping_inputs(swcache, PAIRED_PARAM_DIR, 'igh', None, False, 1, 2, 3) == dg.grouping_inputs(swcache, PAIRED_PARAM_DIR, 'igh', None, False, 4, 5, 6)
    assert dg.grouping_inputs(swcache, PAIRED_PARAM_DIR, 'igh', None, True, 1, 2, 3) != dg.grouping_inputs(swcache, PAIRED_PARAM_DIR, 'igh', None, True, 4, 5, 6)


def test_grouping_inputs_hash_forms(chunked_param_dir, chunked_param_index_dir):
    swcache = '%s/parameters/igh/%s' % (chunked_param_dir, dg.SW_CACHE_FNAME)
    assert [f['xxh3'] for f in dg.grouping_inputs(swcache, None, 'igh', None, False, 1, 2, 3)['files']] == [utils.xxh3_file_hash(swcache)]
    index_pdir = '%s/parameters/igh' % chunked_param_index_dir
    assert [f['xxh3'] for f in dg.grouping_inputs(index_pdir, None, 'igh', None, False, 1, 2, 3)['files']] == [e['xxh3'] for e in read_yaml(dg.sw_cache_index_fname(index_pdir))['sw_caches']]


# ----------------------------------------------------------------------------------------
# temporary refine-only re-run, removed with --rerun-refine


def test_rerun_refine(refine_only_dir, tmp_path):
    outdir = copy_dir(refine_only_dir, tmp_path / 'out')
    before = {l: (stage_mtimes(outdir, l, dg.STAGE_VSEARCH), stage_mtimes(outdir, l, dg.STAGE_REFINE)) for l in LOCI}
    run_partis(paired_partition_args(['--partition-refine', '--rerun-refine'], outdir=outdir), str(tmp_path / 'rerun.log'))
    for locus in LOCI:  # copy_dir keeps mtimes, so a rewritten file is newer
        assert stage_mtimes(outdir, locus, dg.STAGE_VSEARCH) == before[locus][0]
        assert all(new > old for new, old in zip(stage_mtimes(outdir, locus, dg.STAGE_REFINE), before[locus][1]))


def test_overwrite_refused_on_refine_jobs(tmp_path):
    logstr = run_partis_fails(['run-partition-refine-jobs', '--locus', 'igh', '--parameter-dir', '%s/igh' % PAIRED_PARAM_DIR, '--paired-outdir', str(tmp_path), '--overwrite'], str(tmp_path / 'refine.log'))
    assert 'does not apply to \'run-partition-refine-jobs\' (use --rerun-refine)' in logstr


def test_rerun_refine_needs_partition_refine(tmp_path):
    logstr = run_partis_fails(paired_partition_args(['--rerun-refine'], outdir=tmp_path), str(tmp_path / 'rerun.log'))
    assert '--rerun-refine requires --partition-refine' in logstr


def test_merge_fdr_cutoff():
    null = [i / 1000. for i in range(1, 1000)]
    assert prf.merge_fdr_cutoff([1e-10], null[:prf.MERGE_NULL_MIN_SAMPLE - 1]) == (None, 'null sample too small')
    assert prf.merge_fdr_cutoff([1e-10] * 20 + [0.5] * 80, null) == (1e-10, 'fdr')
    assert prf.merge_fdr_cutoff([0.5] * 100, null) == (None, 'no candidate p meets the fdr target')


class _FakeEvidence(object):
    # what sample_null reads, with a p-value that encodes its pair
    def __init__(self, frags):
        self.frags, self.fps = frags, [({}, 1) for _ in frags]
    def pvalue(self, i, j):
        return 10 * min(i, j) + max(i, j)


def test_sample_null_pairs():
    frags = [['a1', 'a2'], ['b'], ['c'], ['d1', 'd2'], ['e1', 'e2']]
    naives = ['AAAAAAAAAA', 'CCCCCCCCCC', 'GGGGGGGGGG', 'AAAAAAAAAA', 'TTTTTTTTTT']
    vdjs = [('v1', 'd', 'j')] * 4 + [('v2', 'd', 'j')]
    out = prf._MergeEvidence.sample_null(_FakeEvidence(frags), 100, naives, vdjs, 0.1)
    # same v and length only, no singleton pair (1, 2), no pair within the threshold (0, 3), each pair once
    assert sorted(out) == [1, 2, 13, 23]


def test_merge_evidence_discounts_entries_every_fragment_carries():
    common, n_frags, n_seqs = (1, 'A'), 20, 40
    shared = [(5, 'T'), (8, 'G'), (12, 'C')]  # fragments 0 and 1 also share these
    muts, frags = {}, []
    for ifrag in range(n_frags):
        own = shared if ifrag < 2 else [(20 + ifrag, 'T')]
        uids = ['f%d-%d' % (ifrag, k) for k in range(2)]
        for uid in uids:
            muts[uid] = dict([common] + own)
        frags.append(uids)
    fps = [prf.get_cluster_fingerprint(f, muts) for f in frags]
    def get_uid_freqs(uid):
        return {k: 0.01 for k in muts[uid].items()}, {}
    ev = prf._MergeEvidence(frags, fps, muts, get_uid_freqs, n_seqs)
    assert ev.group_rate[common] == 1.
    assert ev.pvalue(1, 0) == ev.pvalue(0, 1)
    assert ev.pvalue(0, 1) < prf.MERGE_FDR_TARGET < ev.pvalue(2, 3)
