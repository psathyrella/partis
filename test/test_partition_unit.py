import multiprocessing
import os
import random
import re
import resource
import shutil
import struct
import subprocess
import sys
import time

import numpy as np
import psutil
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


@pytest.mark.parametrize('edited, recorded', [(dg.SW_CACHE_FNAME, dg.SW_CACHE_FNAME), ('sw/%s' % dg.MEAN_MFREQ_FNAME, 'sw/%s' % dg.MEAN_MFREQ_FNAME), ('sw/%s' % utils.mute_freq_sum_fname, 'sw/%s' % utils.mute_freq_sum_fname),
                                              ('sw/all-mean-n-muted.csv', 'sw/all-mean-n-muted.csv'), ('hmm/%s' % dg.MEAN_MFREQ_FNAME, 'hmm')], ids=['sw-cache', 'hfrac-mfreq', 'hfrac-sums', 'hfrac-n-muted', 'hmm-dir'])
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


def hfrac_inputs(pdir, ptype_dir=None):
    return dg.grouping_inputs('%s/igh/%s' % (pdir, dg.SW_CACHE_FNAME), pdir, 'igh', '%s/igh' % pdir, True, 2, 3, 5, ptype_dir=ptype_dir)


def test_sums_file_appears_or_disappears(tmp_path):
    full, no_sums = [copy_dir(PAIRED_PARAM_DIR, tmp_path / name) for name in ['full', 'no-sums']]
    os.remove('%s/igh/sw/%s' % (no_sums, utils.mute_freq_sum_fname))
    with_sums, without_sums = hfrac_inputs(full), hfrac_inputs(no_sums)
    assert len(with_sums['files']) == len(without_sums['files']) + 1
    assert utils.input_diffs(without_sums, hfrac_inputs(no_sums), 'manifest') == []
    assert len(utils.input_diffs(with_sums, without_sums, 'manifest')) == 1
    assert len(utils.input_diffs(without_sums, with_sums, 'manifest')) == 1


def test_hmm_only_layout(tmp_path):
    pdir = copy_dir(PAIRED_PARAM_DIR, tmp_path / 'parameters')
    shutil.rmtree('%s/igh/sw' % pdir)
    mfreq_fnames = [f for f in utils.mean_mfreq_fnames('%s/igh/hmm' % pdir) if os.path.exists(f)]
    assert len(mfreq_fnames) > 0
    assert [f['path'] for f in hfrac_inputs(pdir)['files']][1:] == mfreq_fnames + ['%s/igh/hmm' % pdir]


def test_parameter_type_dir_hashed(tmp_path):
    pdir = copy_dir(PAIRED_PARAM_DIR, tmp_path / 'parameters')
    hmm_dir, multi_dir = ['%s/igh/%s' % (pdir, sdir) for sdir in ['hmm', 'multi-hmm']]
    base = hfrac_inputs(pdir)
    assert hfrac_inputs(pdir, ptype_dir=hmm_dir) == base  # hmm/ is hashed once
    copy_dir(hmm_dir, multi_dir)
    with_multi = hfrac_inputs(pdir, ptype_dir=multi_dir)
    assert [f['path'] for f in with_multi['files']] == [f['path'] for f in base['files']] + [multi_dir]
    with open('%s/%s' % (multi_dir, dg.MEAN_MFREQ_FNAME), 'a') as efile:
        efile.write('\n')
    assert len(utils.input_diffs(with_multi, hfrac_inputs(pdir, ptype_dir=multi_dir), 'manifest')) == 1


def test_multi_hmm_dir_raises_on_rerun(hfrac_partition_dir, tmp_path):
    outdir = copy_dir(hfrac_partition_dir, tmp_path / 'out')
    pdir = copy_dir(PAIRED_PARAM_DIR, tmp_path / 'parameters')
    copy_dir('%s/igh/hmm' % pdir, '%s/igh/multi-hmm' % pdir)
    logstr = run_partis_fails(paired_partition_args(HFRAC_ARGS, param_dir=pdir, outdir=outdir), str(tmp_path / 'rerun.log'))
    assert re.search(r'grouped from different inputs .*: input files .*%s/igh/multi-hmm \(xxh3 .* but manifest has' % pdir, logstr)


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
    def pvalue(self, i, j, cache_pair=True):
        assert not cache_pair, 'null pairs are not cached'
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
    ev = prf._MergeEvidence(frags, fps, muts, lambda uid, pos, base: (0.01, None), n_seqs)
    assert ev.group_rate[common] == 1.
    assert ev.pvalue(1, 0) == ev.pvalue(0, 1)
    assert ev.pvalue(0, 1) < prf.MERGE_FDR_TARGET < ev.pvalue(2, 3)


# ----------------------------------------------------------------------------------------
# the vectorized refine paths against the scalar code they replaced, which must agree exactly


def _dp_reference(wb, shared):
    t_obs = sum(wb[pb][1] for pb in shared if pb in wb)
    if t_obs <= 0:
        return 1.0
    dist, tail = [1.0] + [0.0] * (t_obs - 1), 0.0
    for p, nbin in wb.values():
        if nbin == 0:
            continue
        nxt = [0.0] * t_obs
        for b, mass in enumerate(dist):
            if mass == 0.0:
                continue
            nxt[b] += mass * (1.0 - p)
            if b + nbin >= t_obs:
                tail += mass * p
            else:
                nxt[b + nbin] += mass * p
        dist = nxt
    return min(max(tail, 0.0), 1.0)


def _curve_reference(wb):
    dist = [1.0]
    for p, nbin in wb.values():
        if nbin == 0:
            continue
        nxt = [0.0] * (len(dist) + nbin)
        for b, mass in enumerate(dist):
            if mass == 0.0:
                continue
            nxt[b] += mass * (1.0 - p)
            nxt[b + nbin] += mass * p
        dist = nxt
    survival, running = [0.0] * len(dist), 0.0
    for t in range(len(dist) - 1, -1, -1):
        running += dist[t]
        survival[t] = running
    return survival, len(dist) - 1


def test_dp_matches_scalar_reference():
    rng = random.Random(3)
    for _ in range(500):
        muts = {pos: 'A' for pos in rng.sample(range(400), rng.randint(1, 50))}
        freqs = {(pos, 'A'): 10 ** rng.uniform(-5, 0) for pos in muts}
        wb = prf._weight_bins(muts, freqs, 1000, prf.WEIGHT_GRID_NATS, n_obs={})
        shared = [k for k in wb if rng.random() < 0.3]
        assert prf._tail_dp(wb, prf._t_obs(wb, shared)) == _dp_reference(wb, shared)
        assert prf._build_survival_curve(wb) == _curve_reference(wb)


def test_first_carriers_match_scan():
    rng = random.Random(4)
    for _ in range(50):
        frag = ['u%d' % i for i in range(rng.randint(1, 8))]
        muts = {u: {pos: rng.choice('AC') for pos in rng.sample(range(6), rng.randint(0, 4))} for u in frag if rng.random() > 0.2}
        first = prf.first_carriers(frag, muts)
        assert set(first) == set(e for m in muts.values() for e in m.items())
        for (pos, base), uid in first.items():
            assert uid == next(u for u in frag if u in muts and muts[u].get(pos) == base)


def test_hamming_frac_rows_matches_hamming_frac():
    rng = random.Random(4)
    naives = [None, '', 'NNNN'] + [''.join(rng.choice('ACGTN') for _ in range(rng.randint(1, 12))) for _ in range(60)]
    enc = prf.encode_naives(naives)
    pairs = [(a, b) for a in range(1, len(naives)) for b in range(1, len(naives))]  # row 0 (None) is never drawn
    ia, ib = [a for a, _ in pairs], [b for _, b in pairs]
    assert prf.hamming_frac_rows(enc, ia, ib).tolist() == [prf.hamming_frac(naives[a], naives[b]) for a, b in pairs]


def _naive_threshold_reference(partition, uid_sw_naives):
    hams = []
    for cluster in partition:
        unique_naives = list(set(uid_sw_naives[u] for u in cluster if u in uid_sw_naives))
        if len(cluster) < 3 or len(unique_naives) < 2:
            continue
        hams += [prf.hamming_frac(a, b) for i, a in enumerate(unique_naives) for b in unique_naives[i + 1:]]
    return sorted(hams)[min(int(len(hams) * 0.90), len(hams) - 1)] if hams else prf.NAIVE_THRESHOLD_FALLBACK


@pytest.mark.parametrize('block', [prf.HAMMING_PAIR_BLOCK, 7])
def test_naive_threshold_matches_scalar_reference(monkeypatch, block):
    monkeypatch.setattr(prf, 'HAMMING_PAIR_BLOCK', block)
    rng = random.Random(5)
    for _ in range(100):
        uid_sw_naives, partition = {}, []
        for ic in range(rng.randint(0, 6)):
            cluster = ['c%d-%d' % (ic, k) for k in range(rng.randint(1, 15))]
            for uid in cluster[1:]:  # first uid has no sw naive
                uid_sw_naives[uid] = ''.join(rng.choice('ACGN') for _ in range(rng.randint(0, 8)))
            partition.append(cluster)
        threshold = prf.estimate_naive_threshold(partition, uid_sw_naives)
        assert type(threshold) is float and threshold == _naive_threshold_reference(partition, uid_sw_naives)


def test_naive_threshold_p90_index_on_distinct_distances():
    # tied distances would hide an off-by-one in the order statistic
    naives = ['A' * (2 ** k - 1) + 'C' * 2 ** k + 'A' * (128 - 2 ** (k + 1) + 1) for k in range(7)]  # distance (2^j + 2^k) / 128
    hams = sorted(prf.hamming_frac(a, b) for i, a in enumerate(naives) for b in naives[i + 1:])
    assert len(set(hams)) == len(hams)
    cluster = ['u%d' % i for i in range(len(naives))]
    assert prf.estimate_naive_threshold([cluster], dict(zip(cluster, naives))) == hams[int(len(hams) * 0.90)]


@pytest.mark.parametrize('n', [2, 3, 10, 41])
@pytest.mark.parametrize('block', [1, 5, 1000])
def test_upper_triangle_pair_blocks_match_triu_indices(n, block):
    blocks = list(prf.upper_triangle_pair_blocks(n, block))
    ia, ib = np.triu_indices(n, k=1)
    assert np.array_equal(np.concatenate([a for a, _ in blocks]), ia) and np.array_equal(np.concatenate([b for _, b in blocks]), ib)
    assert all(len(a) < block + n for a, _ in blocks)


def _certified_reference(frag_of_uid, uid_muts, ej_floor):
    uids, n_certified, frag_pairs = sorted(u for u in uid_muts if u in frag_of_uid), 0, set()
    for i, a in enumerate(uids):
        for b in uids[i + 1:]:
            k = sum(1 for pos, base in uid_muts[a].items() if uid_muts[b].get(pos) == base)
            if frag_of_uid[a] == frag_of_uid[b] or k == 0:
                continue
            if k / float(len(set(uid_muts[a]) | set(uid_muts[b]))) >= ej_floor:
                n_certified += 1
                frag_pairs.add(tuple(sorted((frag_of_uid[a], frag_of_uid[b]))))
    return n_certified, frag_pairs


@pytest.mark.parametrize('block_entries', [prf.CROSS_PAIRS_BLOCK_ENTRIES, 7])
def test_certified_fragment_pairs_matches_pairwise(monkeypatch, block_entries):
    monkeypatch.setattr(prf, 'CROSS_PAIRS_BLOCK_ENTRIES', block_entries)
    rng, clusters = random.Random(5), []
    for _ in range(30):
        uid_muts = {'u%d' % i: {pos: rng.choice('ACGT') for pos in rng.sample(range(15), rng.randint(0, 6))} for i in range(rng.randint(0, 25))}
        frag_of_uid = {u: rng.randrange(4) for u in uid_muts if rng.random() > 0.1}  # a uid with no fragment is skipped
        clusters.append((frag_of_uid, uid_muts))
    certified = prf.certified_fragment_pairs(clusters, prf.EJ_SAME_FAMILY_FLOOR)  # all in one call, so clusters share a batch
    assert len(certified) == len(clusters) and sum(n for n, _ in certified) > 0
    for (frag_of_uid, uid_muts), (n_certified, frag_pairs) in zip(clusters, certified):
        assert (n_certified, set(tuple(sorted(p)) for p in frag_pairs)) == _certified_reference(frag_of_uid, uid_muts, prf.EJ_SAME_FAMILY_FLOOR)
        assert [(n_certified, frag_pairs)] == prf.certified_fragment_pairs([(frag_of_uid, uid_muts)], prf.EJ_SAME_FAMILY_FLOOR)


def test_mutations_with_base_rows_match_per_seq():
    rng = random.Random(3)
    naives = [''.join(rng.choice('ACGTN') for _ in range(rng.choice([0, 1, 30, 31]))) for _ in range(40)]
    seqs = [''.join(rng.choice('ACGTN') if rng.random() < 0.3 else c for c in naive) for naive in naives]
    per_seq = [[(i, s) for i, (s, n) in enumerate(zip(seq, naive)) if s != n and prf.N_BASE not in (s, n)] for seq, naive in zip(seqs, naives)]  # base by base
    assert sum(len(m) for m in per_seq) > 0
    assert [list(prf.get_mutations_with_base(s, n).items()) for s, n in zip(seqs, naives)] == per_seq
    for block in (1, 7, 20000):
        assert [list(m.items()) for m in prf.mutations_with_base_rows(seqs, naives, block=block)] == per_seq
    with pytest.raises(Exception, match='different lengths'):
        prf.mutations_with_base_rows(['ACG', 'AC'], ['ACG', 'ACG'])


def _naive_pairs_reference(oh, mask, idx_a, idx_b, threshold, symmetric, chunk_size):
    # every distance of each chunk pair
    for ca in range(0, len(idx_a), chunk_size):
        ra = idx_a[ca:ca + chunk_size]
        oh_a, mask_a = oh[ra], mask[ra]
        cb_start = ca if symmetric else 0
        for cb in range(cb_start, len(idx_b), chunk_size):
            rb = idx_b[cb:cb + chunk_size]
            compared = (mask_a @ mask[rb].T).astype(np.int32)
            matches = (oh_a @ oh[rb].T).astype(np.int32)
            mism = compared - matches
            dists = np.where(compared > 0, mism / np.maximum(compared, 1), 1.0)
            for pa, pb in np.argwhere(dists <= threshold):
                ia, ib = ra[pa], rb[pb]
                if symmetric and ia >= ib:
                    continue
                yield (ia, ib) if ia < ib else (ib, ia)


@pytest.mark.parametrize('threshold', [0., 0.05, 0.2, 1.])
def test_naive_pairs_below_matches_reference(threshold):
    rng = random.Random(11)
    base = ''.join(rng.choice('ACGT') for _ in range(40))
    naives = [''.join(rng.choice('ACGTN') if rng.random() < rng.choice([0.02, 0.3]) else c for c in base) for _ in range(60)]
    naives += naives[:3] + ['N' * 40]  # equal naives, and one that compares no position
    oh, mask = prf.onehot_naives(prf.encode_naives(naives), prf.N_BYTE)
    idx = rng.sample(range(len(naives)), len(naives))
    n_pairs = 0
    for idx_a, idx_b, symmetric in ((idx, idx, True), (idx[:25], idx[25:], False)):
        for chunk_size in (7, 2000):
            pairs = list(prf.naive_pairs_below(oh, mask, idx_a, idx_b, threshold, symmetric, chunk_size))
            assert pairs == list(_naive_pairs_reference(oh, mask, idx_a, idx_b, threshold, symmetric, chunk_size))
            n_pairs += len(pairs)
    assert n_pairs > 0


def _sample_null_reference(ev, n_target, naives, vdjs, naive_threshold):
    buckets = {}
    for i, (vdj, naive) in enumerate(zip(vdjs, naives)):
        if ev.fps[i][1] > 0 and vdj is not None and naive is not None:
            buckets.setdefault((vdj[0], len(naive)), []).append(i)
    keys = [k for k, v in buckets.items() if len(v) >= 2]
    rng = random.Random(prf.MERGE_NULL_SEED)
    out, n_tries, seen = [], 0, set()
    while keys and len(out) < n_target and n_tries < prf.MERGE_NULL_TRIES_PER_PAIR * n_target:
        n_tries += 1
        a, b = rng.sample(buckets[rng.choices(keys, weights=[len(buckets[k]) * (len(buckets[k]) - 1) for k in keys])[0]], 2)
        if len(ev.frags[a]) == 1 and len(ev.frags[b]) == 1 or frozenset((a, b)) in seen:
            continue
        seen.add(frozenset((a, b)))
        if prf.hamming_frac(naives[a], naives[b]) <= naive_threshold:
            continue
        out.append(ev.pvalue(a, b, cache_pair=False))
    return out


@pytest.mark.parametrize('chunk', [prf.MERGE_NULL_DRAW_CHUNK, 7])
def test_sample_null_matches_scalar_reference(monkeypatch, chunk):
    monkeypatch.setattr(prf, 'MERGE_NULL_DRAW_CHUNK', chunk)
    rng = random.Random(6)
    n_frags = 60
    frags = [['f%d-%d' % (i, k) for k in range(rng.choice([1, 1, 2, 3]))] for i in range(n_frags)]
    naives = [''.join(rng.choice('ACGTN') for _ in range(rng.choice([8, 9]))) if rng.random() > 0.05 else None for _ in range(n_frags)]
    vdjs = [(rng.choice(['v1', 'v2', 'v3']), 'd', 'j') for _ in range(n_frags)]
    lens = []
    for n_target in (1, 20, 100, 2000):  # the last wants more distinct pairs than there are, so runs out of tries
        out = prf._MergeEvidence.sample_null(_FakeEvidence(frags), n_target, naives, vdjs, 0.3)
        assert out == _sample_null_reference(_FakeEvidence(frags), n_target, naives, vdjs, 0.3)
        lens.append(len(out))
    assert lens[:3] == [1, 20, 100] and 100 < lens[3] < 2000


def test_sample_null_boundary_and_excluded_fragments():
    # a pair exactly at the threshold is not null, and fragments with no vdj or no mutated member are never drawn
    frags = [['a1', 'a2'], ['b1', 'b2'], ['c1', 'c2'], ['d1', 'd2'], ['e1', 'e2']]
    naives = ['AAAA', 'AAAC', 'AACC', 'CCCC', 'CCCC']
    vdjs = [('v1', 'd', 'j')] * 3 + [None, ('v1', 'd', 'j')]
    ev = _FakeEvidence(frags)
    ev.fps[4] = ({}, 0)
    out = prf._MergeEvidence.sample_null(ev, 100, naives, vdjs, 0.25)
    assert sorted(out) == [2]  # (0, 2) at 0.5; (0, 1) and (1, 2) sit at 0.25


def test_null_pvalue_not_cached():
    frags = [['a1', 'a2'], ['b1', 'b2']]
    muts = {'a1': {5: 'T', 8: 'G'}, 'a2': {5: 'T', 8: 'G'}, 'b1': {5: 'T', 8: 'G'}, 'b2': {5: 'T', 8: 'G'}}
    fps = [prf.get_cluster_fingerprint(f, muts) for f in frags]
    ev = prf._MergeEvidence(frags, fps, muts, lambda uid, pos, base: (0.01, None), 4)
    p = ev.pvalue(0, 1, cache_pair=False)
    assert ev.pair_p == {} and ev.pvalue(1, 0) == p and ev.pair_p == {(0, 1): p}


def test_lone_mutated_member_looked_up_with_no_mutation():
    frags, muts, seen = [['a1', 'a2'], ['b1'], ['c1', 'c2']], {'a1': {}, 'b1': {5: 'T'}, 'c1': {5: 'T'}, 'c2': {5: 'T'}}, []
    fps = [prf.get_cluster_fingerprint(f, muts) for f in frags]
    prf._MergeEvidence(frags, fps, muts, lambda uid, pos, base: (0.01, None), 5, look_up=seen.append)
    assert seen == ['a1', 'b1']


def test_merge_evidence_pvalue_matches_inline_reference():
    # reference builds the raised freqs and the shared mutations inline, as pvalue did before caching them
    rng = random.Random(9)
    frags, muts = [], {}
    for ifrag in range(10):
        base = {pos: rng.choice('AC') for pos in rng.sample(range(15), rng.randint(0, 8))}
        frags.append(['f%d-%d' % (ifrag, i) for i in range(rng.randint(1, 3))])
        for uid in frags[-1]:
            muts[uid] = {**base, rng.randrange(15, 30): 'G'}
    ufreqs = {uid: {k: 10 ** rng.uniform(-4, -0.5) for k in m.items()} for uid, m in muts.items()}
    ev = prf._MergeEvidence(frags, [prf.get_cluster_fingerprint(f, muts) for f in frags], muts, lambda uid, pos, base: (ufreqs[uid][(pos, base)], None), len(muts))
    ref_cache = prf._PerUidPvalCache(len(muts), prf.WEIGHT_GRID_NATS)
    n_below_one = 0
    for i in range(len(frags)):
        for j in range(i + 1, len(frags)):
            (ma, fa, oa), (mb, fb, ob) = ev.cons[i], ev.cons[j]
            fa = {k: max(fa.get(k, prf.NO_GERMLINE_FREQ), ev.group_rate.get(k, 0.0)) for k in ma.items()}
            fb = {k: max(fb.get(k, prf.NO_GERMLINE_FREQ), ev.group_rate.get(k, 0.0)) for k in mb.items()}
            ref = prf.weighted_shared_descent_pvalue(ma, mb, fa, fb, ev.n_seqs, oa, ob, cache=ref_cache, uid_a=i, uid_b=j) if ma and mb else 1.0
            assert ev.pvalue(j, i, cache_pair=(i + j) % 2 == 0) == ref
            n_below_one += ref < 1.0
    assert n_below_one > 0


class _FakeParter(object):
    def __init__(self, workdir):
        self.args = type('args', (), {'workdir': workdir})()


@pytest.mark.parametrize('code', [0, None, 1, 'some error'])
def test_no_sw_pass_only_on_exit_zero(tmp_path, code):
    workdir = tmp_path / 'work'
    workdir.mkdir()
    def step():
        raise SystemExit(code)
    if code == 0:
        assert ha_repartition._run_unless_no_sw_pass(step, _FakeParter(str(workdir))) is False and not workdir.exists()
    else:
        with pytest.raises(SystemExit):
            ha_repartition._run_unless_no_sw_pass(step, _FakeParter(str(workdir)))
        assert workdir.exists()
    assert ha_repartition._run_unless_no_sw_pass(lambda: None, _FakeParter(str(workdir))) is True


def _bundle_specs(sizes):
    return [{'group': {'sequence_count': n}} for n in sizes]


def test_refine_bundles_without_large_groups_match_contiguous_slices():
    for n_groups in range(1, 30):
        for n_bundles in range(1, n_groups + 1):
            size = int(np.ceil(n_groups / float(n_bundles)))
            serial, parallel = prf.refine_bundles(_bundle_specs([5] * n_groups), n_bundles, serial_min_seqs=10)
            assert serial == [] and parallel == [(s, min(size, n_groups - s)) for s in range(0, n_groups, size)]


def test_refine_bundles_large_groups_run_alone():
    rng = random.Random(11)
    assert prf.refine_bundles([], 4, 50) == ([], [])
    for _ in range(300):
        sizes = [rng.choice([1, 5, 50]) for _ in range(rng.randint(1, 25))]
        serial, parallel = prf.refine_bundles(_bundle_specs(sizes), rng.randint(1, 8), serial_min_seqs=50)
        assert serial == [(i, 1) for i, n in enumerate(sizes) if n >= 50]
        covered = sorted(i for start, count in serial + parallel for i in range(start, start + count))
        assert covered == list(range(len(sizes))) and all(sizes[i] < 50 for start, count in parallel for i in range(start, start + count))


def _snap_reference(groups):
    supported = [g for g in groups if len(g[1]) >= 2]
    if len(supported) == 0 or len(supported) == len(groups):
        return groups, 0
    out, n_snapped = dict((nv, list(mem)) for nv, mem in supported), 0
    for nv, mem in groups:
        if len(mem) >= 2:
            continue
        best = min(supported, key=lambda s: (prf.hamming_frac(nv, s[0]), -len(s[1]), s[0]))
        out[best[0]].extend(mem)
        n_snapped += len(mem)
    return list(out.items()), n_snapped


def test_snap_unsupported_matches_scalar_reference():
    rng = random.Random(7)
    for _ in range(200):
        naives = list(set(''.join(rng.choice('ACN') for _ in range(rng.choice([0, 3, 4]))) for _ in range(rng.randint(1, 12))))
        groups = [(nv, ['%s-%d' % (nv, k) for k in range(rng.choice([1, 1, 2, 3]))]) for nv in naives]
        assert prf.snap_unsupported(groups) == _snap_reference(groups)


def test_shared_passed_in_matches_computed():
    rng = random.Random(8)
    for _ in range(200):
        ma, mb = [{pos: rng.choice('AC') for pos in rng.sample(range(30), rng.randint(1, 12))} for _ in range(2)]
        fa, fb = [{k: 10 ** rng.uniform(-4, -0.5) for k in m.items()} for m in (ma, mb)]
        shared = frozenset(ma.items()) & frozenset(mb.items())
        args = (ma, mb, fa, fb, 1000, {}, {})
        assert prf.weighted_shared_descent_pvalue(*args, shared=shared) == prf.weighted_shared_descent_pvalue(*args)


def _run_pool(tmp_path, monkeypatch, groups, n_procs, memory_bytes):
    """map_groups on <groups> of (name, input bytes, peak memory it reports, group whose start it
    waits for before it ends). Returns the names in result order and, per group, the groups
    running when it started (itself among them), from the order of the start and end lines."""
    fake, log = {}, str(tmp_path / 'events')
    monkeypatch.setattr(utils, 'peak_rss', lambda who=resource.RUSAGE_SELF: fake['peak'] if who == resource.RUSAGE_SELF else 0)
    monkeypatch.setattr(utils, 'current_rss', lambda: 0)

    def note(line):
        with open(log, 'a') as f:
            f.write(line + '\n')

    def run(name, peak, wait_for):
        note('start %s' % name)
        fake['peak'] = peak
        deadline = time.time() + 10
        while wait_for is not None and 'start %s\n' % wait_for not in open(log).read() and time.time() < deadline:
            time.sleep(0.01)
        note('end %s' % name)
        return name

    infnames = []
    for name, nbytes, _, _ in groups:
        infnames.append([str(tmp_path / ('%s.in' % name))])
        with open(infnames[-1][0], 'wb') as f:
            f.truncate(nbytes)
    names = utils.map_groups(run, [(name, peak, wait_for) for name, _, peak, wait_for in groups], n_procs, infnames, memory_bytes)
    running, at_start = set(), {}
    for what, name in (l.split() for l in open(log)):
        if what == 'start':
            running.add(name)
            at_start[name] = set(running)
        else:
            running.remove(name)
    return names, at_start


def test_group_pool_stays_within_memory(tmp_path, monkeypatch, capsys):
    # a byte of memory per input byte and 700 in all: the largest alone, then largest first while they fit
    sizes = {'e': 100, 'c': 300, 'top': 1000, 'b': 400, 'over': 800, 'd': 300}
    names, at_start = _run_pool(tmp_path, monkeypatch, [(name, n, n, {'b': 'c', 'c': 'b'}.get(name)) for name, n in sizes.items()], 4, 700)
    assert names == list(sizes)
    alone = {name for name, running in at_start.items() if running == {name}}
    assert {'top', 'over'} <= alone and not any({'top', 'over'} & running for name, running in at_start.items() if name not in ('top', 'over'))
    assert {'b', 'c'} in (at_start['b'], at_start['c'])
    assert all(name in alone or sum(sizes[g] for g in running) <= 700 for name, running in at_start.items())
    out = capsys.readouterr().out
    assert out.count('runs alone') == 1 and 'over.in runs alone' in out


def test_group_pool_rate_comes_from_groups_at_least_as_large(tmp_path, monkeypatch):
    # c takes three times its bytes, which counts for e (as small) and not for b (larger, still running): 500 + 300 fit in 800
    names, at_start = _run_pool(tmp_path, monkeypatch, [('top', 1000, 1000, None), ('b', 500, 500, 'e'), ('c', 100, 300, 'b'), ('e', 100, 100, None)], 2, 800)
    assert names == ['top', 'b', 'c', 'e']
    assert at_start['e'] == {'b', 'e'} and all(len(running) <= 2 for running in at_start.values())


def test_group_pool_takes_sizes_in_place_of_file_bytes(tmp_path, monkeypatch):
    # largest by its file is c, by <sizes> a: a starts first
    monkeypatch.setattr(utils, 'current_rss', lambda: 0)
    order = str(tmp_path / 'order')
    def run(name):
        with open(order, 'a') as f:
            f.write(name + '\n')
        return name
    names = utils.map_groups(run, [('a',), ('b',), ('c',)], 3, _pool_infnames(tmp_path, [1, 2, 3]), 2**50, sizes=[30, 20, 10])
    assert names == ['a', 'b', 'c']  # results stay in group order
    assert open(order).read().split()[0] == 'a'


def test_group_pool_refuses_sizes_of_another_length(tmp_path):
    with pytest.raises(Exception, match='2 sizes for 3 groups'):
        utils.map_groups(lambda i: None, [(0,), (1,), (2,)], 2, _pool_infnames(tmp_path, [1, 1, 1]), sizes=[1, 2])


def test_group_peak_memory_counts_its_subprocess():
    def peak(fcn):
        ctx = multiprocessing.get_context('fork')
        recv, send = ctx.Pipe(duplex=False)
        proc = ctx.Process(target=utils._send_result, args=(send, fcn, ()))
        proc.start()
        assert recv.poll(60)
        _, err, peak_bytes = recv.recv()
        proc.join()
        assert err is None
        return peak_bytes
    held = 200 * 2**20
    with_subproc = peak(lambda: subprocess.check_call([sys.executable, '-c', 'x = b"x" * %d' % held]))
    assert with_subproc - peak(lambda: None) > held / 2


def _pool_infnames(tmp_path, sizes):
    infnames = [[str(tmp_path / ('%d.in' % i))] for i in range(len(sizes))]
    for (fname,), nbytes in zip(infnames, sizes):
        with open(fname, 'wb') as f:
            f.truncate(nbytes)
    return infnames


def test_group_pool_raises_what_a_group_raised(tmp_path):
    def run(i):
        if i == 1:
            raise ValueError('bad group')
        time.sleep(60 if i == 2 else 0)

    start = time.time()
    with pytest.raises(ValueError, match='bad group'):  # group 2 is running then, and is not waited for
        utils.map_groups(run, [(0,), (1,), (2,)], 3, _pool_infnames(tmp_path, [300, 200, 100]), 2**50)
    assert time.time() - start < 30


def test_group_pool_raises_on_a_process_with_no_result(tmp_path):
    with pytest.raises(Exception, match='0.in ended with exit code 3 and no result'):
        utils.map_groups(lambda i: os._exit(3), [(0,), (1,)], 2, _pool_infnames(tmp_path, [2, 1]))


def test_group_pool_raises_on_a_process_ended_partway_through_its_result(tmp_path, monkeypatch):
    def send_part(conn, fcn, args):
        os.write(conn.fileno(), struct.pack('!i', 100) + b'x')  # a 100 byte message, cut after one byte
        os._exit(3)

    monkeypatch.setattr(utils, '_send_result', send_part)
    with pytest.raises(Exception, match='0.in ended with exit code 3 and no result'):
        utils.map_groups(lambda i: None, [(0,), (1,)], 2, _pool_infnames(tmp_path, [2, 1]))


def test_job_memory(monkeypatch):
    monkeypatch.delenv('SLURM_MEM_PER_NODE', raising=False)
    monkeypatch.delenv('SLURM_MEM_PER_CPU', raising=False)
    assert utils.job_memory() == psutil.virtual_memory().total
    monkeypatch.setenv('SLURM_MEM_PER_CPU', '1000')
    monkeypatch.setattr(utils, 'n_available_cpus', lambda: 3)
    assert utils.job_memory() == 3000 * 2**20
    monkeypatch.setenv('SLURM_MEM_PER_NODE', '2048')
    assert utils.job_memory() == 2048 * 2**20
    monkeypatch.setenv('SLURM_MEM_PER_NODE', '0')  # all of the node, so not a bound
    assert utils.job_memory() == 3000 * 2**20
    monkeypatch.setenv('SLURM_MEM_PER_CPU', '0')
    assert utils.job_memory() == psutil.virtual_memory().total
    assert 2**20 < utils.peak_rss() < psutil.virtual_memory().total and 0.05 < utils.current_rss() / utils.peak_rss() < 1.5  # both in bytes
