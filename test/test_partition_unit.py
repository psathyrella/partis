import os
import random
import re

import pytest

from _helpers import HFRAC_ARGS, LOCI, PAIRED_PARAM_DIR, copy_dir, disjoint_dir, group_cmds, group_stage_path, locus_manifest, paired_partition_args, partition_fname, read_log, read_yaml, run_partis, run_partis_fails
from partis import disjointgrouper as dg
from partis import glutils
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
    def get_uid_freqs(uid):
        return {k: 0.01 for k in muts[uid].items()}, {}
    ev = prf._MergeEvidence(frags, fps, muts, get_uid_freqs, n_seqs)
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
        assert prf._conditional_pvalue_from_wb(wb, shared) == _dp_reference(wb, shared)
        assert prf._build_survival_curve(wb) == _curve_reference(wb)


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
    rng = random.Random(5)
    for _ in range(30):
        uid_muts = {'u%d' % i: {pos: rng.choice('ACGT') for pos in rng.sample(range(15), rng.randint(0, 6))} for i in range(rng.randint(0, 25))}
        frag_of_uid = {u: rng.randrange(4) for u in uid_muts if rng.random() > 0.1}  # a uid with no fragment is skipped
        n_certified, frag_pairs = prf.certified_fragment_pairs(frag_of_uid, uid_muts, prf.EJ_SAME_FAMILY_FLOOR)
        assert (n_certified, set(tuple(sorted(p)) for p in frag_pairs)) == _certified_reference(frag_of_uid, uid_muts, prf.EJ_SAME_FAMILY_FLOOR)


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


def test_null_pvalue_not_cached():
    frags = [['a1', 'a2'], ['b1', 'b2']]
    muts = {'a1': {5: 'T', 8: 'G'}, 'a2': {5: 'T', 8: 'G'}, 'b1': {5: 'T', 8: 'G'}, 'b2': {5: 'T', 8: 'G'}}
    fps = [prf.get_cluster_fingerprint(f, muts) for f in frags]
    ev = prf._MergeEvidence(frags, fps, muts, lambda uid: ({k: 0.01 for k in muts[uid].items()}, {}), 4)
    p = ev.pvalue(0, 1, cache_pair=False)
    assert ev.pair_p == {} and ev.pvalue(1, 0) == p and ev.pair_p == {(0, 1): p}


class _FakeParter(object):
    def __init__(self, workdir):
        self.args = type('args', (), {'workdir': workdir})()


@pytest.mark.parametrize('code', [0, None, 1, 'some error'])
def test_passes_sw_only_on_exit_zero(tmp_path, code):
    workdir = tmp_path / 'work'
    workdir.mkdir()
    def step():
        raise SystemExit(code)
    if code == 0:
        assert ha_repartition._passes_sw(step, _FakeParter(str(workdir))) is False and not workdir.exists()
    else:
        with pytest.raises(SystemExit):
            ha_repartition._passes_sw(step, _FakeParter(str(workdir)))
        assert workdir.exists()
    assert ha_repartition._passes_sw(lambda: None, _FakeParter(str(workdir))) is True


def test_write_singletons_result(tmp_path):
    glfo = glutils.read_glfo('%s/igh/hmm/germline-sets' % PAIRED_PARAM_DIR, 'igh')
    fname = str(tmp_path / 'singletons.yaml')
    ha_repartition._write_singletons_result(fname, glfo, ['a', 'b', 'c'])
    _, antns, cpath = utils.read_output(fname, glfo=glfo)
    assert sorted(cpath.best()) == [['a'], ['b'], ['c']] and antns == []


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
