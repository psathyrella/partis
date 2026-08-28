"""Post-partition refinement by locus-appropriate operators.

Pipeline-agnostic module: operates on a partition (list of clusters, each a list
of uids) plus per-sequence annotations and SW naives, and returns a refined
partition. Usable after any partis partition (standard or disjoint grouping).

Which operators run is decided by the locus, resolved before any of them, because
each operator is built on a signal only its locus carries. D-gene presence is the
fallback when no locus is supplied, and a disagreement between the two only warns:

  HEAVY  split_on_naive_identity    exact sw naive proposes, shared mutations veto
         merge_on_naive_similarity  naive hamming proposes, mutation fingerprint
                                    validates, junction guards
  LIGHT  split_on_shared_descent    weighted shared-descent test against a
                                    repertoire-sourced null

Heavy splits then merges; light splits only. The two columns share no operator.

Entry point: refine_partition(); its defaults run with singleton-skip and the
junction guard. Driven by the run-partition-refine-jobs action (and the
integrated --partition-refine flag), which run run_jobs() over the disjoint groups.
"""
import csv
import json
import math
import os
import random
import struct
import tempfile
import time
from collections import defaultdict
import numpy as np


def load_hmm_naives(hmm_annotation_path, sw_naives=None):
    """Load per-sequence naive from HMM annotations, stripping fv_insertion padding.
    Falls back to SW naive if HMM naive length does not match after stripping."""
    with open(hmm_annotation_path) as f:
        data = json.load(f)
    uid_naives = {}
    n_replaced = 0
    n_fallback = 0
    for evt in data['events']:
        fv_len = len(evt.get('fv_insertion', ''))
        core_naive = evt['naive_seq'][fv_len:]
        for uid in evt['unique_ids']:
            if sw_naives is not None and uid in sw_naives:
                if len(core_naive) == len(sw_naives[uid]):
                    uid_naives[uid] = core_naive
                    n_replaced += 1
                else:
                    uid_naives[uid] = sw_naives[uid]
                    n_fallback += 1
            else:
                uid_naives[uid] = core_naive
                n_replaced += 1
    print('  HMM naives: %d replaced, %d fell back to SW (length mismatch)' % (n_replaced, n_fallback))
    return uid_naives


def load_true_partition(simu_path):
    with open(simu_path) as f:
        data = json.load(f)
    families = {}
    events = data['events'] if isinstance(data, dict) else data
    for evt in events:
        if evt.get('invalid'):
            continue
        families.setdefault(evt['reco_id'], []).extend(evt['unique_ids'])
    return list(families.values())


def hamming_frac(s1, s2):
    n = min(len(s1), len(s2))
    if n == 0:
        return 1.0
    compared = 0
    diffs = 0
    for a, b in zip(s1[:n], s2[:n]):
        if a == 'N' or b == 'N':
            continue
        compared += 1
        if a != b:
            diffs += 1
    return diffs / compared if compared > 0 else 1.0


def get_mutations(seq, naive):
    if len(seq) != len(naive):  # zip would silently truncate a frame mismatch
        raise Exception('seq (%d) and naive (%d) different lengths' % (len(seq), len(naive)))
    muts = set()
    for i, (s, n) in enumerate(zip(seq, naive)):
        if s != n and s != 'N' and n != 'N':
            muts.add(i)
    return muts


def jaccard(muts1, muts2):
    if len(muts1) == 0 and len(muts2) == 0:
        return 0.0
    shared = len(muts1 & muts2)
    union = len(muts1 | muts2)
    return shared / union if union > 0 else 0.0


def get_mutations_with_base(seq, naive):
    if len(seq) != len(naive):
        raise Exception('seq (%d) and naive (%d) different lengths' % (len(seq), len(naive)))
    muts = {}
    for i, (s, n) in enumerate(zip(seq, naive)):
        if s != n and s != 'N' and n != 'N':
            muts[i] = s
    return muts


def estimate_naive_threshold(partition, uid_sw_naives):
    within_hammings = []
    for cluster in partition:
        if len(cluster) < 3:
            continue
        cluster_naives = {}
        for uid in cluster:
            if uid in uid_sw_naives:
                cluster_naives[uid] = uid_sw_naives[uid]
        # dedup by sequence, so pairs of identical naives contribute nothing: the quantile below is
        # fitted on the naives that disagree, not on every within-cluster pair
        unique_naives = list(set(cluster_naives.values()))
        if len(unique_naives) < 2:
            continue
        for i in range(len(unique_naives)):
            for j in range(i + 1, len(unique_naives)):
                h = hamming_frac(unique_naives[i], unique_naives[j])
                within_hammings.append(h)
    if len(within_hammings) == 0:
        print('  no within-cluster naive pairs, using default 0.01')
        return 0.01
    within_hammings.sort()
    n = len(within_hammings)
    p90 = within_hammings[min(int(n * 0.90), n - 1)]
    print('  within-cluster SW naive hamming (%d pairs):' % n)
    print('    median: %.4f, p90: %.4f, max: %.4f' % (
        within_hammings[n // 2], p90, within_hammings[-1]))
    print('  adaptive naive threshold (p90): %.4f' % p90)
    return p90


def get_fragment_naive(uids, uid_sw_naives):
    naive_counts = defaultdict(int)
    for uid in uids:
        if uid in uid_sw_naives:
            naive_counts[uid_sw_naives[uid]] += 1
    if not naive_counts:
        return None
    return max(naive_counts, key=naive_counts.get)


def get_fragment_vdj(uids, uid_rearr_features):
    """Modal (v, d, j) gene triple over a fragment's members, or None if unavailable."""
    counts = defaultdict(int)
    for uid in uids:
        feat = uid_rearr_features.get(uid) if uid_rearr_features is not None else None
        if feat is not None and feat.get('vdj') is not None:
            counts[feat['vdj']] += 1
    if len(counts) == 0:
        return None
    return max(counts, key=counts.get)


def get_fragment_junction(uids, uid_rearr_features):
    """Modal (v_3p_del, j_5p_del) over a fragment's members, or None if unavailable."""
    counts = defaultdict(int)
    for uid in uids:
        feat = uid_rearr_features.get(uid) if uid_rearr_features is not None else None
        if feat is None:
            continue
        key = (feat.get('v_3p_del'), feat.get('j_5p_del'))
        if key[0] is None or key[1] is None:
            continue
        counts[key] += 1
    if len(counts) == 0:
        return None
    return max(counts, key=counts.get)


def get_cluster_fingerprint(uids, uid_to_muts_with_base):
    """Build mutation fingerprint for a cluster.

    Returns (fingerprint, n_with_muts), where fingerprint is
    {position: {base: count}} over the cluster's members and n_with_muts is how
    many of them carried any mutations. Positions mutated by multiple members
    with the same base are strong fingerprint positions.
    """
    fingerprint = defaultdict(lambda: defaultdict(int))
    n_with_muts = 0
    for uid in uids:
        muts = uid_to_muts_with_base.get(uid)
        if muts is None:
            continue
        n_with_muts += 1
        for pos, base in muts.items():
            fingerprint[pos][base] += 1
    return fingerprint, n_with_muts


# floor on the strong-position vote count, so a small fragment can't trivially dominate its own fingerprint
FINGERPRINT_MIN_COUNT_FLOOR = 2


def fingerprint_strong_positions(fp, n):
    """(pos, dominant_base) pairs in <fp> where the dominant base clears the strong-position
    vote floor (mutated by at least half of <n> members, or FINGERPRINT_MIN_COUNT_FLOOR)."""
    min_count = max(FINGERPRINT_MIN_COUNT_FLOOR, n // 2)
    out = []
    for pos, bases in fp.items():
        dominant_base = max(bases, key=bases.get)
        if bases[dominant_base] >= min_count:
            out.append((pos, dominant_base))
    return out


def fingerprint_agreement(fp1, n1, fp2, n2, min_fp_positions=0):
    """Measure agreement between two cluster fingerprints.

    For each mutation position in fp1, check if fp2 has the same position
    with the same dominant base. Score = fraction of fp1's positions that
    agree with fp2, symmetrised by taking the max of the two directions, so an
    asymmetric pair passes when the small fragment's positions appear in the
    large one.

    Returns (score, winning_strong_count) where score is in [0, 1] and
    winning_strong_count is the number of strong positions on the side that
    sourced the winning direction (fwd or rev, whichever scored higher).
    High score = strong agreement (likely same family).
    Low score = independent mutations (likely different families).

    Returns score=-1 for insufficient signal: when either cluster has no member
    carrying mutations, or, if min_fp_positions > 0, when both have fewer strong
    positions than the threshold.
    """
    if n1 == 0 or n2 == 0:
        return -1.0, 0

    def directional_agreement(source_fp, source_n, target_fp):
        """Fraction of source's strong positions that appear in target."""
        strong = fingerprint_strong_positions(source_fp, source_n)
        if len(strong) == 0:
            return 0.0, 0
        n_agree = sum(1 for pos, base in strong if pos in target_fp and base in target_fp[pos])
        return n_agree / len(strong), len(strong)

    n_strong_max = max(len(fingerprint_strong_positions(fp1, n1)), len(fingerprint_strong_positions(fp2, n2)))

    # if both clusters have too few strong positions, signal is insufficient
    if min_fp_positions > 0 and n_strong_max < min_fp_positions:
        return -1.0, n_strong_max

    fwd, n_strong_fwd = directional_agreement(fp1, n1, fp2)
    rev, n_strong_rev = directional_agreement(fp2, n2, fp1)

    # asymmetric clusters (one large, one small): take the max, since the small
    # cluster's positions should appear in the large one if they are the same family
    if fwd >= rev:
        return fwd, n_strong_fwd
    return rev, n_strong_rev


def fingerprint_winning_agreement(fp1, n1, fp2, n2):
    """('fwd' or 'rev', agreeing (pos, base) pairs) for whichever direction
    fingerprint_agreement's own comparison would pick as the winner."""
    def directional(source_fp, source_n, target_fp):
        strong = fingerprint_strong_positions(source_fp, source_n)
        if len(strong) == 0:
            return 0.0, []
        agree = [(pos, base) for pos, base in strong if pos in target_fp and base in target_fp[pos]]
        return len(agree) / len(strong), agree

    fwd_score, fwd_agree = directional(fp1, n1, fp2)
    rev_score, rev_agree = directional(fp2, n2, fp1)
    if fwd_score >= rev_score:
        return 'fwd', fwd_agree
    return 'rev', rev_agree


def mutation_carrier(frag_uids, pos, base, uid_to_muts_with_base):
    """uid in <frag_uids> that carries (pos, base), or None."""
    for uid in frag_uids:
        muts = uid_to_muts_with_base.get(uid)
        if muts is not None and muts.get(pos) == base:
            return uid
    return None


MERGE_WEIGHTED_FREQ_FLOOR = 1e-4  # clamp before log10, well below NO_GERMLINE_FREQ
MERGE_WEIGHTED_SCORE_CUTOFF = 2.0


def mute_freq_weighted_score(frag1, fp1, n1, frag2, fp2, n2, uid_to_muts_with_base, get_uid_freqs):
    """Sum of -log10(population mutation frequency) over the winning direction's agreeing
    strong positions: agreement at a rare position counts for more than at a common hotspot.
    get_uid_freqs(uid) -> {(pos, base): freq}, each uid's own per-position mutation
    frequencies from the parameter directory."""
    win_side, agree = fingerprint_winning_agreement(fp1, n1, fp2, n2)
    if len(agree) == 0:
        return 0.0
    source_frag = frag1 if win_side == 'fwd' else frag2
    total = 0.0
    for pos, base in agree:
        carrier = mutation_carrier(source_frag, pos, base, uid_to_muts_with_base)
        freq = get_uid_freqs(carrier).get((pos, base), NO_GERMLINE_FREQ) if carrier else NO_GERMLINE_FREQ
        total += -math.log10(max(freq, MERGE_WEIGHTED_FREQ_FLOOR))
    return total


# shared-mutation fraction at which a cross-fragment pair certifies common descent
EJ_SAME_FAMILY_FLOOR = 0.10


def _union_find(n):
    """(find, union) over n indices, with path halving."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    return find, union


def snap_unsupported(groups):
    """Absorb each one-member naive group into the nearest corroborated (>= 2 member) group.
    Returns (groups, n_members_snapped); no-op unless some group is corroborated and some is not."""
    supported = [g for g in groups if len(g[1]) >= 2]
    if len(supported) == 0 or len(supported) == len(groups):
        return groups, 0
    out = dict((nv, list(mem)) for nv, mem in supported)
    n_snapped = 0
    for nv, mem in groups:
        if len(mem) >= 2:
            continue
        # nearest wins at any distance, and the key reads the frozen list so order cannot matter
        best = min(supported, key=lambda s: (hamming_frac(nv, s[0]), -len(s[1]), s[0]))
        out[best[0]].extend(mem)
        n_snapped += len(mem)
    return list(out.items()), n_snapped


def cross_shared_counts(frag_of_uid, uid_muts):
    """Shared base-specific mutation count for each cross-fragment uid pair sharing at least one.
    Pairs sharing none are absent."""
    inv = defaultdict(list)
    for uid, muts in uid_muts.items():
        for pos, base in muts.items():
            inv[(pos, base)].append(uid)
    pair_k = defaultdict(int)
    for uids in inv.values():
        if len(uids) < 2:
            continue
        for i in range(len(uids)):
            fi = frag_of_uid.get(uids[i])
            if fi is None:
                continue
            for j in range(i + 1, len(uids)):
                fj = frag_of_uid.get(uids[j])
                if fj is None or fi == fj:
                    continue
                a, b = uids[i], uids[j]
                pair_k[(a, b) if a < b else (b, a)] += 1
    return pair_k


# ----------------------------------------------------------------------------
# Convergence-aware veto for the heavy split. Members sharing an inferred
# rearrangement is evidence of common descent only when that rearrangement's total
# insertion length is rare enough that generating it repeatedly by independent
# convergent recombination is implausible, so the veto reads both.
# ----------------------------------------------------------------------------

# members that must share the modal rearrangement for the veto to apply
LENGTH_VETO_MIN_SHARED = 4
LENGTH_SAMPLE_DRAWS = 1000000  # rearrangements drawn to estimate the veto's cutoff
LENGTH_SAMPLE_SEED = 1  # fixed, so one parameter dir always gives one cutoff

# (name, column) for the two gene marginals and the two length tables the cutoff is
# sampled from; d/j gene draws condition len_vd/len_dj so the sampled median reflects
# the parameter dir's real length distribution rather than an unconditioned one
_LENGTH_TABLE_COLUMNS = [('d', 'd_gene'), ('j', 'j_gene'),
                         ('len_vd', 'vd_insertion'), ('len_dj', 'dj_insertion')]


def _length_table_specs():
    """(name, file, [varying column, conditioning columns...]) per table."""
    from partis import utils
    specs = []
    for name, column in _LENGTH_TABLE_COLUMNS:
        cols = [column] + utils.column_dependencies[column]
        specs.append((name, utils.get_parameter_fname(column_and_deps=cols), cols))
    return specs


LENGTH_TABLE_SPECS = _length_table_specs()


def read_length_tables(parameter_dir):
    """Read the d/j gene marginals and the len_vd/len_dj tables from <parameter_dir>/hmm (the
    locus-level parameter dir) and normalise each within its conditioning variable. Returns
    {name: {(value, conditioned-on...): probability}}, or None if <parameter_dir> is None.
    Raises if a dir was passed but a table is missing or unreadable."""
    if parameter_dir is None:
        return None
    tdir = '%s/hmm' % parameter_dir
    if not os.path.isdir(tdir):  # also accept the hmm dir itself
        tdir = parameter_dir
    tables = {}
    for name, fname, cols in LENGTH_TABLE_SPECS:
        counts = defaultdict(dict)
        try:
            with open('%s/%s' % (tdir, fname)) as tfile:
                for row in csv.DictReader(tfile):
                    counts[tuple(row[c] for c in cols[1:])][row[cols[0]]] = float(row['count'])
        except (IOError, OSError, KeyError, ValueError) as terr:
            raise Exception('couldn\'t read length table %s/%s (%s)' % (tdir, fname, terr))
        tables[name] = {}
        for cond, cfo in counts.items():
            tot = sum(cfo.values())
            for val, count in cfo.items():
                tables[name][(val,) + cond] = count / tot if tot > 0 else 0.
    return tables


def get_rearrangement(feat):
    """The full rearrangement (v, d, j, four deletions, two insertion lengths) from one uid's
    rearrangement features, or None if any of them is missing."""
    if feat is None or feat.get('vdj') is None:
        return None
    vals = [feat.get(k) for k in ('v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del', 'len_vd', 'len_dj')]
    if any(v is None for v in vals):
        return None
    return tuple(feat['vdj']) + tuple(vals)


def get_cluster_rearrangement(uids, uid_rearr_features):
    """Modal full rearrangement over a cluster's members and how many carry it, or (None, 0)."""
    counts = defaultdict(int)
    for uid in uids:
        rearr = get_rearrangement(uid_rearr_features.get(uid) if uid_rearr_features is not None else None)
        if rearr is not None:
            counts[rearr] += 1
    if len(counts) == 0:
        return None, 0
    modal = max(sorted(counts), key=counts.get)  # sorted so ties don't depend on dict order
    return modal, counts[modal]


def _sampling_choices(table):
    """{conditioned-on: (values, probabilities)} for weighted sampling from a normalised table."""
    by_cond = defaultdict(list)
    for key, prob in table.items():
        by_cond[key[1:]].append((key[0], prob))
    choices = {}
    for cond, vfo in by_cond.items():
        probs = np.array([p for _, p in vfo])
        if probs.sum() <= 0:  # nothing to draw, and normalising would divide by zero
            continue
        choices[cond] = ([v for v, _ in vfo], probs / probs.sum())  # p= wants an exact sum
    return choices


def _index_groups(picks):
    """{drawn value: index array} over a sampled column, i.e. which draws condition on what."""
    order = np.argsort(picks, kind='stable')
    spicks = picks[order]
    edges = np.flatnonzero(np.concatenate(([True], spicks[1:] != spicks[:-1], [True])))
    return dict((spicks[edges[i]], order[edges[i] : edges[i + 1]]) for i in range(len(edges) - 1))


def _draw_conditioned(rng, choices, groups, genes, n_draws):
    """One draw per index in <groups> ({gene index: draw indices}), from the length
    distribution conditioned on that gene. A gene the table never saw draws length zero."""
    vals = np.zeros(n_draws)
    for igene, idxs in groups.items():
        cond = (genes[igene],)
        if cond not in choices:
            continue
        cvals, cprobs = choices[cond]
        picks = rng.choice(len(cvals), size=len(idxs), p=cprobs)
        vals[idxs] = np.asarray(cvals, dtype=float)[picks]
    return vals


def sample_length_median(tables, n_draws=LENGTH_SAMPLE_DRAWS, seed=LENGTH_SAMPLE_SEED):
    """Median total insertion length (len_vd + len_dj) of rearrangements drawn from <tables>,
    i.e. of the rearrangement distribution the parameter dir itself defines, so it does not
    depend on how the input was clustered. d and j gene are drawn first so len_vd/len_dj come
    from their real gene-conditioned distributions rather than an unconditioned one.
    Deterministic given <seed>."""
    rng = np.random.RandomState(seed)
    choices = dict((name, _sampling_choices(tables[name])) for name, _, _ in LENGTH_TABLE_SPECS)
    genes, groups = {}, {}
    for reg in ('d', 'j'):
        genes[reg], cprobs = choices[reg][()]
        picks = rng.choice(len(genes[reg]), size=n_draws, p=cprobs)
        groups[reg] = _index_groups(picks)
    lens_total = np.zeros(n_draws)
    for name, reg in [('len_vd', 'd'), ('len_dj', 'j')]:
        lens_total += _draw_conditioned(rng, choices[name], groups[reg], genes[reg], n_draws)
    return float(np.median(lens_total))


_length_veto_cache = {}  # parameter dir -> (tables, cutoff), since every group resolves the same dir


def length_veto_inputs(parameter_dir):
    """(tables, length cutoff) for the heavy split's veto, derived from <parameter_dir> alone.
    Both None, with a warning, when no dir was passed; an unreadable dir raises rather than
    warns. Cached, so a locus samples its cutoff once."""
    if parameter_dir in _length_veto_cache:
        return _length_veto_cache[parameter_dir]
    from partis import utils
    tables = read_length_tables(parameter_dir)
    if tables is None:
        print('  %s no parameter dir passed, so the heavy split\'s length veto is off'
              % utils.wrnstr(), flush=True)
        cutoff = None
    else:
        cutoff = sample_length_median(tables)
        print('  length veto: cutoff (model median over %d draws) len_vd + len_dj >= %.1f' % (LENGTH_SAMPLE_DRAWS, cutoff), flush=True)
    _length_veto_cache[parameter_dir] = (tables, cutoff)
    return tables, cutoff


def validate_length_tables(parameter_dir):
    """Read the length veto's tables and discard them, to check the parameter dir is complete."""
    read_length_tables(parameter_dir)


def split_on_naive_identity(partition, uid_sw_naives, uid_muts_sw, min_cluster_size=2,
                            ej_floor=EJ_SAME_FAMILY_FLOOR, uid_rearr_features=None,
                            length_cutoff=None,
                            length_veto_min_shared=LENGTH_VETO_MIN_SHARED):
    """Heavy split: exact sw-naive identity proposes the split, shared mutations veto it.

    Members are grouped by exact per-sequence sw naive, one-member groups are absorbed into the
    nearest corroborated one, then fragments are re-merged when a cross-fragment pair's enhanced
    jaccard reaches ej_floor. uid_muts_sw is mutations against each sequence's own sw naive.

    length_veto_min_shared: keep a cluster whole, unsplit, when this many members share its modal
    rearrangement and that rearrangement's total insertion length (len_vd + len_dj) is at or
    above <length_cutoff>, i.e. long enough that reproducing it by independent convergent
    recombination is implausible. Keyed on the rearrangement rather than on the naive, so
    convergently-similar naives do not trip it. Disabled, with a log line, when no
    <length_cutoff> was derived, i.e. no parameter dir was passed.
    """
    result = []
    ctr = defaultdict(int)
    length_veto_on = (length_cutoff is not None
                      and length_veto_min_shared is not None and length_veto_min_shared > 0)

    def length_vetoed(uids):
        if not length_veto_on:
            return False
        rearr, n_shared = get_cluster_rearrangement(uids, uid_rearr_features)
        if rearr is None or n_shared < length_veto_min_shared:
            return False
        return rearr[7] + rearr[8] >= length_cutoff  # len_vd + len_dj

    for cluster in partition:
        if len(cluster) < min_cluster_size:
            result.append(list(cluster))
            ctr['below_min_size'] += 1
            continue

        if length_vetoed(cluster):
            result.append(list(cluster))
            ctr['length_veto'] += 1
            continue

        uid_list = list(cluster)
        cluster_naives = dict((u, uid_sw_naives[u]) for u in uid_list if u in uid_sw_naives)

        if len(set(cluster_naives.values())) <= 1:  # nothing to propose
            result.append(list(cluster))
            ctr['single_naive'] += 1
            continue
        if len(cluster_naives) < len(uid_list):  # a member with no naive cannot be grouped, so fail closed
            result.append(list(cluster))
            ctr['missing_naive'] += 1
            continue

        groups = defaultdict(list)
        for uid in uid_list:
            groups[cluster_naives[uid]].append(uid)
        proposed, n_snapped = snap_unsupported(list(groups.items()))
        ctr['snapped'] += n_snapped
        if len(proposed) <= 1:
            result.append(list(cluster))
            ctr['no_proposal'] += 1
            continue

        frags = [mem for _, mem in proposed]
        ffind, funion = _union_find(len(frags))
        frag_of_uid = dict((u, i) for i, f in enumerate(frags) for u in f)
        muts_sw = dict((u, uid_muts_sw[u]) for u in uid_list if u in uid_muts_sw)
        npos = dict((u, set(m)) for u, m in muts_sw.items())
        for (a, b), k in cross_shared_counts(frag_of_uid, muts_sw).items():
            if k / float(len(npos[a] | npos[b])) >= ej_floor:  # enhanced jaccard: shared over union
                funion(frag_of_uid[a], frag_of_uid[b])
                ctr['certified_pairs'] += 1

        validated = defaultdict(list)
        for i in range(len(frags)):
            validated[ffind(i)].extend(frags[i])
        groups_out = list(validated.values())
        if len(groups_out) > 1:
            ctr['accepted'] += 1
            result.extend(groups_out)
        else:
            ctr['rejected'] += 1
            result.append(list(cluster))

    # held counts clusters the veto touched, not ones whose outcome it changed
    length_label = ', %d held (length veto)' % ctr['length_veto'] if ctr['length_veto'] > 0 else ''
    print('  naive-identity split: %d accepted, %d rejected (veto), %d skipped (single naive), %d skipped (no proposal), %d rejected (missing naive)%s' % (
        ctr['accepted'], ctr['rejected'], ctr['single_naive'], ctr['no_proposal'], ctr['missing_naive'], length_label), flush=True)
    print('  %d members snapped, %d pairs certified, %d -> %d clusters' % (
        ctr['snapped'], ctr['certified_pairs'], len(partition), len(result)), flush=True)
    return result


# fragment size at which the modal v/d/j call is trusted enough to override a junction mismatch
VDJ_OVERRIDE_MIN_FRAG = 20


def onehot_naives(arr, n_byte):
    """(n, len) uint8 naives -> (one-hot over the non-N symbols present, non-N mask), both float32.

    One-hot inner product counts matching non-N positions, mask inner product counts compared
    positions, so N-aware hamming is a matrix product. Counts are integral in float32.
    """
    symbols = [s for s in np.unique(arr) if s != n_byte]
    n_seqs, seq_len = arr.shape
    oh = np.zeros((n_seqs, seq_len * len(symbols)), dtype=np.float32)
    for si, sym in enumerate(symbols):
        oh[:, si * seq_len:(si + 1) * seq_len] = (arr == sym)
    return oh, (arr != n_byte).astype(np.float32)


def naive_pairs_below(oh, mask, idx_a, idx_b, threshold, symmetric, chunk_size=2000):
    """Yield index pairs whose N-aware naive hamming frac is at or below <threshold>, blocked to
    bound memory. symmetric means idx_a and idx_b are the same set, so only one triangle is walked
    and each pair comes out once; otherwise the two sets must be disjoint."""
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


def merge_on_naive_similarity(split_partition, uid_info, uid_sw_naives,
                              uid_to_muts_with_base, naive_threshold,
                              min_agreement=0.15, min_fp_positions=0,
                              min_weighted_score=MERGE_WEIGHTED_SCORE_CUTOFF,
                              skip_singleton_merge=False, uid_rearr_features=None,
                              junction_guard=True, vdj_override_min=VDJ_OVERRIDE_MIN_FRAG,
                              mute_freq_dir=None, uid_part_antns=None, glfo=None):
    """Heavy merge: incremental naive merge with fingerprint validation.

    For each CDR3 group, find clusters with similar naives (candidates),
    then only merge if their mutation fingerprints agree. This prevents
    false merges of unrelated sequences with similar naives.

    Naive distances are computed only within groups that could survive the junction guard, and as a
    matrix product rather than elementwise. Both leave the partition unchanged; the candidate and
    junction-rejection counters drop, since blocked pairs are no longer examined.

    If skip_singleton_merge=True, skip merge attempts where BOTH clusters
    are singletons. True singletons (survived vsearch + HA) should stay
    singletons. Only attempt merges when at least one side has multiple
    members (i.e., one side has a real cluster-level fingerprint).

    junction_guard: reject a union whose two fragments disagree on (v_3p_del,
    j_5p_del). Rearrangement signal, orthogonal to the naive-distance proposer
    and the fingerprint validator. Inert when uid_rearr_features is None.

    vdj_override_min: rescue a junction-rejected union when both fragments share
    a v/d/j triple and the larger one has at least this many members, i.e. only
    where the modal annotation is trustworthy. 0 disables the override.

    mute_freq_dir, uid_part_antns, glfo: population-level per-position mutation frequency
    table, plus the annotations needed to look each uid up in it. Used to weight agreeing
    positions by rarity when deciding whether to accept a merge.
    """
    missing = [n for n, v in (('mute_freq_dir', mute_freq_dir), ('uid_part_antns', uid_part_antns), ('glfo', glfo)) if v is None]
    if missing:
        raise Exception('merge_on_naive_similarity missing %s' % ', '.join(missing))

    cdr3_frags = defaultdict(list)
    for cluster in split_partition:
        cdr3_len = None
        for uid in cluster:
            if uid in uid_info:
                cdr3_len = uid_info[uid]['cdr3_length']
                break
        cdr3_frags[cdr3_len].append(cluster)

    merged_partition = []
    n_naive_candidates = 0
    n_accepted = 0
    n_rejected_fingerprint = 0
    n_rejected_insufficient = 0
    n_rejected_junction = 0
    n_vdj_override = 0
    n_skipped_singleton = 0
    n_skipped_difflen = 0

    uid_freq_cache = {}

    def get_uid_freqs(uid):
        if uid not in uid_freq_cache:
            freqs, _obs = uid_param_dir_freqs(uid, uid_to_muts_with_base.get(uid, {}), uid_part_antns, glfo, mute_freq_dir)
            uid_freq_cache[uid] = freqs
        return uid_freq_cache[uid]

    for cdr3_len, frags in cdr3_frags.items():
        if len(frags) < 2:
            merged_partition.extend(frags)
            continue

        # precompute fingerprints and naives for all fragments
        frag_naives = [get_fragment_naive(f, uid_sw_naives) for f in frags]
        frag_juncs = ([get_fragment_junction(f, uid_rearr_features) for f in frags]
                      if junction_guard else [None] * len(frags))
        frag_vdjs = ([get_fragment_vdj(f, uid_rearr_features) for f in frags]
                     if junction_guard and vdj_override_min > 0 else [None] * len(frags))
        frag_fps = []
        for f in frags:
            fp, n_w = get_cluster_fingerprint(f, uid_to_muts_with_base)
            frag_fps.append((fp, n_w))

        frag_sizes = [len(f) for f in frags]

        # group fragments by naive to reduce comparisons
        naive_to_frags = defaultdict(list)
        for i, naive in enumerate(frag_naives):
            if naive is not None:
                naive_to_frags[naive].append(i)

        n = len(frags)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        def junction_verdict(i, j):
            # '' to pass, 'block' to reject on a junction mismatch, 'rescue' when a shared
            # vdj on a large fragment overrides that mismatch
            ji, jj = frag_juncs[i], frag_juncs[j]
            if ji is None or jj is None or ji == jj:
                return ''
            if (frag_vdjs[i] is not None and frag_vdjs[i] == frag_vdjs[j]
                    and max(frag_sizes[i], frag_sizes[j]) >= vdj_override_min):
                return 'rescue'
            return 'block'

        def examine(i, j):
            # one naive-similar fragment pair through the junction guard and the fingerprint validator
            nonlocal n_naive_candidates, n_accepted, n_rejected_fingerprint
            nonlocal n_rejected_insufficient, n_rejected_junction, n_vdj_override, n_skipped_singleton
            if find(i) == find(j):
                return
            # skip singleton-singleton merges: both survived vsearch+HA
            # as singletons, likely true singletons not fragments
            if skip_singleton_merge and frag_sizes[i] == 1 and frag_sizes[j] == 1:
                n_skipped_singleton += 1
                return
            n_naive_candidates += 1
            verdict = junction_verdict(i, j)
            if verdict == 'block':
                n_rejected_junction += 1
                return
            if verdict == 'rescue':
                n_vdj_override += 1
            fp_i, n_i = frag_fps[i]
            fp_j, n_j = frag_fps[j]
            agreement, _ = fingerprint_agreement(
                fp_i, n_i, fp_j, n_j, min_fp_positions)
            if agreement < 0:
                n_rejected_insufficient += 1
            elif agreement < min_agreement:
                n_rejected_fingerprint += 1
            else:
                weighted = mute_freq_weighted_score(
                    frags[i], fp_i, n_i, frags[j], fp_j, n_j, uid_to_muts_with_base, get_uid_freqs)
                if weighted >= min_weighted_score:
                    union(i, j)
                    n_accepted += 1
                else:
                    n_rejected_fingerprint += 1

        # within same-naive bucket: always naive-similar, just check fingerprint
        for naive_seq, indices in naive_to_frags.items():
            for ii in range(len(indices)):
                for jj in range(ii + 1, len(indices)):
                    examine(indices[ii], indices[jj])

        # cross-bucket: vectorized pairwise hamming to find near-match naives
        unique_naives = list(naive_to_frags.keys())
        if len(unique_naives) >= 2:
            seq_len = len(unique_naives[0])
            same_len = [n for n in unique_naives if len(n) == seq_len]
            diff_len = [n for n in unique_naives if len(n) != seq_len]

            if len(same_len) >= 2:
                arr = np.frombuffer(''.join(same_len).encode(), dtype=np.uint8).reshape(len(same_len), seq_len)
                oh, mask = onehot_naives(arr, ord('N'))

                def examine_naive_pair(ni_idx, nj_idx):
                    for i in naive_to_frags[same_len[ni_idx]]:
                        for j in naive_to_frags[same_len[nj_idx]]:
                            examine(i, j)

                # only pairs that could survive the junction guard get a distance: same key, either
                # side missing a key, or a shared vdj triple where the rescue is live
                naive_keys, wildcard = [], []
                key_groups, vdj_groups = defaultdict(list), defaultdict(list)
                vdj_live = set()
                if vdj_override_min > 0:
                    for vdj, size in zip(frag_vdjs, frag_sizes):
                        if vdj is not None and size >= vdj_override_min:
                            vdj_live.add(vdj)
                for idx, naive_seq in enumerate(same_len):
                    keys, vdjs = set(), set()
                    for i in naive_to_frags[naive_seq]:
                        if frag_juncs[i] is None:
                            keys = None
                            break
                        keys.add(frag_juncs[i])
                        if frag_vdjs[i] in vdj_live:
                            vdjs.add(frag_vdjs[i])
                    naive_keys.append(keys)
                    if keys is None:
                        wildcard.append(idx)
                    else:
                        for key in keys:
                            key_groups[key].append(idx)
                        for vdj in vdjs:
                            vdj_groups[vdj].append(idx)

                for group in key_groups.values():
                    if len(group) >= 2:
                        for ni_idx, nj_idx in naive_pairs_below(oh, mask, group, group, naive_threshold, True):
                            examine_naive_pair(ni_idx, nj_idx)
                # a keyless naive passes the guard against anything, so it needs all of same_len
                if len(wildcard) > 0:
                    keyed = [idx for idx in range(len(same_len)) if naive_keys[idx] is not None]
                    for ni_idx, nj_idx in naive_pairs_below(oh, mask, wildcard, wildcard, naive_threshold, True):
                        examine_naive_pair(ni_idx, nj_idx)
                    if len(keyed) > 0:
                        for ni_idx, nj_idx in naive_pairs_below(oh, mask, wildcard, keyed, naive_threshold, False):
                            examine_naive_pair(ni_idx, nj_idx)
                # rescue-eligible pairs, minus the same-key ones already done above
                for group in vdj_groups.values():
                    if len(group) >= 2:
                        for ni_idx, nj_idx in naive_pairs_below(oh, mask, group, group, naive_threshold, True):
                            if len(naive_keys[ni_idx] & naive_keys[nj_idx]) > 0:
                                continue
                            examine_naive_pair(ni_idx, nj_idx)

            # handle different-length naives with Python fallback (rare)
            for n1 in diff_len:
                for n2_seq in unique_naives:
                    if n2_seq <= n1:
                        continue
                    if len(n2_seq) != len(n1):  # never compared, so never merged: count so it is not silent
                        n_skipped_difflen += 1
                        continue
                    if hamming_frac(n1, n2_seq) > naive_threshold:
                        continue
                    for i in naive_to_frags[n1]:
                        for j in naive_to_frags[n2_seq]:
                            examine(i, j)

        components = defaultdict(list)
        for i in range(n):
            components[find(i)].extend(frags[i])

        groups = list(components.values())
        merged_partition.extend(groups)

    n_merges = len(split_partition) - len(merged_partition)
    skip_label = ', %d skipped (singleton-singleton)' % n_skipped_singleton if n_skipped_singleton > 0 else ''
    skip_label += ', %d skipped (naive length mismatch)' % n_skipped_difflen if n_skipped_difflen > 0 else ''
    junc_label = ', %d rejected (junction)' % n_rejected_junction if n_rejected_junction > 0 else ''
    junc_label += ', %d rescued (vdj)' % n_vdj_override if n_vdj_override > 0 else ''
    print('  naive-similarity merge: %d naive candidates, %d accepted, %d rejected (fingerprint)%s%s' % (
        n_naive_candidates, n_accepted, n_rejected_fingerprint, junc_label, skip_label))
    print('  %d merges, %d -> %d clusters' % (n_merges, len(split_partition), len(merged_partition)))
    return merged_partition


# ----------------------------------------------------------------------------
# Light link test: weighted shared descent. Each shared mutation is weighted by
# -log of its frequency in the refine input, and the null is conditioned on one
# member's mutation set, so it stays defined for a cluster of any size.
# ----------------------------------------------------------------------------

WEIGHT_GRID_NATS = 0.25  # bin width for the exact tail dp
FREQ_SMOOTH_COUNT = 0.5  # count floor for an unobserved (position, base)
_TINY = 1e-300

# derive_bin_alpha's fdr targets and null-tail fit constants
WEIGHTED_DESCENT_FDR_PRIMARY = 0.05
WEIGHTED_DESCENT_FDR_RETRY = 0.25
WEIGHTED_DESCENT_FIT_FLOOR = 1e-3  # p above this is signal-free bulk, fitted then extrapolated inwards
WEIGHTED_DESCENT_MIN_TAIL_N = 10  # degenerate-fit floor, not a trust threshold
WEIGHTED_DESCENT_SE_Z = 1.0  # rate_hat - z*SE, SE = rate_hat/sqrt(n)
WEIGHTED_DESCENT_RATE_FLOOR = 1e-6
WEIGHTED_DESCENT_CAND = [10 ** (-e / 2.0) for e in range(40, 1, -1)]  # log-spaced candidates, tightest first


def _fit_null_tail_se(sample, z=WEIGHTED_DESCENT_SE_Z, min_n=WEIGHTED_DESCENT_MIN_TAIL_N):
    """Exponential fit to u = -log10(p) over the signal-free bulk, rate pulled down by z*SE.
    Returns (rate, u0, frac), or None below <min_n> tail points or on a degenerate fit."""
    u0 = -math.log10(WEIGHTED_DESCENT_FIT_FLOOR)
    tail = [-math.log10(p) for p in sample if 0 < p < WEIGHTED_DESCENT_FIT_FLOOR]
    n = len(tail)
    if n < min_n:
        return None
    frac = n / float(len(sample))
    mean_excess = sum(u - u0 for u in tail) / n
    if mean_excess <= 0:
        return None
    rate_hat = 1.0 / mean_excess
    se = rate_hat / math.sqrt(n)
    rate = max(rate_hat - z * se, WEIGHTED_DESCENT_RATE_FLOOR)
    return rate, u0, frac


def _expected_null_below(cut, fit, n_total):
    rate, u0, frac = fit
    u = -math.log10(cut)
    if u <= u0:
        return None  # inside the fitted region, so the fit says nothing useful
    return n_total * frac * math.exp(-rate * (u - u0))


def _fdr_cutoff(target, fit, counts, n_total):
    """Loosest candidate in WEIGHTED_DESCENT_CAND whose extrapolated null share stays under <target>."""
    best = None
    for cut, cnt in zip(WEIGHTED_DESCENT_CAND, counts):
        exp_null = _expected_null_below(cut, fit, n_total)
        if exp_null is None or cnt == 0:
            continue
        if min(1.0, exp_null / float(cnt)) <= target:  # keep scanning: loosest passing cutoff wins
            best = cut
    return best


def derive_bin_alpha(sample, counts, n_total):
    """Fit the null tail, then WEIGHTED_DESCENT_FDR_PRIMARY with one retry at WEIGHTED_DESCENT_FDR_RETRY. Returns
    (alpha, diag), alpha None if the fit fails or neither target is reachable."""
    fit = _fit_null_tail_se(sample)
    if fit is None:
        return None, {'fit_rate': None, 'fit_frac': None, 'rule': None}
    diag = {'fit_rate': fit[0], 'fit_frac': fit[2]}
    for target in (WEIGHTED_DESCENT_FDR_PRIMARY, WEIGHTED_DESCENT_FDR_RETRY):
        cutoff = _fdr_cutoff(target, fit, counts, n_total)
        if cutoff is not None:
            return cutoff, dict(diag, rule='fdr:%.3g' % target)
    return None, dict(diag, rule=None)


_PAIR_BLOCK_HEADER = '<I'  # n_pairs in this cluster's block
_PAIR_RECORD = '<iif'  # (i, o): positions in that cluster's own member list; p: float32
_PAIR_RECORD_SIZE = struct.calcsize(_PAIR_RECORD)


def _write_pair_block(spill_f, pairs_here):
    """Append one cluster's block (count header, then that many (i, o, p) records) to the
    bin's spill file. Always writes a block, even count 0, to stay in lockstep with
    clusters_uid_muts for the sequential re-read."""
    spill_f.write(struct.pack(_PAIR_BLOCK_HEADER, len(pairs_here)))
    for i, o, p in pairs_here:
        spill_f.write(struct.pack(_PAIR_RECORD, i, o, p))


def _read_pair_block(spill_f):
    """Read back one written block. Returns {(i, o): p}, or None for an empty
    block (caller falls back to split_by_weighted_descent's own recompute)."""
    n_pairs, = struct.unpack(_PAIR_BLOCK_HEADER, spill_f.read(struct.calcsize(_PAIR_BLOCK_HEADER)))
    if n_pairs == 0:
        return None
    buf = spill_f.read(_PAIR_RECORD_SIZE * n_pairs)
    return {(i, o): p for i, o, p in struct.iter_unpack(_PAIR_RECORD, buf)}


def _scan_bin_full_pairs(clusters_uid_muts, uid_part_antns, glfo, mute_freq_dir, n_seqs, spill_f,
                          counts=None, sample_cap=2000000, seed=1):
    """Full within-cluster pairwise scan across every cluster in one bin, one shared
    _PerUidPvalCache. Every p < 1.0 pair is appended to <spill_f> one cluster-block at a
    time, never held past its own cluster. Returns (sample, cand_counts,
    n_total): sample is a reservoir sample of p for the null-tail fit; cand_counts is the exact
    (unsampled) count of p below each WEIGHTED_DESCENT_CAND candidate."""
    rng = random.Random(seed)
    cache = _PerUidPvalCache(n_seqs, WEIGHT_GRID_NATS)
    freq_cache, obs_cache = {}, {}
    sample, n_seen, n_total = [], 0, 0
    cand_counts = [0] * len(WEIGHTED_DESCENT_CAND)
    for cluster, uid_muts in clusters_uid_muts:
        members = list(cluster)
        n = len(members)
        if n < 2:
            _write_pair_block(spill_f, [])
            continue
        order = sorted(range(n), key=lambda i: -len(uid_muts.get(members[i], {}) or {}))
        for u in members:
            if u not in freq_cache:
                freq_cache[u], obs_cache[u] = uid_param_dir_freqs(
                    u, uid_muts.get(u, {}), uid_part_antns, glfo, mute_freq_dir, counts=counts)
        pairs_here = []
        for idx, i in enumerate(order):
            mi = uid_muts.get(members[i])
            if not mi:
                continue
            ui = members[i]
            for o in order[idx + 1:]:
                mo = uid_muts.get(members[o])
                if not mo:
                    continue
                uo = members[o]
                p = weighted_shared_descent_pvalue(
                    mi, mo, freq_cache[ui], freq_cache[uo], n_seqs,
                    n_obs_a=obs_cache[ui], n_obs_b=obs_cache[uo],
                    cache=cache, uid_a=ui, uid_b=uo)
                if p >= 1.0:
                    continue
                pairs_here.append((i, o, p))
                n_total += 1
                for k, cut in enumerate(WEIGHTED_DESCENT_CAND):
                    if p < cut:
                        cand_counts[k] += 1
                if len(sample) < sample_cap:
                    sample.append(p)
                else:
                    j = rng.randint(0, n_seen)
                    if j < sample_cap:
                        sample[j] = p
                n_seen += 1
        _write_pair_block(spill_f, pairs_here)
    return sample, cand_counts, n_total


def _weight_bins(muts, freqs, n_seqs, grid, n_obs):
    """{(pos, base): (frequency, surprisal in whole grid bins)} for one sequence's mutations.
    Each weight is binned rather than the sums, which is not a bound in either direction: the
    dp returns the exact tail of the rounded statistic, which can sit above or below the tail
    of the unrounded one.

    n_obs: per-(pos, base) observation count from the per-gene mute-freqs table, sets that
    key's floor. A no-signal position is absent from n_obs, but its freq is already
    NO_GERMLINE_FREQ there, well above any floor n_seqs could produce as the .get() default."""
    out = {}
    for pos, base in muts.items():
        p = freqs.get((pos, base), 0.0)
        pfloor = FREQ_SMOOTH_COUNT / max(n_obs.get((pos, base), n_seqs), 1)
        if p < pfloor:
            p = pfloor
        elif p > 1.0:
            p = 1.0
        nbin = int(round(-math.log(p) / grid))
        out[(pos, base)] = (p, max(nbin, 0))
    return out


def _conditional_pvalue_from_wb(wb, shared):
    """Same dp as _conditional_pvalue, given an already-built <wb> instead of rebuilding it."""
    t_obs = sum(wb[pb][1] for pb in shared if pb in wb)
    if t_obs <= 0:
        return 1.0
    dist = [0.0] * t_obs
    dist[0] = 1.0
    tail = 0.0
    for p, nbin in wb.values():
        if nbin == 0:
            continue  # a mutation everyone carries is no evidence either way
        nxt = [0.0] * t_obs
        for b, mass in enumerate(dist):
            if mass == 0.0:
                continue
            nxt[b] += mass * (1.0 - p)
            hit = b + nbin
            if hit >= t_obs:
                tail += mass * p
            else:
                nxt[hit] += mass * p
        dist = nxt
    return min(max(tail, 0.0), 1.0)


def _conditional_pvalue(muts_cond, shared, freqs, n_seqs, grid, n_obs):
    """P(T >= T_obs), where T sums the surprisals of whichever of <muts_cond> another sequence
    carries independently. Exact by dp over binned surprisal, with everything at or above the
    observed value collected into a tail bucket."""
    wb = _weight_bins(muts_cond, freqs, n_seqs, grid, n_obs=n_obs)
    return _conditional_pvalue_from_wb(wb, shared)


def _build_survival_curve(wb):
    """Untruncated dp over <wb>, reduced to the survival function S(t) = P(T >= t)."""
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
    max_t = len(dist) - 1
    survival = [0.0] * (max_t + 1)
    running = 0.0
    for t in range(max_t, -1, -1):
        running += dist[t]
        survival[t] = running
    return survival, max_t


def _survival_lookup(survival, max_t, t_obs):
    if t_obs <= 0:
        return 1.0
    if t_obs > max_t:
        return 0.0
    return min(max(survival[t_obs], 0.0), 1.0)


class _PerUidPvalCache(object):
    """Per-uid weight-bins cache, switching a uid to a full survival curve once its own
    pair count passes break-even."""

    def __init__(self, n_seqs, grid):
        self.n_seqs = n_seqs
        self.grid = grid
        self.wb = {}
        self.curve = {}
        self.use_curve = {}
        self.trunc_cost = {}

    def _get_wb(self, uid, muts, freqs, n_obs):
        wb = self.wb.get(uid)
        if wb is None:
            wb = _weight_bins(muts, freqs, self.n_seqs, self.grid, n_obs=n_obs)
            self.wb[uid] = wb
        return wb

    def conditional_pvalue(self, uid, muts, freqs, n_obs, shared):
        wb = self._get_wb(uid, muts, freqs, n_obs)
        if not wb:
            return 1.0
        if self.use_curve.get(uid):
            survival, max_t = self.curve[uid]
            t_obs = sum(wb[pb][1] for pb in shared if pb in wb)
            return _survival_lookup(survival, max_t, t_obs)
        p = _conditional_pvalue_from_wb(wb, shared)
        t_obs = sum(wb[pb][1] for pb in shared if pb in wb)
        cost = self.trunc_cost.get(uid, 0) + t_obs
        self.trunc_cost[uid] = cost
        max_t = sum(nbin for _, nbin in wb.values() if nbin > 0)
        if cost > max_t:
            self.curve[uid] = _build_survival_curve(wb)
            self.use_curve[uid] = True
        return p


# ----------------------------------------------------------------------------
# Germline mute-freq lookups from the parameter directory: each uid is scored against
# its own V/J gene's per-position table, not one locus-wide statistic.
# ----------------------------------------------------------------------------

_PARAM_MUTE_FREQ_CACHE = {}  # (mute_freq_dir, gene) -> {pos: {base: freq, _N_OBS_KEY: n}} or None
_N_OBS_KEY = 'n_obs'  # not a base, so no collision with the A/C/G/T keys beside it
NO_GERMLINE_FREQ = 0.25  # uniform over the four bases where there is no germline source


def load_param_mute_freq_csv(mute_freq_dir, gene):
    """Per-position mutation frequencies for one germline gene, from
    <mute_freq_dir>/<gene, '*' -> '_star_'>.csv.
    None if the file does not exist or carries no rows (e.g. a light-locus D
    placeholder, header-only). Each row also carries _N_OBS_KEY, the summed *_obs counts
    at that position."""
    key = (mute_freq_dir, gene)
    if key in _PARAM_MUTE_FREQ_CACHE:
        return _PARAM_MUTE_FREQ_CACHE[key]
    import csv
    fname = os.path.join(mute_freq_dir, gene.replace('*', '_star_') + '.csv')
    rows = None
    if os.path.exists(fname):
        rows = {}
        with open(fname) as ffile:
            for row in csv.DictReader(ffile):
                frow = {b: float(row[b]) for b in ('A', 'C', 'G', 'T')}
                frow[_N_OBS_KEY] = sum(float(row.get(b + '_obs') or 0) for b in ('A', 'C', 'G', 'T'))
                rows[int(row['position'])] = frow
        if not rows:
            rows = None
    _PARAM_MUTE_FREQ_CACHE[key] = rows
    return rows


def param_dir_region_bounds(antn, glfo):
    """[start, end) bounds for v/d/j in the partition-frame coordinate <antn> is already
    in, rebuilt from glfo germline lengths since v_gl_seq/d_gl_seq/j_gl_seq are
    add_implicit_info()-only and absent from these no-recompute annotations."""
    def gl_len(region, gene, del_5p, del_3p):
        uneroded = glfo['seqs'].get(region, {}).get(gene)
        if uneroded is None:
            return None
        length = len(uneroded) - del_5p - del_3p
        return length if length >= 0 else None
    try:
        len_v = gl_len('v', antn['v_gene'], antn['v_5p_del'], antn['v_3p_del'])
        len_d = gl_len('d', antn['d_gene'], antn['d_5p_del'], antn['d_3p_del'])
        len_j = gl_len('j', antn['j_gene'], antn['j_5p_del'], antn['j_3p_del'])
    except (KeyError, TypeError):
        return None
    if len_v is None or len_d is None or len_j is None:
        return None
    start_v = len(antn['fv_insertion'])
    end_v = start_v + len_v
    start_d = end_v + len(antn['vd_insertion'])
    end_d = start_d + len_d
    start_j = end_d + len(antn['dj_insertion'])
    end_j = start_j + len_j
    end_jf = end_j + len(antn['jf_insertion'])
    return {'start_v': start_v, 'end_v': end_v, 'start_d': start_d, 'end_d': end_d,
            'start_j': start_j, 'end_j': end_j, 'end_jf': end_jf}


def param_dir_classify_pos(pos, bounds):
    """Region label ('v'/'d'/'j'/'insertion'/'other') for one partition-frame position."""
    if bounds is None:
        return 'other'
    if pos < bounds['start_v']:
        return 'insertion'  # fv_insertion
    if pos < bounds['end_v']:
        return 'v'
    if pos < bounds['start_d']:
        return 'insertion'  # vd_insertion
    if pos < bounds['end_d']:
        return 'd'
    if pos < bounds['start_j']:
        return 'insertion'  # dj_insertion
    if pos < bounds['end_j']:
        return 'j'
    if pos < bounds['end_jf']:
        return 'insertion'  # jf_insertion
    return 'other'


def param_dir_mutation_freq_and_obs(pos, base, antn, glfo, mute_freq_dir):
    """(frequency, n_obs) for one (partition-frame position, base) in the frame of <antn>
    (a partition-frame annotation, matching the frame <pos> is already in). n_obs is the
    observation count behind that germline position.

    (None, None) where there is no germline source: D-region (light loci have no real D)
    and every insertion (fv/vd/dj/jf) have no germline base by definition."""
    bounds = param_dir_region_bounds(antn, glfo)
    region = param_dir_classify_pos(pos, bounds)
    if region == 'v':
        offset = max(0, len(antn['fv_insertion']) - antn['v_5p_del'])
        gl_pos, gene = pos - offset, antn['v_gene']
    elif region == 'j':
        gl_pos, gene = pos - bounds['start_j'] + antn['j_5p_del'], antn['j_gene']
    elif region == 'd':
        gl_pos, gene = pos - bounds['start_d'] + antn['d_5p_del'], antn['d_gene']
    else:
        return None, None
    rows = load_param_mute_freq_csv(mute_freq_dir, gene)
    if rows is None or gl_pos not in rows:
        return None, None
    n_obs = rows[gl_pos].get(_N_OBS_KEY)
    return rows[gl_pos].get(base), (n_obs if n_obs else None)


def param_dir_mutation_freq(pos, base, antn, glfo, mute_freq_dir):
    """Germline mutation frequency for one (partition-frame position, base)."""
    return param_dir_mutation_freq_and_obs(pos, base, antn, glfo, mute_freq_dir)[0]


def uid_param_dir_freqs(uid, muts, uid_part_antns, glfo, mute_freq_dir, counts=None):
    """({(pos, base): freq}, {(pos, base): n_obs}) for one uid's own mutations, sourced
    from the parameter directory. NO_GERMLINE_FREQ wherever there is no germline source or
    no partition-frame annotation; those (pos, base) are absent from the n_obs map.

    counts, if passed, is incremented in place so the caller can report how often the
    fallback fires and why: 'no_antn' (uid has no partition-frame annotation) and
    'gene_missing_from_glfo' (the uid's v/d/j gene call is not in <glfo>, so
    param_dir_region_bounds can't be built and every one of the uid's positions falls
    back to NO_GERMLINE_FREQ, not the smoothing floor, since there is no per-gene
    mute-freqs signal at all)."""
    antn = uid_part_antns.get(uid)
    if antn is None:
        if counts is not None:
            counts['no_antn'] = counts.get('no_antn', 0) + 1
        return {(pos, base): NO_GERMLINE_FREQ for pos, base in muts.items()}, {}
    if counts is not None and param_dir_region_bounds(antn, glfo) is None:
        counts['gene_missing_from_glfo'] = counts.get('gene_missing_from_glfo', 0) + 1
    out, obs = {}, {}
    for pos, base in muts.items():
        f, n_obs = param_dir_mutation_freq_and_obs(pos, base, antn, glfo, mute_freq_dir)
        out[(pos, base)] = f if f is not None else NO_GERMLINE_FREQ
        if n_obs is not None:
            obs[(pos, base)] = n_obs
    return out, obs


def weighted_shared_descent_pvalue(muts_a, muts_b, freqs_a, freqs_b, n_seqs, n_obs_a, n_obs_b,
                                    grid=WEIGHT_GRID_NATS, cache=None, uid_a=None, uid_b=None):
    """Probability that two unrelated sequences would share mutations this improbable, with
    each sequence scored against its own per-gene mute-freqs table (freqs_a, freqs_b) rather
    than one shared locus-wide table, since the germline frequency of a position depends on which
    V/J gene that sequence used. Returns 1.0 when nothing is shared. Pass <cache>/<uid_a>/
    <uid_b> to route through a _PerUidPvalCache instead of recomputing from scratch."""
    shared = [(pos, base) for pos, base in muts_a.items() if muts_b.get(pos) == base]
    if not shared:
        return 1.0
    if cache is not None:
        pa = cache.conditional_pvalue(uid_a, muts_a, freqs_a, n_obs_a, shared)
        pb = cache.conditional_pvalue(uid_b, muts_b, freqs_b, n_obs_b, shared)
    else:
        pa = _conditional_pvalue(muts_a, shared, freqs_a, n_seqs, grid, n_obs=n_obs_a)
        pb = _conditional_pvalue(muts_b, shared, freqs_b, n_seqs, grid, n_obs=n_obs_b)
    return math.exp(0.5 * (math.log(max(pa, _TINY)) + math.log(max(pb, _TINY))))


def split_by_weighted_descent(cluster, uid_muts, uid_part_antns, glfo, mute_freq_dir, n_seqs,
                               alpha, counts=None, pair_pvals=None):
    """Split one cluster by weighted shared descent, assigning members to non-transitive
    greedy centroids, scoring each pair against the two sequences' own per-gene mute-freqs
    frequencies rather than a single locus-wide table. Returns a list of sub-clusters
    (lists of uids). counts: optional dict, incremented in place per fallback reason when
    a uid's own mutation frequency can't be sourced from the parameter directory.

    pair_pvals: optional {(i, o): p}, precomputed p-values keyed by this same <cluster>'s
    member-list positions. When given, this is a pure dict-lookup pass with no further
    p-value computation."""
    members = list(cluster)
    n = len(members)
    if n <= 1:
        return [list(members)]
    order = sorted(range(n), key=lambda i: -len(uid_muts.get(members[i], {}) or {}))
    if pair_pvals is None:
        cache = {u: uid_param_dir_freqs(u, uid_muts.get(u, {}), uid_part_antns, glfo, mute_freq_dir, counts=counts) for u in members}
        freq_cache = {u: fo[0] for u, fo in cache.items()}
        obs_cache = {u: fo[1] for u, fo in cache.items()}
        pval_cache = _PerUidPvalCache(n_seqs, WEIGHT_GRID_NATS)

    def pval(i, mi, o, mo, ui, uo):
        if pair_pvals is not None:
            return pair_pvals.get((i, o), pair_pvals.get((o, i), 1.0))
        return weighted_shared_descent_pvalue(
            mi, mo, freq_cache[ui], freq_cache[uo], n_seqs,
            n_obs_a=obs_cache[ui], n_obs_b=obs_cache[uo],
            cache=pval_cache, uid_a=ui, uid_b=uo)

    assigned = [False] * n
    clusters = []
    for i in order:
        if assigned[i]:
            continue
        assigned[i] = True
        sub = [members[i]]
        mi = uid_muts.get(members[i])
        ui = members[i]
        for o in order:
            if assigned[o]:
                continue
            mo = uid_muts.get(members[o])
            if not mi or not mo:  # no mutations is no evidence either way
                continue
            uo = members[o]
            if pval(i, mi, o, mo, ui, uo) < alpha:
                sub.append(members[o])
                assigned[o] = True
        clusters.append(sub)
    return clusters


def _partition_has_real_d(uid_rearr_features):
    """True if any uid carries a real (non-placeholder) D gene (heavy chain).
    Light loci (igk/igl) have no real D, so this is False for them, which is how
    refine_partition picks between the heavy locus and light loci operators when
    light_chain is not passed explicitly."""
    if not uid_rearr_features:
        return False
    for feat in uid_rearr_features.values():
        if not isinstance(feat, dict):
            continue
        d_gene = feat.get('vdj', ('', '', ''))[1]
        if d_gene and 'x-x' not in d_gene and 'Dx' not in d_gene:
            return True
    return False


def split_on_shared_descent(partition, uid_info, uid_sw_naives, uid_part_antns, glfo, mute_freq_dir,
                            n_seqs, alpha=None):
    """Light split: split over-merged clusters by weighted shared descent. Every
    cluster of size >= 2 is passed to the proposer. Light chain only.

    uid_part_antns, glfo, mute_freq_dir: source each uid's own V/J germline mute-freqs
    table. Positions the table cannot cover (insertions, D) fall back to
    NO_GERMLINE_FREQ, not a smoothing floor.

    alpha: link threshold. None (default): fit one per-bin threshold from the pair
    p-value null tail, no split at all for this bin if that fails. Pass an explicit
    float to use one fixed threshold for every cluster instead."""
    result = []
    n_resplit = n_skipped = n_input_seqs = 0
    param_dir_counts = {}

    clusters_uid_muts = []
    for cluster in partition:
        if len(cluster) < 2:
            result.append(cluster)
            n_skipped += 1
            continue
        # per-cell mutations vs the SW naive, for the proposer
        uid_muts = {}
        for uid in cluster:
            if uid in uid_sw_naives and uid in uid_info:
                uid_muts[uid] = get_mutations_with_base(uid_info[uid]['seq'], uid_sw_naives[uid])
        n_input_seqs += len(cluster)
        clusters_uid_muts.append((cluster, uid_muts))

    bin_alpha, spill_path = alpha, None
    if alpha is None:
        scratch_dir = os.environ.get('SLURM_TMPDIR') or os.environ.get('TMPDIR') or '/tmp'
        spill_fd, spill_path = tempfile.mkstemp(prefix='refine-pairpvals-', dir=scratch_dir)
        with os.fdopen(spill_fd, 'wb') as spill_f:
            sample, cand_counts, n_scanned = _scan_bin_full_pairs(
                clusters_uid_muts, uid_part_antns, glfo, mute_freq_dir, n_seqs, spill_f,
                counts=param_dir_counts)
        bin_alpha, fit_diag = derive_bin_alpha(sample, cand_counts, n_scanned)
        if fit_diag['fit_rate'] is None:
            print('  shared-descent null fit: FAILED (fewer than %d tail p-values in %d pairs scanned)' % (
                WEIGHTED_DESCENT_MIN_TAIL_N, n_scanned), flush=True)
        else:
            print('  shared-descent null fit: rate=%.4g frac=%.4g (%d pairs scanned)' % (
                fit_diag['fit_rate'], fit_diag['fit_frac'], n_scanned), flush=True)
        if bin_alpha is None:
            print('  shared-descent: no cutoff met fdr:%.3g or the fdr:%.3g retry -- NO ACTION, '
                  'all %d clusters in this bin left unchanged' % (
                      WEIGHTED_DESCENT_FDR_PRIMARY, WEIGHTED_DESCENT_FDR_RETRY, len(clusters_uid_muts)), flush=True)
        else:
            print('  shared-descent: cutoff %.4g derived via %s' % (bin_alpha, fit_diag['rule']), flush=True)

    try:
        spill_r = open(spill_path, 'rb') if bin_alpha is not None and spill_path is not None else None
        try:
            for cluster, uid_muts in clusters_uid_muts:
                if bin_alpha is None:
                    result.append(cluster)  # no defensible cutoff for this bin: leave every cluster as-is
                    continue
                pair_pvals = _read_pair_block(spill_r) if spill_r is not None else None
                pieces = split_by_weighted_descent(
                    list(cluster), uid_muts, uid_part_antns, glfo, mute_freq_dir, n_seqs, bin_alpha,
                    counts=param_dir_counts, pair_pvals=pair_pvals)
                result.extend(pieces)
                if len(pieces) > 1:
                    n_resplit += 1
        finally:
            if spill_r is not None:
                spill_r.close()
    finally:
        if spill_path is not None:
            os.remove(spill_path)

    n_result_singletons = sum(1 for c in result if len(c) == 1)
    alpha_label = '%.3g' % bin_alpha if bin_alpha is not None else 'none (no split)'
    print('  shared-descent split [alpha=%s]: %d re-split (%d seqs processed), %d skipped, %d -> %d clusters (%d singletons)' % (
        alpha_label, n_resplit, n_input_seqs, n_skipped,
        len(partition), len(result), n_result_singletons), flush=True)
    n_gene_missing = param_dir_counts.get('gene_missing_from_glfo', 0)
    n_no_antn = param_dir_counts.get('no_antn', 0)
    if n_gene_missing > 0 or n_no_antn > 0:
        from partis import utils
        print('  %s per-gene mute-freqs: %d uid-lookups had a v/d/j gene call missing from glfo '
              '(no region bounds, every position fell back to NO_GERMLINE_FREQ, not the smoothing floor, '
              'since there was no per-gene mute-freqs signal), %d uids had no partition-frame annotation' % (
                  utils.wrnstr(), n_gene_missing, n_no_antn), flush=True)
    return result


def calc_metrics(true_partition, inf_partition):
    uid_to_fam = {}
    for i, fam in enumerate(true_partition):
        for u in fam:
            uid_to_fam[u] = i

    family_sizes = defaultdict(int)
    for fam_id in uid_to_fam.values():
        family_sizes[fam_id] += 1

    total_weight = 0
    purity_sum = 0
    for cluster in inf_partition:
        fams = defaultdict(int)
        for uid in cluster:
            if uid in uid_to_fam:
                fams[uid_to_fam[uid]] += 1
        if not fams:
            continue
        dominant = max(fams.values())
        size = sum(fams.values())
        purity_sum += size * (dominant / size)
        total_weight += size
    purity = purity_sum / total_weight if total_weight > 0 else 0

    fam_clusters = defaultdict(lambda: defaultdict(int))
    uid_to_clust = {}
    for i, cluster in enumerate(inf_partition):
        for uid in cluster:
            uid_to_clust[uid] = i

    for uid, fam_id in uid_to_fam.items():
        if uid in uid_to_clust:
            fam_clusters[fam_id][uid_to_clust[uid]] += 1

    comp_sum = 0
    comp_weight = 0
    for fam_id, clust_counts in fam_clusters.items():
        fam_size = family_sizes[fam_id]
        if fam_size < 2:
            continue
        dominant = max(clust_counts.values())
        comp_sum += fam_size * (dominant / fam_size)
        comp_weight += fam_size
    completeness = comp_sum / comp_weight if comp_weight > 0 else 0

    return purity, completeness


def refine_partition(partition, uid_info, uid_sw_naives, uid_rearr_features=None,
                     naive_threshold=None,
                     min_agreement=0.15, min_fp_positions=0,
                     min_weighted_score=MERGE_WEIGHTED_SCORE_CUTOFF,
                     skip_singleton_merge=True,
                     min_cluster_size=2, light_chain=None, alpha=None,
                     parameter_dir=None, length_veto_min_shared=LENGTH_VETO_MIN_SHARED,
                     verbose=True, random_seed=None,
                     mute_freq_dir=None, uid_part_antns=None, glfo=None):
    """Run refinement and return the refined partition.

    partition: list of clusters, each a list of uids.
    uid_info: uid -> {'seq', 'naive', 'cdr3_length'} (partition annotations).
    uid_sw_naives: uid -> per-sequence SW naive_seq.
    uid_rearr_features: uid -> {'vdj': (v, d, j), 'v_3p_del', 'j_5p_del', 'd_5p_del',
    'd_3p_del', 'len_vd', 'len_dj'}.

    Which operators run forks on chain, since each is built on the signal its locus
    provides: heavy splits on naive identity then merges on naive similarity, light
    splits on shared descent and nothing else. light_chain: if None, inferred from
    D-gene presence in uid_rearr_features.
    alpha: light shared-descent link threshold. None (default): derive per-bin, see
    split_on_shared_descent.

    parameter_dir: the locus-level parameter dir, and the only input the heavy locus's
    length veto takes: its length cutoff is derived from it here, so any caller passing the
    same dir refines the same way. Absent, the veto warns and stays off; present but
    incomplete, it raises instead. It reaches only the heavy locus.

    mute_freq_dir, uid_part_antns, glfo: required for both chains. Population-level
    per-position mutation frequency table, plus the annotations needed to look each uid
    up in it.

    random_seed: seeds the global RNG. Refinement reads no RNG, so this changes nothing.
    """
    if random_seed is not None:
        import random
        random.seed(random_seed)
    partition = [list(c) for c in partition]
    all_uids = set(uid for c in partition for uid in c)
    # sw naives have to already be in the partition frame (read_refine_inputs does this), since a
    # frame mismatch shifts every position silently
    badfo = [(u, len(uid_info[u]['seq']), len(uid_sw_naives[u])) for u in all_uids
             if u in uid_info and u in uid_sw_naives and len(uid_info[u]['seq']) != len(uid_sw_naives[u])]
    if len(badfo) > 0:
        raise Exception('%d uids whose sw naive is not in the partition frame (e.g. %s); pass naives through pad_sw_naive()' % (len(badfo), badfo[:3]))
    n_no_naive = sum(1 for uid in all_uids if uid not in uid_sw_naives)
    if n_no_naive > 0:  # validation keys off naives, so a split fails closed and keeps these fragments intact
        print('  warning: %d/%d input seqs have no sw naive (refine validates on naives; their fragments are kept intact)' % (n_no_naive, len(all_uids)), flush=True)
    uid_to_muts_with_base, uid_to_muts_sw = {}, {}
    for uid, info in uid_info.items():
        uid_to_muts_with_base[uid] = get_mutations_with_base(info['seq'], info['naive'])
        if uid in uid_sw_naives:  # mutations against each sequence's own sw naive
            uid_to_muts_sw[uid] = get_mutations_with_base(info['seq'], uid_sw_naives[uid])

    # resolved before any refinement, so the d-gene test gates all of it
    has_d = _partition_has_real_d(uid_rearr_features)
    if light_chain is None:
        if not uid_rearr_features:  # otherwise the light path is taken silently
            print('  warning: no rearrangement features, so assuming light chain', flush=True)
        light_chain = not has_d
    elif light_chain == has_d:  # cross check the caller's locus, since a stray d call means one of them is wrong
        from partis import utils
        print('  %s locus says %s chain but the d genes say %s: using the locus, but check the annotations'
              % (utils.wrnstr(), 'light' if light_chain else 'heavy', 'heavy' if has_d else 'light'), flush=True)

    if light_chain:
        missing = [n for n, v in (('mute_freq_dir', mute_freq_dir), ('uid_part_antns', uid_part_antns), ('glfo', glfo)) if v is None]
        if missing:
            raise Exception('light-chain refine missing %s' % ', '.join(missing))
        if verbose:
            alpha_label = '%.3g' % alpha if alpha is not None else 'auto'
            print('\n=== light: shared-descent split (alpha=%s) ===' % alpha_label, flush=True)
        tstart = time.time()
        out = split_on_shared_descent(partition, uid_info, uid_sw_naives, uid_part_antns, glfo,
                                      mute_freq_dir, len(uid_to_muts_sw), alpha=alpha)
        print('  timing: shared-descent split %.2f s' % (time.time() - tstart), flush=True)
        return out

    missing = [n for n, v in (('mute_freq_dir', mute_freq_dir), ('uid_part_antns', uid_part_antns), ('glfo', glfo)) if v is None]
    if missing:
        raise Exception('heavy-chain refine missing %s' % ', '.join(missing))

    naive_thresh = (naive_threshold if naive_threshold is not None
                    else estimate_naive_threshold(partition, uid_sw_naives))

    _, length_cutoff = (
        length_veto_inputs(parameter_dir) if length_veto_min_shared else (None, None))
    if verbose:
        length_str = ('' if length_cutoff is None
                      else ', length veto at %d shared and len_vd + len_dj >= %.1f' % (
                          length_veto_min_shared, length_cutoff))
        print('\n=== heavy: naive-identity split (EJ veto >= %.2f%s) ===' % (EJ_SAME_FAMILY_FLOOR, length_str), flush=True)
    tstart = time.time()
    split_partition = split_on_naive_identity(
        partition, uid_sw_naives, uid_to_muts_sw, min_cluster_size,
        uid_rearr_features=uid_rearr_features, length_cutoff=length_cutoff,
        length_veto_min_shared=length_veto_min_shared)
    tsplit = time.time()
    print('  timing: naive-identity split %.2f s' % (tsplit - tstart), flush=True)

    if verbose:
        print('\n=== heavy: naive-similarity merge (naive <= %.4f, min_agreement %.2f) ===' % (
            naive_thresh, min_agreement), flush=True)
    out = merge_on_naive_similarity(
        split_partition, uid_info, uid_sw_naives, uid_to_muts_with_base, naive_thresh,
        min_agreement=min_agreement, min_fp_positions=min_fp_positions, min_weighted_score=min_weighted_score,
        skip_singleton_merge=skip_singleton_merge, uid_rearr_features=uid_rearr_features,
        mute_freq_dir=mute_freq_dir, uid_part_antns=uid_part_antns, glfo=glfo)
    print('  timing: naive-similarity merge %.2f s' % (time.time() - tsplit), flush=True)
    return out


# ----------------------------------------------------------------------------
# File I/O layer (reads already-generated partition + sw-cache, writes full
# partis output). Kept separate from the pure algorithm above; imports
# partis.utils lazily so refine_partition() stays dependency-light. Works as a
# standalone per-group job: no SW, vsearch, bcrham, or HMM is re-run.
# ----------------------------------------------------------------------------

def sw_frame_offset(part_antn, sw_antn):
    """Offset of the sw annotation's frame within the partition annotation's frame. Each
    annotation records where its own (indel-reversed) frame puts the v codon, and the partition
    frame is the sw frame with N padding on each side, so the difference is the left pad."""
    return part_antn['codon_positions']['v'] - sw_antn['codon_positions']['v']


def pad_sw_naive(sw_seq, sw_naive, part_seq, ioff):
    """Re-pad a per-sequence sw naive into the partition frame at offset <ioff> (from
    sw_frame_offset()). The partition annotation is padded (N) relative to the sw annotation, so
    comparing the two frames position-wise is a frame shift; everything downstream keys mutations
    by position across cluster members, so they all have to be in the partition frame. Returns
    None if the sw sequence does not sit at <ioff> in <part_seq>, i.e. the two annotations
    disagree on coordinates and the caller should drop this uid's naive."""
    nright = len(part_seq) - ioff - len(sw_seq)
    if ioff < 0 or nright < 0:
        return None
    # partis masks the odd base in the padded frame, e.g. a 1-base jf_insertion, so N matches anything
    if any(p != 'N' and p != s for p, s in zip(part_seq[ioff : ioff + len(sw_seq)], sw_seq)):
        return None
    padded = 'N' * ioff + sw_naive + 'N' * nright
    if len(padded) != len(part_seq):  # shm indel: different coord systems, padding cannot align them
        return None
    return padded


def get_antn_key(antn, key, label):
    """The annotation's <key>, with no default. A default would flow onward as a real annotation
    saying something false, e.g. an empty gene name into the heavy junction guard or a zero
    cdr3_length into the heavy merge's length bins. Failed queries carry none of these keys, so skip
    those before asking (partis writes them as invalid stubs of unique_ids and input_seqs)."""
    if key not in antn:
        raise Exception('no \'%s\' in %s annotation for %s (implicit info has to be added when reading)' % (key, label, ':'.join(antn['unique_ids'][:3])))
    return antn[key]


def get_ir_seqs(antn, label):
    """The annotation's indel-reversed seqs, i.e. the ones that line up with naive_seq. There is
    deliberately no input_seqs fallback: input_seqs are not indel-reversed, so substituting them
    shifts every mutation position on any sequence carrying an shm indel."""
    return get_antn_key(antn, 'seqs', label)


def read_refine_inputs(partition_fname, sw_cache_fname):
    """Read an existing partition file and sw-cache (no recompute) and build the
    inputs refine_partition() needs, plus the per-sequence sw_info and the
    partition annotations used to write full output. Returns a dict with keys:
    glfo, part_glfo, partition, uid_info, uid_sw_naives, uid_rearr_features,
    sw_info, uid_part_antns."""
    from partis import utils
    from partis import disjointgrouper
    # refine reads only stored keys
    part_glfo, part_antns, cpath = utils.read_output(partition_fname, dont_add_implicit_info=True)
    disjointgrouper.check_stage_file_complete(partition_fname, part_antns, cpath)  # a truncated input silently shrinks the refined output
    partition = [list(c) for c in (cpath.best() if cpath is not None else [])]
    uid_info, uid_part_antns = {}, {}
    n_failed, n_short = 0, 0
    for antn in part_antns:
        failed = antn.get('invalid', False)  # sw/hmm failures, written as stubs with no annotation
        if failed:
            n_failed += 1
        else:
            naive = get_antn_key(antn, 'naive_seq', 'partition')
            cdr3 = get_antn_key(antn, 'cdr3_length', 'partition')
            iseqs = get_ir_seqs(antn, 'partition')
        for i, uid in enumerate(antn['unique_ids']):
            uid_part_antns[uid] = antn  # write path synthesizes from these, so output stays in the input's padded frame
            if not failed:
                if i < len(iseqs):
                    uid_info[uid] = {'seq': iseqs[i], 'naive': naive, 'cdr3_length': cdr3}
                else:  # fewer seqs than uids, i.e. a malformed annotation
                    n_short += 1
    if n_failed > 0:
        print('  skipped %d failed queries (no annotation) in %s' % (n_failed, partition_fname), flush=True)
    if n_short > 0:
        print('  %s dropped %d uids whose annotation had fewer seqs than unique_ids in %s' % (utils.wrnstr(), n_short, partition_fname), flush=True)
    # they carry no partition-frame seq to pad an sw naive against, and are dropped from the output
    # anyway, so drop them here rather than letting them into grouping decisions
    n_before = sum(len(c) for c in partition)
    partition = [[u for u in c if u in uid_info] for c in partition]
    partition = [c for c in partition if len(c) > 0]
    n_unannotated = n_before - sum(len(c) for c in partition)
    if n_unannotated > 0:
        print('  %s dropped %d partitioned uids with no usable annotation before refining %s' % (utils.wrnstr(), n_unannotated, partition_fname), flush=True)

    sw_glfo, sw_antns, _ = utils.read_output(sw_cache_fname, dont_add_implicit_info=True)
    sw_info, uid_sw_naives, uid_rearr_features = {}, {}, {}
    n_unalignable = 0
    for antn in sw_antns:  # the sw cache has no failed queries in it (waterer.write_cachefile)
        sw_seqs = get_ir_seqs(antn, 'sw cache')
        sw_naive = get_antn_key(antn, 'naive_seq', 'sw cache')
        vdj = tuple(get_antn_key(antn, '%s_gene' % r, 'sw cache') for r in utils.regions)
        # junction boundaries for the heavy merge guard
        v_3p_del, j_5p_del = get_antn_key(antn, 'v_3p_del', 'sw cache'), get_antn_key(antn, 'j_5p_del', 'sw cache')
        # the rest of the rearrangement, for the heavy split's length veto
        d_5p_del, d_3p_del = get_antn_key(antn, 'd_5p_del', 'sw cache'), get_antn_key(antn, 'd_3p_del', 'sw cache')
        len_vd, len_dj = [len(get_antn_key(antn, '%s_insertion' % b, 'sw cache')) for b in ('vd', 'dj')]
        for i, uid in enumerate(antn['unique_ids']):
            sw_info[uid] = antn
            uid_rearr_features[uid] = {'vdj': vdj, 'v_3p_del': v_3p_del, 'j_5p_del': j_5p_del,
                                       'd_5p_del': d_5p_del, 'd_3p_del': d_3p_del,
                                       'len_vd': len_vd, 'len_dj': len_dj}
            if uid not in uid_info or i >= len(sw_seqs):  # nothing to pad against, and an
                continue                                  # unpadded sw naive is in the wrong frame
            naive = pad_sw_naive(sw_seqs[i], sw_naive, uid_info[uid]['seq'],
                                 sw_frame_offset(uid_part_antns[uid], antn))
            if naive is None:  # frames disagree: drop, refine fails closed on missing naives
                n_unalignable += 1
                continue
            uid_sw_naives[uid] = naive
    if n_unalignable > 0:
        print('  %s could not align %d sw naives to the partition frame (dropped)' % (utils.wrnstr(), n_unalignable), flush=True)
    return {'glfo': sw_glfo, 'part_glfo': part_glfo, 'partition': partition, 'uid_info': uid_info,
            'uid_sw_naives': uid_sw_naives, 'uid_rearr_features': uid_rearr_features,
            'sw_info': sw_info, 'uid_part_antns': uid_part_antns}


def write_full_output(outfname, glfo, refined_partition, ant_info, label='refine'):
    """Write a full partis output file (germline-info + annotations + partition)
    for a refined partition. Each cluster's annotation is synthesized from the
    per-uid annotations in <ant_info> via the same no-recompute path partis uses for
    --fast (synthesize_multi_seq_line_from_reco_info). Both callers pass the partition
    step's annotations: synthesis takes the family keys from the first uid and the seqs
    from all of them, so every uid must be in one padded frame (the sw cache is unpadded,
    each sequence in its own germline frame). A cluster whose multi-sequence synthesis
    fails is still emitted, as singletons, rather than dropped, so the written partition
    and annotation list always match, since paired clustering requires an annotation for every
    partition cluster. Uids whose annotation is invalid are dropped up front so one does
    not shatter its cluster.
    label: names the calling stage in the drop message, since ha-repartition calls this too.
    Returns counts: n_written, n_dropped, n_fallback, n_fallback_uids."""
    import os
    from partis import utils
    from partis import clusterpath
    from partis import indelutils

    def _annotate(uids):
        antn = utils.synthesize_multi_seq_line_from_reco_info(uids, ant_info, warn=False)
        utils.remove_all_implicit_info(antn)
        utils.add_implicit_info(glfo, antn, reset_indel_genes=True)
        return antn

    def _fast_singleton(uid):
        # slice the stored annotation instead of rebuilding its implicit info, same result for one uid
        src = ant_info[uid]
        iseq = src['unique_ids'].index(uid)
        if indelutils.has_indels_line(src, iseq):  # these need _annotate's reset_indel_genes
            return None
        antn = utils.synthesize_single_seq_line(src, iseq)
        antn['indelfos'] = [indelutils.get_empty_indel()]
        for key in utils.special_indel_columns_for_output:  # writer takes key order from indelfos
            antn.pop(key, None)
        return antn

    def _pad_to_uniform_length(antns):
        # pad naive_seqs to uniform length per cdr3 class (a no-op when the input annotations
        # are already uniformly padded, kept as a safety net for inputs that are not)
        maxfo = {}  # cdr3_length: [max_gl_cpos, max_gl_cpos_to_j_end]
        for antn in antns:
            cpos, seqlen = antn['codon_positions']['v'], len(antn['seqs'][0])
            gl_cpos = glfo['cyst-positions'][antn['v_gene']] + max(0, len(antn['fv_insertion']) - antn['v_5p_del'])
            gl_cpos_to_j_end = seqlen - cpos + antn['j_3p_del']
            cdr3 = antn['cdr3_length']
            if cdr3 not in maxfo:
                maxfo[cdr3] = [gl_cpos, gl_cpos_to_j_end]
            else:
                maxfo[cdr3][0] = max(maxfo[cdr3][0], gl_cpos)
                maxfo[cdr3][1] = max(maxfo[cdr3][1], gl_cpos_to_j_end)
        for antn in antns:
            cpos, seqlen = antn['codon_positions']['v'], len(antn['seqs'][0])
            padleft = maxfo[antn['cdr3_length']][0] - cpos
            padright = maxfo[antn['cdr3_length']][1] - (seqlen - cpos)
            if padleft < 0 or padright < 0:  # should never be negative (max minus value)
                raise Exception('bad padding %d %d for cluster %s' % (padleft, padright, antn['unique_ids'][0]))
            if padleft != 0 or padright != 0:
                utils.re_pad_atn(padleft, padright, antn, glfo)

    annotation_list, out_partition = [], []
    n_dropped, n_fallback, n_fallback_uids = 0, 0, 0
    first_err = None
    for cluster in refined_partition:
        known = [uid for uid in cluster if uid in ant_info]
        n_dropped += len(cluster) - len(known)  # no annotation at all, so nothing to synthesize from
        good = [uid for uid in known if not ant_info[uid].get('invalid', False)]
        n_dropped += len(known) - len(good)
        if len(good) == 0:
            continue
        try:
            antn = _fast_singleton(good[0]) if len(good) == 1 else None
            annotation_list.append(antn if antn is not None else _annotate(good))
            out_partition.append(list(good))
        except Exception as e:  # fall back to singletons on synthesis failure
            first_err = first_err if first_err is not None else repr(e)
            n_fallback += 1
            n_fallback_uids += len(good)
            for uid in good:
                try:
                    annotation_list.append(_annotate([uid]))
                    out_partition.append([uid])
                except Exception as e2:
                    first_err = first_err if first_err is not None else repr(e2)
                    n_dropped += 1
    _pad_to_uniform_length(annotation_list)
    cpath = clusterpath.ClusterPath(partition=out_partition)
    partition_lines = cpath.get_partition_lines()
    outdir = os.path.dirname(outfname)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
    utils.write_annotations(outfname, glfo, annotation_list, utils.annotation_headers,
                            partition_lines=partition_lines)
    if n_dropped or n_fallback:
        msgs = []
        if n_dropped > 0:
            msgs.append('dropped %d unannotatable uids (invalid annotation)' % n_dropped)
        if n_fallback > 0:
            msgs.append('%d clusters (%d uids) fell back to singletons (synthesis failed)' % (n_fallback, n_fallback_uids))
        print('  %s %s output: %s%s' % (utils.wrnstr(), label, '; '.join(msgs),
                                        ('; first error: %s' % first_err) if first_err else ''), flush=True)
    return {'n_written': len(annotation_list), 'n_dropped': n_dropped,
            'n_fallback': n_fallback, 'n_fallback_uids': n_fallback_uids}


# ----------------------------------------------------------------------------
# Locus-level orchestration over the disjoint-groups manifest (mirrors the
# ha_repartition module). The unit is the group (refine is pure-python, no per-cluster
# fan-out), so there is no separate prepare step; assemble is assemble_groups.
# ----------------------------------------------------------------------------
def bundle_marker_fname(disjoint_dir, job_start):
    # done-marker for the group slice starting at <job_start>, written by the slice and looked
    # for by whatever dispatched it
    return '%s/refine-bundle-%d.done' % (disjoint_dir, job_start)


def group_specs(disjoint_dir, groups, locus):
    """Per-group refine I/O paths. Refine runs on the HA re-partition if it exists,
    else the vsearch partition; output is partition-refine-<locus>.yaml. Only groups
    whose input partition and sw-cache exist are returned."""
    from partis import disjointgrouper
    specs = []
    n_harep = n_vsearch = 0
    for group in groups:
        fasta_dir = os.path.dirname(group['fasta_path'])
        harep_p = '%s/%s/%s' % (disjoint_dir, fasta_dir, disjointgrouper.stage_fname(disjointgrouper.STAGE_HAREP, locus))
        vsearch_p = '%s/%s/%s' % (disjoint_dir, fasta_dir, disjointgrouper.stage_fname(disjointgrouper.STAGE_VSEARCH, locus))
        sw = '%s/%s/%s' % (disjoint_dir, fasta_dir, disjointgrouper.group_sw_cache_fname(locus))
        inp = harep_p if os.path.exists(harep_p) else vsearch_p  # existence is the ha-repartition completion signal
        if not (os.path.exists(inp) and os.path.exists(sw)):
            continue
        if inp == harep_p:
            n_harep += 1
        else:
            n_vsearch += 1
        refined_rel = '%s/%s' % (fasta_dir, disjointgrouper.stage_fname(disjointgrouper.STAGE_REFINE, locus))
        specs.append({'group': group, 'input': inp, 'sw_cache': sw,
                      'refined_out': '%s/%s' % (disjoint_dir, refined_rel),
                      'refined_rel': refined_rel})
    if n_vsearch > 0:  # any vsearch input means ha-repartition did not run
        from partis import utils
        print('  %s refine input: %d groups from ha-repartition, %d from vsearch (run ha-repartition first to refine its output)'
              % (utils.wrnstr(), n_harep, n_vsearch), flush=True)
    return specs


def estimate_locuswide_threshold(specs):
    """Estimate the naive threshold over the full-locus partition (all cdr3 groups in <specs>
    merged). A global statistic, so estimate once here and pass into each per-group refine
    rather than letting refine_partition estimate it per-group."""
    partition, uid_sw_naives = [], {}
    for spec in specs:
        inp = read_refine_inputs(spec['input'], spec['sw_cache'])
        partition.extend(inp['partition'])
        uid_sw_naives.update(inp['uid_sw_naives'])
    return estimate_naive_threshold(partition, uid_sw_naives)


def locuswide_threshold_fname(disjoint_dir, locus):
    return '%s/naive-threshold-%s.json' % (disjoint_dir, locus)


def locuswide_threshold(disjoint_dir, specs, locus, overwrite=False):
    """Locus-wide naive threshold, cached beside the manifest, None where nothing reads it.

    Invalidated by a change in the spec count or in any input's mtime.
    """
    from . import utils
    if len(specs) == 0 or not utils.has_d_gene(locus):
        return None
    mtimes = [os.path.getmtime(spec['input']) for spec in specs]
    signal = {'n_specs': len(specs), 'input_mtime_max': max(mtimes), 'input_mtime_sum': sum(mtimes)}
    fname = locuswide_threshold_fname(disjoint_dir, locus)
    if not overwrite and os.path.exists(fname):
        with open(fname) as tfile:
            cfo = json.load(tfile)
        if all(cfo.get(key) == val for key, val in signal.items()):
            return float(cfo['threshold'])
        print('  locus-wide threshold cache is stale, re-estimating: %s' % fname, flush=True)
    cfo = dict(signal, threshold=repr(estimate_locuswide_threshold(specs)))  # repr so the float round trips exactly
    tmpfname = '%s.tmp.%d' % (fname, os.getpid())  # atomic, so concurrent slices race without corrupting
    with open(tmpfname, 'w') as tfile:
        json.dump(cfo, tfile)
    os.replace(tmpfname, fname)
    return float(cfo['threshold'])


def validate_mute_freq_tables(mfdir):
    """Check the per-gene mute-freqs tables in <mfdir> are there. Any locus can have them; the
    shared-descent split is the only consumer today. Raises rather than falling back: without them
    the split would silently revert to estimating frequencies over its own input, which is the
    scope dependence reading them removes."""
    import glob
    if not os.path.isdir(mfdir) or len(glob.glob('%s/*.csv' % mfdir)) == 0:
        raise Exception('no per-gene mutation frequency tables in %s. a parameter dir merged before '
                        'mute-freqs merging existed will not have them, so re-merge into a fresh dir.' % mfdir)


def run_jobs(specs, naive_threshold=None, overwrite=False, locus=None, parameter_dir=None,
             length_veto_min_shared=LENGTH_VETO_MIN_SHARED, mute_freq_dir=None):
    """Run refinement on a list of group specs (from group_specs), writing each group's
    refined partition, with the production defaults (singleton-skip, junction guard, vdj
    override) that the standalone CLI and integrated pipeline both use. Groups whose
    refined output already exists are skipped unless <overwrite>. The naive threshold
    defaults to a locus-wide estimate over <specs>; when running a slice, pass one
    estimated over the full group list. Passing <locus> skips that estimate on a locus with
    no D gene, where no operator reads it, and pins the heavy/light fork for every group.
    <parameter_dir> is the locus-level parameter dir, and is passed straight through: refine
    derives the heavy split's length veto from it.

    mute_freq_dir: required, both chains read it."""
    from argparse import Namespace
    from partis import utils
    oargs = Namespace(overwrite=overwrite)
    heavy = locus is None or utils.has_d_gene(locus)
    if naive_threshold is None and heavy:
        naive_threshold = estimate_locuswide_threshold(specs)
    light_chain = None if locus is None else not heavy  # locus wins over the per-group d-gene test
    totals, n_run = defaultdict(int), 0
    tlocus = time.time()
    for spec in specs:
        # says what it skipped, and a zero-length output is removed and re-made rather than counted as done
        if utils.output_exists(oargs, spec['refined_out'], outlabel='refine', offset=4):
            continue
        tgroup = time.time()
        inp = read_refine_inputs(spec['input'], spec['sw_cache'])
        refined = refine_partition(
            inp['partition'], inp['uid_info'], inp['uid_sw_naives'],
            uid_rearr_features=inp['uid_rearr_features'],
            naive_threshold=naive_threshold, light_chain=light_chain,
            parameter_dir=parameter_dir, length_veto_min_shared=length_veto_min_shared,
            skip_singleton_merge=True, min_agreement=0.15, verbose=False,
            mute_freq_dir=mute_freq_dir, uid_part_antns=inp['uid_part_antns'],
            # gene calls on uid_part_antns are partition-frame (from spec['input']), so the glfo
            # that resolves them has to be part_glfo, not the sw glfo inp['glfo'] (sw_cache_fname)
            glfo=inp['part_glfo'])
        cfo = write_full_output(spec['refined_out'], inp['part_glfo'], refined, inp['uid_part_antns'])
        print('  timing: group %s total %.2f s' % (os.path.dirname(spec['refined_rel']), time.time() - tgroup), flush=True)
        n_run += 1
        for key, val in cfo.items():
            totals[key] += val
    print('  timing: %d groups total %.2f s' % (n_run, time.time() - tlocus), flush=True)
    if totals['n_dropped'] > 0 or totals['n_fallback'] > 0:  # per-group counts are easy to miss, so total them
        print('  %s refine output over %d groups: dropped %d unannotatable uids, %d clusters (%d uids) fell back to singletons'
              % (utils.wrnstr(), n_run, totals['n_dropped'], totals['n_fallback'], totals['n_fallback_uids']), flush=True)
