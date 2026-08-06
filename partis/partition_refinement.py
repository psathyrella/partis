"""Post-partition refinement by locus-appropriate operators.

Pipeline-agnostic module: operates on a partition (list of clusters, each a list
of uids) plus per-sequence annotations and SW naives, and returns a refined
partition. Usable after any partis partition (standard or disjoint grouping).

Which operators run is decided by D-gene presence, resolved before any of them,
because each operator is built on a signal only its locus carries:

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
import json
import math
import os
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


def fingerprint_agreement(fp1, n1, fp2, n2, min_fp_positions=0):
    """Measure agreement between two cluster fingerprints.

    For each mutation position in fp1, check if fp2 has the same position
    with the same dominant base. Score = fraction of fp1's positions that
    agree with fp2 (symmetric: take the average of both directions).

    Returns (score, n_strong_max) where score is in [0, 1] and n_strong_max
    is the number of strong positions in the larger fingerprint.
    High score = strong agreement (likely same family).
    Low score = independent mutations (likely different families).

    Returns score=-1 for insufficient signal: when either cluster has no member
    carrying mutations, or, if min_fp_positions > 0, when both have fewer strong
    positions than the threshold.
    """
    if n1 == 0 or n2 == 0:
        return -1.0, 0

    def count_strong(fp, n):
        min_count = max(1, n // 2)
        count = 0
        for pos, bases in fp.items():
            dominant_base = max(bases, key=bases.get)
            if bases[dominant_base] >= min_count:
                count += 1
        return count

    def directional_agreement(source_fp, source_n, target_fp, target_n):
        """Fraction of source's strong positions that appear in target."""
        if source_n == 0:
            return 0.0, 0
        # "strong" positions: mutated by at least half the source members
        # for singletons: all positions are strong (count=1, n=1, 1/1 >= 0.5)
        min_count = max(1, source_n // 2)
        strong_positions = []
        for pos, bases in source_fp.items():
            dominant_base = max(bases, key=bases.get)
            if bases[dominant_base] >= min_count:
                strong_positions.append((pos, dominant_base))

        if len(strong_positions) == 0:
            return 0.0, 0

        n_agree = 0
        for pos, base in strong_positions:
            if pos in target_fp and base in target_fp[pos]:
                n_agree += 1

        return n_agree / len(strong_positions), len(strong_positions)

    n_strong_1 = count_strong(fp1, n1)
    n_strong_2 = count_strong(fp2, n2)
    n_strong_max = max(n_strong_1, n_strong_2)

    # if both clusters have too few strong positions, signal is insufficient
    if min_fp_positions > 0 and n_strong_max < min_fp_positions:
        return -1.0, n_strong_max

    fwd, _ = directional_agreement(fp1, n1, fp2, n2)
    rev, _ = directional_agreement(fp2, n2, fp1, n1)

    # asymmetric clusters (one large, one small): take the max, since the small
    # cluster's positions should appear in the large one if they are the same family
    return max(fwd, rev), n_strong_max


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


def split_on_naive_identity(partition, uid_sw_naives, uid_muts_sw, min_cluster_size=2,
                            ej_floor=EJ_SAME_FAMILY_FLOOR):
    """Heavy split: exact sw-naive identity proposes the split, shared mutations veto it.

    Members are grouped by exact per-sequence sw naive, one-member groups are absorbed into the
    nearest corroborated one, then fragments are re-merged when a cross-fragment pair's enhanced
    jaccard reaches ej_floor. uid_muts_sw is mutations against each sequence's own sw naive.
    """
    result = []
    ctr = defaultdict(int)

    for cluster in partition:
        if len(cluster) < min_cluster_size:
            result.append(list(cluster))
            ctr['below_min_size'] += 1
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

    print('  naive-identity split: %d accepted, %d rejected (veto), %d skipped (single naive), %d rejected (missing naive)' % (
        ctr['accepted'], ctr['rejected'], ctr['single_naive'], ctr['missing_naive']), flush=True)
    print('  %d members snapped, %d pairs certified, %d -> %d clusters' % (
        ctr['snapped'], ctr['certified_pairs'], len(partition), len(result)), flush=True)
    return result


# fragment size at which the modal v/d/j call is trusted enough to override a junction mismatch
VDJ_OVERRIDE_MIN_FRAG = 20


def merge_on_naive_similarity(split_partition, uid_info, uid_sw_naives,
                              uid_to_muts_with_base, naive_threshold,
                              min_agreement=0.15, min_fp_positions=0,
                              skip_singleton_merge=False, uid_rearr_features=None,
                              junction_guard=True, vdj_override_min=VDJ_OVERRIDE_MIN_FRAG):
    """Heavy merge: incremental naive merge with fingerprint validation.

    For each CDR3 group, find clusters with similar naives (candidates),
    then only merge if their mutation fingerprints agree. This prevents
    false merges of unrelated sequences with similar naives.

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
    """
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

        # within same-naive bucket: always naive-similar, just check fingerprint
        for naive_seq, indices in naive_to_frags.items():
            for ii in range(len(indices)):
                for jj in range(ii + 1, len(indices)):
                    i, j = indices[ii], indices[jj]
                    if find(i) == find(j):
                        continue
                    # skip singleton-singleton merges: both survived vsearch+HA
                    # as singletons, likely true singletons not fragments
                    if skip_singleton_merge and frag_sizes[i] == 1 and frag_sizes[j] == 1:
                        n_skipped_singleton += 1
                        continue
                    n_naive_candidates += 1
                    verdict = junction_verdict(i, j)
                    if verdict == 'block':
                        n_rejected_junction += 1
                        continue
                    if verdict == 'rescue':
                        n_vdj_override += 1
                    fp_i, n_i = frag_fps[i]
                    fp_j, n_j = frag_fps[j]
                    agreement, n_strong = fingerprint_agreement(
                        fp_i, n_i, fp_j, n_j, min_fp_positions)
                    if agreement < 0:
                        n_rejected_insufficient += 1
                    elif agreement >= min_agreement:
                        union(i, j)
                        n_accepted += 1
                    else:
                        n_rejected_fingerprint += 1

        # cross-bucket: vectorized pairwise hamming to find near-match naives
        unique_naives = list(naive_to_frags.keys())
        if len(unique_naives) >= 2:
            seq_len = len(unique_naives[0])
            same_len = [n for n in unique_naives if len(n) == seq_len]
            diff_len = [n for n in unique_naives if len(n) != seq_len]

            if len(same_len) >= 2:
                arr = np.frombuffer(''.join(same_len).encode(), dtype=np.uint8).reshape(len(same_len), seq_len)
                n_byte = ord('N')
                # chunked pairwise to bound memory when there are many unique naives
                chunk_size = 2000
                for ci in range(0, len(same_len), chunk_size):
                    ci_end = min(ci + chunk_size, len(same_len))
                    chunk_i = arr[ci:ci_end]
                    for cj in range(ci, len(same_len), chunk_size):
                        cj_end = min(cj + chunk_size, len(same_len))
                        chunk_j = arr[cj:cj_end]
                        # N-aware hamming (matches hamming_frac): drop positions where either naive is N
                        valid = (chunk_i[:, None, :] != n_byte) & (chunk_j[None, :, :] != n_byte)
                        mism = ((chunk_i[:, None, :] != chunk_j[None, :, :]) & valid).sum(axis=2)
                        compared = valid.sum(axis=2)
                        dists = np.where(compared > 0, mism / np.maximum(compared, 1), 1.0)
                        pairs = np.argwhere(dists <= naive_threshold)
                        for pi, pj in pairs:
                            ni_idx = ci + pi
                            nj_idx = cj + pj
                            if ni_idx >= nj_idx:
                                continue
                            for i in naive_to_frags[same_len[ni_idx]]:
                                for j in naive_to_frags[same_len[nj_idx]]:
                                    if find(i) == find(j):
                                        continue
                                    if skip_singleton_merge and frag_sizes[i] == 1 and frag_sizes[j] == 1:
                                        n_skipped_singleton += 1
                                        continue
                                    n_naive_candidates += 1
                                    verdict = junction_verdict(i, j)
                                    if verdict == 'block':
                                        n_rejected_junction += 1
                                        continue
                                    if verdict == 'rescue':
                                        n_vdj_override += 1
                                    fp_i, n_i = frag_fps[i]
                                    fp_j, n_j = frag_fps[j]
                                    agreement, n_strong = fingerprint_agreement(
                                        fp_i, n_i, fp_j, n_j, min_fp_positions)
                                    if agreement < 0:
                                        n_rejected_insufficient += 1
                                    elif agreement >= min_agreement:
                                        union(i, j)
                                        n_accepted += 1
                                    else:
                                        n_rejected_fingerprint += 1

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
                            if find(i) == find(j):
                                continue
                            if skip_singleton_merge and frag_sizes[i] == 1 and frag_sizes[j] == 1:
                                n_skipped_singleton += 1
                                continue
                            n_naive_candidates += 1
                            verdict = junction_verdict(i, j)
                            if verdict == 'block':
                                n_rejected_junction += 1
                                continue
                            if verdict == 'rescue':
                                n_vdj_override += 1
                            fp_i, n_i = frag_fps[i]
                            fp_j, n_j = frag_fps[j]
                            agreement, n_strong = fingerprint_agreement(
                                fp_i, n_i, fp_j, n_j, min_fp_positions)
                            if agreement < 0:
                                n_rejected_insufficient += 1
                            elif agreement >= min_agreement:
                                union(i, j)
                                n_accepted += 1
                            else:
                                n_rejected_fingerprint += 1

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

WEIGHTED_DESCENT_ALPHA = 0.005  # link threshold (weighted-surprisal tail cutoff)
WEIGHT_GRID_NATS = 0.25  # bin width for the exact tail dp
FREQ_SMOOTH_COUNT = 0.5  # count floor for an unobserved (position, base)
_TINY = 1e-300


def repertoire_mutation_freqs(uid_muts, n_seqs=None):
    """Per-(position, base) mutation frequency over every sequence in <uid_muts>, and the
    denominator used. Pass a whole refine input rather than one cluster: the null must not be
    estimated from the sequences being compared."""
    counts = defaultdict(int)
    for muts in uid_muts.values():
        for pos, base in muts.items():
            counts[(pos, base)] += 1
    n = len(uid_muts) if n_seqs is None else n_seqs
    n = max(int(n), 1)
    return {pb: c / n for pb, c in counts.items()}, n


def _weight_bins(muts, freqs, n_seqs, grid):
    """{(pos, base): (frequency, surprisal in whole grid bins)} for one sequence's mutations.
    Each weight is binned, not the sums, which keeps the tail below exact."""
    floor = FREQ_SMOOTH_COUNT / n_seqs
    out = {}
    for pos, base in muts.items():
        p = freqs.get((pos, base), 0.0)
        if p < floor:
            p = floor
        elif p > 1.0:
            p = 1.0
        nbin = int(round(-math.log(p) / grid))
        out[(pos, base)] = (p, max(nbin, 0))
    return out


def _conditional_pvalue(muts_cond, shared, freqs, n_seqs, grid):
    """P(T >= T_obs), where T sums the surprisals of whichever of <muts_cond> another sequence
    carries independently. Exact by dp over binned surprisal, with everything at or above the
    observed value collected into a tail bucket."""
    wb = _weight_bins(muts_cond, freqs, n_seqs, grid)
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


def weighted_shared_descent_pvalue(muts_a, muts_b, freqs, n_seqs, grid=WEIGHT_GRID_NATS):
    """Probability that two unrelated sequences would share mutations this improbable. Both
    mutation dicts are against each sequence's own sw naive. Returns 1.0 when nothing is shared.
    Symmetric: each direction conditions on one of the two mutation sets and the pair of
    p-values is combined as a geometric mean."""
    shared = [(pos, base) for pos, base in muts_a.items() if muts_b.get(pos) == base]
    if not shared:
        return 1.0
    pa = _conditional_pvalue(muts_a, shared, freqs, n_seqs, grid)
    pb = _conditional_pvalue(muts_b, shared, freqs, n_seqs, grid)
    return math.exp(0.5 * (math.log(max(pa, _TINY)) + math.log(max(pb, _TINY))))


def split_by_weighted_descent(cluster, uid_muts, freqs, n_seqs, alpha=WEIGHTED_DESCENT_ALPHA):
    """Split one cluster by weighted shared descent, assigning members to non-transitive greedy
    centroids. Returns a list of sub-clusters (lists of uids)."""
    members = list(cluster)
    n = len(members)
    if n <= 1:
        return [list(members)]
    order = sorted(range(n), key=lambda i: -len(uid_muts.get(members[i], {}) or {}))
    assigned = [False] * n
    clusters = []
    for i in order:
        if assigned[i]:
            continue
        assigned[i] = True
        sub = [members[i]]
        mi = uid_muts.get(members[i])
        for o in order:
            if assigned[o]:
                continue
            mo = uid_muts.get(members[o])
            if not mi or not mo:  # no mutations is no evidence either way
                continue
            if weighted_shared_descent_pvalue(mi, mo, freqs, n_seqs) < alpha:
                sub.append(members[o])
                assigned[o] = True
        clusters.append(sub)
    return clusters


def _partition_has_real_d(uid_rearr_features):
    """True if any uid carries a real (non-placeholder) D gene (heavy chain).
    Light loci (igk/igl) have no real D, so this is False for them, which is how
    refine_partition picks the heavy vs light fork when light_chain is not passed
    explicitly."""
    if not uid_rearr_features:
        return False
    for feat in uid_rearr_features.values():
        if not isinstance(feat, dict):
            continue
        d_gene = feat.get('vdj', ('', '', ''))[1]
        if d_gene and 'x-x' not in d_gene and 'Dx' not in d_gene:
            return True
    return False


def split_on_shared_descent(partition, uid_info, uid_sw_naives, freqs, n_seqs,
                            alpha=WEIGHTED_DESCENT_ALPHA):
    """Light split: split over-merged clusters by weighted shared descent
    (split_by_weighted_descent). Every cluster of size >= 2 is passed to the
    proposer. Light chain only (see refine_partition).

    freqs, n_seqs: from repertoire_mutation_freqs over the whole refine input.
    alpha: link threshold for the weighted shared-descent test
    (default WEIGHTED_DESCENT_ALPHA).
    """
    result = []
    n_resplit = n_skipped = n_input_seqs = 0

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
        pieces = split_by_weighted_descent(list(cluster), uid_muts, freqs, n_seqs, alpha)
        result.extend(pieces)
        if len(pieces) > 1:
            n_resplit += 1

    n_result_singletons = sum(1 for c in result if len(c) == 1)
    print('  shared-descent split [alpha=%.3g]: %d re-split (%d seqs processed), %d skipped, %d -> %d clusters (%d singletons)' % (
        alpha, n_resplit, n_input_seqs, n_skipped,
        len(partition), len(result), n_result_singletons), flush=True)
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
                     min_agreement=0.15, min_fp_positions=0, skip_singleton_merge=True,
                     min_cluster_size=2, light_chain=None, alpha=WEIGHTED_DESCENT_ALPHA,
                     verbose=True, random_seed=None):
    """Run refinement and return the refined partition.

    partition: list of clusters, each a list of uids.
    uid_info: uid -> {'seq', 'naive', 'cdr3_length'} (partition annotations).
    uid_sw_naives: uid -> per-sequence SW naive_seq.
    uid_rearr_features: uid -> {'vdj': (v, d, j), 'v_3p_del', 'j_5p_del'}.

    Which operators run forks on chain, since each is built on the signal its locus
    provides: heavy splits on naive identity then merges on naive similarity, light
    splits on shared descent and nothing else. light_chain: if None, inferred from
    D-gene presence in uid_rearr_features.
    alpha: link threshold for the light shared-descent test.

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
    if light_chain is None:
        if not uid_rearr_features:  # otherwise the light path is taken silently
            print('  warning: no rearrangement features, so assuming light chain', flush=True)
        light_chain = not _partition_has_real_d(uid_rearr_features)

    if light_chain:
        if verbose:
            print('\n=== light: shared-descent split (alpha=%.3g) ===' % alpha, flush=True)
        tstart = time.time()
        wd_freqs, wd_n_seqs = repertoire_mutation_freqs(uid_to_muts_sw)
        out = split_on_shared_descent(partition, uid_info, uid_sw_naives, wd_freqs, wd_n_seqs,
                                      alpha=alpha)
        print('  timing: shared-descent split %.2f s' % (time.time() - tstart), flush=True)
        return out

    naive_thresh = (naive_threshold if naive_threshold is not None
                    else estimate_naive_threshold(partition, uid_sw_naives))

    if verbose:
        print('\n=== heavy: naive-identity split (EJ veto >= %.2f) ===' % EJ_SAME_FAMILY_FLOOR, flush=True)
    tstart = time.time()
    split_partition = split_on_naive_identity(
        partition, uid_sw_naives, uid_to_muts_sw, min_cluster_size)
    tsplit = time.time()
    print('  timing: naive-identity split %.2f s' % (tsplit - tstart), flush=True)

    if verbose:
        print('\n=== heavy: naive-similarity merge (naive <= %.4f, min_agreement %.2f) ===' % (
            naive_thresh, min_agreement), flush=True)
    out = merge_on_naive_similarity(
        split_partition, uid_info, uid_sw_naives, uid_to_muts_with_base, naive_thresh,
        min_agreement, min_fp_positions, skip_singleton_merge=skip_singleton_merge,
        uid_rearr_features=uid_rearr_features)
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
    part_glfo, part_antns, cpath = utils.read_output(partition_fname)
    disjointgrouper.check_stage_file_complete(partition_fname, part_antns, cpath)  # a truncated input silently shrinks the refined output
    partition = [list(c) for c in (cpath.best() if cpath is not None else [])]
    uid_info, uid_part_antns = {}, {}
    n_failed = 0
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
            if not failed and i < len(iseqs):
                uid_info[uid] = {'seq': iseqs[i], 'naive': naive, 'cdr3_length': cdr3}
    if n_failed > 0:
        print('  skipped %d failed queries (no annotation) in %s' % (n_failed, partition_fname), flush=True)

    sw_glfo, sw_antns, _ = utils.read_output(sw_cache_fname)
    sw_info, uid_sw_naives, uid_rearr_features = {}, {}, {}
    n_unalignable = 0
    for antn in sw_antns:  # the sw cache has no failed queries in it (waterer.write_cachefile)
        sw_seqs = get_ir_seqs(antn, 'sw cache')
        sw_naive = get_antn_key(antn, 'naive_seq', 'sw cache')
        vdj = tuple(get_antn_key(antn, '%s_gene' % r, 'sw cache') for r in utils.regions)
        # junction boundaries for the heavy merge guard
        v_3p_del, j_5p_del = get_antn_key(antn, 'v_3p_del', 'sw cache'), get_antn_key(antn, 'j_5p_del', 'sw cache')
        for i, uid in enumerate(antn['unique_ids']):
            sw_info[uid] = antn
            uid_rearr_features[uid] = {'vdj': vdj, 'v_3p_del': v_3p_del, 'j_5p_del': j_5p_del}
            naive = sw_naive
            if uid in uid_info and i < len(sw_seqs):  # re-pad into the partition frame
                naive = pad_sw_naive(sw_seqs[i], naive, uid_info[uid]['seq'],
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

    def _annotate(uids):
        antn = utils.synthesize_multi_seq_line_from_reco_info(uids, ant_info, warn=False)
        utils.remove_all_implicit_info(antn)
        utils.add_implicit_info(glfo, antn, reset_indel_genes=True)
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
        cluster = [uid for uid in cluster if uid in ant_info]
        good = [uid for uid in cluster if not ant_info[uid].get('invalid', False)]
        n_dropped += len(cluster) - len(good)
        if len(good) == 0:
            continue
        try:
            annotation_list.append(_annotate(good))
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


def run_jobs(specs, naive_threshold=None, overwrite=False, locus=None):
    """Run refinement on a list of group specs (from group_specs), writing each group's
    refined partition, with the production defaults (singleton-skip, junction guard, vdj
    override) that the standalone CLI and integrated pipeline both use. Groups whose
    refined output already exists are skipped unless <overwrite>. The naive threshold
    defaults to a locus-wide estimate over <specs>; when running a slice, pass one
    estimated over the full group list. Passing <locus> skips that estimate on a locus with
    no D gene, where no operator reads it."""
    from argparse import Namespace
    from partis import utils
    oargs = Namespace(overwrite=overwrite)
    if naive_threshold is None and (locus is None or utils.has_d_gene(locus)):
        naive_threshold = estimate_locuswide_threshold(specs)
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
            naive_threshold=naive_threshold,
            skip_singleton_merge=True, min_agreement=0.15, verbose=False)
        cfo = write_full_output(spec['refined_out'], inp['part_glfo'], refined, inp['uid_part_antns'])
        print('  timing: group %s total %.2f s' % (os.path.dirname(spec['refined_rel']), time.time() - tgroup), flush=True)
        n_run += 1
        for key, val in cfo.items():
            totals[key] += val
    print('  timing: %d groups total %.2f s' % (n_run, time.time() - tlocus), flush=True)
    if totals['n_dropped'] > 0 or totals['n_fallback'] > 0:  # per-group counts are easy to miss, so total them
        print('  %s refine output over %d groups: dropped %d unannotatable uids, %d clusters (%d uids) fell back to singletons'
              % (utils.wrnstr(), n_run, totals['n_dropped'], totals['n_fallback'], totals['n_fallback_uids']), flush=True)
