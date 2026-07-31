"""Post-partition refinement (three-step split/merge).

Pipeline-agnostic module: operates on a partition (list of clusters, each a list
of uids) plus per-sequence annotations and SW naives, and returns a refined
partition. Usable after any partis partition (standard or disjoint grouping).

  Step 1: Validated split (Jaccard proposes, per-sequence naive validates)
  Step 2: Incremental naive merge (merge only if mutation fingerprints agree)
  Step 3: Shared-descent split of over-merged clusters (guards fork on chain)

Entry point: refine_partition(); its defaults run step 3 with singleton-skip and
the rearrangement guard. Driven by the run-partition-refine-jobs action (and the
integrated --partition-refine flag), which run run_jobs() over the disjoint groups.
"""
import json
import math
import os
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


def compute_naive_frequencies(uid_sw_naives):
    naive_freq = defaultdict(int)
    for uid, naive in uid_sw_naives.items():
        naive_freq[naive] += 1
    return naive_freq


def get_cluster_fingerprint(uids, uid_to_muts_with_base):
    """Build mutation fingerprint for a cluster.

    Returns dict of {position: {base: count}} representing the mutation
    profile of the cluster. Positions mutated by multiple members with
    the same base are strong fingerprint positions.
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

    If min_fp_positions > 0, returns score=-1 when both clusters have
    fewer strong positions than the threshold (insufficient signal).
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

    # for asymmetric clusters (one large, one small), take the max
    # the small cluster's positions should appear in the large cluster
    # if they're from the same family
    return max(fwd, rev), n_strong_max


def step1_validated_split(partition, uid_to_muts, uid_sw_naives,
                          jaccard_threshold, naive_threshold, min_cluster_size=2):
    """Step 1: Propose splits with Jaccard, validate with naive."""
    result = []
    n_proposed = 0
    n_accepted = 0
    n_rejected = 0
    n_single_naive_skip = 0

    for cluster in partition:
        if len(cluster) < min_cluster_size:
            result.append(list(cluster))
            continue

        uid_list = list(cluster)
        n = len(uid_list)

        cluster_naives = {}
        for uid in uid_list:
            if uid in uid_sw_naives:
                cluster_naives[uid] = uid_sw_naives[uid]

        # optimization: if all sequences share the same naive (within threshold),
        # any proposed split will be rejected by naive validation. Skip Jaccard.
        unique_naives = list(set(cluster_naives.values()))
        if len(unique_naives) <= 1:
            result.append(list(cluster))
            n_single_naive_skip += 1
            continue
        # check if all unique naives are within threshold of each other
        all_similar = True
        for i in range(len(unique_naives)):
            for j in range(i + 1, len(unique_naives)):
                if hamming_frac(unique_naives[i], unique_naives[j]) > naive_threshold:
                    all_similar = False
                    break
            if not all_similar:
                break
        if all_similar:
            result.append(list(cluster))
            n_single_naive_skip += 1
            continue

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

        for i in range(n):
            for j in range(i + 1, n):
                m1 = uid_to_muts.get(uid_list[i])
                m2 = uid_to_muts.get(uid_list[j])
                if m1 is not None and m2 is not None:
                    if jaccard(m1, m2) >= jaccard_threshold:
                        union(i, j)

        components = defaultdict(list)
        for i in range(n):
            components[find(i)].append(uid_list[i])

        proposed_fragments = list(components.values())

        if len(proposed_fragments) <= 1:
            result.append(list(cluster))
            continue

        n_proposed += 1

        frag_naives = [get_fragment_naive(f, uid_sw_naives) for f in proposed_fragments]
        if any(fn is None for fn in frag_naives):  # a fragment with no sw naive cannot be validated; fail closed (keep the cluster intact rather than accept an unvalidated split)
            n_rejected += 1
            result.append(list(cluster))
            continue

        nf = len(proposed_fragments)
        frag_parent = list(range(nf))

        def frag_find(x):
            while frag_parent[x] != x:
                frag_parent[x] = frag_parent[frag_parent[x]]
                x = frag_parent[x]
            return x

        def frag_union(a, b):
            ra, rb = frag_find(a), frag_find(b)
            if ra != rb:
                frag_parent[ra] = rb

        for i in range(nf):
            for j in range(i + 1, nf):
                if hamming_frac(frag_naives[i], frag_naives[j]) <= naive_threshold:
                    frag_union(i, j)

        validated = defaultdict(list)
        for i in range(nf):
            validated[frag_find(i)].extend(proposed_fragments[i])

        validated_groups = list(validated.values())

        if len(validated_groups) > 1:
            n_accepted += 1
            result.extend(validated_groups)
        else:
            n_rejected += 1
            result.append(list(cluster))

    print('  step 1: %d proposed, %d accepted, %d rejected, %d skipped (single naive)' % (
        n_proposed, n_accepted, n_rejected, n_single_naive_skip), flush=True)
    print('  %d -> %d clusters' % (len(partition), len(result)), flush=True)
    return result


def step2_incremental_merge(split_partition, uid_info, uid_sw_naives,
                            uid_to_muts_with_base, naive_threshold,
                            min_agreement=0.15, min_fp_positions=0,
                            skip_singleton_merge=False):
    """Step 2: Incremental naive merge with fingerprint validation.

    For each CDR3 group, find clusters with similar naives (candidates),
    then only merge if their mutation fingerprints agree. This prevents
    false merges of unrelated sequences with similar naives.

    If skip_singleton_merge=True, skip merge attempts where BOTH clusters
    are singletons. True singletons (survived vsearch + HA) should stay
    singletons. Only attempt merges when at least one side has multiple
    members (i.e., one side has a real cluster-level fingerprint).
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
    n_skipped_singleton = 0
    n_skipped_difflen = 0

    for cdr3_len, frags in cdr3_frags.items():
        if len(frags) < 2:
            merged_partition.extend(frags)
            continue

        # precompute fingerprints and naives for all fragments
        frag_naives = [get_fragment_naive(f, uid_sw_naives) for f in frags]
        frag_fps = []
        for f in frags:
            fp, n_w = get_cluster_fingerprint(f, uid_to_muts_with_base)
            frag_fps.append((fp, n_w))

        # precompute fragment sizes for singleton skip
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
            # convert to numpy byte array for vectorized hamming
            seq_len = len(unique_naives[0])
            same_len = [n for n in unique_naives if len(n) == seq_len]
            diff_len = [n for n in unique_naives if len(n) != seq_len]

            # handle same-length naives with numpy
            if len(same_len) >= 2:
                arr = np.frombuffer(''.join(same_len).encode(), dtype=np.uint8).reshape(len(same_len), seq_len)
                n_byte = ord('N')
                # chunked pairwise to avoid huge memory for very large U
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
                    if len(n2_seq) != len(n1):  # never compared, so never merged: count so it isn't silent
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
    print('  step 2: %d naive candidates, %d accepted, %d rejected (fingerprint)%s' % (
        n_naive_candidates, n_accepted, n_rejected_fingerprint, skip_label))
    print('  %d merges, %d -> %d clusters' % (n_merges, len(split_partition), len(merged_partition)))
    return merged_partition


# ----------------------------------------------------------------------------
# Step-3 shared-descent split (non-transitive greedy centroid). For a cluster, a
# per-position base-specific mutation propensity g(pos,base) is measured from its
# members; the expected number of shared specific mutations for an unrelated pair
# is mu = sum over (pos,base) of g(pos,base)^2. Two cells are linked only when
# their observed shared specific-mutation count exceeds that chance null (Poisson
# upper tail < alpha), i.e. they share more mutations than convergence would
# produce, so the threshold self-calibrates per cluster.
# ----------------------------------------------------------------------------

# positions mutated in >= HOTSPOT_FRAC of the cluster (or in more than
# HOTSPOT_ABS members) are treated as uninformative (likely shared by chance)
# and excluded from the pair graph; this also bounds the pairwise cost to rare
# mutations. NOTE the abs bound isn't only a cost bound: once n > 2*HOTSPOT_ABS
# it excludes more than the frac rule does, so it tightens with cluster size.
HOTSPOT_FRAC = 0.5
HOTSPOT_ABS = 500
SHARED_DESCENT_ALPHA = 0.01  # link threshold (Poisson upper-tail cutoff) for the shared-descent test


def _poisson_upper_tail(k, mu):
    """P(Poisson(mu) >= k) = 1 - CDF(k-1, mu), via a stable running sum of the
    pmf terms exp(-mu) * mu^i / i! for i in 0..k-1 (no scipy)."""
    if k <= 0:
        return 1.0
    term = math.exp(-mu)  # i = 0
    cdf = term
    for i in range(1, k):
        term *= mu / i
        cdf += term
    cdf = min(cdf, 1.0)
    return max(0.0, 1.0 - cdf)


def shared_descent_pair_stats(cluster, uid_muts):
    """For one cluster, estimate the base-specific per-position mutation propensity
    g(pos,base) = (# members whose mutation at pos == base) / cluster_size, the
    null expected shared-count mu = sum over (pos,base) of g^2, and the observed
    shared specific-mutation count per member-index pair. Hotspot/consensus
    positions are excluded (informative positions only). Returns
    (members, pair_k, mu)."""
    members = list(cluster)
    n = len(members)
    inv = defaultdict(list)  # (pos, base): member indices
    for idx, uid in enumerate(members):
        for pos, base in uid_muts.get(uid, {}).items():
            inv[(pos, base)].append(idx)
    mu = 0.0
    pair_k = defaultdict(int)
    for idxs in inv.values():
        c = len(idxs)
        if c < 2:
            continue
        if c >= HOTSPOT_FRAC * n or c > HOTSPOT_ABS:
            continue  # skip uninformative (hotspot/consensus) position
        g = c / n
        mu += g * g
        for i in range(c):
            ai = idxs[i]
            for j in range(i + 1, c):
                bj = idxs[j]
                key = (ai, bj) if ai < bj else (bj, ai)
                pair_k[key] += 1
    return members, pair_k, mu


def split_by_shared_descent(cluster, uid_muts, alpha=SHARED_DESCENT_ALPHA):
    """Split one cluster by shared descent in a non-transitive greedy-centroid
    assignment: sort members by mutation count (descending); the first unassigned
    member seeds a centroid; every other unassigned member joins it iff the pair
    shares more specific mutations than the cluster chance null (Poisson upper tail
    < alpha); assigned members are removed (non-transitive -> no chaining).
    Unlinked members become singletons. Returns a list of sub-clusters (lists of
    uids)."""
    members, pair_k, mu = shared_descent_pair_stats(cluster, uid_muts)
    n = len(members)
    if n <= 1:
        return [list(members)]
    order = sorted(range(n), key=lambda i: -len(uid_muts.get(members[i], {})))
    assigned = [False] * n
    clusters = []
    for i in order:
        if assigned[i]:
            continue
        assigned[i] = True
        sub = [members[i]]
        for o in order:
            if assigned[o]:
                continue
            key = (i, o) if i < o else (o, i)
            if _poisson_upper_tail(pair_k.get(key, 0), mu) < alpha:
                sub.append(members[o])
                assigned[o] = True
        clusters.append(sub)
    return clusters


def _partition_has_real_d(uid_rearr_features):
    """True if any uid carries a real (non-placeholder) D gene (heavy chain).
    Light loci (igk/igl) have no real D, so this is False for them, which is how
    step 3 decides the heavy (guarded) vs light (relaxed) fork when light_chain is
    not passed explicitly."""
    if not uid_rearr_features:
        return False
    for feat in uid_rearr_features.values():
        if not isinstance(feat, dict):
            continue
        d_gene = feat.get('vdj', ('', '', ''))[1]
        if d_gene and 'x-x' not in d_gene and 'Dx' not in d_gene:
            return True
    return False


def step3_split(merged_partition, uid_info, uid_sw_naives,
                naive_freq=None, naive_freq_threshold=10,
                max_expansion_ratio=2.0,
                uid_rearr_features=None, rearrangement_guard=False,
                light_chain=None, alpha=SHARED_DESCENT_ALPHA):
    """Step 3: split over-merged clusters by shared descent
    (split_by_shared_descent). The guards fork on chain:

      HEAVY (light_chain False): a cluster is split only if its dominant SW naive
        occurs >= naive_freq_threshold times AND it fails both the D-gene guard
        (rearrangement_guard: real D present and >= 80% shared V+D+J -> keep) and
        the concentration-ratio guard (cluster_size / dominant-naive-count <=
        max_expansion_ratio -> keep). Otherwise kept intact.

      LIGHT (light_chain True): the naive-freq trigger and the concentration-ratio
        guard are skipped (they do not discriminate when every cell shares one
        light naive); every cluster of size >= 2 is passed to the proposer.

    light_chain: if None, inferred from the absence of real D genes in
    uid_rearr_features (heavy has real D; light does not).
    alpha: link threshold (Poisson upper-tail cutoff) for the shared-descent test
    (default SHARED_DESCENT_ALPHA).
    """
    if light_chain is None:
        light_chain = not _partition_has_real_d(uid_rearr_features)

    result = []
    n_resplit = n_skipped = n_expansion_skip = n_must_collision = 0
    n_rearr_guard_skip = n_input_seqs = 0

    for cluster in merged_partition:
        if len(cluster) < 2:
            result.append(cluster)
            n_skipped += 1
            continue

        # per-cell mutations vs the SW naive, for the proposer
        uid_muts = {}
        for uid in cluster:
            if uid in uid_sw_naives and uid in uid_info:
                uid_muts[uid] = get_mutations_with_base(uid_info[uid]['seq'], uid_sw_naives[uid])

        if not light_chain:
            # HEAVY guard block: only genuine collisions reach the proposer.
            rep_naive = None
            if naive_freq is not None and naive_freq_threshold is not None:
                rep_naive = get_fragment_naive(cluster, uid_sw_naives)
            if rep_naive is None or naive_freq.get(rep_naive, 0) < naive_freq_threshold:
                result.append(cluster)
                n_skipped += 1
                continue
            # D-gene guard: real D present and >= 80% share V+D+J, so keep (deep tree)
            if rearrangement_guard and uid_rearr_features is not None:
                vdj_counts = {}
                has_real_d = False
                n_with_feat = 0
                for uid in cluster:
                    feat = uid_rearr_features.get(uid)
                    if feat is not None:
                        vdj = feat['vdj']
                        vdj_counts[vdj] = vdj_counts.get(vdj, 0) + 1
                        n_with_feat += 1
                        d_gene = vdj[1]
                        if d_gene and 'x-x' not in d_gene and 'Dx' not in d_gene:
                            has_real_d = True
                if has_real_d and n_with_feat > 0 and max(vdj_counts.values()) / n_with_feat >= 0.80:
                    result.append(cluster)
                    n_rearr_guard_skip += 1
                    continue
            # concentration-ratio guard: star-like expansion, so keep
            if rep_naive is not None and max_expansion_ratio > 0:
                naive_count = sum(1 for uid in cluster
                                  if uid in uid_sw_naives and uid_sw_naives[uid] == rep_naive)
                ratio = len(cluster) / naive_count if naive_count > 0 else 999
                if ratio <= max_expansion_ratio:
                    result.append(cluster)
                    n_expansion_skip += 1
                    continue
            n_must_collision += 1

        # shared-descent split (light: every cluster; heavy: those that passed guards)
        n_input_seqs += len(cluster)
        pieces = split_by_shared_descent(list(cluster), uid_muts, alpha)
        result.extend(pieces)
        if len(pieces) > 1:
            n_resplit += 1

    n_result_singletons = sum(1 for c in result if len(c) == 1)
    chain = 'light (relaxed)' if light_chain else 'heavy (guarded)'
    coll_label = ', %d must-collision' % n_must_collision if n_must_collision > 0 else ''
    exp_label = ', %d expansion-skip' % n_expansion_skip if n_expansion_skip > 0 else ''
    rearr_label = ', %d rearr-guard-skip' % n_rearr_guard_skip if (not light_chain and rearrangement_guard) else ''
    print('  step 3 [%s, alpha=%.3g]: %d re-split (%d seqs processed), %d skipped%s%s%s, %d -> %d clusters (%d singletons)' % (
        chain, alpha, n_resplit, n_input_seqs, n_skipped, coll_label, exp_label, rearr_label,
        len(merged_partition), len(result), n_result_singletons), flush=True)
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


def estimate_jaccard_threshold(partition, uid_to_muts, seed=0):
    """Adaptive step-1 Jaccard threshold: p25 of sampled within-cluster pairwise
    Jaccard (fallback 0.10 if no signal). Fixed-seed local RNG keeps it deterministic."""
    import random
    rng = random.Random(seed)
    all_jaccards = []
    for cluster in partition:
        if len(cluster) < 3:
            continue
        uid_list = list(cluster)
        n = len(uid_list)
        if n * (n - 1) // 2 <= 100:
            pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        else:
            pairs = [tuple(rng.sample(range(n), 2)) for _ in range(100)]
        for i, j in pairs:
            m1 = uid_to_muts.get(uid_list[i])
            m2 = uid_to_muts.get(uid_list[j])
            if m1 is not None and m2 is not None:
                all_jaccards.append(jaccard(m1, m2))
    if all_jaccards:
        all_jaccards.sort()
        return all_jaccards[int(len(all_jaccards) * 0.25)]
    return 0.10


def refine_partition(partition, uid_info, uid_sw_naives, uid_rearr_features=None,
                     jaccard_threshold=None, naive_threshold=None,
                     min_agreement=0.15, min_fp_positions=0, skip_singleton_merge=True,
                     naive_freq_threshold=10, max_expansion_ratio=2.0, min_cluster_size=2,
                     rearrangement_guard=True, light_chain=None, alpha=SHARED_DESCENT_ALPHA,
                     verbose=True, random_seed=None):
    """Run the three-step refinement and return the refined partition.

    partition: list of clusters, each a list of uids.
    uid_info: uid -> {'seq', 'naive', 'cdr3_length'} (partition annotations).
    uid_sw_naives: uid -> per-sequence SW naive_seq.
    uid_rearr_features: uid -> {'vdj': (v, d, j)} (required if rearrangement_guard).

    Step 3 splits over-merged clusters by shared descent (step3_split), with guards
    that fork on chain: heavy (D-gene present) keeps the D-gene + concentration
    guards; light relaxes them. light_chain: if None, inferred from D-gene presence.
    alpha: link threshold (Poisson upper-tail cutoff) for the step-3 shared-descent test.

    random_seed: seeds the global RNG. NOTE estimate_jaccard_threshold uses its own
    fixed-seed RNG, so refinement is already deterministic without this.
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
        raise Exception('%d uids whose sw naive isn\'t in the partition frame (e.g. %s); pass naives through pad_sw_naive()' % (len(badfo), badfo[:3]))
    n_no_naive = sum(1 for uid in all_uids if uid not in uid_sw_naives)
    if n_no_naive > 0:  # surface incomplete sw-naive coverage: validation keys off naives, so these seqs' fragments are kept intact (step 1 fails closed)
        print('  warning: %d/%d input seqs have no sw naive (refine validates on naives; their fragments are kept intact)' % (n_no_naive, len(all_uids)), flush=True)
    uid_to_muts, uid_to_muts_with_base = {}, {}
    for uid, info in uid_info.items():
        uid_to_muts[uid] = get_mutations(info['seq'], info['naive'])
        uid_to_muts_with_base[uid] = get_mutations_with_base(info['seq'], info['naive'])

    naive_thresh = (naive_threshold if naive_threshold is not None
                    else estimate_naive_threshold(partition, uid_sw_naives))
    jaccard_thresh = (jaccard_threshold if jaccard_threshold is not None
                      else estimate_jaccard_threshold(partition, uid_to_muts))

    if verbose:
        print('\n=== Step 1: Validated split (Jaccard >= %.4f, naive <= %.4f) ===' % (
            jaccard_thresh, naive_thresh), flush=True)
    split_partition = step1_validated_split(
        partition, uid_to_muts, uid_sw_naives, jaccard_thresh, naive_thresh, min_cluster_size)

    if verbose:
        print('\n=== Step 2: Incremental merge (naive <= %.4f, min_agreement %.2f) ===' % (
            naive_thresh, min_agreement), flush=True)
    merged_partition = step2_incremental_merge(
        split_partition, uid_info, uid_sw_naives, uid_to_muts_with_base, naive_thresh,
        min_agreement, min_fp_positions, skip_singleton_merge=skip_singleton_merge)

    naive_freq = compute_naive_frequencies(uid_sw_naives)
    if light_chain is None:
        light_chain = not _partition_has_real_d(uid_rearr_features)
    if verbose:
        chain = 'light (relaxed)' if light_chain else 'heavy (guarded)'
        print('\n=== Step 3: shared-descent split (%s, alpha=%.3g, naive_freq >= %d) ===' % (
            chain, alpha, naive_freq_threshold), flush=True)
    final_partition = step3_split(
        merged_partition, uid_info, uid_sw_naives,
        naive_freq=naive_freq, naive_freq_threshold=naive_freq_threshold,
        max_expansion_ratio=max_expansion_ratio,
        uid_rearr_features=uid_rearr_features, rearrangement_guard=rearrangement_guard,
        light_chain=light_chain, alpha=alpha)

    return final_partition


# ----------------------------------------------------------------------------
# File I/O layer (reads already-generated partition + sw-cache, writes full
# partis output). Kept separate from the pure algorithm above; imports
# partis.utils lazily so refine_partition() stays dependency-light. Works as a
# standalone per-group job: no SW, vsearch, bcrham, or HMM is re-run.
# ----------------------------------------------------------------------------

def find_sw_offset(part_seq, sw_seq):
    """Offset of <sw_seq> within <part_seq>, or -1. Falls back to an N-tolerant scan (N in
    part_seq matches anything) since partis masks the odd base in the padded frame, e.g. a
    1-base jf_insertion at the 3' end, which defeats an exact find."""
    ioff = part_seq.find(sw_seq)
    if ioff >= 0:
        return ioff
    for ioff in range(len(part_seq) - len(sw_seq) + 1):
        if all(p == 'N' or p == s for p, s in zip(part_seq[ioff:ioff + len(sw_seq)], sw_seq)):
            return ioff
    return -1


def pad_sw_naive(sw_seq, sw_naive, part_seq):
    """Re-pad a per-sequence sw naive into the partition frame. The partition annotation is
    padded (N) relative to the sw annotation, so comparing the two frames position-wise is a
    frame shift; everything downstream keys mutations by position across cluster members, so
    they all have to be in the partition frame. Returns None if <sw_seq> can't be located in
    <part_seq>, i.e. the offset is unknown and the caller should drop this uid's naive."""
    if len(part_seq) == len(sw_naive):
        return sw_naive
    ioff = find_sw_offset(part_seq, sw_seq)
    if ioff < 0:
        return None
    nright = len(part_seq) - ioff - len(sw_seq)
    if nright < 0:
        return None
    padded = 'N' * ioff + sw_naive + 'N' * nright
    if len(padded) != len(part_seq):  # shm indel: different coord systems, padding can't align them
        return None
    return padded


def read_refine_inputs(partition_fname, sw_cache_fname):
    """Read an existing partition file and sw-cache (no recompute) and build the
    inputs refine_partition() needs, plus the per-sequence sw_info and the
    partition annotations used to write full output. Returns a dict with keys:
    glfo, part_glfo, partition, uid_info, uid_sw_naives, uid_rearr_features,
    sw_info, uid_part_antns."""
    from partis import utils
    part_glfo, part_antns, cpath = utils.read_output(partition_fname)
    partition = [list(c) for c in (cpath.best() if cpath is not None else [])]
    uid_info, uid_part_antns = {}, {}
    for antn in part_antns:
        naive = antn.get('naive_seq')
        cdr3 = antn.get('cdr3_length', 0)
        iseqs = antn.get('seqs', antn.get('input_seqs', []))  # indel-reversed: these are what line up with naive_seq
        for i, uid in enumerate(antn['unique_ids']):
            uid_part_antns[uid] = antn  # write path synthesizes from these, so output stays in the input's padded frame
            if naive is not None and i < len(iseqs):
                uid_info[uid] = {'seq': iseqs[i], 'naive': naive, 'cdr3_length': cdr3}

    sw_glfo, sw_antns, _ = utils.read_output(sw_cache_fname)
    sw_info, uid_sw_naives, uid_rearr_features = {}, {}, {}
    n_unalignable = 0
    for antn in sw_antns:
        sw_seqs = antn.get('seqs', antn.get('input_seqs', []))
        for i, uid in enumerate(antn['unique_ids']):
            sw_info[uid] = antn
            uid_rearr_features[uid] = {'vdj': (antn.get('v_gene', ''),
                                               antn.get('d_gene', ''),
                                               antn.get('j_gene', ''))}
            naive = antn['naive_seq']
            if uid in uid_info and i < len(sw_seqs):  # re-pad into the partition frame
                naive = pad_sw_naive(sw_seqs[i], naive, uid_info[uid]['seq'])
            if naive is None:  # unknown offset: drop, refine fails closed on missing naives
                n_unalignable += 1
                continue
            uid_sw_naives[uid] = naive
    if n_unalignable > 0:
        print('  %s couldn\'t align %d sw naives to the partition frame (dropped)' % (utils.wrnstr(), n_unalignable), flush=True)
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
    and annotation list always match -- paired clustering requires an annotation for every
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
        # drop invalid uids so one bad annotation does not shatter the cluster
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
def group_specs(disjoint_dir, groups, locus):
    """Per-group refine I/O paths. Refine runs on the HA re-partition if it exists,
    else the vsearch partition; output is partition-refine-<locus>.yaml. Only groups
    whose input partition and sw-cache exist are returned."""
    specs = []
    n_harep = n_vsearch = 0
    for group in groups:
        fasta_dir = os.path.dirname(group['fasta_path'])
        harep_p = '%s/%s/ha-repartition-%s.yaml' % (disjoint_dir, fasta_dir, locus)
        vsearch_p = '%s/%s/partition-%s.yaml' % (disjoint_dir, fasta_dir, locus)
        sw = '%s/%s/sw-cache-%s.yaml' % (disjoint_dir, fasta_dir, locus)
        inp = harep_p if os.path.exists(harep_p) else vsearch_p
        if not (os.path.exists(inp) and os.path.exists(sw)):
            continue
        if inp == harep_p:
            n_harep += 1
        else:
            n_vsearch += 1
        specs.append({'group': group, 'input': inp, 'sw_cache': sw,
                      'refined_out': '%s/%s/partition-refine-%s.yaml' % (disjoint_dir, fasta_dir, locus),
                      'refined_rel': '%s/partition-refine-%s.yaml' % (fasta_dir, locus)})
    if n_harep > 0 and n_vsearch > 0:  # some ha-repartition jobs didn't finish
        from partis import utils
        print('  %s refine input is a MIX: %d groups from ha-repartition, %d from vsearch' % (utils.wrnstr(), n_harep, n_vsearch), flush=True)
    return specs


def estimate_locuswide_thresholds(specs):
    """Estimate the step-1 naive/jaccard thresholds over the full-locus partition (all cdr3
    groups in <specs> merged). Global statistics, so estimate once here and pass into each
    per-group refine rather than letting refine_partition estimate them per-group."""
    partition, uid_sw_naives, uid_to_muts = [], {}, {}
    for spec in specs:
        inp = read_refine_inputs(spec['input'], spec['sw_cache'])
        partition.extend(inp['partition'])
        uid_sw_naives.update(inp['uid_sw_naives'])
        for uid, info in inp['uid_info'].items():
            uid_to_muts[uid] = get_mutations(info['seq'], info['naive'])
    return estimate_naive_threshold(partition, uid_sw_naives), estimate_jaccard_threshold(partition, uid_to_muts)


def run_jobs(specs, naive_threshold=None, jaccard_threshold=None):
    """Run refinement on a list of group specs (from group_specs), writing each group's
    refined partition. Production defaults (shared-descent step 3, singleton-skip, rearrangement
    guard) -- the same config the standalone CLI and integrated step use. Groups whose
    refined output already exists are skipped. naive/jaccard thresholds default to a
    locus-wide estimate over <specs>; when running a slice, pass thresholds estimated over
    the full group list."""
    from partis import utils
    if naive_threshold is None or jaccard_threshold is None:
        naive_threshold, jaccard_threshold = estimate_locuswide_thresholds(specs)
    totals, n_run = defaultdict(int), 0
    for spec in specs:
        if os.path.exists(spec['refined_out']):
            continue
        inp = read_refine_inputs(spec['input'], spec['sw_cache'])
        refined = refine_partition(
            inp['partition'], inp['uid_info'], inp['uid_sw_naives'],
            uid_rearr_features=inp['uid_rearr_features'],
            naive_threshold=naive_threshold, jaccard_threshold=jaccard_threshold,
            skip_singleton_merge=True, min_agreement=0.15,
            naive_freq_threshold=10, rearrangement_guard=True, verbose=False)
        cfo = write_full_output(spec['refined_out'], inp['part_glfo'], refined, inp['uid_part_antns'])
        n_run += 1
        for key, val in cfo.items():
            totals[key] += val
    if totals['n_dropped'] > 0 or totals['n_fallback'] > 0:  # per-group counts are easy to miss, so total them
        print('  %s refine output over %d groups: dropped %d unannotatable uids, %d clusters (%d uids) fell back to singletons'
              % (utils.wrnstr(), n_run, totals['n_dropped'], totals['n_fallback'], totals['n_fallback_uids']), flush=True)
