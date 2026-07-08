"""Post-partition refinement (three-step split/merge).

Pipeline-agnostic module: operates on a partition (list of clusters, each a list
of uids) plus per-sequence annotations and SW naives, and returns a refined
partition. Usable after any partis partition (standard or disjoint grouping).

  Step 1: Validated split (Jaccard proposes, per-sequence naive validates)
  Step 2: Incremental naive merge (merge only if mutation fingerprints agree)
  Step 3: EJ centroid split on collision groups (naive_freq >= threshold)

Entry point: refine_partition(). Its defaults are the validated production
config (relative-mode step 3, singleton-skip, rearrangement guard). See
refinement-pipeline.md for the derivation of every threshold. The standalone
CLI lives in bin/partition-refinement.py, which imports this module.
"""
import json
from collections import defaultdict
import numpy as np


def load_partition(path):
    with open(path) as f:
        data = json.load(f)
    return data


def load_sw_naives(sw_cache_path):
    with open(sw_cache_path) as f:
        data = json.load(f)
    uid_naives = {}
    for evt in data['events']:
        if evt.get('invalid'):
            continue
        naive = evt.get('naive_seq')
        if naive is None:
            continue
        for uid in evt['unique_ids']:
            uid_naives[uid] = naive
    return uid_naives


def load_rearrangement_features(sw_cache_path):
    """Load per-sequence rearrangement features from sw-cache.

    Returns dict: uid -> {'vdj': (v_gene, d_gene, j_gene)}.
    The VDJ tuple is used for the step-3 D-gene majority guard. (Boundary
    features v_3p_del/j_5p_del/dj_ins were explored for a rearrangement-first
    split but reverted -- guard only. See refinement-pipeline.md.)
    """
    with open(sw_cache_path) as f:
        data = json.load(f)
    uid_features = {}
    for evt in data['events']:
        if evt.get('invalid'):
            continue
        feat = {
            'vdj': (evt.get('v_gene', ''), evt.get('d_gene', ''), evt.get('j_gene', '')),
        }
        for uid in evt['unique_ids']:
            uid_features[uid] = feat
    return uid_features


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


def load_annotations(data):
    uid_info = {}
    for antn in data['events']:
        if antn.get('invalid'):
            continue
        naive = antn.get('naive_seq')
        cdr3_len = antn.get('cdr3_length', 0)
        input_seqs = antn.get('input_seqs', [])
        for i, uid in enumerate(antn['unique_ids']):
            seq = input_seqs[i] if i < len(input_seqs) else None
            if naive is not None and seq is not None:
                uid_info[uid] = {'seq': seq, 'naive': naive, 'cdr3_length': cdr3_len}
    return uid_info


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
    muts = {}
    for i, (s, n) in enumerate(zip(seq, naive)):
        if s != n and s != 'N' and n != 'N':
            muts[i] = s
    return muts


def enhanced_jaccard(muts1, muts2):
    all_pos = set(muts1.keys()) | set(muts2.keys())
    if len(all_pos) == 0:
        return 0.0
    shared = sum(1 for p in muts1 if p in muts2 and muts1[p] == muts2[p])
    return shared / len(all_pos)


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
            if frag_naives[i] is None:
                continue
            for j in range(i + 1, nf):
                if frag_naives[j] is None:
                    continue
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

        local_min_agreement = min_agreement

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
                    elif agreement >= local_min_agreement:
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
                # chunked pairwise to avoid huge memory for very large U
                chunk_size = 2000
                for ci in range(0, len(same_len), chunk_size):
                    ci_end = min(ci + chunk_size, len(same_len))
                    chunk_i = arr[ci:ci_end]
                    for cj in range(ci, len(same_len), chunk_size):
                        cj_end = min(cj + chunk_size, len(same_len))
                        chunk_j = arr[cj:cj_end]
                        dists = (chunk_i[:, None, :] != chunk_j[None, :, :]).sum(axis=2) / seq_len
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
                                    elif agreement >= local_min_agreement:
                                        union(i, j)
                                        n_accepted += 1
                                    else:
                                        n_rejected_fingerprint += 1

            # handle different-length naives with Python fallback (rare)
            for n1 in diff_len:
                for n2_seq in unique_naives:
                    if n2_seq <= n1 or len(n2_seq) != len(n1):
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
                            elif agreement >= local_min_agreement:
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
    print('  step 2: %d naive candidates, %d accepted, %d rejected (fingerprint)%s' % (
        n_naive_candidates, n_accepted, n_rejected_fingerprint, skip_label))
    print('  %d merges, %d -> %d clusters' % (n_merges, len(split_partition), len(merged_partition)))
    return merged_partition


def ej_centroid_cluster(uid_list, uid_muts, threshold, rescue_threshold=None):
    """Greedy centroid clustering using enhanced Jaccard.

    If rescue_threshold is set, after the first pass any non-centroid
    singletons are reassigned to the best-matching centroid if
    EJ >= rescue_threshold. This recovers star-tree family members
    that had low EJ to the centroid but nonzero shared ancestry.
    """
    if len(uid_list) <= 1:
        return [list(uid_list)]

    uid_nmuts = [(uid, len(uid_muts.get(uid, {}))) for uid in uid_list]
    uid_nmuts.sort(key=lambda x: -x[1])

    assigned = set()
    clusters = []
    centroid_uids = []

    for centroid_uid, _ in uid_nmuts:
        if centroid_uid in assigned:
            continue

        centroid_muts = uid_muts.get(centroid_uid)
        if centroid_muts is None:
            clusters.append([centroid_uid])
            assigned.add(centroid_uid)
            continue

        cluster = [centroid_uid]
        assigned.add(centroid_uid)
        centroid_uids.append(centroid_uid)

        for other_uid, _ in uid_nmuts:
            if other_uid in assigned:
                continue
            other_muts = uid_muts.get(other_uid)
            if other_muts is None:
                continue
            if enhanced_jaccard(centroid_muts, other_muts) >= threshold:
                cluster.append(other_uid)
                assigned.add(other_uid)

        clusters.append(cluster)

    for uid in uid_list:
        if uid not in assigned:
            clusters.append([uid])

    # rescue pass: reassign non-centroid singletons to best-matching centroid
    if rescue_threshold is not None and rescue_threshold < threshold and len(centroid_uids) > 1:
        # build centroid index
        centroid_to_idx = {}
        for idx, c in enumerate(clusters):
            if c[0] in centroid_uids:
                centroid_to_idx[c[0]] = idx

        # build cluster fingerprints for each centroid cluster
        centroid_fps = {}
        for cent_uid in centroid_uids:
            cidx = centroid_to_idx[cent_uid]
            fp, n_w = get_cluster_fingerprint(clusters[cidx], uid_muts)
            centroid_fps[cent_uid] = (fp, n_w)

        rescued = 0
        non_centroid_clusters = []
        for idx, c in enumerate(clusters):
            if c[0] in centroid_uids:
                continue  # centroid clusters handled separately
            if len(c) == 1:
                uid = c[0]
                uid_m = uid_muts.get(uid)
                if uid_m is None:
                    non_centroid_clusters.append(c)
                    continue
                # build singleton fingerprint
                sing_fp = defaultdict(lambda: defaultdict(int))
                for pos, base in uid_m.items():
                    sing_fp[pos][base] = 1
                sing_n = 1

                best_agreement = 0.0
                best_centroid = None
                for cent_uid in centroid_uids:
                    cfp, cn = centroid_fps[cent_uid]
                    agreement, _ = fingerprint_agreement(
                        sing_fp, sing_n, cfp, cn, min_fp_positions=0)
                    if agreement > best_agreement:
                        best_agreement = agreement
                        best_centroid = cent_uid
                if best_agreement >= rescue_threshold and best_centroid is not None:
                    clusters[centroid_to_idx[best_centroid]].append(uid)
                    # update the centroid fingerprint with the rescued member
                    for pos, base in uid_m.items():
                        centroid_fps[best_centroid][0][pos][base] += 1
                    centroid_fps[best_centroid] = (centroid_fps[best_centroid][0],
                                                   centroid_fps[best_centroid][1] + 1)
                    rescued += 1
                else:
                    non_centroid_clusters.append(c)
            else:
                non_centroid_clusters.append(c)

        # rebuild: centroid clusters (possibly grown by rescue) + non-rescued
        final = [c for c in clusters if c[0] in centroid_uids]
        final.extend(non_centroid_clusters)
        return final

    return clusters


def estimate_group_ej_noise(uid_list, uid_muts, max_sample=500):
    """Sample pairwise EJ within a group to estimate the noise floor.

    Returns (median, p75) of sampled pairwise EJ.
    - p75: used by 'adaptive' mode as the threshold (above which only
      ~25% of random pairs fall)
    - median: used by 'relative' mode as the noise baseline
    """
    import random
    uids_with_muts = [u for u in uid_list if u in uid_muts and len(uid_muts[u]) > 0]
    n = len(uids_with_muts)
    if n < 2:
        return 0.0, 0.0

    ejs = []
    n_pairs = n * (n - 1) // 2
    if n_pairs <= max_sample:
        for i in range(n):
            for j in range(i + 1, n):
                ejs.append(enhanced_jaccard(uid_muts[uids_with_muts[i]],
                                           uid_muts[uids_with_muts[j]]))
    else:
        for _ in range(max_sample):
            i, j = random.sample(range(n), 2)
            ejs.append(enhanced_jaccard(uid_muts[uids_with_muts[i]],
                                       uid_muts[uids_with_muts[j]]))

    if not ejs:
        return 0.0, 0.0
    ejs.sort()
    median = ejs[len(ejs) // 2]
    p75 = ejs[int(len(ejs) * 0.75)]
    return median, p75


def step3_ej_centroid(merged_partition, uid_info, uid_sw_naives,
                      threshold=0.10, naive_freq=None, naive_freq_threshold=None,
                      rescue_threshold=None, step3_mode='fixed',
                      min_mutations=0, ej_margin=0.05,
                      max_expansion_ratio=2.0, max_family_size=None,
                      uid_rearr_features=None, rearrangement_guard=False):
    """Step 3: EJ centroid split on collision groups.

    step3_mode controls how the EJ threshold is determined:
      'fixed':    use the fixed threshold for all groups (original behavior)
      'gate':     skip groups where avg mutations < min_mutations
      'adaptive': set threshold = max(fixed, p75 of sampled pairwise EJ)
      'relative': set threshold = max(fixed, group_median_EJ + ej_margin)

    When rearrangement_guard is True and uid_rearr_features is provided:
      - V+D+J majority guard: if >= 80% of cluster members agree on V+D+J
        gene assignment, skip splitting (rearrangement-consistent, one family).
        Replaces the concentration ratio guard.

    When rearrangement_guard is False, uses the original two-level expansion
    guard (concentration ratio + max_family_size).
    """
    result = []
    n_resplit = 0
    n_skipped = 0
    n_gated = 0
    n_expansion_skip = 0
    n_must_collision = 0
    n_rearr_guard_skip = 0
    n_input_seqs = 0

    for cluster in merged_partition:
        if len(cluster) < 2:
            result.append(cluster)
            n_skipped += 1
            continue

        # check naive frequency
        should_split = False
        rep_naive = None
        if naive_freq is not None and naive_freq_threshold is not None:
            rep_naive = get_fragment_naive(cluster, uid_sw_naives)
            if rep_naive is not None and naive_freq.get(rep_naive, 0) >= naive_freq_threshold:
                should_split = True

        if not should_split:
            result.append(cluster)
            n_skipped += 1
            continue

        if rearrangement_guard and uid_rearr_features is not None:
            # D-gene guard: additional protection for IGH deep trees.
            # Only triggers when real D genes are present AND V+D+J
            # majority >= 80%. On light chain (no real D gene), this
            # does not trigger and falls through to the expansion guard.
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
            if has_real_d and n_with_feat > 0:
                top_count = max(vdj_counts.values())
                if top_count / n_with_feat >= 0.80:
                    result.append(cluster)
                    n_rearr_guard_skip += 1
                    continue

        # two-level expansion guard
        rep_freq = naive_freq.get(rep_naive, 0) if rep_naive else 0

        # level 1: if naive_freq > max_family_size, must be collision
        if max_family_size is not None and rep_freq > max_family_size:
            n_must_collision += 1
            # definitely collision, proceed to splitting
        else:
            # level 2: concentration ratio check (fixed 2.0). Adaptive p50
            # from data was tested but reverted -- 2.0 is stable across
            # datasets (see refinement-pipeline.md).
            ratio_threshold = max_expansion_ratio
            if rep_naive is not None and ratio_threshold > 0:
                naive_count_in_cluster = sum(1 for uid in cluster
                                            if uid in uid_sw_naives and uid_sw_naives[uid] == rep_naive)
                ratio = len(cluster) / naive_count_in_cluster if naive_count_in_cluster > 0 else 999
                if ratio <= ratio_threshold:
                    result.append(cluster)
                    n_expansion_skip += 1
                    continue

        # EJ centroid split (runs on clusters that passed through the guard)
        uid_muts = {}
        for uid in cluster:
            if uid in uid_sw_naives and uid in uid_info:
                naive = uid_sw_naives[uid]
                seq = uid_info[uid]['seq']
                uid_muts[uid] = get_mutations_with_base(seq, naive)

        # signal gate
        if step3_mode == 'gate' and min_mutations > 0:
            avg_muts = sum(len(m) for m in uid_muts.values()) / max(len(uid_muts), 1)
            if avg_muts < min_mutations:
                result.append(cluster)
                n_gated += 1
                continue

        # adaptive/relative threshold
        if step3_mode in ('adaptive', 'relative'):
            median_ej, p75_ej = estimate_group_ej_noise(list(cluster), uid_muts)
            if step3_mode == 'adaptive':
                local_threshold = max(threshold, p75_ej)
            else:
                local_threshold = max(threshold, median_ej + ej_margin)
        else:
            local_threshold = threshold

        n_input_seqs += len(cluster)
        pieces = ej_centroid_cluster(list(cluster), uid_muts, local_threshold, rescue_threshold)
        result.extend(pieces)
        if len(pieces) > 1:
            n_resplit += 1

    n_result_singletons = sum(1 for c in result if len(c) == 1)
    mode_label = ' [%s]' % step3_mode if step3_mode != 'fixed' else ''
    gate_label = ', %d gated (low signal)' % n_gated if n_gated > 0 else ''
    exp_label = ', %d expansion-skip' % n_expansion_skip if n_expansion_skip > 0 else ''
    coll_label = ', %d must-collision' % n_must_collision if n_must_collision > 0 else ''
    rearr_label = ''
    if rearrangement_guard:
        rearr_label = ', %d rearr-guard-skip' % n_rearr_guard_skip
    print('  step 3: %d re-split (%d seqs processed), %d skipped%s%s%s%s, %d -> %d clusters (%d singletons)%s' % (
        n_resplit, n_input_seqs, n_skipped, gate_label, coll_label, exp_label, rearr_label, len(merged_partition), len(result), n_result_singletons, mode_label), flush=True)
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


def estimate_jaccard_threshold(partition, uid_to_muts):
    """Adaptive step-1 Jaccard threshold: p25 of sampled within-cluster pairwise
    Jaccard (fallback 0.10 if no signal)."""
    import random
    all_jaccards = []
    for cluster in partition:
        if len(cluster) < 3:
            continue
        uid_list = list(cluster)
        n = len(uid_list)
        if n * (n - 1) // 2 <= 100:
            pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        else:
            pairs = [tuple(random.sample(range(n), 2)) for _ in range(100)]
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
                     max_family_size=None, jaccard_threshold=None, naive_threshold=None,
                     min_agreement=0.15, min_fp_positions=0, skip_singleton_merge=True,
                     step3_threshold=0.10, step3_mode='relative', min_mutations=20,
                     ej_margin=0.05, rescue_threshold=None, naive_freq_threshold=10,
                     max_expansion_ratio=2.0, min_cluster_size=2,
                     rearrangement_guard=True, verbose=True):
    """Run the three-step refinement and return the refined partition.

    partition: list of clusters, each a list of uids.
    uid_info: uid -> {'seq', 'naive', 'cdr3_length'} (partition annotations).
    uid_sw_naives: uid -> per-sequence SW naive_seq.
    uid_rearr_features: uid -> {'vdj': (v, d, j)} (required if rearrangement_guard).
    max_family_size: from cluster_size.csv, for the step-3 level-1 expansion guard.

    Defaults are the validated production config (relative-mode step 3,
    singleton-skip, rearrangement guard). See refinement-pipeline.md.
    """
    partition = [list(c) for c in partition]
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
    if verbose:
        n_high = sum(1 for f in naive_freq.values() if f >= naive_freq_threshold)
        print('\n=== Step 3: EJ centroid (threshold %.2f, naive_freq >= %d, %d qualifying naives, mode=%s) ===' % (
            step3_threshold, naive_freq_threshold, n_high, step3_mode), flush=True)
    final_partition = step3_ej_centroid(
        merged_partition, uid_info, uid_sw_naives, step3_threshold, naive_freq,
        naive_freq_threshold, rescue_threshold, step3_mode, min_mutations, ej_margin,
        max_expansion_ratio, max_family_size,
        uid_rearr_features=uid_rearr_features, rearrangement_guard=rearrangement_guard)

    return final_partition
