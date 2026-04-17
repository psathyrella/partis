from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import sys
import glob
import yaml
import collections

from . import utils

MANIFEST_FNAME = 'manifest.yaml'

# ----------------------------------------------------------------------------------------
def group_sequences_by_cdr3_length(annotation_list):
    # group uids and their input sequences by cdr3 length from sw annotation list
    # also extracts naive_seq for hfrac sub-grouping
    # returns {cdr3_length : [{'name': uid, 'seq': seq, 'naive_seq': naive_seq}, ...]}, n_failed
    seqfo_map = {}  # uid -> {'name': uid, 'seq': seq, 'cdr3_length': int, 'naive_seq': str}
    n_failed = 0
    for line in annotation_list:
        if 'cdr3_length' not in line or line['cdr3_length'] is None:
            n_failed += len(line['unique_ids'])
            continue
        naive_seq = line.get('naive_seq', '')
        for uid, seq in zip(line['unique_ids'], line['input_seqs']):
            seqfo_map[uid] = {'name' : uid, 'seq' : seq, 'cdr3_length' : line['cdr3_length'], 'naive_seq' : naive_seq}
    if n_failed > 0:
        print('  %s %d sequences had no cdr3_length and were excluded from grouping' % (utils.color('yellow', 'warning'), n_failed))
    if len(seqfo_map) == 0:
        return collections.OrderedDict(), n_failed
    uid_groups = utils.group_seqs_by_value(list(seqfo_map.keys()), lambda u: seqfo_map[u]['cdr3_length'], return_values=True)
    groups = collections.OrderedDict()
    for c3len, uids in sorted(uid_groups, key=lambda x: x[0]):
        groups[c3len] = [seqfo_map[u] for u in uids]
    return groups, n_failed

# ----------------------------------------------------------------------------------------
def _read_vsearch_uc_with_centroids(cluster_file):
    # parse vsearch UC output preserving centroid info per cluster
    # returns list of (centroid_uid, [member_uids]) pairs in cluster order
    # vsearch writes an 'S' row first for each cluster (centroid), then 'H' rows (hits)
    import csv
    clusters = {}  # cluster_id -> {'centroid': uid, 'members': [uids]}
    with open(cluster_file) as f:
        reader = csv.reader(f, delimiter=str('\t'))
        for row in reader:
            if len(row) < 9:
                continue
            rtype = row[0]
            if rtype == 'C':
                continue
            cid = int(row[1])
            uid = row[8]
            if cid not in clusters:
                clusters[cid] = {'centroid': None, 'members': []}
            if rtype == 'S':
                clusters[cid]['centroid'] = uid
            clusters[cid]['members'].append(uid)
    return [(clusters[cid]['centroid'], clusters[cid]['members']) for cid in sorted(clusters)]

# ----------------------------------------------------------------------------------------
def _build_round2_cc_cmd(centroid_naives, round2_threshold, workdir):
    # build a vsearch --allpairs_global command for round 2 CC on centroid naive sequences
    # returns (cmdfo dict, pairs_outfname) or (None, None) if fewer than 2 centroids
    if len(centroid_naives) < 2:
        return None, None
    utils.prep_dir(workdir)
    infname = '%s/centroids.fa' % workdir
    outfname = '%s/pairs.tsv' % workdir
    with open(infname, 'w') as f:
        for name, seq in centroid_naives.items():
            f.write('>%s\n%s\n' % (name, seq))
    vsearch_binary = '%s/bin/vsearch-2.4.3-%s-x86_64' % (utils.get_partis_dir(), utils.get_platform_binstr())
    # --maxaccepts 0 --maxrejects 0 so vsearch considers all pairs (not just first match)
    # this is a local override for hfrac CC round 2 only, does not affect other vsearch calls in partis
    cmd = ('%s --allpairs_global %s --id %s --userout %s --userfields query+target+id '
           '--gapopen 1000I/2E --match 2 --mismatch -4 --threads 1 --quiet '
           '--maxaccepts 0 --maxrejects 0') % (
        vsearch_binary, infname, str(1.0 - round2_threshold), outfname)
    return {'cmd_str': cmd, 'outfname': outfname, 'workdir': workdir}, outfname

# ----------------------------------------------------------------------------------------
def _cc_merge_from_pairs(pairs_file, centroid_uids):
    # parse vsearch pairs output, build adjacency graph, find connected components via BFS
    # returns list of component centroid lists (each component = a merged bin's centroids)
    adj = collections.defaultdict(set)
    if os.path.exists(pairs_file):
        with open(pairs_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    a, b = parts[0], parts[1]
                    adj[a].add(b)
                    adj[b].add(a)
    for c in centroid_uids:
        if c not in adj:
            adj[c] = set()
    visited = set()
    components = []
    for c in centroid_uids:
        if c in visited:
            continue
        comp = []
        queue = collections.deque([c])
        visited.add(c)
        while queue:
            node = queue.popleft()
            comp.append(node)
            for nbr in adj[node]:
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append(nbr)
        components.append(comp)
    return components

# ----------------------------------------------------------------------------------------
def _merge_r1_subgroups_by_components(r1_results, components, max_bin_size=None):
    # r1_results: list of (centroid_uid, [member_uids]) pairs from round 1
    # components: list of lists of centroid uids (one list per component)
    # max_bin_size: if set, bin-pack round-1 sub-groups so no merged bin exceeds this size
    #   (20% tolerance). round-1 sub-groups remain indivisible so family safety is preserved.
    # returns: list of merged member uid lists (one list per output bin)
    centroid_to_comp = {}
    for cidx, comp in enumerate(components):
        for c in comp:
            centroid_to_comp[c] = cidx
    # group round-1 sub-groups by component (preserving sub-group boundaries for bin-packing)
    bins_r1 = [[] for _ in range(len(components))]
    for centroid, members in r1_results:
        cidx = centroid_to_comp.get(centroid, 0)
        bins_r1[cidx].append(members)
    bins_r1 = [b for b in bins_r1 if b]

    if max_bin_size is None or max_bin_size <= 0:
        return [[m for sg in b for m in sg] for b in bins_r1]

    # split oversize bins via First Fit Decreasing bin-packing on round-1 sub-groups
    tolerance = 1.2  # allow up to 1.2x target before splitting
    cap = max_bin_size * tolerance
    out_bins = []
    for r1_subs in bins_r1:
        bin_size = sum(len(s) for s in r1_subs)
        if bin_size <= cap or len(r1_subs) == 1:
            # within tolerance or cannot split (single indivisible sub-group over cap)
            out_bins.append([m for sg in r1_subs for m in sg])
            continue
        sorted_subs = sorted(r1_subs, key=len, reverse=True)
        packed = [[]]
        for sg in sorted_subs:
            placed = False
            for pb in packed:
                if sum(len(s) for s in pb) + len(sg) <= cap:
                    pb.append(sg)
                    placed = True
                    break
            if not placed:
                packed.append([sg])
        for pb in packed:
            if pb:
                out_bins.append([m for sg in pb for m in sg])
    return out_bins

# ----------------------------------------------------------------------------------------
def subgroup_by_naive_hamming(seqfos, hi_bound, workdir, min_group_size=100):
    # split a CDR3 group into sub-groups by clustering naive sequences with vsearch
    # uses the hi hamming bound as the identity threshold so sequences with similar
    # naive sequences end up in the same sub-group
    # note: vsearch greedy centroid clustering may split families in rare edge cases
    # with very high SHM (>15%)
    # returns list of sub-group seqfo lists
    if len(seqfos) < min_group_size:
        return [seqfos]
    naive_seqdict = {}
    for sfo in seqfos:
        if sfo.get('naive_seq', ''):
            naive_seqdict[sfo['name']] = sfo['naive_seq']
    if len(naive_seqdict) == 0:
        return [seqfos]
    partition = utils.run_vsearch('cluster', naive_seqdict, workdir, hi_bound, no_indels=True)
    uid_to_cluster = {}
    for iclust, cluster in enumerate(partition):
        for uid in cluster:
            uid_to_cluster[uid] = iclust
    n_clusters = len(partition)
    sub_groups = [[] for _ in range(n_clusters)]
    for sfo in seqfos:
        iclust = uid_to_cluster.get(sfo['name'], 0)
        sub_groups[iclust].append(sfo)
    sub_groups = [sg for sg in sub_groups if len(sg) > 0]
    return sub_groups

# ----------------------------------------------------------------------------------------
def write_group_fastas(groups, outdir, locus):
    # write per-group fasta files, returns list of group info dicts for the manifest
    group_infos = []
    for gid, (c3len, seqfos) in enumerate(sorted(groups.items())):
        group_dir = '%s/groups/cdr3-%d' % (outdir, c3len)
        fasta_path = '%s/%s.fa' % (group_dir, locus)
        utils.write_fasta(fasta_path, seqfos)
        rel_fasta_path = 'groups/cdr3-%d/%s.fa' % (c3len, locus)
        group_infos.append({
            'group_id' : gid,
            'cdr3_length' : c3len,
            'locus' : locus,
            'sequence_count' : len(seqfos),
            'fasta_path' : rel_fasta_path,
            'partition_path' : None,
        })
    return group_infos

# ----------------------------------------------------------------------------------------
def write_group_sw_caches(groups, glfo, annotation_list, outdir, locus):
    # write per-group sw-cache subsets so partition subprocesses don't have to read the huge full sw cache
    uid_to_c3len = {}
    for c3len, seqfos in groups.items():
        for sfo in seqfos:
            uid_to_c3len[sfo['name']] = c3len
    antns_by_c3len = collections.defaultdict(list)
    for line in annotation_list:
        assert len(line['unique_ids']) == 1  # sw cache always has single-sequence annotations
        uid = line['unique_ids'][0]
        if uid in uid_to_c3len:
            antns_by_c3len[uid_to_c3len[uid]].append(line)
    for c3len in sorted(groups):
        group_dir = '%s/groups/cdr3-%d' % (outdir, c3len)
        sw_cache_path = '%s/sw-cache-%s.yaml' % (group_dir, locus)
        utils.write_annotations(sw_cache_path, glfo, antns_by_c3len.get(c3len, []), utils.sw_cache_headers)

# ----------------------------------------------------------------------------------------
def _apply_hfrac_and_write(groups, hi_bound, outdir, locus, glfo, annotation_list, merge_factor=3.0, max_bin_size=15000):
    # apply hfrac sub-grouping within each CDR3 group, write per-sub-group outputs
    # round 1: vsearch greedy clustering at hi_bound produces safe sub-groups per CDR3 group
    # round 2 (if merge_factor > 0): vsearch --allpairs_global on round 1 centroids at
    #   merge_factor * hi_bound, then connected components merge via BFS to reduce
    #   sub-group count and heal round 1 greedy edge cases (family splits)
    # the greedy-round-1 code below always runs; round 2 CC is conditional on merge_factor > 0
    import shutil

    uid_to_antn = {}
    for line in annotation_list:
        assert len(line['unique_ids']) == 1
        uid_to_antn[line['unique_ids'][0]] = line

    min_group_size = 100
    vsearch_binary = '%s/bin/vsearch-2.4.3-%s-x86_64' % (utils.get_partis_dir(), utils.get_platform_binstr())

    # round 1, step 1: write naive FASTAs and build vsearch commands
    cmdfos = []
    vsearch_groups = {}
    small_groups = set()
    for c3len, seqfos in sorted(groups.items()):
        if len(seqfos) < min_group_size:
            small_groups.add(c3len)
            continue
        naive_seqdict = {sfo['name']: sfo['naive_seq'] for sfo in seqfos if sfo.get('naive_seq', '')}
        if len(naive_seqdict) == 0:
            small_groups.add(c3len)
            continue
        workdir = '%s/groups/cdr3-%d/_vsearch_work' % (outdir, c3len)
        utils.prep_dir(workdir)
        infname = workdir + '/input.fa'
        with open(infname, 'w') as f:
            for name, seq in naive_seqdict.items():
                f.write('>%s\n%s\n' % (name, seq))
        outfname = workdir + '/vsearch-clusters.txt'
        cmd = '%s --cluster_fast %s --id %s --uc %s --gapopen 1000I/2E --match 2 --mismatch -4 --threads 1' % (
            vsearch_binary, infname, str(1. - hi_bound), outfname)
        cmdfos.append({'cmd_str': cmd, 'outfname': outfname, 'workdir': workdir})
        vsearch_groups[c3len] = workdir

    # round 1, step 2: run vsearch in parallel
    if len(cmdfos) > 0:
        n_procs = min(8, len(cmdfos))
        print('        running %d vsearch hfrac jobs (%d concurrent)' % (len(cmdfos), n_procs))
        utils.run_cmds(cmdfos, n_max_procs=n_procs)

    # round 1, step 3: parse UC files, extract centroids and members per CDR3 group
    r1_by_c3len = {}  # c3len -> list of (centroid_uid, [member_uids])
    for c3len in sorted(groups):
        if c3len in small_groups:
            continue
        cluster_file = '%s/vsearch-clusters.txt' % vsearch_groups[c3len]
        r1_by_c3len[c3len] = _read_vsearch_uc_with_centroids(cluster_file)

    # round 2 (optional): CC merge via vsearch --allpairs_global on centroids
    comps_by_c3len = {}  # c3len -> list of components (each a list of centroid uids)
    if merge_factor > 0:
        r2_cmdfos = []
        r2_workdirs = {}
        r2_pairs_files = {}
        round2_threshold = min(merge_factor * hi_bound, 0.49)
        for c3len, r1_results in r1_by_c3len.items():
            if len(r1_results) < 2:
                comps_by_c3len[c3len] = [[c] for c, _ in r1_results]
                continue
            naive_by_uid = {sfo['name']: sfo['naive_seq'] for sfo in groups[c3len] if sfo.get('naive_seq', '')}
            centroid_naives = {c: naive_by_uid[c] for c, _ in r1_results if c in naive_by_uid}
            r2_workdir = '%s/groups/cdr3-%d/_round2_cc' % (outdir, c3len)
            cmdfo, pairs_fname = _build_round2_cc_cmd(centroid_naives, round2_threshold, r2_workdir)
            if cmdfo is None:
                comps_by_c3len[c3len] = [[c] for c, _ in r1_results]
                continue
            r2_cmdfos.append(cmdfo)
            r2_workdirs[c3len] = r2_workdir
            r2_pairs_files[c3len] = pairs_fname
        if len(r2_cmdfos) > 0:
            n_procs2 = min(8, len(r2_cmdfos))
            print('        running %d round 2 CC jobs (%d concurrent, %.2fx hi_bound = %.4f threshold)' % (
                len(r2_cmdfos), n_procs2, merge_factor, round2_threshold))
            utils.run_cmds(r2_cmdfos, n_max_procs=n_procs2)
        for c3len, pairs_fname in r2_pairs_files.items():
            centroid_uids = [c for c, _ in r1_by_c3len[c3len]]
            comps_by_c3len[c3len] = _cc_merge_from_pairs(pairs_fname, centroid_uids)

    # step 4: build sub_groups_list per CDR3 group (merged by CC if round 2 enabled, else raw round 1)
    all_group_infos = []
    flattened_groups = collections.OrderedDict()
    for c3len, seqfos in sorted(groups.items()):
        if c3len in small_groups:
            sub_groups_list = [seqfos]
        else:
            r1_results = r1_by_c3len[c3len]
            if merge_factor > 0 and c3len in comps_by_c3len:
                merged_bins_uids = _merge_r1_subgroups_by_components(r1_results, comps_by_c3len[c3len], max_bin_size=max_bin_size)
            else:
                merged_bins_uids = [members for _, members in r1_results]
            uid_to_sfo = {sfo['name']: sfo for sfo in seqfos}
            sub_groups_list = []
            for bin_uids in merged_bins_uids:
                sg = [uid_to_sfo[u] for u in bin_uids if u in uid_to_sfo]
                if sg:
                    sub_groups_list.append(sg)
            shutil.rmtree(vsearch_groups[c3len])
            if merge_factor > 0 and c3len in r2_workdirs:
                shutil.rmtree(r2_workdirs[c3len], ignore_errors=True)
        for isub, sub_seqfos in enumerate(sub_groups_list):
            sub_dir = '%s/groups/cdr3-%d/sub-groups/sub-%03d' % (outdir, c3len, isub)
            fasta_path = '%s/%s.fa' % (sub_dir, locus)
            utils.write_fasta(fasta_path, sub_seqfos)
            sub_uids = set(sfo['name'] for sfo in sub_seqfos)
            sub_antns = [uid_to_antn[uid] for uid in sub_uids if uid in uid_to_antn]
            sw_cache_path = '%s/sw-cache-%s.yaml' % (sub_dir, locus)
            utils.write_annotations(sw_cache_path, glfo, sub_antns, utils.sw_cache_headers)
            unique_naive = len(set(sfo.get('naive_seq', '') for sfo in sub_seqfos if sfo.get('naive_seq', '')))
            rel_fasta = 'groups/cdr3-%d/sub-groups/sub-%03d/%s.fa' % (c3len, isub, locus)
            all_group_infos.append({
                'group_id' : len(all_group_infos),
                'cdr3_length' : c3len,
                'sub_group_id' : isub,
                'locus' : locus,
                'sequence_count' : len(sub_seqfos),
                'unique_naive_count' : unique_naive,
                'fasta_path' : rel_fasta,
                'partition_path' : None,
            })
            flattened_groups[(c3len, isub)] = sub_seqfos
        if len(sub_groups_list) > 1:
            print('        cdr3-%d: %d seqs -> %d sub-groups (sizes: %s)' % (c3len, len(seqfos), len(sub_groups_list), ' '.join(str(len(sg)) for sg in sub_groups_list)))
    return flattened_groups, all_group_infos

# ----------------------------------------------------------------------------------------
def _apply_hfrac_two_pass(groups, hi_bound, outdir, locus, glfo, merge_factor=3.0, max_bin_size=15000):
    # memory-efficient hfrac for multi-cache path: reads one CDR3 group SW cache at a time
    # round 1 pass 1: write naive FASTAs, build vsearch commands
    # round 1 dispatch: run all vsearch in parallel
    # round 2 (if merge_factor > 0): CC merge on round 1 centroids via --allpairs_global
    # pass 2: re-read SW caches, parse results, write sub-group outputs
    import shutil

    min_group_size = 100
    vsearch_binary = '%s/bin/vsearch-2.4.3-%s-x86_64' % (utils.get_partis_dir(), utils.get_platform_binstr())

    # round 1 pass 1: write naive FASTAs for groups needing splitting
    cmdfos = []
    vsearch_groups = {}
    small_groups = set()
    for c3len, seqfos in sorted(groups.items()):
        if len(seqfos) < min_group_size:
            small_groups.add(c3len)
            continue
        naive_seqdict = {sfo['name']: sfo['naive_seq'] for sfo in seqfos if sfo.get('naive_seq', '')}
        if len(naive_seqdict) == 0:
            small_groups.add(c3len)
            continue
        workdir = '%s/groups/cdr3-%d/_vsearch_work' % (outdir, c3len)
        utils.prep_dir(workdir)
        infname = workdir + '/input.fa'
        with open(infname, 'w') as f:
            for name, seq in naive_seqdict.items():
                f.write('>%s\n%s\n' % (name, seq))
        outfname = workdir + '/vsearch-clusters.txt'
        cmd = '%s --cluster_fast %s --id %s --uc %s --gapopen 1000I/2E --match 2 --mismatch -4 --threads 1' % (
            vsearch_binary, infname, str(1. - hi_bound), outfname)
        cmdfos.append({'cmd_str': cmd, 'outfname': outfname, 'workdir': workdir})
        vsearch_groups[c3len] = workdir

    # round 1 dispatch: run all vsearch in parallel
    if len(cmdfos) > 0:
        n_procs = min(8, len(cmdfos))
        print('        running %d vsearch hfrac jobs (%d concurrent)' % (len(cmdfos), n_procs))
        utils.run_cmds(cmdfos, n_max_procs=n_procs)

    # round 1 parse: extract centroids and members per CDR3 group
    r1_by_c3len = {}
    for c3len in sorted(groups):
        if c3len in small_groups:
            continue
        cluster_file = '%s/vsearch-clusters.txt' % vsearch_groups[c3len]
        r1_by_c3len[c3len] = _read_vsearch_uc_with_centroids(cluster_file)

    # round 2 (optional): CC merge on centroids via --allpairs_global
    comps_by_c3len = {}
    r2_workdirs = {}
    if merge_factor > 0:
        r2_cmdfos = []
        r2_pairs_files = {}
        round2_threshold = min(merge_factor * hi_bound, 0.49)
        for c3len, r1_results in r1_by_c3len.items():
            if len(r1_results) < 2:
                comps_by_c3len[c3len] = [[c] for c, _ in r1_results]
                continue
            naive_by_uid = {sfo['name']: sfo['naive_seq'] for sfo in groups[c3len] if sfo.get('naive_seq', '')}
            centroid_naives = {c: naive_by_uid[c] for c, _ in r1_results if c in naive_by_uid}
            r2_workdir = '%s/groups/cdr3-%d/_round2_cc' % (outdir, c3len)
            cmdfo, pairs_fname = _build_round2_cc_cmd(centroid_naives, round2_threshold, r2_workdir)
            if cmdfo is None:
                comps_by_c3len[c3len] = [[c] for c, _ in r1_results]
                continue
            r2_cmdfos.append(cmdfo)
            r2_workdirs[c3len] = r2_workdir
            r2_pairs_files[c3len] = pairs_fname
        if len(r2_cmdfos) > 0:
            n_procs2 = min(8, len(r2_cmdfos))
            print('        running %d round 2 CC jobs (%d concurrent, %.2fx hi_bound = %.4f threshold)' % (
                len(r2_cmdfos), n_procs2, merge_factor, round2_threshold))
            utils.run_cmds(r2_cmdfos, n_max_procs=n_procs2)
        for c3len, pairs_fname in r2_pairs_files.items():
            centroid_uids = [c for c, _ in r1_by_c3len[c3len]]
            comps_by_c3len[c3len] = _cc_merge_from_pairs(pairs_fname, centroid_uids)

    # pass 2: re-read SW caches, build sub_groups_list (from CC components if round 2, else round 1), write outputs
    all_group_infos = []
    flattened_groups = collections.OrderedDict()
    for c3len, seqfos in sorted(groups.items()):
        if c3len in small_groups:
            sub_groups_list = [seqfos]
        else:
            r1_results = r1_by_c3len[c3len]
            if merge_factor > 0 and c3len in comps_by_c3len:
                merged_bins_uids = _merge_r1_subgroups_by_components(r1_results, comps_by_c3len[c3len], max_bin_size=max_bin_size)
            else:
                merged_bins_uids = [members for _, members in r1_results]
            uid_to_sfo = {sfo['name']: sfo for sfo in seqfos}
            sub_groups_list = []
            for bin_uids in merged_bins_uids:
                sg = [uid_to_sfo[u] for u in bin_uids if u in uid_to_sfo]
                if sg:
                    sub_groups_list.append(sg)
            shutil.rmtree(vsearch_groups[c3len])
            if merge_factor > 0 and c3len in r2_workdirs:
                shutil.rmtree(r2_workdirs[c3len], ignore_errors=True)

        # re-read SW cache for this CDR3 group to get annotations AND the merged glfo
        # (multi-chunk path: merge_yamls reconciled glfos across chunks; using the outer glfo
        # from the first chunk would lose novel alleles from later chunks)
        swc_path = '%s/groups/cdr3-%d/sw-cache-%s.yaml' % (outdir, c3len, locus)
        uid_to_antn = {}
        group_glfo = glfo
        if os.path.exists(swc_path):
            group_glfo, antn_list, _ = utils.read_yaml_output(swc_path, dont_add_implicit_info=True)
            for line in antn_list:
                if len(line['unique_ids']) == 1:
                    uid_to_antn[line['unique_ids'][0]] = line
            del antn_list

        for isub, sub_seqfos in enumerate(sub_groups_list):
            sub_dir = '%s/groups/cdr3-%d/sub-groups/sub-%03d' % (outdir, c3len, isub)
            fasta_path = '%s/%s.fa' % (sub_dir, locus)
            utils.write_fasta(fasta_path, sub_seqfos)
            sub_uids = set(sfo['name'] for sfo in sub_seqfos)
            sub_antns = [uid_to_antn[uid] for uid in sub_uids if uid in uid_to_antn]
            sw_cache_path = '%s/sw-cache-%s.yaml' % (sub_dir, locus)
            utils.write_annotations(sw_cache_path, group_glfo, sub_antns, utils.sw_cache_headers)
            unique_naive = len(set(sfo.get('naive_seq', '') for sfo in sub_seqfos if sfo.get('naive_seq', '')))
            rel_fasta = 'groups/cdr3-%d/sub-groups/sub-%03d/%s.fa' % (c3len, isub, locus)
            all_group_infos.append({
                'group_id' : len(all_group_infos),
                'cdr3_length' : c3len,
                'sub_group_id' : isub,
                'locus' : locus,
                'sequence_count' : len(sub_seqfos),
                'unique_naive_count' : unique_naive,
                'fasta_path' : rel_fasta,
                'partition_path' : None,
            })
            flattened_groups[(c3len, isub)] = sub_seqfos
        if len(sub_groups_list) > 1:
            print('        cdr3-%d: %d seqs -> %d sub-groups (sizes: %s)' % (c3len, len(seqfos), len(sub_groups_list), ' '.join(str(len(sg)) for sg in sub_groups_list)))
    return flattened_groups, all_group_infos

# ----------------------------------------------------------------------------------------
def write_manifest(group_infos, outdir, locus, total_input, n_failed, parameter_dir=None, hfrac=False):
    # write manifest yaml to outdir
    manifest = {
        'grouping-info' : {
            'method' : 'cdr3-length+hfrac' if hfrac else 'cdr3-length',
            'locus' : locus,
            'total_input_sequences' : total_input,
            'total_grouped_sequences' : total_input - n_failed,
            'failed_sequences' : n_failed,
            'parameter_dir' : parameter_dir,
        },
        'groups' : group_infos,
        'assembly' : {
            'status' : 'pending',
            'merged_output_path' : None,
            'validation' : {
                'uids_unique' : None,
                'sequence_count_preserved' : None,
                # TODO add gene_lists_consistent check (verify germline gene lists are compatible across groups)
            },
        },
    }
    manifest_path = '%s/%s' % (outdir, MANIFEST_FNAME)
    utils.mkdir(manifest_path, isfile=True)
    with open(manifest_path, 'w') as mfile:
        yaml.dump(manifest, mfile, width=400, default_flow_style=False)
    print('      wrote manifest to %s' % manifest_path)
    return manifest

# ----------------------------------------------------------------------------------------
def read_manifest(manifest_path):
    # read and validate manifest yaml
    if not os.path.exists(manifest_path):
        raise Exception('manifest file does not exist: %s' % manifest_path)
    with open(manifest_path) as mfile:
        manifest = yaml.safe_load(mfile)
    for required_key in ['grouping-info', 'groups', 'assembly']:
        if required_key not in manifest:
            raise Exception('missing required key \'%s\' in manifest %s' % (required_key, manifest_path))
    for ginfo in manifest['groups']:
        for required_key in ['group_id', 'cdr3_length', 'locus', 'sequence_count', 'fasta_path']:
            if required_key not in ginfo:
                raise Exception('missing required key \'%s\' in group entry in manifest %s' % (required_key, manifest_path))
    return manifest

# ----------------------------------------------------------------------------------------
def validate_sequence_count(manifest):
    # verify that sum of group sequence counts equals total_grouped_sequences
    total_grouped = manifest['grouping-info']['total_grouped_sequences']
    group_sum = sum(g['sequence_count'] for g in manifest['groups'])
    if group_sum != total_grouped:
        raise Exception('sequence count mismatch: sum of group counts %d does not equal total_grouped_sequences %d' % (group_sum, total_grouped))
    total_input = manifest['grouping-info']['total_input_sequences']
    n_failed = manifest['grouping-info']['failed_sequences']
    if total_grouped + n_failed != total_input:
        raise Exception('sequence count mismatch: total_grouped %d + failed %d does not equal total_input %d' % (total_grouped, n_failed, total_input))
    print('      sequence count validated: %d grouped + %d failed = %d total' % (total_grouped, n_failed, total_input))

# ----------------------------------------------------------------------------------------
def get_partition_paths(manifest, manifest_dir):
    # collect and verify partition file paths for a single locus
    # if partition_path is set in manifest, use it directly
    # if partition_path is None, try to discover the partition file in the group dir
    # (supports standalone partition jobs that do not update the manifest)
    paths = []
    skipped_groups = []
    missing_files = []
    for ginfo in manifest['groups']:
        ppath = ginfo.get('partition_path')
        if ppath is None:
            # check for partition file in the same directory as the fasta
            fasta_dir = os.path.dirname(ginfo['fasta_path'])
            default_ppath = '%s/partition-%s.yaml' % (fasta_dir, ginfo['locus'])
            if os.path.exists('%s/%s' % (manifest_dir, default_ppath)):
                ppath = default_ppath
            else:
                skipped_groups.append(ginfo['group_id'])
                continue
        full_ppath = '%s/%s' % (manifest_dir, ppath)
        if not os.path.exists(full_ppath):
            missing_files.append(ginfo['group_id'])
            continue
        if os.path.getsize(full_ppath) == 0:
            raise Exception('partition file is empty for group %d: %s' % (ginfo['group_id'], full_ppath))
        paths.append(full_ppath)
    if len(skipped_groups) > 0:
        print('      skipping %d groups with no partition output (e.g. too small): %s' % (len(skipped_groups), skipped_groups))
    if len(missing_files) > 0:
        raise Exception('partition files missing for %d groups (partition_path set but file not found): %s' % (len(missing_files), missing_files))
    return paths

# ----------------------------------------------------------------------------------------
def validate_assembly(manifest, manifest_dir):
    # validate uid uniqueness and sequence counts by reading partitioned groups
    all_uids = set()
    total_seqs = 0
    skipped = [g for g in manifest['groups'] if g.get('partition_path') is None]
    skipped_seqs = sum(g['sequence_count'] for g in skipped)
    for ppath in get_partition_paths(manifest, manifest_dir):
        _, annotation_list, _ = utils.read_yaml_output(ppath, dont_add_implicit_info=True)
        for line in annotation_list:
            for uid in line['unique_ids']:
                if uid in all_uids:
                    raise Exception('duplicate uid %s found across groups' % uid)
                all_uids.add(uid)
        total_seqs += sum(len(line['unique_ids']) for line in annotation_list)
    expected = manifest['grouping-info']['total_grouped_sequences'] - skipped_seqs
    if total_seqs != expected:
        raise Exception('sequence count mismatch after assembly: found %d in partition files, expected %d (total %d minus %d skipped)' % (total_seqs, expected, manifest['grouping-info']['total_grouped_sequences'], skipped_seqs))
    print('      assembly validation passed: %d sequences from %d groups (%d sequences in %d groups skipped)' % (total_seqs, len(manifest['groups']) - len(skipped), skipped_seqs, len(skipped)))

# ----------------------------------------------------------------------------------------
def resolve_sw_cache_paths(sw_cache_paths):
    # resolve <sw_cache_paths> to a list: accepts a single path string, a list of paths, or a directory
    # for directory input, checks two expected patterns:
    #   paired/chunked layout: {dir}/*/parameters/*/sw-cache*.yaml
    #   flat unpaired layout:  {dir}/*/parameters/sw-cache*.yaml
    if isinstance(sw_cache_paths, str):
        if os.path.isdir(sw_cache_paths):
            paths = sorted(glob.glob('%s/*/parameters/*/sw-cache*.yaml' % sw_cache_paths))
            if len(paths) == 0:
                paths = sorted(glob.glob('%s/*/parameters/sw-cache*.yaml' % sw_cache_paths))
            if len(paths) == 0:
                paths = sorted(glob.glob('%s/*/sw-cache*.yaml' % sw_cache_paths))
            if len(paths) == 0:
                paths = sorted(glob.glob('%s/sw-cache*.yaml' % sw_cache_paths))
            if len(paths) == 0:
                raise Exception('no sw-cache*.yaml files found in %s (checked */parameters/*/sw-cache*.yaml, */parameters/sw-cache*.yaml, */sw-cache*.yaml, sw-cache*.yaml)' % sw_cache_paths)
            return paths
        else:
            return [sw_cache_paths]
    return list(sw_cache_paths)

# ----------------------------------------------------------------------------------------
def create_cdr3_groups(locus, sw_cache_paths, outdir, parameter_dir, hfrac=False, hfrac_merge_factor=3.0, hfrac_max_bin_size=15000):
    # read sw cache(s) for a single locus, group sequences by CDR3 length,
    # optionally sub-group by naive hamming fraction (--hfrac),
    # write per-group (or per-sub-group) fastas and sw-cache subsets, write manifest.
    # <sw_cache_paths>: single path string, list of paths, or directory (for chunked cache-parameters at scale).
    # For multiple caches, processes one chunk at a time to limit peak memory:
    #   - per-group FASTAs are written after all chunks are grouped (seqfos are lightweight)
    #   - per-group sw-cache fragments are written per chunk, then merged and cleaned up
    sw_cache_paths = resolve_sw_cache_paths(sw_cache_paths)
    multi_cache = len(sw_cache_paths) > 1

    # compute hi hamming bound for hfrac sub-grouping
    hi_bound = None
    if hfrac:
        # get_mean_mfreq expects the dir containing all-mean-mute-freqs.csv
        # try nested layout ({pdir}/{locus}/sw/) then flat layout ({pdir}/sw/)
        mfreq_dir = None
        for candidate in ['%s/%s/sw' % (parameter_dir, locus), '%s/%s/hmm' % (parameter_dir, locus),
                          '%s/sw' % parameter_dir, '%s/hmm' % parameter_dir]:
            if os.path.exists('%s/all-mean-mute-freqs.csv' % candidate):
                mfreq_dir = candidate
                break
        if mfreq_dir is None:
            raise Exception('could not find all-mean-mute-freqs.csv in %s (checked {locus}/sw, {locus}/hmm, sw, hmm)' % parameter_dir)
        _, hi_bound = utils.get_naive_hamming_bounds('likelihood', mfreq_dir)
        print('      hfrac sub-grouping enabled (hi bound: %.4f)' % hi_bound)

    if not multi_cache:
        # single sw cache: read once, process everything in memory (existing behavior)
        print('      reading sw cache for %s from %s' % (locus, sw_cache_paths[0]))
        glfo, annotation_list, _ = utils.read_yaml_output(sw_cache_paths[0], dont_add_implicit_info=True)
        groups, n_failed = group_sequences_by_cdr3_length(annotation_list)
        n_seqs = sum(len(seqfos) for seqfos in groups.values()) + n_failed
        if hfrac:
            groups, group_infos = _apply_hfrac_and_write(groups, hi_bound, outdir, locus, glfo, annotation_list, merge_factor=hfrac_merge_factor, max_bin_size=hfrac_max_bin_size)
        else:
            group_infos = write_group_fastas(groups, outdir, locus)
            write_group_sw_caches(groups, glfo, annotation_list, outdir, locus)
    else:
        # multiple sw caches: process one chunk at a time
        print('      processing %d sw cache files for %s' % (len(sw_cache_paths), locus))
        glfo = None
        all_groups = collections.OrderedDict()  # cdr3_length -> [seqfos] (lightweight: uid + seq only)
        n_failed = 0
        n_seqs = 0
        chunk_fragments = collections.defaultdict(list)  # cdr3_length -> list of fragment file paths

        for ichunk, swpath in enumerate(sw_cache_paths):
            print('      chunk %d/%d: %s' % (ichunk + 1, len(sw_cache_paths), swpath))
            tglfo, tantn_list, _ = utils.read_yaml_output(swpath, dont_add_implicit_info=True)
            if glfo is None:
                glfo = tglfo
            chunk_groups, chunk_failed = group_sequences_by_cdr3_length(tantn_list)
            n_failed += chunk_failed
            n_seqs += sum(len(seqfos) for seqfos in chunk_groups.values()) + chunk_failed

            # accumulate seqfos for FASTA writing, write sw-cache fragment per group
            for c3len, seqfos in chunk_groups.items():
                all_groups.setdefault(c3len, []).extend(seqfos)
                group_dir = '%s/groups/cdr3-%d' % (outdir, c3len)
                frag_path = '%s/sw-cache-chunk%03d.yaml' % (group_dir, ichunk)
                uid_set = set(sfo['name'] for sfo in seqfos)
                chunk_antns = [line for line in tantn_list if len(line['unique_ids']) == 1 and line['unique_ids'][0] in uid_set]
                utils.mkdir(frag_path, isfile=True)
                utils.write_annotations(frag_path, tglfo, chunk_antns, utils.sw_cache_headers)
                chunk_fragments[c3len].append(frag_path)

            del tantn_list  # free chunk annotations

        # write per-group FASTAs
        groups = collections.OrderedDict(sorted(all_groups.items()))
        group_infos = write_group_fastas(groups, outdir, locus)

        # merge per-chunk sw-cache fragments into final per-group files, then clean up
        for c3len in sorted(groups):
            final_swc = '%s/groups/cdr3-%d/sw-cache-%s.yaml' % (outdir, c3len, locus)
            frags = chunk_fragments.get(c3len, [])
            if len(frags) == 1:
                os.rename(frags[0], final_swc)
            elif len(frags) > 1:
                utils.merge_yamls(final_swc, frags, utils.sw_cache_headers, dont_write_git_info=True)
                for frag in frags:
                    os.remove(frag)
        # NOTE if multiple chunks inferred different novel alleles, glfo from the first chunk is used.
        # For proper germline reconciliation across chunks, merge parameter dirs before grouping.

        # apply hfrac after merging: two-pass approach for memory efficiency
        # pass 1: read each CDR3 sw cache, write naive FASTAs (lightweight)
        # then dispatch all vsearch jobs in parallel
        # pass 2: read each CDR3 sw cache again, parse vsearch results, write sub-group outputs
        if hfrac:
            _, group_infos = _apply_hfrac_two_pass(groups, hi_bound, outdir, locus, glfo, merge_factor=hfrac_merge_factor, max_bin_size=hfrac_max_bin_size)

    n_cdr3_groups = len(set(g['cdr3_length'] for g in group_infos)) if len(group_infos) > 0 else 0
    print('      %s: %d sequences in %d cdr3 length groups (%d failed)' % (locus, n_seqs - n_failed, n_cdr3_groups, n_failed))
    if hfrac and len(group_infos) > 0:
        print('      hfrac: %d sub-groups total' % len(group_infos))
    if len(group_infos) > 0 and not hfrac:
        cw = max(len(str(g['cdr3_length'])) for g in group_infos)
        print('        cdr3 lengths : %s' % '  '.join('%*d' % (cw, g['cdr3_length']) for g in group_infos))
        print('        N seqs       : %s' % '  '.join('%*d' % (cw, g['sequence_count']) for g in group_infos))
    manifest = write_manifest(group_infos, outdir, locus, n_seqs, n_failed, parameter_dir=parameter_dir, hfrac=hfrac)
    validate_sequence_count(manifest)
    return manifest

# ----------------------------------------------------------------------------------------
def assemble_groups(locus, disjoint_dir, outfname):
    # validate and concatenate per-group partition results for a single locus
    manifest_path = '%s/%s' % (disjoint_dir, MANIFEST_FNAME)
    print('    assembling groups from %s' % manifest_path)
    manifest = read_manifest(manifest_path)
    disjoint_dir = os.path.abspath(disjoint_dir)

    validate_assembly(manifest, disjoint_dir)
    manifest['assembly']['validation']['uids_unique'] = True
    manifest['assembly']['validation']['sequence_count_preserved'] = True

    utils.mkdir(outfname, isfile=True)
    yaml_list = get_partition_paths(manifest, disjoint_dir)
    headers = list(utils.annotation_headers)
    print('      merging %d partition files for %s:' % (len(yaml_list), locus))
    utils.merge_yamls(outfname, yaml_list, headers, best_partition_only=True, dont_write_git_info=True, debug=True)

    manifest['assembly']['status'] = 'merged'
    manifest['assembly']['merged_output_path'] = outfname
    with open(manifest_path, 'w') as mfile:
        yaml.dump(manifest, mfile, width=400, default_flow_style=False)
    print('      updated manifest')
