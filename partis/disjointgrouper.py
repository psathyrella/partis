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
    partition = utils.run_vsearch('cluster', naive_seqdict, workdir, hi_bound, no_indels=True, maxaccepts=0, maxrejects=0)
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
def _apply_hfrac_and_write(groups, hi_bound, outdir, locus, glfo, annotation_list):
    # apply hfrac sub-grouping within each CDR3 group, write per-sub-group outputs
    # step 1: write naive FASTAs for groups needing splitting
    # step 2: run vsearch in parallel via utils.run_cmds()
    # step 3: parse results, write per-sub-group FASTAs and SW caches
    import shutil

    uid_to_antn = {}
    for line in annotation_list:
        assert len(line['unique_ids']) == 1
        uid_to_antn[line['unique_ids'][0]] = line

    min_group_size = 100
    vsearch_binary = '%s/bin/vsearch-2.4.3-%s-x86_64' % (utils.get_partis_dir(), utils.get_platform_binstr())

    # step 1: write naive FASTAs and build vsearch commands
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
        cmd = '%s --cluster_fast %s --id %s --uc %s --maxaccepts 0 --maxrejects 0 --gapopen 1000I/2E --match 2 --mismatch -4 --threads 1' % (
            vsearch_binary, infname, str(1. - hi_bound), outfname)
        cmdfos.append({'cmd_str': cmd, 'outfname': outfname, 'workdir': workdir})
        vsearch_groups[c3len] = workdir

    # step 2: run vsearch in parallel
    if len(cmdfos) > 0:
        n_procs = min(8, len(cmdfos))
        print('        running %d vsearch hfrac jobs (%d concurrent)' % (len(cmdfos), n_procs))
        utils.run_cmds(cmdfos, n_max_procs=n_procs)

    # step 3: parse results and write outputs
    all_group_infos = []
    flattened_groups = collections.OrderedDict()
    for c3len, seqfos in sorted(groups.items()):
        if c3len in small_groups:
            sub_groups_list = [seqfos]
        else:
            cluster_file = '%s/vsearch-clusters.txt' % vsearch_groups[c3len]
            partition = utils.read_vsearch_cluster_file(cluster_file)
            uid_to_cluster = {}
            for iclust, cluster in enumerate(partition):
                for uid in cluster:
                    uid_to_cluster[uid] = iclust
            sub_groups_list = [[] for _ in range(len(partition))]
            for sfo in seqfos:
                iclust = uid_to_cluster.get(sfo['name'], 0)
                sub_groups_list[iclust].append(sfo)
            sub_groups_list = [sg for sg in sub_groups_list if len(sg) > 0]
            shutil.rmtree(vsearch_groups[c3len])
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
def create_cdr3_groups(locus, sw_cache_paths, outdir, parameter_dir, hfrac=False):
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
            groups, group_infos = _apply_hfrac_and_write(groups, hi_bound, outdir, locus, glfo, annotation_list)
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

        # apply hfrac after merging: read all per-CDR3-group sw caches and pass to
        # _apply_hfrac_and_write at once so vsearch jobs can run in parallel across groups
        if hfrac:
            all_antn_list = []
            for c3len in sorted(groups):
                swc_path = '%s/groups/cdr3-%d/sw-cache-%s.yaml' % (outdir, c3len, locus)
                _, antn_list, _ = utils.read_yaml_output(swc_path, dont_add_implicit_info=True)
                all_antn_list.extend(antn_list)
            _, group_infos = _apply_hfrac_and_write(groups, hi_bound, outdir, locus, glfo, all_antn_list)
            del all_antn_list

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
