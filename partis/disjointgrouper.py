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
    # returns {cdr3_length : [{'name': uid, 'seq': seq}, ...]}, n_failed
    seqfo_map = {}  # uid -> {'name': uid, 'seq': seq, 'cdr3_length': int}
    n_failed = 0
    for line in annotation_list:
        if 'cdr3_length' not in line or line['cdr3_length'] is None:
            n_failed += len(line['unique_ids'])
            continue
        for uid, seq in zip(line['unique_ids'], line['input_seqs']):
            seqfo_map[uid] = {'name' : uid, 'seq' : seq, 'cdr3_length' : line['cdr3_length']}
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
def write_manifest(group_infos, outdir, locus, total_input, n_failed, parameter_dir=None):
    # write manifest yaml to outdir
    manifest = {
        'grouping-info' : {
            'method' : 'cdr3-length',
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
    paths = []
    missing_partitions = []
    for ginfo in manifest['groups']:
        ppath = ginfo.get('partition_path')
        if ppath is None:
            missing_partitions.append(ginfo['group_id'])
            continue
        full_ppath = '%s/%s' % (manifest_dir, ppath)
        if not os.path.exists(full_ppath):
            missing_partitions.append(ginfo['group_id'])
            continue
        if os.path.getsize(full_ppath) == 0:
            raise Exception('partition file is empty for group %d: %s' % (ginfo['group_id'], full_ppath))
        paths.append(full_ppath)
    if len(missing_partitions) > 0:
        raise Exception('partition files missing for %d groups: %s' % (len(missing_partitions), missing_partitions))
    return paths

# ----------------------------------------------------------------------------------------
def validate_assembly(manifest, manifest_dir):
    # validate uid uniqueness and sequence counts by reading each group one at a time
    all_uids = set()
    total_seqs = 0
    for ppath in get_partition_paths(manifest, manifest_dir):
        _, annotation_list, _ = utils.read_yaml_output(ppath, dont_add_implicit_info=True)
        for line in annotation_list:
            for uid in line['unique_ids']:
                if uid in all_uids:
                    raise Exception('duplicate uid %s found across groups' % uid)
                all_uids.add(uid)
        total_seqs += sum(len(line['unique_ids']) for line in annotation_list)
    expected = manifest['grouping-info']['total_grouped_sequences']
    if total_seqs != expected:
        raise Exception('sequence count mismatch after assembly: found %d uids in partition files, expected %d' % (total_seqs, expected))
    print('      assembly validation passed: %d unique sequences across %d groups' % (total_seqs, len(manifest['groups'])))

# ----------------------------------------------------------------------------------------
def resolve_sw_cache_paths(sw_cache_paths):
    # resolve <sw_cache_paths> to a list: accepts a single path string, a list of paths, or a directory
    # (recursively globs for sw-cache*.yaml)
    if isinstance(sw_cache_paths, str):
        if os.path.isdir(sw_cache_paths):
            paths = sorted(glob.glob('%s/**/sw-cache*.yaml' % sw_cache_paths, recursive=True))
            if len(paths) == 0:
                raise Exception('no sw-cache*.yaml files found in %s or its subdirectories' % sw_cache_paths)
            return paths
        else:
            return [sw_cache_paths]
    return list(sw_cache_paths)

# ----------------------------------------------------------------------------------------
def create_cdr3_groups(locus, sw_cache_paths, outdir, parameter_dir):
    # read sw cache(s) for a single locus, group sequences by CDR3 length,
    # write per-group fastas and sw-cache subsets, write manifest.
    # <sw_cache_paths>: single path string, list of paths, or directory (for chunked cache-parameters at scale).
    # For multiple caches, processes one chunk at a time to limit peak memory:
    #   - per-group FASTAs are written after all chunks are grouped (seqfos are lightweight)
    #   - per-group sw-cache fragments are written per chunk, then merged and cleaned up
    sw_cache_paths = resolve_sw_cache_paths(sw_cache_paths)
    multi_cache = len(sw_cache_paths) > 1

    if not multi_cache:
        # single sw cache: read once, process everything in memory (existing behavior)
        print('      reading sw cache for %s from %s' % (locus, sw_cache_paths[0]))
        glfo, annotation_list, _ = utils.read_yaml_output(sw_cache_paths[0], dont_add_implicit_info=True)
        groups, n_failed = group_sequences_by_cdr3_length(annotation_list)
        n_seqs = sum(len(seqfos) for seqfos in groups.values()) + n_failed
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

    print('      %s: %d sequences in %d cdr3 length groups (%d failed)' % (locus, n_seqs - n_failed, len(groups), n_failed))
    if len(group_infos) > 0:
        cw = max(len(str(g['cdr3_length'])) for g in group_infos)
        print('        cdr3 lengths : %s' % '  '.join('%*d' % (cw, g['cdr3_length']) for g in group_infos))
        print('        N seqs       : %s' % '  '.join('%*d' % (cw, g['sequence_count']) for g in group_infos))
    manifest = write_manifest(group_infos, outdir, locus, n_seqs, n_failed, parameter_dir=parameter_dir)
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
