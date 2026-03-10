from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import sys
import json
import yaml
import collections

from . import utils
from . import glutils
from . import paircluster

# ----------------------------------------------------------------------------------------
def group_sequences_by_cdr3_length(annotation_list):
    # group uids and their input sequences by cdr3 length from sw annotation list
    # each annotation has unique_ids (list) and input_seqs (list, parallel), plus cdr3_length (int)
    groups = collections.OrderedDict()  # {cdr3_length : [{'name': uid, 'seq': seq}, ...]}
    n_failed = 0
    for line in annotation_list:
        if 'cdr3_length' not in line or line['cdr3_length'] is None:
            n_failed += len(line['unique_ids'])
            continue
        c3len = line['cdr3_length']
        if c3len not in groups:
            groups[c3len] = []
        for uid, seq in zip(line['unique_ids'], line['input_seqs']):
            groups[c3len].append({'name' : uid, 'seq' : seq})
    if n_failed > 0:
        print('  %s %d sequences had no cdr3_length and were excluded from grouping' % (utils.color('yellow', 'warning'), n_failed))
    return groups, n_failed

# ----------------------------------------------------------------------------------------
def write_group_fastas(groups, outdir, locus):
    # write per-group fasta files, one per cdr3 length, streaming writes
    # returns list of group info dicts for the manifest
    group_infos = []
    for gid, (c3len, seqfos) in enumerate(sorted(groups.items())):
        group_dir = '%s/groups/cdr3-%d' % (outdir, c3len)
        fasta_path = '%s/%s.fa' % (group_dir, locus)
        utils.mkdir(fasta_path, isfile=True)
        with open(fasta_path, 'w') as ffile:
            for sfo in seqfos:
                ffile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))
        rel_fasta_path = 'groups/cdr3-%d/%s.fa' % (c3len, locus)
        group_infos.append({
            'group_id' : gid,
            'cdr3_length' : c3len,
            'locus' : locus,
            'sequence_count' : len(seqfos),
            'fasta_path' : rel_fasta_path,
            'partition_path' : None,
        })
        print('    group %d: cdr3 length %d, %d sequences -> %s' % (gid, c3len, len(seqfos), rel_fasta_path))
    return group_infos

# ----------------------------------------------------------------------------------------
def write_manifest(group_infos, outdir, loci, total_input, n_failed, parameter_dir=None):
    # write manifest yaml to outdir
    manifest = {
        'version-info' : {'partis-yaml' : 0.2},
        'grouping-info' : {
            'method' : 'cdr3-length',
            'loci' : loci,
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
                'gene_lists_consistent' : None,
                'uids_unique' : None,
                'sequence_count_preserved' : None,
            },
        },
    }
    manifest_path = '%s/manifest.yaml' % outdir
    utils.mkdir(manifest_path, isfile=True)
    with open(manifest_path, 'w') as mfile:
        yaml.dump(manifest, mfile, width=400, default_flow_style=False)
    print('    wrote manifest to %s' % manifest_path)
    return manifest

# ----------------------------------------------------------------------------------------
def read_manifest(manifest_path):
    # read and validate manifest yaml
    if not os.path.exists(manifest_path):
        raise Exception('manifest file does not exist: %s' % manifest_path)
    with open(manifest_path) as mfile:
        manifest = yaml.safe_load(mfile)
    for required_key in ['version-info', 'grouping-info', 'groups', 'assembly']:
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
    print('    sequence count validated: %d grouped + %d failed = %d total' % (total_grouped, n_failed, total_input))

# ----------------------------------------------------------------------------------------
def get_sw_cache_path(parameter_dir, locus):
    locus_path = '%s/%s/sw-cache.yaml' % (parameter_dir, locus)
    if os.path.exists(locus_path):
        return locus_path
    flat_path = '%s/sw-cache.yaml' % parameter_dir
    if os.path.exists(flat_path):
        return flat_path
    raise Exception('sw cache not found at %s or %s' % (locus_path, flat_path))

# ----------------------------------------------------------------------------------------
def get_loci(args):
    if args.paired_loci:
        return utils.sub_loci(args.ig_or_tr)
    else:
        return [args.locus]

# ----------------------------------------------------------------------------------------
def run_disjoint_group(args):
    if args.disjoint_dir is None:
        raise Exception('--disjoint-dir must be set for disjoint-group')
    if args.parameter_dir is None:
        raise Exception('--parameter-dir must be set for disjoint-group (auto-caching not yet implemented)')  # TODO auto-trigger sw annotation

    outdir = args.disjoint_dir
    loci = get_loci(args)
    print('  running disjoint-group on %s with parameter dir %s' % (' '.join(loci), args.parameter_dir))

    all_group_infos = []
    total_input = 0
    total_failed = 0
    for ltmp in loci:
        sw_cache_path = get_sw_cache_path(args.parameter_dir, ltmp)
        print('    reading sw cache for %s from %s' % (ltmp, sw_cache_path))
        glfo, annotation_list, _ = utils.read_yaml_output(sw_cache_path, dont_add_implicit_info=True)
        groups, n_failed = group_sequences_by_cdr3_length(annotation_list)
        n_seqs = sum(len(seqfos) for seqfos in groups.values()) + n_failed
        total_input += n_seqs
        total_failed += n_failed
        print('    %s: %d sequences in %d cdr3 length groups (%d failed)' % (ltmp, n_seqs - n_failed, len(groups), n_failed))
        group_infos = write_group_fastas(groups, outdir, ltmp)
        all_group_infos.extend(group_infos)

    manifest = write_manifest(all_group_infos, outdir, loci, total_input, total_failed, parameter_dir=args.parameter_dir)
    validate_sequence_count(manifest)

# ----------------------------------------------------------------------------------------
def validate_assembly(manifest, manifest_dir):
    # check all partition files exist and are non-empty, verify uid uniqueness across groups
    all_uids = set()
    total_seqs = 0
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
        glfo, annotation_list, _ = utils.read_yaml_output(full_ppath, dont_add_implicit_info=True)
        group_uids = set()
        for line in annotation_list:
            for uid in line['unique_ids']:
                if uid in all_uids:
                    raise Exception('duplicate uid %s found across groups' % uid)
                group_uids.add(uid)
                all_uids.add(uid)
        total_seqs += len(group_uids)
    if len(missing_partitions) > 0:
        raise Exception('partition files missing for %d groups: %s' % (len(missing_partitions), missing_partitions))
    expected = manifest['grouping-info']['total_grouped_sequences']
    if total_seqs != expected:
        raise Exception('sequence count mismatch after assembly: found %d uids in partition files, expected %d' % (total_seqs, expected))
    print('    assembly validation passed: %d unique sequences across %d groups' % (total_seqs, len(manifest['groups'])))

# ----------------------------------------------------------------------------------------
def run_assemble_groups(args):
    if args.disjoint_dir is None:
        raise Exception('--disjoint-dir must be set for assemble-groups')

    manifest_path = '%s/manifest.yaml' % args.disjoint_dir
    print('  running assemble-groups from %s' % manifest_path)
    manifest = read_manifest(manifest_path)
    manifest_dir = os.path.dirname(os.path.abspath(manifest_path))

    validate_assembly(manifest, manifest_dir)

    # update manifest with validation results
    manifest['assembly']['validation']['uids_unique'] = True
    manifest['assembly']['validation']['sequence_count_preserved'] = True
    manifest['assembly']['status'] = 'validated'
    with open(manifest_path, 'w') as mfile:
        yaml.dump(manifest, mfile, width=400, default_flow_style=False)
    print('    updated manifest with validation results')

    if args.merge_output:
        print('    --merge-output not yet implemented')  # TODO implement merge via merge_paired_yamls()
