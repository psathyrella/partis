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
def run_disjoint_group(args):
    print('  running disjoint-group')  # placeholder

# ----------------------------------------------------------------------------------------
def run_assemble_groups(args):
    print('  running assemble-groups')  # placeholder
