from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import sys
import json
import yaml
import glob
import shutil
import collections

from . import utils
from . import glutils
from . import paircluster

MANIFEST_FNAME = 'manifest.yaml'

# ----------------------------------------------------------------------------------------
def group_sequences_by_cdr3_length(annotation_list):
    # group uids and their input sequences by cdr3 length from sw annotation list
    # returns {cdr3_length : [{'name': uid, 'seq': seq}, ...]}, n_failed
    # uses utils.group_seqs_by_value() for the grouping, then collects sequences per group
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
    # write per-group fasta files using utils.write_fasta(), returns list of group info dicts for the manifest
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
        print('    group %d: cdr3 length %d, %d sequences -> %s' % (gid, c3len, len(seqfos), rel_fasta_path))
    return group_infos

# ----------------------------------------------------------------------------------------
def write_manifest(group_infos, outdir, loci, total_input, n_failed, parameter_dir=None):
    # write manifest yaml to outdir
    manifest = {
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
    print('    wrote manifest to %s' % manifest_path)
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
    print('    sequence count validated: %d grouped + %d failed = %d total' % (total_grouped, n_failed, total_input))

# ----------------------------------------------------------------------------------------
def get_sw_cache_path(parameter_dir, locus, paired):
    # look in the expected location based on paired vs unpaired, rather than searching both
    search_dir = '%s/%s' % (parameter_dir, locus) if paired else parameter_dir
    fnames = glob.glob(search_dir + '/sw-cache*.yaml')
    if len(fnames) == 0:
        return None
    if len(fnames) > 1:
        # multiple cache files from different inputs (partis does not clean up old sw-cache-{hash}.yaml files
        # when re-running cache-parameters with different input). Use most recent by mtime, since sw/hmm
        # parameter dirs are overwritten on re-run so the newest cache file matches the current parameters.
        fnames.sort(key=os.path.getmtime, reverse=True)
        print('  %s found %d sw cache files in %s, using most recent: %s' % (utils.color('yellow', 'warning'), len(fnames), search_dir, fnames[0]))
    return fnames[0]

# ----------------------------------------------------------------------------------------
def run_disjoint_group(args):
    if args.disjoint_dir is None:
        base = utils.non_none([getattr(args, 'paired_outdir', None), getattr(args, 'workdir', None)])
        if base is None:
            raise Exception('--disjoint-dir, --paired-outdir, or --workdir must be set for disjoint-group')
        args.disjoint_dir = '%s/disjoint-groups' % base
        print('  note: --disjoint-dir not set, using %s' % args.disjoint_dir)
    has_input = args.infname is not None or getattr(args, 'paired_indir', None) is not None
    if not has_input and args.parameter_dir is None:
        raise Exception('--infname (or --paired-indir with --paired-loci) or --parameter-dir must be set for disjoint-group')

    outdir = args.disjoint_dir
    loci = utils.sub_loci(args.ig_or_tr) if args.paired_loci else [args.locus]

    if args.parameter_dir is None:
        # match getpdir() in run_all_loci(): prefer {paired_outdir}/parameters, fall back to _output/ convention
        if getattr(args, 'paired_outdir', None) is not None:
            args.parameter_dir = '%s/parameters' % args.paired_outdir
        else:
            instr = args.paired_indir if args.paired_loci and getattr(args, 'paired_indir', None) is not None else args.infname
            args.parameter_dir = '_output/%s' % utils.getprefix(instr).replace('/', '_')
        print('  note: --parameter-dir not set, so using default: %s' % args.parameter_dir)

    # auto-trigger sw annotation if sw cache does not exist
    auto_cache_parameters(args, loci)

    print('  running disjoint-group on %s with parameter dir %s' % (' '.join(loci), args.parameter_dir))

    all_group_infos = []
    total_input = 0
    total_failed = 0
    for ltmp in loci:
        sw_cache_path = get_sw_cache_path(args.parameter_dir, ltmp, args.paired_loci)
        if sw_cache_path is None:
            search_dir = '%s/%s' % (args.parameter_dir, ltmp) if args.paired_loci else args.parameter_dir
            raise Exception('sw cache not found in %s/ (run cache-parameters --only-smith-waterman first)' % search_dir)
        print('    reading sw cache for %s from %s' % (ltmp, sw_cache_path))
        _, annotation_list, _ = utils.read_yaml_output(sw_cache_path, dont_add_implicit_info=True)
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
def auto_cache_parameters(args, loci):
    # run cache-parameters --only-smith-waterman for each locus that does not have an sw cache
    # TODO run_step() in run_all_loci() handles full arg forwarding via sys.argv deep copy,
    # but it is a nested function with closure dependencies and not callable from here.
    # This only forwards --n-procs, --is-simu, and --dry-run. See CR-01 in tasks-code-review.md.
    partis_cmd = 'partis'
    if not shutil.which('partis'):
        partis_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        partis_path = '%s/bin/partis' % partis_dir
        if os.path.exists(partis_path):
            partis_cmd = partis_path
        else:
            raise Exception('could not find partis binary in PATH or at %s' % partis_path)
    for ltmp in loci:
        if get_sw_cache_path(args.parameter_dir, ltmp, args.paired_loci) is not None:
            print('    %s: sw cache already exists in %s' % (ltmp, args.parameter_dir))
            continue
        if args.paired_loci:
            if args.paired_indir is not None:
                yfn, ffn = [paircluster.paired_fn(args.paired_indir, ltmp, suffix=sx) for sx in ['.yaml', '.fa']]
                if os.path.exists(yfn) and os.path.exists(ffn):
                    raise Exception('both %s and %s exist, not sure which to use' % (yfn, ffn))
                infname = yfn if os.path.exists(yfn) else ffn
            else:
                raise Exception('--paired-loci requires --paired-indir for auto-caching')
            pdir = '%s/%s' % (args.parameter_dir, ltmp)
        else:
            infname = args.infname
            pdir = args.parameter_dir
        print('    %s: running cache-parameters --only-smith-waterman (input: %s, parameter-dir: %s)' % (ltmp, infname, pdir))
        cmd = '%s cache-parameters --infname %s --parameter-dir %s --locus %s --only-smith-waterman' % (partis_cmd, infname, pdir, ltmp)
        if hasattr(args, 'n_procs') and args.n_procs is not None:
            cmd += ' --n-procs %d' % args.n_procs
        if hasattr(args, 'is_simu') and args.is_simu:
            cmd += ' --is-simu'
        utils.simplerun(cmd, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def get_partition_paths_by_locus(manifest, manifest_dir):
    # collect and verify partition file paths grouped by locus
    paths_by_locus = collections.OrderedDict()
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
        ltmp = ginfo['locus']
        if ltmp not in paths_by_locus:
            paths_by_locus[ltmp] = []
        paths_by_locus[ltmp].append(full_ppath)
    if len(missing_partitions) > 0:
        raise Exception('partition files missing for %d groups: %s' % (len(missing_partitions), missing_partitions))
    return paths_by_locus

# ----------------------------------------------------------------------------------------
def validate_assembly(manifest, manifest_dir):
    # validate uid uniqueness and sequence counts by reading each group one at a time
    all_uids = set()
    total_seqs = 0
    paths_by_locus = get_partition_paths_by_locus(manifest, manifest_dir)
    for ltmp, yaml_list in paths_by_locus.items():
        for ppath in yaml_list:
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
    print('    assembly validation passed: %d unique sequences across %d groups' % (total_seqs, len(manifest['groups'])))

# ----------------------------------------------------------------------------------------
def assemble_merged_output(manifest, disjoint_dir):
    # merge per-group partition yamls into single per-locus output using merge_yamls with best_partition_only
    assembled_dir = '%s/assembled' % disjoint_dir
    utils.mkdir(assembled_dir)
    paths_by_locus = get_partition_paths_by_locus(manifest, disjoint_dir)
    headers = list(utils.annotation_headers)
    for ltmp, yaml_list in paths_by_locus.items():
        outfname = '%s/partition-%s.yaml' % (assembled_dir, ltmp)
        print('    merging %d partition files for %s -> %s' % (len(yaml_list), ltmp, outfname))
        utils.merge_yamls(outfname, yaml_list, headers, best_partition_only=True, dont_write_git_info=True, debug=True)  # debug=True is intentional: prints per-group sequence/cluster counts during assembly
    manifest['assembly']['merged_output_path'] = 'assembled/'

# ----------------------------------------------------------------------------------------
def run_assemble_groups(args):
    if args.disjoint_dir is None:
        base = utils.non_none([getattr(args, 'paired_outdir', None), getattr(args, 'workdir', None)])
        if base is None:
            raise Exception('--disjoint-dir, --paired-outdir, or --workdir must be set for assemble-groups')
        args.disjoint_dir = '%s/disjoint-groups' % base
        print('  note: --disjoint-dir not set, using %s' % args.disjoint_dir)

    manifest_path = '%s/%s' % (args.disjoint_dir, MANIFEST_FNAME)
    print('  running assemble-groups from %s' % manifest_path)
    manifest = read_manifest(manifest_path)
    disjoint_dir = os.path.abspath(args.disjoint_dir)

    validate_assembly(manifest, disjoint_dir)
    manifest['assembly']['validation']['uids_unique'] = True
    manifest['assembly']['validation']['sequence_count_preserved'] = True

    if args.no_merge_output:
        manifest['assembly']['status'] = 'validated'
        print('    --no-merge-output: skipping merged output write (per-group files remain separate)')
    else:
        assemble_merged_output(manifest, disjoint_dir)
        manifest['assembly']['status'] = 'merged'

    with open(manifest_path, 'w') as mfile:
        yaml.dump(manifest, mfile, width=400, default_flow_style=False)
    print('    updated manifest')
