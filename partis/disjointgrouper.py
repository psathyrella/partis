from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import sys
import json
import yaml
import collections
import shutil

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
    import glob
    # try per-locus subdir first (paired data), then flat (unpaired)
    for search_dir in ['%s/%s' % (parameter_dir, locus), parameter_dir]:
        fnames = glob.glob(search_dir + '/sw-cache*.yaml')
        if len(fnames) > 0:
            return fnames[0]
    raise Exception('sw cache not found in %s/%s/ or %s/' % (parameter_dir, locus, parameter_dir))

# ----------------------------------------------------------------------------------------
def get_loci(args):
    if args.paired_loci:
        return utils.sub_loci(args.ig_or_tr)
    else:
        return [args.locus]

# ----------------------------------------------------------------------------------------
def find_partis_cmd():
    # find the partis binary, same logic as find_cmd() in bin/partis
    if shutil.which('partis'):
        return 'partis'
    partis_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    partis_path = '%s/bin/partis' % partis_dir
    if os.path.exists(partis_path):
        return partis_path
    raise Exception('could not find partis binary in PATH or at %s' % partis_path)

# ----------------------------------------------------------------------------------------
def get_infname_for_locus(args, locus):
    # get the input file path for a given locus
    if args.paired_loci:
        if args.paired_indir is not None:
            # try both .yaml and .fa suffixes (same logic as getifn() in run_all_loci)
            yfn, ffn = [paircluster.paired_fn(args.paired_indir, locus, suffix=sx) for sx in ['.yaml', '.fa']]
            if os.path.exists(yfn) and os.path.exists(ffn):
                raise Exception('both %s and %s exist, not sure which to use' % (yfn, ffn))
            return yfn if os.path.exists(yfn) else ffn
        elif args.infname is not None:
            raise Exception('--paired-loci with --infname requires locus splitting (use --paired-indir with pre-split files, or run split-loci.py first)')
        else:
            raise Exception('--paired-indir or --infname must be set for auto-caching with --paired-loci')
    else:
        return args.infname

# ----------------------------------------------------------------------------------------
def get_parameter_dir_for_locus(args, locus):
    # get the per-locus parameter dir (paired data uses subdirs per locus)
    if args.paired_loci:
        return '%s/%s' % (args.parameter_dir, locus)
    else:
        return args.parameter_dir

# ----------------------------------------------------------------------------------------
def auto_cache_parameters(args, loci):
    # run cache-parameters --only-smith-waterman for each locus that does not have an sw cache
    partis_cmd = find_partis_cmd()
    for ltmp in loci:
        try:
            get_sw_cache_path(args.parameter_dir, ltmp)
            print('    %s: sw cache already exists in %s' % (ltmp, args.parameter_dir))
            continue
        except Exception:
            pass  # no sw cache found, need to run cache-parameters
        infname = get_infname_for_locus(args, ltmp)
        pdir = get_parameter_dir_for_locus(args, ltmp)
        print('    %s: running cache-parameters --only-smith-waterman (input: %s, parameter-dir: %s)' % (ltmp, infname, pdir))
        cmd = '%s cache-parameters --infname %s --parameter-dir %s --locus %s --only-smith-waterman' % (partis_cmd, infname, pdir, ltmp)
        if hasattr(args, 'n_procs') and args.n_procs is not None:
            cmd += ' --n-procs %d' % args.n_procs
        if hasattr(args, 'is_simu') and args.is_simu:
            cmd += ' --is-simu'
        utils.simplerun(cmd)

# ----------------------------------------------------------------------------------------
def run_disjoint_group(args):
    if args.disjoint_dir is None:
        raise Exception('--disjoint-dir must be set for disjoint-group')
    has_input = args.infname is not None or (hasattr(args, 'paired_indir') and args.paired_indir is not None)
    if not has_input and args.parameter_dir is None:
        raise Exception('--infname (or --paired-indir with --paired-loci) or --parameter-dir must be set for disjoint-group')

    outdir = args.disjoint_dir
    loci = get_loci(args)

    # set default parameter dir if not provided (same convention as run_partitiondriver in bin/partis)
    if args.parameter_dir is None:
        instr = args.paired_indir if args.paired_loci and hasattr(args, 'paired_indir') and args.paired_indir is not None else args.infname
        args.parameter_dir = '_output/%s' % utils.getprefix(instr).replace('/', '_')
        print('  note: --parameter-dir not set, so using default: %s' % args.parameter_dir)

    # auto-trigger sw annotation if sw cache does not exist
    auto_cache_parameters(args, loci)

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
    # (does not load all annotations into memory at once, so scales to large datasets)
    all_uids = set()
    total_seqs = 0
    paths_by_locus = get_partition_paths_by_locus(manifest, manifest_dir)
    for ltmp, yaml_list in paths_by_locus.items():
        for ppath in yaml_list:
            glfo, annotation_list, _ = utils.read_yaml_output(ppath, dont_add_implicit_info=True)
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
def assemble_merged_output(manifest, manifest_dir, disjoint_dir):
    # merge per-group partition yamls into single per-locus output
    # cannot use utils.merge_yamls() because per-group partition files have different numbers
    # of partition steps (from hierarchical agglomeration), and merge_yamls() asserts they match.
    # instead, read each file, extract the best partition and annotations, reconcile germlines,
    # and write a combined output.
    from .clusterpath import ClusterPath
    assembled_dir = '%s/assembled' % disjoint_dir
    utils.mkdir(assembled_dir)
    paths_by_locus = get_partition_paths_by_locus(manifest, manifest_dir)
    for ltmp, yaml_list in paths_by_locus.items():
        outfname = '%s/partition-%s.yaml' % (assembled_dir, ltmp)
        print('    merging %d partition files for %s -> %s' % (len(yaml_list), ltmp, outfname))
        merged_annotation_list = []
        merged_partition = []
        merged_glfo = None
        for infname in yaml_list:
            glfo, annotation_list, cpath = utils.read_yaml_output(infname, dont_add_implicit_info=True)
            print('        %d sequences in %d clusters from %s' % (sum(len(l['unique_ids']) for l in annotation_list), len(annotation_list), infname))
            # reconcile germline info
            if merged_glfo is None:
                merged_glfo = glfo
            elif glfo is not None:
                merged_glfo, name_mapping = glutils.get_merged_glfo(glfo, merged_glfo)
                utils.update_gene_names_in_annotation_list(merged_annotation_list, name_mapping)
            merged_annotation_list += annotation_list
            # take only the best partition from each group
            if cpath is not None and cpath.i_best is not None:
                merged_partition += cpath.partitions[cpath.i_best]
        # build a single-entry cluster path with the combined best partitions
        merged_cpath = ClusterPath()
        merged_cpath.add_partition(merged_partition, logprob=0., n_procs=1)
        utils.write_annotations(outfname, merged_glfo, merged_annotation_list, utils.annotation_headers,
                                partition_lines=merged_cpath.get_partition_lines(), dont_write_git_info=True)
    manifest['assembly']['merged_output_path'] = 'assembled/'

# ----------------------------------------------------------------------------------------
def run_assemble_groups(args):
    if args.disjoint_dir is None:
        raise Exception('--disjoint-dir must be set for assemble-groups')

    manifest_path = '%s/manifest.yaml' % args.disjoint_dir
    print('  running assemble-groups from %s' % manifest_path)
    manifest = read_manifest(manifest_path)
    manifest_dir = os.path.dirname(os.path.abspath(manifest_path))

    validate_assembly(manifest, manifest_dir)
    manifest['assembly']['validation']['uids_unique'] = True
    manifest['assembly']['validation']['sequence_count_preserved'] = True

    if args.no_merge_output:
        manifest['assembly']['status'] = 'validated'
        print('    --no-merge-output: skipping merged output write (per-group files remain separate)')
    else:
        assemble_merged_output(manifest, manifest_dir, args.disjoint_dir)
        manifest['assembly']['status'] = 'merged'

    with open(manifest_path, 'w') as mfile:
        yaml.dump(manifest, mfile, width=400, default_flow_style=False)
    print('    updated manifest')
