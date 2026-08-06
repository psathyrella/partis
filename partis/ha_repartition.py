"""HA re-partition (vsearch base, per-cluster keep-or-split HA).

Pipeline-agnostic, file-based module with a prepare, run, assemble shape, so it works
both in-process (driven by run_cmds) and as external array jobs.

  prepare(): for every vsearch cluster >= min_cluster_size, write a per-cluster FASTA
             and a task list (external array orchestration).
  run_jobs_in_process(): run HA (`partis partition --no-indels`, no --fast, so bcrham
             decides keep-or-split) on each cluster, reusing one pre-loaded glfo.
  assemble(): merge the per-cluster HA results over the vsearch partition (re-partitioned
             clusters replaced by their HA sub-clusters, others kept) and write a full
             partis output file (annotations synthesized from the sw-cache, reusing
             partition_refinement.write_full_output, no recompute).

Per-cluster HA is a single keep-or-split decision with no between-cluster comparisons,
so there is no max-cluster-size guard; min_cluster_size defaults to 3 (singletons and
pairs cannot be improved by HA).
"""
import json
import os


def read_fasta(path):
    seqs, uid, parts = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if uid is not None:
                    seqs[uid] = ''.join(parts)
                uid = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
    if uid is not None:
        seqs[uid] = ''.join(parts)
    return seqs


def _load_best_partition(fname):
    # best partition is partitions[-1] (for keep-or-split HA runs partis writes it at the
    # merged end), not data['events'] (the most-fragmented step's annotations).
    with open(fname) as f:
        data = json.load(f)
    return [list(c) for c in data['partitions'][-1]['partition']]


def cluster_id(idx, cluster):
    return 'cluster-%06d-size-%d' % (idx, len(cluster))


def prepare(partition_fname, fasta_fname, workdir, min_cluster_size=3):
    """Write a per-cluster FASTA for every vsearch cluster >= min_cluster_size, plus a
    task_list.txt (infname, outfname, n_seqs per line) for external array runs. Returns
    the list of job dicts {cluster_id, infname, outfname, n_seqs}."""
    clusters = _load_best_partition(partition_fname)
    seqs = read_fasta(fasta_fname)
    cdir = os.path.join(workdir, 'clusters')
    os.makedirs(cdir, exist_ok=True)
    jobs = []
    for idx, cluster in enumerate(clusters):
        if len(cluster) < min_cluster_size:
            continue
        cid = cluster_id(idx, cluster)
        d = os.path.join(cdir, cid)
        os.makedirs(d, exist_ok=True)
        infa = os.path.join(d, 'input.fa')
        n_written = 0
        with open(infa, 'w') as f:
            for uid in cluster:
                if uid in seqs:
                    f.write('>%s\n%s\n' % (uid, seqs[uid]))
                    n_written += 1
        if n_written != len(cluster):
            print('    warning: ha-repartition prepare wrote %d/%d cluster seqs to %s (rest missing from input fasta)' % (n_written, len(cluster), infa))
        jobs.append({'cluster_id': cid, 'infname': infa,
                     'outfname': os.path.join(d, 'partition.yaml'), 'n_seqs': len(cluster)})
    with open(os.path.join(workdir, 'task_list.txt'), 'w') as f:
        for j in jobs:
            f.write('%s\t%s\t%d\n' % (j['infname'], j['outfname'], j['n_seqs']))
    return jobs


def assemble(partition_fname, workdir, sw_cache_fname, out_fname, min_cluster_size=3):
    """Merge per-cluster HA results over the vsearch partition and write full partis
    output. A cluster >= min_cluster_size with an HA result whose sub-clusters account
    for all its uids is replaced by those sub-clusters; otherwise it is kept as-is
    (conservative: never drop uids). Returns (repartitioned, n_split, n_kept)."""
    from partis import utils, partition_refinement, disjointgrouper
    clusters = _load_best_partition(partition_fname)
    cdir = os.path.join(workdir, 'clusters')
    repartitioned, n_split, n_kept = [], 0, 0
    for idx, cluster in enumerate(clusters):
        result = os.path.join(cdir, cluster_id(idx, cluster), 'partition.yaml')
        # existence is the per-cluster HA job's completion signal; a result that does not cover the
        # cluster is rejected below, so a half-written one keeps the cluster whole rather than losing uids
        if len(cluster) >= min_cluster_size and os.path.exists(result):
            subs = _load_best_partition(result)
            if subs and set().union(*subs) == set(cluster):
                repartitioned.extend(subs)
                n_split += 1 if len(subs) > 1 else 0
                n_kept += 1 if len(subs) == 1 else 0
                continue
        repartitioned.append(cluster)
        n_kept += 1
    # write full output: synthesize each cluster's annotation from the persistent
    # partition-step annotations (full-length HMM naive_seq/input_seqs), not the sw-cache
    # (which stores SW-trimmed, variable-length seqs). partition-refine consumes these
    # downstream and is accurate only on the full-length convention it was validated on.
    glfo, part_antns, part_cpath = utils.read_output(partition_fname)
    disjointgrouper.check_stage_file_complete(partition_fname, part_antns, part_cpath)  # free here: the file is already read
    ant_info = {}
    for antn in part_antns:
        for uid in antn['unique_ids']:
            ant_info[uid] = antn
    partition_refinement.write_full_output(out_fname, glfo, repartitioned, ant_info, label='ha-repartition')
    return repartitioned, n_split, n_kept


# ----------------------------------------------------------------------------
# Locus-level orchestration over the disjoint-groups manifest (shared by the
# integrated run_cmds path and external SLURM-array jobs).
# ----------------------------------------------------------------------------
MIN_CLUSTER_SIZE = 3


def _group_spec(disjoint_dir, group, locus):
    """Resolve the per-group HA re-partition I/O paths from a manifest group entry. Paths
    are derived from the group's fasta_path so this works with only the manifest on disk."""
    from partis import disjointgrouper
    fasta_dir = os.path.dirname(group['fasta_path'])
    harep_rel = '%s/%s' % (fasta_dir, disjointgrouper.stage_fname(disjointgrouper.STAGE_HAREP, locus))
    return {'group': group, 'fasta_dir': fasta_dir,
            'vsearch': '%s/%s/%s' % (disjoint_dir, fasta_dir, disjointgrouper.stage_fname(disjointgrouper.STAGE_VSEARCH, locus)),
            'fasta': '%s/%s' % (disjoint_dir, group['fasta_path']),
            'sw_cache': '%s/%s/%s' % (disjoint_dir, fasta_dir, disjointgrouper.group_sw_cache_fname(locus)),
            'workdir': '%s/%s/ha-repartition' % (disjoint_dir, fasta_dir),
            'harep_out_rel': harep_rel,
            'harep_out': '%s/%s' % (disjoint_dir, harep_rel)}


def group_specs(disjoint_dir, groups, locus):
    """Per-group HA re-partition path bundles for every group whose vsearch partition and
    sw-cache exist (the groups HA re-partition can run on)."""
    specs = []
    for group in groups:
        spec = _group_spec(disjoint_dir, group, locus)
        if os.path.exists(spec['vsearch']) and os.path.exists(spec['sw_cache']):
            specs.append(spec)
    return specs


def task_list_fname(disjoint_dir, locus):
    return '%s/ha-repartition-task-list-%s.txt' % (disjoint_dir, locus)


def bundle_marker_fname(disjoint_dir, job_start):
    # done-marker for the task-list slice starting at <job_start>, written by the slice and
    # looked for by whatever dispatched it
    return '%s/ha-repartition-bundle-%d.done' % (disjoint_dir, job_start)


def prepare_all(disjoint_dir, groups, locus, min_cluster_size=MIN_CLUSTER_SIZE):
    """Write per-cluster HA inputs for every group, plus a locus-level task list (one
    line per cluster: cluster_id, infname, outfname, n_seqs). Groups whose HA re-partition
    already exists emit no jobs. Returns (specs, jobs)."""
    specs = group_specs(disjoint_dir, groups, locus)
    all_jobs = []
    for spec in specs:
        if os.path.exists(spec['harep_out']):  # existence is ha-repartition's completion signal, so a truncated output is never re-made
            continue
        spec['jobs'] = prepare(spec['vsearch'], spec['fasta'], spec['workdir'],
                               min_cluster_size=min_cluster_size)
        all_jobs.extend(spec['jobs'])
    with open(task_list_fname(disjoint_dir, locus), 'w') as f:
        for j in all_jobs:
            f.write('%s\t%s\t%s\t%d\n' % (j['cluster_id'], j['infname'], j['outfname'], j['n_seqs']))
    return specs, all_jobs


def read_task_list(path):
    """Read a locus task list (written by prepare_all) into job dicts."""
    jobs = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            cid, infname, outfname, n_seqs = line.split('\t')
            jobs.append({'cluster_id': cid, 'infname': infname,
                         'outfname': outfname, 'n_seqs': int(n_seqs)})
    return jobs


def run_jobs_in_process(jobs, args, glfo, timing_csv=None):
    """Run HA re-partition on a list of cluster jobs within a single process, reusing one
    pre-loaded glfo. Each cluster is partitioned exactly as a standalone `partis
    partition` would (keep-or-split), writing the same per-cluster output; the only
    difference from a per-cluster subprocess is that Python/import/glfo startup is paid
    once for the whole batch instead of once per cluster (the dominant per-cluster cost
    at these sizes). Returns (cluster_id, n_seqs, wall_sec) timing rows."""
    import copy
    import time
    import csv as _csv
    from partis import seqfileopener
    from partis.partitiondriver import PartitionDriver
    base_workdir = args.workdir
    rows = []
    for i, job in enumerate(jobs):
        if os.path.exists(job['outfname']):  # existence is this cluster's completion signal
            continue
        a = copy.copy(args)
        a.infname = job['infname']
        a.outfname = job['outfname']
        a.workdir = '%s/ha-work-%d' % (base_workdir, i)
        input_info, reco_info, _ = seqfileopener.read_sequence_file(job['infname'], True, args=a, quiet=True)
        t0 = time.time()
        parter = PartitionDriver(a, glfo, input_info, None, reco_info)
        parter.run(['partition'])
        parter.clean()
        rows.append((job['cluster_id'], job['n_seqs'], round(time.time() - t0, 3)))
    if timing_csv and rows:
        write_header = not os.path.exists(timing_csv)
        with open(timing_csv, 'a') as f:
            w = _csv.writer(f)
            if write_header:
                w.writerow(['cluster_id', 'n_seqs', 'wall_sec'])
            w.writerows(rows)
    return rows


def assemble_all(specs, min_cluster_size=MIN_CLUSTER_SIZE):
    """Assemble each group's HA re-partition (writing ha-repartition-<locus>.yaml
    where missing). Returns [(group, harep_out_rel), ...] so the caller can repoint
    downstream inputs to the re-partitioned partitions."""
    out = []
    for spec in specs:
        if not os.path.exists(spec['harep_out']):  # existence is ha-repartition's completion signal
            assemble(spec['vsearch'], spec['workdir'], spec['sw_cache'], spec['harep_out'],
                     min_cluster_size=min_cluster_size)
        out.append((spec['group'], spec['harep_out_rel']))
    return out
