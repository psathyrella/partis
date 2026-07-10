"""Hybrid HA correction (vsearch + per-cluster HA).

Pipeline-agnostic, file-based module with a prepare -> run -> assemble shape, so it
works both in-memory (driven by run_cmds) and as external array jobs.

  prepare(): for every vsearch cluster >= min_cluster_size, write a per-cluster FASTA
             and a task list (external array orchestration).
  ha_command(): the per-cluster HA command (`partis partition --no-indels`, no --fast,
             so it runs HA/bcrham on the isolated cluster: keep-or-split).
  assemble(): merge the per-cluster HA results over the vsearch partition (corrected
             clusters replaced by their HA sub-clusters, others kept) and write a full
             partis output file (annotations synthesized from the sw-cache, reusing
             refine.write_full_output -- no recompute).

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


def _load_clusters(partition_fname):
    with open(partition_fname) as f:
        data = json.load(f)
    return [list(c) for c in data['partitions'][-1]['partition']]


def _load_ha_subclusters(result_fname):
    with open(result_fname) as f:
        data = json.load(f)
    return [list(evt['unique_ids']) for evt in data.get('events', []) if not evt.get('invalid')]


def cluster_id(idx, cluster):
    return 'cluster-%06d-size-%d' % (idx, len(cluster))


def prepare(partition_fname, fasta_fname, workdir, min_cluster_size=3):
    """Write a per-cluster FASTA for every vsearch cluster >= min_cluster_size, plus a
    task_list.txt (infname, outfname, n_seqs per line) for external array runs. Returns
    the list of job dicts {cluster_id, infname, outfname, n_seqs}."""
    clusters = _load_clusters(partition_fname)
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
        with open(infa, 'w') as f:
            for uid in cluster:
                if uid in seqs:
                    f.write('>%s\n%s\n' % (uid, seqs[uid]))
        jobs.append({'cluster_id': cid, 'infname': infa,
                     'outfname': os.path.join(d, 'partition.yaml'), 'n_seqs': len(cluster)})
    with open(os.path.join(workdir, 'task_list.txt'), 'w') as f:
        for j in jobs:
            f.write('%s\t%s\t%d\n' % (j['infname'], j['outfname'], j['n_seqs']))
    return jobs


def _ha_args(job, partis_bin, parameter_dir, locus):
    """Per-cluster HA arg list. No --fast, so partis runs HA (bcrham) on the isolated
    cluster and decides keep-or-split."""
    return [partis_bin, 'partition', '--infname', job['infname'],
            '--parameter-dir', parameter_dir, '--outfname', job['outfname'],
            '--locus', locus, '--n-procs', '1', '--no-indels', '--random-seed', '0',
            '--refuse-to-cache-parameters']


def ha_command(job, partis_bin, parameter_dir, locus):
    return ' '.join(_ha_args(job, partis_bin, parameter_dir, locus))


def run_batch(jobs, partis_bin, parameter_dir, locus, timing_csv=None):
    """Run HA on a list of cluster jobs, one partis subprocess per cluster, timing each.
    Used by both internal (run_cmds over batches) and external (array) modes. When
    timing_csv is set, appends (cluster_id, n_seqs, wall_sec) rows to build the
    per-cluster HA runtime-vs-size curve. Returns the timing rows."""
    import subprocess
    import time
    import csv as _csv
    rows = []
    for job in jobs:
        if os.path.exists(job['outfname']):
            continue
        t0 = time.time()
        subprocess.run(_ha_args(job, partis_bin, parameter_dir, locus), check=False)
        rows.append((job['cluster_id'], job['n_seqs'], round(time.time() - t0, 3)))
    if timing_csv and rows:
        write_header = not os.path.exists(timing_csv)
        with open(timing_csv, 'a') as f:
            w = _csv.writer(f)
            if write_header:
                w.writerow(['cluster_id', 'n_seqs', 'wall_sec'])
            w.writerows(rows)
    return rows


def assemble(partition_fname, workdir, sw_cache_fname, out_fname, min_cluster_size=3):
    """Merge per-cluster HA results over the vsearch partition and write full partis
    output. A cluster >= min_cluster_size with an HA result whose sub-clusters account
    for all its uids is replaced by those sub-clusters; otherwise it is kept as-is
    (conservative: never drop uids). Returns (hybrid_partition, n_split, n_kept)."""
    from partis import utils, refine
    clusters = _load_clusters(partition_fname)
    cdir = os.path.join(workdir, 'clusters')
    hybrid, n_split, n_kept = [], 0, 0
    for idx, cluster in enumerate(clusters):
        result = os.path.join(cdir, cluster_id(idx, cluster), 'partition.yaml')
        if len(cluster) >= min_cluster_size and os.path.exists(result):
            subs = _load_ha_subclusters(result)
            if subs and sum(len(s) for s in subs) == len(cluster):
                hybrid.extend(subs)
                n_split += 1 if len(subs) > 1 else 0
                n_kept += 1 if len(subs) == 1 else 0
                continue
        hybrid.append(cluster)
        n_kept += 1
    # write full output: synthesize each cluster's annotation from the sw-cache
    glfo, sw_antns, _ = utils.read_output(sw_cache_fname)
    sw_info = {}
    for antn in sw_antns:
        for uid in antn['unique_ids']:
            sw_info[uid] = antn
    refine.write_full_output(out_fname, glfo, hybrid, sw_info)
    return hybrid, n_split, n_kept


# ----------------------------------------------------------------------------
# Locus-level orchestration over the disjoint-groups manifest. These drive the
# per-group prepare/assemble above, so the integrated (run_cmds) and external
# (SLURM array) paths share one code path.
# ----------------------------------------------------------------------------
MIN_CLUSTER_SIZE = 3


def _group_spec(disjoint_dir, group, locus):
    """Resolve the per-group hybrid I/O paths from a manifest group entry. Paths are
    derived from the group's fasta_path so this works with only the manifest on disk."""
    fasta_dir = os.path.dirname(group['fasta_path'])
    return {'group': group, 'fasta_dir': fasta_dir,
            'vsearch': '%s/%s/partition-%s.yaml' % (disjoint_dir, fasta_dir, locus),
            'fasta': '%s/%s' % (disjoint_dir, group['fasta_path']),
            'sw_cache': '%s/%s/sw-cache-%s.yaml' % (disjoint_dir, fasta_dir, locus),
            'workdir': '%s/%s/hybrid' % (disjoint_dir, fasta_dir),
            'hybrid_out_rel': '%s/hybrid-partition-%s.yaml' % (fasta_dir, locus),
            'hybrid_out': '%s/%s/hybrid-partition-%s.yaml' % (disjoint_dir, fasta_dir, locus)}


def group_specs(disjoint_dir, groups, locus):
    """Per-group hybrid path bundles for every group whose vsearch partition and
    sw-cache exist (the groups hybrid can run on)."""
    specs = []
    for group in groups:
        spec = _group_spec(disjoint_dir, group, locus)
        if os.path.exists(spec['vsearch']) and os.path.exists(spec['sw_cache']):
            specs.append(spec)
    return specs


def task_list_fname(disjoint_dir, locus):
    return '%s/hybrid-task-list-%s.txt' % (disjoint_dir, locus)


def prepare_all(disjoint_dir, groups, locus, min_cluster_size=MIN_CLUSTER_SIZE):
    """Write per-cluster HA inputs for every group, plus a locus-level task list (one
    line per cluster: cluster_id, infname, outfname, n_seqs). Groups whose hybrid
    partition already exists emit no jobs. Returns (specs, jobs)."""
    specs = group_specs(disjoint_dir, groups, locus)
    all_jobs = []
    for spec in specs:
        if os.path.exists(spec['hybrid_out']):
            continue
        spec['jobs'] = prepare(spec['vsearch'], spec['fasta'], spec['workdir'],
                               min_cluster_size=min_cluster_size)
        all_jobs.extend(spec['jobs'])
    with open(task_list_fname(disjoint_dir, locus), 'w') as f:
        for j in all_jobs:
            f.write('%s\t%s\t%s\t%d\n' % (j['cluster_id'], j['infname'], j['outfname'], j['n_seqs']))
    return specs, all_jobs


def read_task_list(path):
    """Read a locus task list (written by prepare_all) into run_batch job dicts."""
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
    """Run hybrid HA on a list of cluster jobs within a single process, reusing one
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
        if os.path.exists(job['outfname']):
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
    """Assemble each group's hybrid partition (writing hybrid-partition-<locus>.yaml
    where missing). Returns [(group, hybrid_out_rel), ...] so the caller can repoint
    downstream inputs to the hybrid partitions."""
    out = []
    for spec in specs:
        if not os.path.exists(spec['hybrid_out']):
            assemble(spec['vsearch'], spec['workdir'], spec['sw_cache'], spec['hybrid_out'],
                     min_cluster_size=min_cluster_size)
        out.append((spec['group'], spec['hybrid_out_rel']))
    return out
