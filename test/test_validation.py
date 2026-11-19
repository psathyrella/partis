import argparse
import glob
import json
import os
import random
import re
import runpy
import shutil
import sys

import pytest

from _helpers import BACKEND, LOCI, N_PROCS, PAIRED_PARAM_DIR, REPO_DIR, PAIRED_SIMU_DIR, UNPAIRED_IGH_FNAME, copy_dir, disjoint_dir, drop_ha_results, fixture_log, ha_group_dir, ha_result_files, locus_manifest, manifest_fname, merge_subsets_args, paired_partition_args, read_ha_results, read_log, run_partis, run_partis_fails, write_ha_results
from partis import disjointgrouper as dg
from partis import ha_repartition
from partis import clusterpath
from partis import waterer
from partis.processargs import single_locus_actions


def paired_disjoint_args(tmp_path, extra_args, **kwargs):
    return paired_partition_args(extra_args, outdir=tmp_path / 'out', **kwargs)


# ----------------------------------------------------------------------------------------
# argument checks that fire before any work
@pytest.mark.parametrize('action', single_locus_actions)
def test_single_locus_action_needs_locus(tmp_path, action):
    log = run_partis_fails([action, '--parameter-dir', PAIRED_PARAM_DIR, '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert '\'%s\' needs --locus set explicitly' % action in log


def test_single_locus_action_locus_with_equals(tmp_path):
    log = run_partis_fails(['create-disjoint-groups', '--locus=igh', '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert 'needs --parameter-dir (or --sw-cachefname)' in log  # past the --locus check


@pytest.mark.parametrize('hfarg', ['--hfrac-min-seqs', '--hfrac-max-bin-size', '--hfrac-merge-factor'])
def test_negative_hfrac_arg(tmp_path, hfarg):
    log = run_partis_fails(paired_disjoint_args(tmp_path, ['--hfrac', hfarg, '-1']), str(tmp_path / 'partis.log'))
    assert '%s must not be negative, but got -1' % hfarg in log


def test_seed_with_disjoint_groups(tmp_path):
    log = run_partis_fails(['partition', '--infname', UNPAIRED_IGH_FNAME, '--seed-unique-id', 'x', '--disjoint-groups', '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert '--seed-unique-id is not yet supported with \'partition --disjoint-groups\'' in log


def test_more_subsets_than_seqs(tmp_path):
    outdir = tmp_path / 'out'
    log = run_partis_fails(['cache-parameters', '--infname', UNPAIRED_IGH_FNAME, '--parameter-dir', str(outdir), '--n-subsets', '100'], str(tmp_path / 'partis.log'))
    assert 'subsets are empty: --n-subsets 100 is too large' in log
    assert glob.glob('%s/**/input-seqs.fa' % outdir, recursive=True) == []


# ----------------------------------------------------------------------------------------
# create-disjoint-groups sw cache input
def test_group_without_parameter_dir(tmp_path):
    log = run_partis_fails(['create-disjoint-groups', '--locus', 'igh', '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert 'needs --parameter-dir (or --sw-cachefname)' in log


def test_group_from_nonexistent_sw_cache(tmp_path):
    swfn = tmp_path / 'sw-cache.yaml'
    log = run_partis_fails(['create-disjoint-groups', '--locus', 'igh', '--sw-cachefname', str(swfn), '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert 'sw cache file does not exist: %s' % swfn in log


def test_group_from_moved_index(tmp_path, chunked_param_index_dir):
    # copying the locus dir away breaks the index's relative paths
    pdir = tmp_path / 'moved'
    copy_dir(chunked_param_index_dir / 'parameters' / 'igh', pdir / 'igh')
    log = run_partis_fails(['create-disjoint-groups', '--locus', 'igh', '--parameter-dir', str(pdir), '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert 'sw cache files listed in %s/igh/%s do not exist' % (pdir, dg.SW_CACHE_INDEX_FNAME) in log


def test_group_from_dir_with_no_sw_cache(tmp_path):
    pdir = tmp_path / 'params'
    os.makedirs(str(pdir / 'igh'))
    log = run_partis_fails(['create-disjoint-groups', '--locus', 'igh', '--parameter-dir', str(pdir), '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert 'no sw cache to group from in %s/igh: neither %s nor %s exists' % (pdir, dg.SW_CACHE_FNAME, dg.SW_CACHE_INDEX_FNAME) in log


def test_group_from_single_locus_cache_name(tmp_path):
    pdir = tmp_path / 'params'
    os.makedirs(str(pdir / 'igh'))
    (pdir / 'igh' / 'sw-cache-abc123.yaml').touch()
    (pdir / 'igh' / dg.group_sw_cache_fname('igh')).touch()  # not a hashed name, so not suggested
    log = run_partis_fails(['create-disjoint-groups', '--locus', 'igh', '--parameter-dir', str(pdir), '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert 'but found %s/igh/sw-cache-abc123.yaml (pass it with --sw-cachefname)' % pdir in log


def test_hfrac_with_no_group_above_min_seqs(tmp_path):
    log = tmp_path / 'partis.log'
    run_partis(['create-disjoint-groups', '--locus', 'igh', '--parameter-dir', PAIRED_PARAM_DIR, '--paired-outdir', str(tmp_path / 'out'), '--hfrac', '--hfrac-min-seqs', '100000'], str(log))
    assert re.search(r'--hfrac has no effect: no cdr3 group qualifies for hfrac vsearch \(--hfrac-min-seqs 100000, largest group [0-9]+ seqs\)', read_log(str(log)))


def test_backfill_with_no_subset_caches(tmp_path, chunked_param_index_dir):
    # merged hmm and no merged sw cache, with the per-subset caches gone
    basedir = copy_dir(chunked_param_index_dir, tmp_path / 'params')
    for swfn in glob.glob('%s/parameters/*/%s' % (basedir, dg.SW_CACHE_INDEX_FNAME)) + glob.glob('%s/**/%s' % (basedir, dg.SW_CACHE_FNAME), recursive=True):
        os.remove(swfn)
    log = tmp_path / 'merge.log'
    run_partis(merge_subsets_args(basedir), str(log))
    assert 'subset-merged input exists but has no merged sw cache, and there are no per-subset caches to merge into it' in read_log(str(log))


# ----------------------------------------------------------------------------------------
# integrated route, before grouping
def test_paired_merge_with_index_only_params(tmp_path, chunked_param_index_dir):
    log = run_partis_fails(paired_disjoint_args(tmp_path, [], param_dir=chunked_param_index_dir / 'parameters'), str(tmp_path / 'partis.log'))
    assert 'only the per-subset index %s: the paired merge after \'partition --disjoint-groups\' needs one' % dg.SW_CACHE_INDEX_FNAME in log
    assert not os.path.exists(disjoint_dir(tmp_path / 'out', 'igh'))


# a light locus breaks, so the check has to run before igh is grouped
def test_ha_without_germline_sets(tmp_path):
    pdir = copy_dir(PAIRED_PARAM_DIR, tmp_path / 'params')
    for gldir in glob.glob('%s/igk/*/germline-sets' % pdir):
        shutil.rmtree(gldir)
    log = run_partis_fails(paired_disjoint_args(tmp_path, ['--ha-repartition'], param_dir=pdir), str(tmp_path / 'partis.log'))
    assert re.search(r'HA re-partition needs germline sets in \S+/igk/[a-z-]+/germline-sets/igk, which does not exist', log)
    assert not os.path.exists(manifest_fname(tmp_path / 'out', 'igh'))


def test_refine_without_mute_freq_tables(tmp_path):
    pdir = copy_dir(PAIRED_PARAM_DIR, tmp_path / 'params')
    shutil.rmtree(str(pdir / 'igk' / 'hmm' / 'mute-freqs'))
    log = run_partis_fails(paired_disjoint_args(tmp_path, ['--partition-refine'], param_dir=pdir), str(tmp_path / 'partis.log'))
    assert 'no per-gene mutation frequency tables' in log
    assert not os.path.exists(manifest_fname(tmp_path / 'out', 'igh'))


def test_dry_run_writes_no_groups(tmp_path):
    log = tmp_path / 'partis.log'
    run_partis(paired_disjoint_args(tmp_path, ['--dry-run']), str(log))
    for locus in LOCI:
        assert '--dry-run: would group %s with parameter dir' % locus in read_log(str(log))
        assert not os.path.exists(manifest_fname(tmp_path / 'out', locus))


def test_skipped_loci_unpaired(unpaired_auto_partition_dir):
    # unpaired input leaves igk and igl empty by design
    unpaired_log = read_log(str(fixture_log(unpaired_auto_partition_dir)))
    for locus in ['igk', 'igl']:
        assert '(expected with unpaired input) %s:' % locus in unpaired_log


def test_skipped_loci_paired(tmp_path):
    indir = copy_dir(PAIRED_SIMU_DIR, tmp_path / 'simu')
    os.remove(str(indir / 'igk.yaml'))
    pdir = copy_dir(PAIRED_PARAM_DIR, tmp_path / 'params')
    shutil.rmtree(str(pdir / 'igl'))
    log = tmp_path / 'partis.log'
    run_partis(paired_disjoint_args(tmp_path, ['--dry-run'], param_dir=pdir, indir=indir), str(log))
    logstr = read_log(str(log))
    assert re.search(r'warning igk: input file missing, skipping disjoint partition', logstr)
    assert re.search(r'warning igl: parameter dir missing, skipping disjoint partition', logstr)


# ----------------------------------------------------------------------------------------
# standalone HA jobs on a copy of an integrated HA run
def ha_job_args(outdir, pdir=os.path.join(PAIRED_PARAM_DIR, 'igh'), job_start=0):
    return ['run-ha-repartition-jobs', '--locus', 'igh', '--parameter-dir', str(pdir), '--paired-outdir', str(outdir), '--job-start', str(job_start), '--job-count', '1']


def ha_copy_with_first_job_to_do(ha_only_dir, tmp_path):
    # the task list holds absolute paths, so point the copy's at the copy
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    task_list = ha_repartition.task_list_fname(disjoint_dir(outdir, 'igh'), 'igh')
    with open(task_list) as tfile:
        tlines = tfile.read().replace(str(ha_only_dir), str(outdir))
    with open(task_list, 'w') as tfile:
        tfile.write(tlines)
    job = ha_repartition.read_task_list(task_list)[0]
    drop_ha_results(ha_group_dir(job), [job['cluster_id']])
    assert ha_repartition.jobs_todo([job]) == [job]
    return outdir, job


def test_ha_jobs_with_parent_parameter_dir(tmp_path, ha_only_dir):
    outdir, job = ha_copy_with_first_job_to_do(ha_only_dir, tmp_path)
    run_partis(ha_job_args(outdir, pdir=PAIRED_PARAM_DIR), str(tmp_path / 'partis.log'))
    assert ha_repartition.jobs_todo([job]) == []


def test_ha_jobs_with_sw_cache_index_only(tmp_path, ha_only_dir):
    # a parameter dir merged with --no-merged-sw-cache has the index: a cluster runs its own sw all the same
    outdir, job = ha_copy_with_first_job_to_do(ha_only_dir, tmp_path)
    pdir = copy_dir(os.path.join(PAIRED_PARAM_DIR, 'igh'), tmp_path / 'igh')
    open(dg.sw_cache_index_fname(str(pdir)), 'w').close()
    run_partis(ha_job_args(outdir, pdir=pdir), str(tmp_path / 'partis.log'))
    assert ha_repartition.jobs_todo([job]) == []


def test_ha_jobs_without_parameter_dir(tmp_path):
    log = run_partis_fails(['run-ha-repartition-jobs', '--locus', 'igh', '--paired-outdir', str(tmp_path / 'out'), '--job-start', '0'], str(tmp_path / 'partis.log'))
    assert 'HA re-partition needs --parameter-dir' in log


def test_ha_plan_with_no_bundles(tmp_path):
    log = run_partis_fails(['create-ha-repartition-jobs', '--locus', 'igh', '--paired-outdir', str(tmp_path / 'out'), '--ha-max-bundles', '0'], str(tmp_path / 'partis.log'))
    assert '--ha-max-bundles must be positive, but got 0' in log


def test_ha_jobs_without_task_list(tmp_path):
    log = run_partis_fails(ha_job_args(tmp_path / 'out'), str(tmp_path / 'partis.log'))
    assert 'HA task list %s does not exist (run create-ha-repartition-jobs first)' % ha_repartition.task_list_fname(disjoint_dir(tmp_path / 'out', 'igh'), 'igh') in log


def test_ha_jobs_with_malformed_task_list(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    task_list = ha_repartition.task_list_fname(disjoint_dir(outdir, 'igh'), 'igh')
    with open(task_list) as tfile:
        n_lines = len(tfile.readlines())
    with open(task_list, 'a') as tfile:
        tfile.write('bad\tline\n')
    log = run_partis_fails(ha_job_args(outdir), str(tmp_path / 'partis.log'))
    assert 'malformed line %d in HA task list %s' % (n_lines + 1, task_list) in log


def test_ha_jobs_past_end(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    log = tmp_path / 'partis.log'
    run_partis(ha_job_args(outdir, job_start=100000), str(log))
    assert re.search(r'warning --job-start 100000 is past the end of the [0-9]+ jobs in', read_log(str(log)))


def test_ha_assemble_with_short_result(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    specs = ha_repartition.group_specs(ddir, manifest['groups'], 'igh')
    for spec in specs:  # so assemble redoes every group
        os.remove(spec['harep_out'])
    results = [(f, res) for f in ha_result_files(ddir) for res in read_ha_results(f)]
    assert len(results) >= 2
    (_, first), (rfname, second) = results[:2]
    write_ha_results(rfname, [dict(res, partition=first['partition']) if res == second else res for res in read_ha_results(rfname)])  # the second cluster's result now holds the first one's uids
    log = run_partis_fails(['assemble-ha-repartition', '--locus', 'igh', '--paired-outdir', str(outdir)], str(tmp_path / 'partis.log'))
    assert re.search(r'1 HA results in \S+ are empty or hold uids from outside their cluster', log)


def test_ha_assemble_with_missing_uid(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    specs = ha_repartition.group_specs(ddir, manifest['groups'], 'igh')
    for spec in specs:  # so assemble redoes every group
        os.remove(spec['harep_out'])
    rfname = ha_result_files(ddir)[0]
    results = read_ha_results(rfname)
    big = max([c for res in results for c in res['partition']], key=len)
    assert len(big) > 1
    dropped = big.pop()  # remove one uid from the result
    write_ha_results(rfname, results)
    log = tmp_path / 'partis.log'
    run_partis(['assemble-ha-repartition', '--locus', 'igh', '--paired-outdir', str(outdir)], str(log))
    counts = [int(n) for n in re.findall(r'warning HA results in \S+: ([0-9]+) miss some of their cluster\'s uids, kept 1 of them as singletons', read_log(str(log)))]
    assert sum(counts) == 1
    out_clusters = []
    for spec in specs:
        with open(spec['harep_out']) as f:
            out_clusters += json.load(f)['partitions'][-1]['partition']
    assert [dropped] in out_clusters


def test_ha_assemble_singletons_for_missing_most_uids(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    specs = ha_repartition.group_specs(ddir, manifest['groups'], 'igh')
    for spec in specs:  # so assemble redoes every group
        os.remove(spec['harep_out'])
    rfname = ha_result_files(ddir)[0]
    results = read_ha_results(rfname)
    all_uids = [u for c in results[0]['partition'] for u in c]
    assert len(all_uids) >= 3
    results[0]['partition'] = [all_uids[:1]]  # result keeps one uid of its cluster
    write_ha_results(rfname, results)
    log = tmp_path / 'partis.log'
    run_partis(['assemble-ha-repartition', '--locus', 'igh', '--paired-outdir', str(outdir)], str(log))
    assert re.search(r'warning HA results in \S+: 1 miss some of their cluster\'s uids, kept [0-9]+ of them as singletons', read_log(str(log)))
    out_clusters = []
    for spec in specs:
        with open(spec['harep_out']) as f:
            out_clusters += json.load(f)['partitions'][-1]['partition']
    assert all([u] in out_clusters for u in all_uids)


def test_ha_prepare_plan_markers(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igh'):
        os.remove(spec['harep_out'])
    for fname in [ha_repartition.task_list_fname(ddir, 'igh'), ha_repartition.bundle_list_fname(ddir, 'igh')] + glob.glob(ha_repartition.bundle_marker_fname(ddir, '*', '*')):
        if os.path.exists(fname):
            os.remove(fname)
    stale = ha_repartition.bundle_marker_fname(ddir, 'oldplan', 12345)
    open(stale, 'w').close()
    ha_repartition.prepare_all(ddir, manifest['groups'], 'igh', 4, write_plan=False)
    assert os.path.exists(stale)
    assert not os.path.exists(ha_repartition.bundle_list_fname(ddir, 'igh')) and not os.path.exists(ha_repartition.task_list_fname(ddir, 'igh'))
    ha_repartition.prepare_all(ddir, manifest['groups'], 'igh', 4)
    assert not os.path.exists(stale)
    assert os.path.exists(ha_repartition.bundle_list_fname(ddir, 'igh'))


def test_ha_prepare_in_parallel_matches_serial(tmp_path, ha_only_dir, capsys):
    def prepared(name, n_procs):
        outdir = copy_dir(ha_only_dir, tmp_path / name)
        ddir, manifest = locus_manifest(outdir, 'igh')
        for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igh'):
            os.remove(spec['harep_out'])
            shutil.rmtree(spec['workdir'], ignore_errors=True)
        _, jobs, bundles = ha_repartition.prepare_all(ddir, manifest['groups'], 'igh', 4, n_procs=n_procs)
        files = {os.path.relpath(f, ddir): open(f).read().replace(str(outdir), '') for f in glob.glob('%s/**/ha-repartition/**' % ddir, recursive=True) if os.path.isfile(f)}
        return [dict(j, infname=os.path.relpath(j['infname'], ddir), outdir=os.path.relpath(j['outdir'], ddir)) for j in jobs], bundles, files, capsys.readouterr().out.replace(str(outdir), '')
    serial, parallel = prepared('serial', 1), prepared('parallel', 3)
    assert len(serial[0]) > 1 and len(serial[2]) > 1 and all(os.path.basename(f) == 'task_list.txt' for f in serial[2])
    assert serial == parallel


def test_ha_assemble_in_parallel_matches_serial(tmp_path, ha_only_dir, capsys):
    def assembled(name, n_procs, memory_bytes=None):
        outdir = copy_dir(ha_only_dir, tmp_path / name)
        ddir, manifest = locus_manifest(outdir, 'igh')
        specs = ha_repartition.group_specs(ddir, manifest['groups'], 'igh')
        for spec in specs:
            os.remove(spec['harep_out'])
        rfname = ha_result_files(ddir)[0]
        results = read_ha_results(rfname)
        max([c for res in results for c in res['partition']], key=len).pop()  # a result that misses a uid, so its group prints a warning
        write_ha_results(rfname, results)
        rels = ha_repartition.assemble_all(specs, n_procs=n_procs, memory_bytes=memory_bytes)
        return [rel for _, rel in rels], [open(s['harep_out']).read().replace(str(outdir), '') for s in specs], capsys.readouterr().out.replace(str(outdir), '')
    serial, parallel = assembled('serial', 1), assembled('parallel', 3)
    assert len(serial[1]) > 1 and 'miss some of their cluster' in serial[2]
    assert serial == parallel
    no_room = assembled('no-room', 3, memory_bytes=1)  # every group alone, with a warning from the second on
    assert no_room[:2] == serial[:2] and no_room[2].count('runs alone') == len(serial[1]) - 1


def test_ha_plan_bundles():
    big, n_multi = ha_repartition.HA_MULTIPROC_MIN_COLLAPSED, ha_repartition.HA_MULTIPROC_N_PROCS
    jobs = [{'cluster_id': str(i), 'infname': 'g%d.fa' % (i % 2), 'n_seqs': u, 'n_collapsed': u, 'n_procs': n_multi if u >= big else 1} for i, u in enumerate([10, big + 100, 40, big + 400, 20, 30])]
    ordered, bundles = ha_repartition.plan_bundles(jobs, 2)
    assert [j['n_collapsed'] for j in ordered[:2]] == [big + 400, big + 100]
    assert bundles[:2] == [(0, 1, n_multi), (1, 1, n_multi)]
    assert all(n_procs == 1 for _, _, n_procs in bundles[2:])
    assert sum(count for _, count, _ in bundles) == len(jobs)
    assert sorted(j['cluster_id'] for j in ordered) == sorted(j['cluster_id'] for j in jobs)
    assert all(len(set(j['infname'] for j in ordered[start : start + count])) == 1 for start, count, _ in bundles)
    # a group with nearly all the cost gets every bundle asked for, the other group one
    jobs = [{'cluster_id': str(i), 'infname': 'g%d.fa' % (i >= 8), 'n_seqs': 50 if i < 8 else 3, 'n_collapsed': 50 if i < 8 else 3, 'n_procs': 1} for i in range(10)]
    ordered, bundles = ha_repartition.plan_bundles(jobs, 4)
    assert [(b.job_count, ordered[b.job_start]['infname']) for b in bundles] == [(2, 'g0.fa')] * 4 + [(2, 'g1.fa')]
    # a group gets no more bundles than its cost over that of the costliest job, which is then alone
    jobs = [{'cluster_id': str(i), 'infname': 'g0.fa', 'n_seqs': u, 'n_collapsed': u, 'n_procs': 1} for i, u in enumerate([100] + [10] * 20)]
    ordered, bundles = ha_repartition.plan_bundles(jobs, 16)
    assert [b.job_count for b in bundles] == [1, 20] and ordered[0]['n_collapsed'] == 100
    # many seqs with few distinct naive seqs cost more than fewer seqs that are all distinct
    jobs = [{'cluster_id': str(i), 'infname': 'g0.fa', 'n_seqs': n, 'n_collapsed': u, 'n_procs': 1} for i, (n, u) in enumerate([(45, 40)] * 6 + [(5000, 5)])]
    ordered, bundles = ha_repartition.plan_bundles(jobs, 16)
    assert [b.job_count for b in bundles] == [1, 6] and ordered[0]['n_seqs'] == 5000


def test_ha_sw_naives_match_read_output(ha_only_dir):
    from partis import utils
    ddir, manifest = locus_manifest(ha_only_dir, 'igh')
    for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igh')[:3]:
        _, antns, _ = utils.read_output(spec['sw_cache'], dont_add_implicit_info=True)
        naives = {uid: antn['naive_seq'] for antn in antns for uid in antn['unique_ids']}
        assert len(naives) > 0 and ha_repartition.read_sw_naives(spec['sw_cache']) == naives


def test_ha_cluster_seqs_from_group_inputs(ha_only_dir, capsys):
    ddir, _ = locus_manifest(ha_only_dir, 'igh')
    job = ha_repartition.read_task_list(ha_repartition.task_list_fname(ddir, 'igh'))[0]
    idx = ha_repartition.cluster_index(job['cluster_id'])
    cluster = ha_repartition._load_best_partition(ha_repartition.group_partition_fname(job['infname'], 'igh'))[idx]
    fasta = ha_repartition.read_fasta(job['infname'])
    assert len(cluster) == job['n_seqs'] and ha_repartition.cluster_seqs(job, 'igh', {}) == [(u, fasta[u]) for u in cluster]
    inputs = {'fasta': job['infname'], 'clusters': [cluster] * (idx + 1), 'seqs': {u: fasta[u] for u in cluster[1:]}}  # a group read whose fasta lacks one seq
    assert ha_repartition.cluster_seqs(job, 'igh', inputs) == [(u, fasta[u]) for u in cluster[1:]]
    assert 'found %d/%d seqs of %s' % (len(cluster) - 1, len(cluster), job['cluster_id']) in capsys.readouterr().out
    with pytest.raises(Exception, match='is not a cluster of'):
        ha_repartition.cluster_seqs(dict(job, cluster_id='cluster-%06d-size-%d' % (idx, job['n_seqs'] + 1)), 'igh', {})
    with pytest.raises(Exception, match='needs a re-run of create-ha-repartition-jobs'):  # an old task list names a per-cluster fasta
        ha_repartition.cluster_seqs(dict(job, infname=os.path.join(str(ha_only_dir), 'none', 'input.fa')), 'igh', {})


def test_ha_prepare_multiproc_tier(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igh'):
        os.remove(spec['harep_out'])
    drop_ha_results(ddir)  # finished clusters skip the confirmation
    # confirm half the candidates, so both outcomes of the confirmation are covered
    confirmed = []
    def count_collapsed(candidates):
        confirmed.extend(j['cluster_id'] for j in candidates[::2])
        return [2 if j['cluster_id'] in confirmed else 1 for j in candidates]
    _, jobs, bundles = ha_repartition.prepare_all(ddir, manifest['groups'], 'igh', 4, count_collapsed=count_collapsed, multiproc_min_collapsed=2)
    multi = [j for j in jobs if j['n_procs'] > 1]
    assert len(multi) > 0 and sorted(j['cluster_id'] for j in multi) == sorted(confirmed)
    assert jobs[:len(multi)] == multi
    assert bundles == ha_repartition.read_bundle_list(ha_repartition.bundle_list_fname(ddir, 'igh'), len(jobs))
    assert ha_repartition.read_task_list(ha_repartition.task_list_fname(ddir, 'igh')) == jobs
    lddir, lmanifest = locus_manifest(outdir, 'igk')
    for spec in ha_repartition.group_specs(lddir, lmanifest['groups'], 'igk'):
        os.remove(spec['harep_out'])
    _, light_jobs, _ = ha_repartition.prepare_all(lddir, lmanifest['groups'], 'igk', 4, multiproc_min_collapsed=2)
    assert len(light_jobs) > 0 and all(j['n_procs'] == 1 for j in light_jobs)


def test_ha_confirm_skips_finished_candidates(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igh'):
        os.remove(spec['harep_out'])
    asked = []
    def count_collapsed(candidates):
        asked.extend((j['outdir'], j['cluster_id']) for j in candidates)
        return [2] * len(candidates)
    _, jobs, _ = ha_repartition.prepare_all(ddir, manifest['groups'], 'igh', 4, count_collapsed=count_collapsed, multiproc_min_collapsed=2)
    multi = [j for j in jobs if j['n_procs'] > 1]
    assert len(multi) >= 2 and ha_repartition.jobs_todo(multi) == [] and asked == []  # the fixture ran every cluster
    drop_ha_results(ha_group_dir(multi[0]), [multi[0]['cluster_id']])
    ha_repartition.prepare_all(ddir, manifest['groups'], 'igh', 4, count_collapsed=count_collapsed, multiproc_min_collapsed=2)
    assert asked == [(multi[0]['outdir'], multi[0]['cluster_id'])]


def test_ha_bundle_list_must_cover_task_list(tmp_path):
    path = str(tmp_path / 'bundles.txt')
    with open(path, 'w') as f:
        f.write('0\t2\t8\n2\t3\t1\n')
    assert ha_repartition.read_bundle_list(path, 5) == [ha_repartition.Bundle(0, 2, 8), ha_repartition.Bundle(2, 3, 1)]
    for n_jobs, text in [(6, '0\t2\t8\n2\t3\t1\n'), (5, '0\t2\t8\n3\t3\t1\n'), (5, '2\t3\t1\n0\t2\t8\n')]:
        with open(path, 'w') as f:
            f.write(text)
        with pytest.raises(Exception, match='does not cover'):
            ha_repartition.read_bundle_list(path, n_jobs)


def test_ha_jobs_refuse_more_procs_than_cpus(tmp_path, monkeypatch):
    from partis import utils
    monkeypatch.setattr(utils, 'n_available_cpus', lambda: 2)
    job = {'cluster_id': 'c', 'infname': 'x.fa', 'outdir': str(tmp_path), 'n_seqs': 3, 'n_collapsed': 3, 'n_procs': 8}
    with pytest.raises(Exception, match='need 8 procs but only 2 cpus'):
        ha_repartition.run_jobs_in_process([job], argparse.Namespace(workdir=str(tmp_path)), None)
    write_ha_results(ha_repartition.results_fname(str(tmp_path), 0), [{'cluster_id': 'c', 'partition': [['x', 'y', 'z']]}])  # a finished job needs no cpus
    assert ha_repartition.run_jobs_in_process([job], argparse.Namespace(workdir=str(tmp_path)), None) == []


def test_ha_result_lines(tmp_path):
    outdir = str(tmp_path)
    jobs = [{'cluster_id': c, 'outdir': outdir} for c in 'ab']
    first = ha_repartition.results_fname(outdir, 0)
    write_ha_results(first, [{'cluster_id': 'a', 'partition': [['x', 'y']]}])
    with open(first, 'a') as f:
        f.write('{"cluster_id": "b", "partition": [["z"')  # the cut last line of a killed slice
    assert ha_repartition.read_results(outdir) == {'a': [['x', 'y']]}
    assert ha_repartition.jobs_todo(jobs) == jobs[1:]
    with ha_repartition._open_results(first) as out:  # the re-run starts on a new line
        out.write('%s\n' % json.dumps({'cluster_id': 'b', 'partition': [['z']]}))
    assert ha_repartition.read_results(outdir) == {'a': [['x', 'y']], 'b': [['z']]}
    write_ha_results(ha_repartition.results_fname(outdir, 5), [{'cluster_id': 'a', 'partition': [['x', 'y']]}])  # the same result from another slice
    assert ha_repartition.jobs_todo(jobs) == []
    write_ha_results(ha_repartition.results_fname(outdir, 7), [{'cluster_id': 'a', 'partition': [['x'], ['y']]}])
    with pytest.raises(Exception, match='HA results of a differ'):
        ha_repartition.read_results(outdir)


def test_ha_assemble_removes_intermediates(tmp_path, ha_only_dir):
    def assembled(name, extra_args):
        outdir = copy_dir(ha_only_dir, tmp_path / name)
        ddir, manifest = locus_manifest(outdir, 'igh')
        specs = ha_repartition.group_specs(ddir, manifest['groups'], 'igh')
        os.remove(specs[0]['harep_out'])  # one group to assemble, the rest already are
        run_partis(['assemble-ha-repartition', '--locus', 'igh', '--paired-outdir', str(outdir)] + extra_args, str(tmp_path / ('%s.log' % name)))
        assert all(os.path.exists(s['harep_out']) for s in specs)
        left = ha_result_files(ddir) + glob.glob('%s/**/ha-repartition' % ddir, recursive=True) + glob.glob(ha_repartition.bundle_marker_fname(ddir, '*', '*'))
        left += [f for f in [ha_repartition.task_list_fname(ddir, 'igh'), ha_repartition.bundle_list_fname(ddir, 'igh')] if os.path.exists(f)]
        return len(left), os.path.exists('%s/ha-repartition-timing-igh.csv' % ddir)
    n_kept, timing_kept = assembled('kept', ['--keep-intermediates'])
    assert n_kept > 4 and timing_kept
    assert assembled('removed', []) == (0, True)


def test_ha_intermediates_stay_for_a_group_not_assembled(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    specs = ha_repartition.group_specs(ddir, manifest['groups'], 'igh')
    os.remove(specs[0]['harep_out'])
    ha_repartition.remove_intermediates(ddir, specs, 'igh')
    assert os.path.exists(os.path.join(specs[0]['workdir'], 'task_list.txt'))
    assert not any(os.path.exists(s['workdir']) for s in specs[1:])
    assert os.path.exists(ha_repartition.task_list_fname(ddir, 'igh')) and os.path.exists(ha_repartition.bundle_list_fname(ddir, 'igh'))


def test_ha_slice_appends_a_group_results_to_one_file(tmp_path, ha_only_dir):
    # one bundle per group, run as one slice: the results of the fixture's run, one file per group
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igk')
    specs = ha_repartition.group_specs(ddir, manifest['groups'], 'igk')
    expected = {s['workdir']: ha_repartition.read_results(s['workdir']) for s in specs}
    assert max(len(res) for res in expected.values()) > 1  # a group with more than one job
    for spec in specs:
        os.remove(spec['harep_out'])
    drop_ha_results(ddir)
    pdir = os.path.join(PAIRED_PARAM_DIR, 'igk')
    run_partis(['create-ha-repartition-jobs', '--locus', 'igk', '--parameter-dir', pdir, '--paired-outdir', str(outdir), '--ha-max-bundles', '1'], str(tmp_path / 'create.log'))
    run_partis(['run-ha-repartition-jobs', '--locus', 'igk', '--parameter-dir', pdir, '--paired-outdir', str(outdir), '--job-start', '0'], str(tmp_path / 'run.log'))
    for workdir, results in expected.items():
        assert ha_repartition.read_results(workdir) == results
        if len(results) > 0:
            assert sorted(res['cluster_id'] for res in read_ha_results(ha_repartition.results_fname(workdir, 0))) == sorted(results)


def _run_partis_in_process(args):
    sys.argv = ['partis'] + args + ['--n-procs', str(N_PROCS), '--random-seed', '1', '--dont-write-git-info'] + (['--zig'] if BACKEND == 'zig' else [])
    os.environ.setdefault('PYTHONHASHSEED', '0')  # otherwise bin/partis re-execs itself, replacing the pytest process
    runpy.run_path(os.path.join(REPO_DIR, 'bin', 'partis'), run_name='__main__')


def test_ha_sw_cache_cleanup_stays_in_the_workdir(tmp_path, ha_only_dir, monkeypatch):
    # sw on a cluster removes no cache outside its workdir: not the --sw-cachefname file, nor one in the parameter dir
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igk')
    for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igk'):
        os.remove(spec['harep_out'])
    drop_ha_results(ddir)
    swfn, workdir, cleaned = tmp_path / 'sw-cache.yaml', tmp_path / 'wd', []
    swfn.write_text('not read\n')
    real_clean = waterer.Waterer.clean_cache
    def clean_cache(self, cache_path):
        cleaned.append(cache_path)
        return real_clean(self, cache_path)
    monkeypatch.setattr(waterer.Waterer, 'clean_cache', clean_cache)
    args = ['--locus', 'igk', '--parameter-dir', os.path.join(PAIRED_PARAM_DIR, 'igk'), '--sw-cachefname', str(swfn), '--paired-outdir', str(outdir)]
    _run_partis_in_process(['create-ha-repartition-jobs'] + args + ['--ha-max-bundles', '1'])
    _run_partis_in_process(['run-ha-repartition-jobs'] + args + ['--job-start', '0', '--workdir', str(workdir)])
    assert len(cleaned) > 0 and all(os.path.realpath(path).startswith('%s/' % os.path.realpath(workdir)) for path in cleaned)
    assert swfn.exists()


def test_ha_collapsed_count_matches_partition(tmp_path, ha_only_dir, monkeypatch):
    # the scheduling count must match the initial partition that partition builds
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igh'):
        os.remove(spec['harep_out'])
    drop_ha_results(ddir)
    real_prepare, real_count = ha_repartition.prepare_all, ha_repartition.count_ha_collapsed
    counts, memory = {}, []
    def count(jobs, args, glfo, n_procs, memory_bytes=None):
        memory.append(memory_bytes)
        res = real_count(jobs, args, glfo, n_procs, memory_bytes=memory_bytes)
        inputs = {}
        assert res == [ha_repartition._ha_collapsed_count(i, ha_repartition.cluster_seqs(j, args.locus, inputs), args, glfo) for i, j in enumerate(jobs)]  # the pool against this process
        counts.update(((j['outdir'], j['cluster_id']), n) for j, n in zip(jobs, res))
        return res
    monkeypatch.setattr(ha_repartition, 'prepare_all', lambda *a, **kw: real_prepare(*a, **dict(kw, multiproc_min_collapsed=2)))
    monkeypatch.setattr(ha_repartition, 'count_ha_collapsed', count)
    pdir = os.path.join(PAIRED_PARAM_DIR, 'igh')
    _run_partis_in_process(['create-ha-repartition-jobs', '--locus', 'igh', '--parameter-dir', pdir, '--paired-outdir', str(outdir), '--pool-memory-gb', '4'])
    jobs = ha_repartition.read_task_list(ha_repartition.task_list_fname(ddir, 'igh'))
    checked = [j for j in jobs if (j['outdir'], j['cluster_id']) in counts][:3]
    assert len(checked) > 0
    assert set(memory) == {4 * 2**30}  # --pool-memory-gb reaches the count
    ha_repartition.write_task_list(ddir, 'igh', [dict(j, n_procs=1) for j in jobs])  # one proc, as the count runs
    initial = []
    real_add = clusterpath.ClusterPath.add_partition
    def add_partition(self, partition, *a, **kw):
        initial.append(len(partition))
        return real_add(self, partition, *a, **kw)
    monkeypatch.setattr(clusterpath.ClusterPath, 'add_partition', add_partition)
    for job in checked:
        del initial[:]
        _run_partis_in_process(['run-ha-repartition-jobs', '--locus', 'igh', '--parameter-dir', pdir, '--paired-outdir', str(outdir),
                                '--job-start', str(jobs.index(job)), '--job-count', '1'])
        assert len(initial) > 0 and initial[0] == counts[(job['outdir'], job['cluster_id'])], job['cluster_id']


def test_ha_cluster_with_no_sw_pass_writes_singletons(tmp_path, ha_only_dir, monkeypatch):
    # a cluster where no query passes sw gets singletons and a collapsed count of 0, and the bundle goes on
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igh'):
        os.remove(spec['harep_out'])
    drop_ha_results(ddir)
    pdir = os.path.join(PAIRED_PARAM_DIR, 'igh')
    _run_partis_in_process(['create-ha-repartition-jobs', '--locus', 'igh', '--parameter-dir', pdir, '--paired-outdir', str(outdir)])
    jobs = ha_repartition.read_task_list(ha_repartition.task_list_fname(ddir, 'igh'))
    assert len(jobs) >= 2
    bad = jobs[0]
    uids = ha_repartition._load_best_partition(ha_repartition.group_partition_fname(bad['infname'], 'igh'))[ha_repartition.cluster_index(bad['cluster_id'])]
    seqs = ha_repartition.read_fasta(bad['infname'])
    rng = random.Random(1)
    for uid in uids:  # the cluster's seqs in its group's fasta
        seqs[uid] = ''.join(rng.choice('ACGT') for _ in range(40))
    with open(bad['infname'], 'w') as f:
        f.writelines('>%s\n%s\n' % us for us in seqs.items())
    real_run = ha_repartition.run_jobs_in_process
    counts = []
    def run(run_jobs, args, glfo, timing_csv=None, job_start=0):
        counts.extend(ha_repartition.count_ha_collapsed(run_jobs[:2], args, glfo, 2))
        return real_run(run_jobs, args, glfo, timing_csv=timing_csv, job_start=job_start)
    monkeypatch.setattr(ha_repartition, 'run_jobs_in_process', run)
    workdir = tmp_path / 'wd'
    _run_partis_in_process(['run-ha-repartition-jobs', '--locus', 'igh', '--parameter-dir', pdir, '--paired-outdir', str(outdir),
                            '--job-start', '0', '--workdir', str(workdir)])
    assert counts[0] == 0 and counts[1] > 0
    assert ha_repartition.jobs_todo(jobs) == []
    assert sorted(ha_repartition.read_results(bad['outdir'])[bad['cluster_id']]) == sorted([u] for u in uids)
    assert glob.glob('%s/ha-*' % workdir) == []
