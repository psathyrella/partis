import glob
import json
import os
import re
import shutil

import pytest

from _helpers import LOCI, PAIRED_PARAM_DIR, PAIRED_SIMU_DIR, UNPAIRED_IGH_FNAME, copy_dir, disjoint_dir, fixture_log, locus_manifest, manifest_fname, merge_subsets_args, paired_partition_args, read_log, run_partis, run_partis_fails
from partis import disjointgrouper as dg
from partis import ha_repartition
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


def test_ha_jobs_with_parent_parameter_dir(tmp_path, ha_only_dir):
    # the task list holds absolute paths, so point the copy's at the copy
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    task_list = ha_repartition.task_list_fname(disjoint_dir(outdir, 'igh'), 'igh')
    with open(task_list) as tfile:
        tlines = tfile.read().replace(str(ha_only_dir), str(outdir))
    with open(task_list, 'w') as tfile:
        tfile.write(tlines)
    job = ha_repartition.read_task_list(task_list)[0]
    os.remove(job['outfname'])
    run_partis(ha_job_args(outdir, pdir=PAIRED_PARAM_DIR), str(tmp_path / 'partis.log'))
    assert os.path.exists(job['outfname'])


def test_ha_jobs_without_parameter_dir(tmp_path):
    log = run_partis_fails(['run-ha-repartition-jobs', '--locus', 'igh', '--paired-outdir', str(tmp_path / 'out'), '--job-start', '0'], str(tmp_path / 'partis.log'))
    assert 'HA re-partition needs --parameter-dir' in log


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
    results = sorted(glob.glob('%s/**/clusters/*/partition.yaml' % ddir, recursive=True))
    assert len(results) >= 2
    shutil.copy(results[0], results[1])  # cluster 1's result now holds cluster 0's uids
    log = tmp_path / 'partis.log'
    run_partis(['assemble-ha-repartition', '--locus', 'igh', '--paired-outdir', str(outdir)], str(log))
    counts = [int(n) for n in re.findall(r'warning HA results in \S+: ([0-9]+) not covering their cluster, so those clusters are kept whole: cluster-', read_log(str(log)))]
    assert sum(counts) == 1


def test_ha_assemble_with_missing_uid(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    specs = ha_repartition.group_specs(ddir, manifest['groups'], 'igh')
    for spec in specs:  # so assemble redoes every group
        os.remove(spec['harep_out'])
    results = sorted(glob.glob('%s/**/clusters/*/partition.yaml' % ddir, recursive=True))
    with open(results[0]) as f:
        rjs = json.load(f)
    last = rjs['partitions'][-1]['partition']
    big = max(last, key=len)
    assert len(big) > 1
    dropped = big.pop()  # remove one uid from the result
    with open(results[0], 'w') as f:
        json.dump(rjs, f)
    log = tmp_path / 'partis.log'
    run_partis(['assemble-ha-repartition', '--locus', 'igh', '--paired-outdir', str(outdir)], str(log))
    counts = [int(n) for n in re.findall(r'warning HA results in \S+: ([0-9]+) with missing uids, dropped 1 uids that failed HA annotation', read_log(str(log)))]
    assert sum(counts) == 1
    assert 'not covering their cluster' not in read_log(str(log))
    out_uids = set()
    for spec in specs:
        with open(spec['harep_out']) as f:
            out_uids |= set(u for c in json.load(f)['partitions'][-1]['partition'] for u in c)
    assert dropped not in out_uids


def test_ha_plan_bundles():
    big, n_multi = ha_repartition.HA_MULTIPROC_MIN_COLLAPSED, ha_repartition.HA_MULTIPROC_N_PROCS
    jobs = [{'cluster_id': str(i), 'n_collapsed': u, 'n_procs': n_multi if u >= big else 1} for i, u in enumerate([10, big + 100, 40, big + 400, 20, 30])]
    ordered, bundles = ha_repartition.plan_bundles(jobs, 2)
    assert [j['n_collapsed'] for j in ordered[:2]] == [big + 400, big + 100]
    assert bundles[:2] == [(0, 1, n_multi), (1, 1, n_multi)]
    assert all(n_procs == 1 for _, _, n_procs in bundles[2:])
    assert sum(count for _, count, _ in bundles) == len(jobs)
    assert sorted(j['cluster_id'] for j in ordered) == sorted(j['cluster_id'] for j in jobs)


def test_ha_prepare_multiproc_tier(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igh'):
        os.remove(spec['harep_out'])
    # confirm half the candidates, so both outcomes of the confirmation are covered
    confirmed = []
    def count_collapsed(candidates):
        confirmed.extend(j['cluster_id'] for j in candidates[::2])
        return [2 if j['cluster_id'] in confirmed else 1 for j in candidates]
    _, jobs, bundles = ha_repartition.prepare_all(ddir, manifest['groups'], 'igh', 4, count_collapsed=count_collapsed, multiproc_min_collapsed=2)
    multi = [j for j in jobs if j['n_procs'] > 1]
    assert len(multi) > 0 and sorted(j['cluster_id'] for j in multi) == sorted(confirmed)
    assert jobs[:len(multi)] == multi
    assert bundles == ha_repartition.read_bundle_list(ha_repartition.bundle_list_fname(ddir, 'igh'))
    assert ha_repartition.read_task_list(ha_repartition.task_list_fname(ddir, 'igh')) == jobs
    lddir, lmanifest = locus_manifest(outdir, 'igk')
    for spec in ha_repartition.group_specs(lddir, lmanifest['groups'], 'igk'):
        os.remove(spec['harep_out'])
    _, light_jobs, _ = ha_repartition.prepare_all(lddir, lmanifest['groups'], 'igk', 4, multiproc_min_collapsed=2)
    assert len(light_jobs) > 0 and all(j['n_procs'] == 1 for j in light_jobs)


@pytest.mark.parametrize('k, n, f', [(0, 5, .1), (1, 5, .1), (3, 10, .05), (5, 10, .3), (20, 40, .2), (2, 3, .9), (7, 7, .5), (60, 200, .1)])
def test_binomial_sf_matches_scipy(k, n, f):
    from scipy.stats import binom
    from partis.partition_refinement import binomial_sf
    want = 1. if k <= n * f else binom.sf(k - 1, n, f)  # binomial_sf returns 1 at or below the mean
    assert binomial_sf(k, n, f) == pytest.approx(want, rel=1e-9, abs=1e-300)
