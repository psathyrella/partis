import glob
import os
import shutil

import pytest

from _helpers import HFRAC_ARGS, LOCI, N_SUBSETS, PAIRED_PARAM_DIR, PAIRED_SIMU_DIR, UNPAIRED_IGH_FNAME, check_subsets_complete, copy_dir, group_stage_path, locus_manifest, paired_partition_args, partition_fname, read_log, read_yaml, run_partis, run_partis_fails, write_yaml
from partis import utils
from partis import disjointgrouper as dg
from partis import ha_repartition


def paired_chunked_args(basedir):
    return ['cache-parameters', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--n-subsets', str(N_SUBSETS), '--paired-outdir', str(basedir)]


def single_chunked_args(basedir):
    return ['cache-parameters', '--infname', UNPAIRED_IGH_FNAME, '--n-subsets', str(N_SUBSETS), '--parameter-dir', str(basedir)]


def external_subset_args(form, sdir):
    # one subset's job, the way an external array template runs it
    infname = '%s/%s' % (sdir, utils.SUBSET_INPUT_FNAME)
    if form == 'paired':
        return ['cache-parameters', '--paired-loci', '--infname', infname, '--input-metafnames', '%s/meta.yaml' % sdir, '--paired-outdir', sdir]
    return ['cache-parameters', '--infname', infname, '--parameter-dir', '%s/parameters' % sdir]


def remove_markers(sdir):
    for fname in [utils.SUBSET_COMPLETE_FNAME, utils.LEGACY_SUBSET_COMPLETE_FNAME]:
        if os.path.exists('%s/%s' % (sdir, fname)):
            os.remove('%s/%s' % (sdir, fname))


# ----------------------------------------------------------------------------------------
# subset completion markers

@pytest.mark.parametrize('form', ['paired', 'single'])
def test_external_cache_parameters_subset_marks_itself(tmp_path, form):
    basedir = tmp_path / 'out'
    cmd_args = paired_chunked_args(basedir) if form == 'paired' else single_chunked_args(basedir)
    run_partis(cmd_args + ['--write-subsets-only'], str(tmp_path / 'split.log'))
    sdir = utils.parameter_subset_dir(str(basedir), 0)
    run_partis(external_subset_args(form, sdir), str(tmp_path / 'subset-0.log'))
    assert utils.subset_is_marked_complete(sdir)
    run_partis(cmd_args, str(tmp_path / 'resume.log'))
    assert 'running 1 subset jobs' in read_log(tmp_path / 'resume.log')  # only the subset that was not run externally
    check_subsets_complete(str(basedir))


def test_external_subset_partition_job_marks_itself(tmp_path):
    outdir = tmp_path / 'out'
    run_partis(['subset-partition', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--disjoint-groups', '--n-subsets', str(N_SUBSETS), '--write-subsets-only', '--paired-outdir', str(outdir)], str(tmp_path / 'split.log'))
    sdir = str(outdir / 'isub-0')
    run_partis(['partition', '--paired-loci', '--disjoint-groups', '--infname', '%s/%s' % (sdir, utils.SUBSET_INPUT_FNAME), '--input-metafnames', '%s/meta.yaml' % sdir, '--paired-outdir', sdir], str(tmp_path / 'subset-0.log'))
    assert utils.subset_is_marked_complete(sdir)
    assert not utils.subset_is_marked_complete(str(outdir / 'isub-1'))


def test_unmarked_subset_partition_reruns(tmp_path, subset_hfrac_dir):
    outdir = copy_dir(subset_hfrac_dir, tmp_path / 'out')
    remove_markers(str(outdir / 'isub-0'))
    log = tmp_path / 'partis.log'
    run_partis(['subset-partition', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--disjoint-groups'] + HFRAC_ARGS + ['--n-subsets', str(N_SUBSETS), '--paired-outdir', str(outdir)], str(log))
    assert '%s has partition output but no completion marker' % (outdir / 'isub-0') in read_log(log)
    assert 'running 1 subset jobs' in read_log(log)
    assert utils.subset_is_marked_complete(str(outdir / 'isub-0'))


def test_legacy_subset_index_refused(tmp_path, chunked_param_dir):
    basedir = copy_dir(chunked_param_dir, tmp_path / 'out')
    index = read_yaml(utils.subset_index_fname(str(basedir)))
    del index['completion_markers']
    write_yaml(utils.subset_index_fname(str(basedir)), index)
    for isub in range(N_SUBSETS):
        remove_markers(utils.parameter_subset_dir(str(basedir), isub))
    shutil.rmtree(str(basedir / 'parameters'))  # so the run gets as far as the subsets
    log = run_partis_fails(paired_chunked_args(basedir), str(tmp_path / 'partis.log'))
    assert 'predates completion markers, so whether its job finished is unknown' in log


def test_merge_refuses_unmarked_subset(tmp_path, chunked_param_dir):
    basedir = copy_dir(chunked_param_dir, tmp_path / 'out')
    for locus in LOCI:  # unmerged, as after externally run subset jobs
        shutil.rmtree(str(basedir / 'parameters' / locus))
    remove_markers(utils.parameter_subset_dir(str(basedir), 1))
    log = run_partis_fails(['merge-parameter-subsets', '--paired-loci', '--paired-outdir', str(basedir)], str(tmp_path / 'merge.log'))
    assert 'merge-parameter-subsets: 1 of %d subsets have no completion marker' % N_SUBSETS in log
    assert utils.parameter_subset_dir(str(basedir), 1) in log


# ----------------------------------------------------------------------------------------
# earlier stages missing at HA and refine

def remove_vsearch_partition(outdir, locus='igh'):
    ddir, manifest = locus_manifest(outdir, locus)
    group = max(manifest['groups'], key=lambda g: g['sequence_count'])
    os.remove(group_stage_path(ddir, group, dg.STAGE_VSEARCH, locus))
    return group


def test_ha_jobs_need_vsearch_partition(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    group = remove_vsearch_partition(outdir)
    log = run_partis_fails(['create-ha-repartition-jobs', '--locus', 'igh', '--paired-outdir', str(outdir)], str(tmp_path / 'partis.log'))
    assert 'are missing their vsearch partition or group sw cache, so the earlier stage did not finish' in log
    assert dg.group_str(group) in log


def test_refine_jobs_need_input_partition(tmp_path, refine_only_dir):
    outdir = copy_dir(refine_only_dir, tmp_path / 'out')
    group = remove_vsearch_partition(outdir)
    log = run_partis_fails(['run-partition-refine-jobs', '--locus', 'igh', '--parameter-dir', os.path.join(PAIRED_PARAM_DIR, 'igh'), '--paired-outdir', str(outdir), '--job-start', '0', '--job-count', '1'], str(tmp_path / 'partis.log'))
    assert 'are missing their input partition (ha-repartition or vsearch) or group sw cache' in log
    assert dg.group_str(group) in log


def test_ha_assemble_refuses_missing_result(tmp_path, ha_only_dir):
    outdir = copy_dir(ha_only_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    for spec in ha_repartition.group_specs(ddir, manifest['groups'], 'igh'):  # so assemble redoes every group
        os.remove(spec['harep_out'])
    results = sorted(glob.glob('%s/**/clusters/*/partition.yaml' % ddir, recursive=True))
    os.remove(results[0])
    log = run_partis_fails(['assemble-ha-repartition', '--locus', 'igh', '--paired-outdir', str(outdir)], str(tmp_path / 'partis.log'))
    assert '1 HA results missing in' in log
    assert os.path.basename(os.path.dirname(results[0])) in log  # the cluster id


# ----------------------------------------------------------------------------------------
# group partitions skipped on re-run

def drop_first_cluster(fname):
    yinfo = read_yaml(fname)
    for ptn in yinfo['partitions']:
        ptn['partition'] = ptn['partition'][1:]
    write_yaml(fname, yinfo)


@pytest.mark.parametrize('damage', ['empty', 'dropped-cluster'])
def test_rerun_refuses_incomplete_group_partition(tmp_path, hfrac_partition_dir, damage):
    outdir = copy_dir(hfrac_partition_dir, tmp_path / 'out')
    ddir, manifest = locus_manifest(outdir, 'igh')
    group = max(manifest['groups'], key=lambda g: g['sequence_count'])
    ppath = group_stage_path(ddir, group, dg.STAGE_VSEARCH, 'igh')
    if damage == 'empty':
        open(ppath, 'w').close()
    else:
        drop_first_cluster(ppath)
    log = run_partis_fails(paired_partition_args(HFRAC_ARGS + ['--paired-outdir', str(outdir)]), str(tmp_path / 'partis.log'))
    assert 'existing partition for %s is incomplete' % dg.group_str(group) in log


# ----------------------------------------------------------------------------------------
# single file and multifile dir side by side

def add_stale_single_file(outdir, locus='igh'):
    sc_fname = partition_fname(outdir, locus, single_chain=True)
    assert os.path.isdir(dg.multifile_dir_path(sc_fname)) and not os.path.exists(sc_fname)
    ddir, manifest = locus_manifest(outdir, locus)
    shutil.copy(group_stage_path(ddir, manifest['groups'][0], dg.STAGE_VSEARCH, locus), sc_fname)
    return sc_fname


def test_read_refuses_both_output_forms(tmp_path, plain_multifile_partition_dir):
    sc_fname = add_stale_single_file(copy_dir(plain_multifile_partition_dir, tmp_path / 'out'))
    with pytest.raises(Exception, match='exists both as a single file and as a multifile dir'):
        utils.read_output(sc_fname, skip_annotations=True)


def test_assemble_refuses_stale_single_file(tmp_path, plain_multifile_partition_dir):
    outdir = copy_dir(plain_multifile_partition_dir, tmp_path / 'out')
    sc_fname = add_stale_single_file(outdir)
    log = run_partis_fails(['assemble-groups', '--locus', 'igh', '--paired-outdir', str(outdir), '--outfname', sc_fname, '--multifile-min-seqs', '1', '--multifile-max-seqs-per-file', '10'], str(tmp_path / 'partis.log'))
    assert 'single file output from an earlier run is in the way of the multifile output for igh' in log
