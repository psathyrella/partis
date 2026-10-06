import glob
import os
import re
import shutil

import pytest

from _helpers import LOCI, MIN_ANNOTATED_FRAC, N_PROCS, N_SUBSETS, PAIRED_DATA_DIR, PAIRED_SIMU_DIR, UNPAIRED_IGH_FNAME, check_subsets_complete, copy_dir, fasta_uids, fixture_log, input_uids, merge_subsets_args, read_log, run_partis, run_partis_fails, sw_cache_uids
from partis import utils
from partis.disjointgrouper import SW_CACHE_FNAME, SW_CACHE_INDEX_FNAME

# ----------------------------------------------------------------------------------------
# chunked cache-parameters, --write-subsets-only and merge-parameter-subsets


def test_single_locus_chunked(single_locus_chunked_dir):
    check_subsets_complete(single_locus_chunked_dir)
    uids = sw_cache_uids(single_locus_chunked_dir)
    assert uids <= fasta_uids(UNPAIRED_IGH_FNAME)
    assert len(uids) >= MIN_ANNOTATED_FRAC * len(fasta_uids(UNPAIRED_IGH_FNAME))
    assert (single_locus_chunked_dir / 'hmm' / 'all-mean-mute-freqs.csv').is_file()


@pytest.mark.parametrize('locus', LOCI)
def test_unpaired_chunked(unpaired_chunked_dir, locus):
    check_subsets_complete(unpaired_chunked_dir)
    uids = sw_cache_uids(unpaired_chunked_dir / 'parameters' / locus)
    assert len(uids) > 0
    assert uids <= fasta_uids(os.path.join(PAIRED_DATA_DIR, '%s.fa' % locus))


def test_subset_job_procs_within_budget(unpaired_chunked_dir):
    match = re.search(r'running %d subset jobs \((\d+) concurrent, (\d+) procs each\)' % N_SUBSETS, read_log(fixture_log(unpaired_chunked_dir)))
    assert match is not None
    n_jobs, n_procs = int(match.group(1)), int(match.group(2))
    assert n_jobs * n_procs <= N_PROCS
    assert n_procs == N_PROCS  # per-job target is above the budget, so one job gets all of it


@pytest.mark.parametrize('form', ['paired', 'single'])
def test_write_subsets_only(tmp_path, form):
    if form == 'paired':
        args, basedir = ['--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--paired-outdir', str(tmp_path / 'out')], tmp_path / 'out'
        n_input = sum(len(input_uids(l)) for l in LOCI)
    else:
        args, basedir = ['--infname', UNPAIRED_IGH_FNAME, '--parameter-dir', str(tmp_path / 'params')], tmp_path / 'params'
        n_input = len(fasta_uids(UNPAIRED_IGH_FNAME))
    run_partis(['cache-parameters', '--n-subsets', str(N_SUBSETS), '--write-subsets-only'] + args, str(tmp_path / 'partis.log'))
    assert utils.read_subset_index(str(basedir))['n_subsets'] == N_SUBSETS
    n_written = 0
    for isub in range(N_SUBSETS):
        sdir = utils.parameter_subset_dir(str(basedir), isub)
        n_written += len(fasta_uids('%s/input-seqs.fa' % sdir))
        assert not utils.subset_is_marked_complete(sdir)
        assert not os.path.exists('%s/parameters' % sdir)
    assert n_written == n_input


def test_merge_action_paired(tmp_path, chunked_param_dir):
    # a copy whose subsets are done but never merged, as after externally run subset jobs
    basedir = copy_dir(chunked_param_dir, tmp_path / 'out')
    for locus in LOCI:
        shutil.rmtree(str(basedir / 'parameters' / locus))
    run_partis(merge_subsets_args(basedir), str(tmp_path / 'merge.log'))
    for locus in LOCI:
        assert sw_cache_uids(basedir / 'parameters' / locus) == input_uids(locus)
        assert os.path.isdir(basedir / 'parameters' / locus / 'hmm' / 'germline-sets')


def test_merge_action_single_locus(tmp_path, single_locus_chunked_dir):
    basedir = copy_dir(single_locus_chunked_dir, tmp_path / 'params')
    expected = sw_cache_uids(basedir)
    for path in basedir.iterdir():
        if path.name == 'parameter-subsets':
            continue
        if path.is_dir():
            shutil.rmtree(str(path))
        else:
            os.remove(str(path))
    run_partis(['merge-parameter-subsets', '--parameter-dir', str(basedir), '--locus', 'igh'], str(tmp_path / 'merge.log'))
    assert sw_cache_uids(basedir) == expected
    assert os.path.isdir(basedir / 'hmm' / 'germline-sets')


def test_merged_hmms_rebuilt_from_counts(chunked_param_dir):
    for locus in LOCI:
        merged = glob.glob('%s/parameters/%s/hmm/hmms/*.yaml' % (chunked_param_dir, locus))
        assert len(merged) > 0
        assert not any(os.path.islink(f) for f in merged)
        subset_genes = set()
        for isub in range(N_SUBSETS):
            subset_genes |= set(os.path.basename(f) for f in glob.glob('%s/parameters/%s/hmm/hmms/*.yaml' % (utils.parameter_subset_dir(str(chunked_param_dir), isub), locus)))
        assert set(os.path.basename(f) for f in merged) >= subset_genes


@pytest.mark.parametrize('to_index', [True, False])
def test_overwrite_remerge_switches_cache_form(tmp_path, chunked_param_dir, chunked_param_index_dir, to_index):
    basedir = copy_dir(chunked_param_dir if to_index else chunked_param_index_dir, tmp_path / 'out')
    run_partis(merge_subsets_args(basedir, extra_args=['--overwrite'] + (['--no-merged-sw-cache'] if to_index else [])), str(tmp_path / 'merge.log'))
    for locus in LOCI:
        pdir = basedir / 'parameters' / locus
        assert (pdir / SW_CACHE_INDEX_FNAME).exists() == to_index
        assert (pdir / SW_CACHE_FNAME).exists() != to_index
        if not to_index:
            assert sw_cache_uids(pdir) == input_uids(locus)


def test_symlinked_hmm_refused_on_remerge(tmp_path, chunked_param_dir):
    basedir = copy_dir(chunked_param_dir, tmp_path / 'out')
    merged = sorted(glob.glob('%s/parameters/igh/hmm/hmms/*.yaml' % basedir))[0]
    subset_model = '%s/parameters/igh/hmm/hmms/%s' % (utils.parameter_subset_dir(str(basedir), 0), os.path.basename(merged))
    os.remove(merged)
    os.symlink(subset_model, merged)
    log = run_partis_fails(merge_subsets_args(basedir), str(tmp_path / 'merge.log'))
    assert 'are symlinks from a pre-rebuild merge' in log


def test_index_size_warning(chunked_param_index_dir):
    log = read_log(fixture_log(chunked_param_index_dir))
    for locus in LOCI:
        match = re.search(r'warning %s: (\d+) sequences, small enough for a merged sw cache' % locus[2], log)  # loci print as their last letter
        assert match is not None
        assert int(match.group(1)) == len(input_uids(locus))


@pytest.mark.parametrize('action', ['annotate', 'partition'])
def test_index_only_param_dir_refused(tmp_path, chunked_param_index_dir, action):
    log = run_partis_fails([action, '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--parameter-dir', str(chunked_param_index_dir / 'parameters'), '--paired-outdir', str(tmp_path / 'out')], str(tmp_path / 'partis.log'))
    assert 'only the per-subset index %s' % SW_CACHE_INDEX_FNAME in log


@pytest.mark.parametrize('extra_args, errstr', [
    (['--n-subsets', '0'], '--n-subsets must be positive'),
    (['--write-subsets-only'], '--write-subsets-only requires --n-subsets'),
    (['--no-merged-sw-cache'], '--no-merged-sw-cache requires --n-subsets'),
])
def test_chunking_args_validated(tmp_path, extra_args, errstr):
    log = run_partis_fails(['cache-parameters', '--paired-loci', '--paired-indir', PAIRED_SIMU_DIR, '--paired-outdir', str(tmp_path / 'out')] + extra_args, str(tmp_path / 'partis.log'))
    assert errstr in log
