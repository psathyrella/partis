import os

import pytest

from _helpers import N_SUBSETS, PAIRED_DATA_DIR, PAIRED_SIMU_DIR, copy_dir, edit_yaml, merge_subsets_args, paired_chunked_args, read_log, run_partis, run_partis_fails
from partis import utils
from partis import disjointgrouper as dg

# ----------------------------------------------------------------------------------------
# multi-cache grouping goes through the hash-checked sw cache index


def igh_pdir(basedir):
    return basedir / 'parameters' / 'igh'


def subset_caches(basedir):
    return dg.resolve_sw_cache_paths(str(igh_pdir(basedir)))[0]


def group_from(tmp_path, sw_cachefname):
    return run_partis_fails(['create-disjoint-groups', '--locus', 'igh', '--sw-cachefname', str(sw_cachefname), '--paired-outdir', str(tmp_path / 'groups')], str(tmp_path / 'create.log'))


def test_dir_without_index_refused(tmp_path, chunked_param_index_dir):
    basedir = copy_dir(chunked_param_index_dir, tmp_path / 'out')
    os.remove(dg.sw_cache_index_fname(str(igh_pdir(basedir))))
    log = group_from(tmp_path, basedir)  # holds parameter-subsets/subset-*/parameters/igh/sw-cache.yaml
    assert '%s is a dir with no %s' % (basedir, dg.SW_CACHE_INDEX_FNAME) in log


def test_colon_list_refused(tmp_path, chunked_param_index_dir):
    log = group_from(tmp_path, ':'.join(subset_caches(chunked_param_index_dir)))
    assert 'colon-separated lists are not accepted' in log


def test_one_subset_index_hash_checked(tmp_path, chunked_param_index_dir):
    basedir = copy_dir(chunked_param_index_dir, tmp_path / 'out')
    def keep_first(index):  # the shape --n-subsets 1 --no-merged-sw-cache writes
        index['n_subsets'], index['sw_caches'] = 1, index['sw_caches'][:1]
    edit_yaml(dg.sw_cache_index_fname(str(igh_pdir(basedir))), keep_first)
    swfn = subset_caches(basedir)[0]
    assert not os.path.islink(swfn)
    with open(swfn, 'a') as swfile:
        swfile.write('\n')
    log = group_from(tmp_path, igh_pdir(basedir))
    assert 'sw cache %s does not match the index' % swfn in log


# ----------------------------------------------------------------------------------------
# subset-index.yaml: hashes required, checked at merge, input recorded and compared on re-run


def merge_fails(tmp_path, basedir):
    return run_partis_fails(merge_subsets_args(basedir), str(tmp_path / 'merge.log'))


def test_subset_index_without_hashes_refused(tmp_path, chunked_param_dir):
    basedir = copy_dir(chunked_param_dir, tmp_path / 'out')
    def drop_hashes(index):
        for sfo in index['subsets']:
            del sfo['xxh3']
    edit_yaml(utils.subset_index_fname(str(basedir)), drop_hashes)
    assert '0 of %d subsets have an xxh3 hash, expected all of them' % N_SUBSETS in merge_fails(tmp_path, basedir)


def test_changed_subset_input_refused_at_merge(tmp_path, chunked_param_dir):
    basedir = copy_dir(chunked_param_dir, tmp_path / 'out')
    infname = '%s/%s' % (utils.parameter_subset_dir(str(basedir), 0), utils.SUBSET_INPUT_FNAME)
    with open(infname, 'a') as infile:
        infile.write('>extra-seq\nACGT\n')
    assert 'subset input %s does not match the index' % infname in merge_fails(tmp_path, basedir)


def split_args(indir, outdir, extra_args=None):
    return paired_chunked_args(outdir, indir=indir, extra_args=['--write-subsets-only'] + (extra_args or []))


def change_input(fname):
    if fname.suffix == '.fa':
        utils.write_fasta(str(fname), utils.read_fastx(str(fname))[:-1])
    else:
        edit_yaml(fname, lambda meta: meta.pop(sorted(meta)[0]))


@pytest.mark.parametrize('layout, changed_fname', [('simu', 'all-seqs.fa'), ('simu', 'meta.yaml'), ('data', 'igh.fa')])
def test_changed_input_refused_on_rerun(tmp_path, layout, changed_fname):
    indir = copy_dir(PAIRED_SIMU_DIR if layout == 'simu' else PAIRED_DATA_DIR, tmp_path / 'indir')
    if layout == 'data':
        os.remove(str(indir / 'all-seqs.fa'))  # so the split concatenates the per-locus files
    outdir = tmp_path / 'out'
    run_partis(split_args(indir, outdir), str(tmp_path / 'split.log'))
    change_input(indir / changed_fname)
    log = run_partis_fails(split_args(indir, outdir), str(tmp_path / 'rerun.log'))
    assert 'subset index %s does not match this run\'s inputs' % utils.subset_index_fname(str(outdir)) in log
    assert str(indir / changed_fname) in log


@pytest.mark.parametrize('extra_args, errstr', [(['--n-max-queries', '50'], '--n-max-queries 50 but index has -1'),
                                               (['--queries-to-include', 'a:b'], "--queries-to-include ['a', 'b'] but index has None")])
def test_changed_split_arg_refused_on_rerun(tmp_path, extra_args, errstr):
    outdir = tmp_path / 'out'
    run_partis(split_args(PAIRED_SIMU_DIR, outdir), str(tmp_path / 'split.log'))
    log = run_partis_fails(split_args(PAIRED_SIMU_DIR, outdir, extra_args=extra_args), str(tmp_path / 'rerun.log'))
    assert errstr in log


def test_unchanged_rerun_accepted(tmp_path):
    outdir = tmp_path / 'out'
    run_partis(split_args(PAIRED_SIMU_DIR, outdir), str(tmp_path / 'split.log'))
    run_partis(split_args(PAIRED_SIMU_DIR, outdir), str(tmp_path / 'rerun.log'))
    assert 'subset input files exist' in read_log(tmp_path / 'rerun.log')


def test_index_without_split_inputs_refused(tmp_path):
    outdir = tmp_path / 'out'
    run_partis(split_args(PAIRED_SIMU_DIR, outdir), str(tmp_path / 'split.log'))
    edit_yaml(utils.subset_index_fname(str(outdir)), lambda index: index.pop('split_inputs'))
    log = run_partis_fails(split_args(PAIRED_SIMU_DIR, outdir), str(tmp_path / 'rerun.log'))
    assert 'is missing required key \'split_inputs\'' in log
