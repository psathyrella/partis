import collections
import os

import pytest

from _helpers import HFRAC_ARGS, PAIRED_PARAM_DIR, copy_dir, disjoint_dir, edit_yaml, fasta_uids, input_uids, locus_manifest, output_partition_uids, partition_uids, read_yaml, run_partis
from partis import disjointgrouper as dg

# ----------------------------------------------------------------------------------------
# create-disjoint-groups and assemble-groups, run standalone


def create_groups(outdir, locus, pdir, extra_args=None):
    logfname = '%s/create-%s.log' % (outdir, locus)
    run_partis(['create-disjoint-groups', '--locus', locus, '--parameter-dir', str(pdir), '--paired-outdir', str(outdir)] + (extra_args or []), logfname)
    return locus_manifest(outdir, locus)[1]


def grouped_uids(outdir, locus, manifest):
    ddir = disjoint_dir(outdir, locus)
    return set().union(*[fasta_uids('%s/%s' % (ddir, g['fasta_path'])) for g in manifest['groups']])


@pytest.mark.parametrize('pdir_level', ['parent', 'locus'])
def test_create_groups_plain(tmp_path, pdir_level):
    pdir = PAIRED_PARAM_DIR if pdir_level == 'parent' else os.path.join(PAIRED_PARAM_DIR, 'igh')
    manifest = create_groups(tmp_path, 'igh', pdir)
    assert manifest['grouping-info']['method'] == 'cdr3-length'
    assert all('/sub-groups/' not in g['fasta_path'] for g in manifest['groups'])
    assert len(set(g['cdr3_length'] for g in manifest['groups'])) == len(manifest['groups'])  # one group per cdr3 length
    assert grouped_uids(tmp_path, 'igh', manifest) == input_uids('igh')


def test_create_groups_hfrac(tmp_path):
    manifest = create_groups(tmp_path, 'igh', PAIRED_PARAM_DIR, extra_args=HFRAC_ARGS)
    assert manifest['grouping-info']['method'] == 'cdr3-length+hfrac'
    assert all('/sub-groups/sub-' in g['fasta_path'] for g in manifest['groups'])
    assert max(collections.Counter(g['cdr3_length'] for g in manifest['groups']).values()) > 1
    assert grouped_uids(tmp_path, 'igh', manifest) == input_uids('igh')


def test_create_groups_from_chunked_merged_cache(tmp_path, chunked_param_dir):
    manifest = create_groups(tmp_path, 'igh', chunked_param_dir / 'parameters')
    assert grouped_uids(tmp_path, 'igh', manifest) == input_uids('igh')


@pytest.mark.parametrize('use_hfrac', [False, True])
def test_create_groups_from_index(tmp_path, chunked_param_index_dir, use_hfrac):
    manifest = create_groups(tmp_path, 'igh', chunked_param_index_dir / 'parameters', extra_args=HFRAC_ARGS if use_hfrac else None)
    assert grouped_uids(tmp_path, 'igh', manifest) == input_uids('igh')
    assert manifest['grouping-info']['method'] == ('cdr3-length+hfrac' if use_hfrac else 'cdr3-length')


# ----------------------------------------------------------------------------------------
# assemble-groups on a copy of an integrated run's disjoint dir

def assemble(outdir, outfname, extra_args=None, locus='igh'):
    run_partis(['assemble-groups', '--locus', locus, '--paired-outdir', str(outdir), '--outfname', str(outfname)] + (extra_args or []), '%s.log' % outfname)


def test_assemble_single_file(tmp_path, hfrac_partition_dir):
    outdir = copy_dir(hfrac_partition_dir, tmp_path / 'out')
    outfname = tmp_path / 'assembled-igh.yaml'
    assemble(outdir, outfname)
    assert partition_uids(outfname) == input_uids('igh')
    _, manifest = locus_manifest(outdir, 'igh')
    assert manifest['assembly']['status'] == 'merged'
    assert manifest['assembly']['merged_output_path'] == str(outfname)
    assert all(manifest['assembly']['validation'].values())


def test_assemble_multifile(tmp_path, hfrac_partition_dir):
    outdir = copy_dir(hfrac_partition_dir, tmp_path / 'out')
    outfname = tmp_path / 'assembled-igh.yaml'
    assemble(outdir, outfname, extra_args=['--multifile-min-seqs', '1', '--multifile-max-seqs-per-file', '3'])
    assert not outfname.exists()
    mfdir = dg.multifile_dir_path(str(outfname))
    index = dg.read_multifile_index('%s/%s' % (mfdir, dg.MULTIFILE_INDEX_FNAME))
    files_per_cdr3 = collections.Counter(f['cdr3_length'] for f in index['files'])
    assert max(files_per_cdr3.values()) > 1  # sub-groups of one cdr3 length packed into more than one file
    assert output_partition_uids(outfname) == input_uids('igh')  # readers resolve the multifile dir from the outfname


def test_assemble_rerun_is_idempotent(tmp_path, hfrac_partition_dir):
    outdir = copy_dir(hfrac_partition_dir, tmp_path / 'out')
    outfname = tmp_path / 'assembled-igh.yaml'
    assemble(outdir, outfname)
    first = outfname.read_text()
    assemble(outdir, outfname)
    assert outfname.read_text() == first


# ----------------------------------------------------------------------------------------
# assemble error paths, in-process on a mutated copy


def group_ppath(ddir, ginfo):
    return '%s/%s' % (ddir, dg.resolve_partition_path(ginfo, ddir)[0])


def drop_manifest_tail(ddir, manifest):
    manifest['groups'] = manifest['groups'][:-1]


def remove_group_partition(ddir, manifest):
    os.remove(group_ppath(ddir, manifest['groups'][0]))


def empty_group_partition(ddir, manifest):
    open(group_ppath(ddir, manifest['groups'][0]), 'w').close()


def rewrite_group(ddir, ginfo, fcn):
    edit_yaml(group_ppath(ddir, ginfo), fcn)


def read_group(ddir, ginfo):
    return read_yaml(group_ppath(ddir, ginfo))


def change_a_gene_seq(ddir, manifest):
    g0, g1 = [read_group(ddir, g)['germline-info']['seqs']['v'] for g in manifest['groups'][:2]]
    gene = sorted(set(g0) & set(g1))[0]
    def fcn(yinfo):
        seq = yinfo['germline-info']['seqs']['v'][gene]
        yinfo['germline-info']['seqs']['v'][gene] = ('T' if seq[0] != 'T' else 'A') + seq[1:]
    rewrite_group(ddir, manifest['groups'][1], fcn)


def duplicate_a_uid(ddir, manifest):
    # rename a uid in one group to a uid from another, in both its annotation and its partition
    uid = read_group(ddir, manifest['groups'][0])['events'][0]['unique_ids'][0]
    def fcn(yinfo):
        old = yinfo['events'][0]['unique_ids'][0]
        for event in yinfo['events']:
            event['unique_ids'] = [uid if u == old else u for u in event['unique_ids']]
        for pinfo in yinfo['partitions']:
            pinfo['partition'] = [[uid if u == old else u for u in c] for c in pinfo['partition']]
    rewrite_group(ddir, manifest['groups'][1], fcn)


def lower_manifest_total(ddir, manifest):
    manifest['grouping-info']['total_grouped_sequences'] -= 1
    manifest['grouping-info']['total_input_sequences'] -= 1
    manifest['groups'][0]['sequence_count'] -= 1


ASSEMBLE_ERRORS = [
    ('missing-manifest', None, 'manifest file does not exist'),
    ('truncated-manifest', drop_manifest_tail, 'sequence count mismatch'),
    ('group-without-partition', remove_group_partition, 'have a manifest partition_path whose file is gone'),
    ('empty-partition', empty_group_partition, 'partition file is empty'),
    ('gene-with-two-seqs', change_a_gene_seq, 'has different sequences in groups'),
    ('duplicate-uid', duplicate_a_uid, 'duplicate uid'),
    ('count-above-manifest', lower_manifest_total, 'sequence count exceeds expected'),
]


@pytest.mark.parametrize('name, mutate, errstr', ASSEMBLE_ERRORS, ids=[e[0] for e in ASSEMBLE_ERRORS])
def test_assemble_error_paths(tmp_path, plain_multifile_partition_dir, name, mutate, errstr):
    outdir = copy_dir(plain_multifile_partition_dir, tmp_path / 'out')
    ddir = disjoint_dir(outdir, 'igh')
    mfname = '%s/%s' % (ddir, dg.MANIFEST_FNAME)
    if mutate is None:
        os.remove(mfname)
    else:
        edit_yaml(mfname, lambda manifest: mutate(ddir, manifest))  # unvalidated, since some mutations break the manifest on purpose
    with pytest.raises(Exception, match=errstr):
        dg.assemble_groups('igh', ddir, str(tmp_path / 'assembled-igh.yaml'))


# ----------------------------------------------------------------------------------------
# pure helpers

def test_uid_group_mapping_covers_every_grouped_uid(hfrac_partition_dir):
    ddir, manifest = locus_manifest(hfrac_partition_dir, 'igh')
    mapping = dg.build_uid_group_mapping(manifest, ddir)
    assert set(mapping) == input_uids('igh')
    assert set(mapping.values()) == set(g['group_id'] for g in manifest['groups'])


def test_group_by_cdr3_length():
    lines = [{'unique_ids': ['a', 'b'], 'input_seqs': ['AA', 'CC'], 'cdr3_length': 45, 'naive_seq': 'NN'},
             {'unique_ids': ['c'], 'input_seqs': ['GG'], 'cdr3_length': 33, 'naive_seq': 'MM'},
             {'unique_ids': ['d', 'e'], 'input_seqs': ['TT', 'TT'], 'cdr3_length': None}]
    groups, n_failed = dg.group_sequences_by_cdr3_length(lines)
    assert n_failed == 2
    assert list(groups) == [33, 45]
    assert [s['name'] for s in groups[45]] == ['a', 'b']
    assert groups[33][0]['naive_seq'] == 'MM'


def manifest_groups(counts):
    return [{'group_id': i, 'cdr3_length': 30 + i, 'locus': 'igh', 'sequence_count': n, 'fasta_path': 'groups/cdr3-%d/igh.fa' % (30 + i)} for i, n in enumerate(counts)]


def test_manifest_round_trip(tmp_path):
    written = dg.write_manifest(manifest_groups([3, 4]), str(tmp_path), 'igh', 9, 2, parameter_dir='x', hfrac=True)
    read = dg.read_manifest(str(tmp_path / dg.MANIFEST_FNAME))
    assert read == written
    assert read['grouping-info']['method'] == 'cdr3-length+hfrac'


@pytest.mark.parametrize('total_grouped, total_input, errstr', [(8, 9, 'sum of the 2 group counts'), (7, 10, 'does not equal total_input')])
def test_validate_sequence_count_raises(tmp_path, total_grouped, total_input, errstr):
    manifest = dg.write_manifest(manifest_groups([3, 4]), str(tmp_path), 'igh', 9, 2)
    manifest['grouping-info'].update({'total_grouped_sequences': total_grouped, 'total_input_sequences': total_input})
    with pytest.raises(Exception, match=errstr):
        dg.validate_sequence_count(manifest)


def test_read_manifest_missing_group_key(tmp_path):
    groups = manifest_groups([3])
    del groups[0]['fasta_path']
    dg.write_manifest(groups, str(tmp_path), 'igh', 3, 0)
    with pytest.raises(Exception, match='missing required key \'fasta_path\''):
        dg.read_manifest(str(tmp_path / dg.MANIFEST_FNAME))


def test_stage_names():
    assert dg.stage_fname(dg.STAGE_REFINE, 'igk') == 'partition-refine-igk.yaml'
    assert dg.stage_from_path('a/b/partition-refine-igh.yaml') == dg.STAGE_REFINE
    assert dg.stage_from_path('a/b/ha-repartition-igh.yaml') == dg.STAGE_HAREP
    assert dg.stage_from_path('a/b/partition-igh.yaml') == dg.STAGE_VSEARCH
    assert dg.stage_from_path('a/b/other.yaml') is None


def test_pack_multifile_output_never_spans_cdr3_groups():
    gpaths = [({'group_id': i, 'cdr3_length': c3}, 'p%d' % i) for i, c3 in enumerate([30, 30, 30, 42, 42])]
    counts = {i: {'sequence_count': n} for i, n in enumerate([2, 2, 5, 1, 1])}
    fspecs = dg.pack_multifile_output(gpaths, counts, max_seqs_per_file=4)
    assert [(c3, [p for _, p in glist]) for c3, glist in fspecs] == [(30, ['p0', 'p1']), (30, ['p2']), (42, ['p3', 'p4'])]


def test_resolve_sw_cache_paths(tmp_path):
    swfn = tmp_path / dg.SW_CACHE_FNAME
    swfn.write_text('')
    assert dg.resolve_sw_cache_paths(str(swfn)) == ([str(swfn)], [None])


def test_read_vsearch_uc_with_centroids(tmp_path):
    rows = [['S', '0', '10', '*', '*', '*', '*', '*', 'a', '*'],
            ['H', '0', '10', '99', '+', '0', '0', '=', 'b', 'a'],
            ['S', '1', '10', '*', '*', '*', '*', '*', 'c', '*'],
            ['C', '0', '2', '*', '*', '*', '*', '*', 'a', '*']]
    ucfname = tmp_path / 'clusters.uc'
    ucfname.write_text(''.join('\t'.join(r) + '\n' for r in rows))
    assert dg._read_vsearch_uc_with_centroids(str(ucfname)) == [('a', ['a', 'b']), ('c', ['c'])]


def test_merge_r1_subgroups_by_components():
    r1 = [('a', ['a', 'a2']), ('b', ['b']), ('c', ['c', 'c2', 'c3', 'c4', 'c5', 'c6']), ('d', ['d'])]
    components = [['a', 'b'], ['c'], ['d']]
    assert dg._merge_r1_subgroups_by_components(r1, components, max_bin_size=0) == [['a', 'a2', 'b'], ['c', 'c2', 'c3', 'c4', 'c5', 'c6'], ['d']]
    # big component above the cap stays whole, small ones are bundled up to the cap
    bins = dg._merge_r1_subgroups_by_components(r1, components, max_bin_size=4, min_group_size=5)
    assert sorted(map(sorted, bins)) == sorted([sorted(['c', 'c2', 'c3', 'c4', 'c5', 'c6']), ['a', 'a2', 'b', 'd']])
