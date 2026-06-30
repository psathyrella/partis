#!/usr/bin/env python3
"""Diff two partis paired-loci output dirs for byte-equality on partition + annotations.

Usage: zig-compare.py <baseline_dir> <new_dir>

Companion to bin/partis-30k-parity-test.sh (issue #375 cpp-mirror gate).
Exit 0 if identical, 1 otherwise.
"""
import hashlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from partis import utils
from partis.clusterpath import ClusterPath


def best_partition(path):
    cpath = ClusterPath()
    cpath.readfile(path)
    return sorted(tuple(sorted(c)) for c in cpath.partitions[cpath.i_best])


FIELDS = ('naive_seq', 'v_gene', 'd_gene', 'j_gene',
          'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del',
          'vd_insertion', 'dj_insertion', 'cdr3_length')

_MD5_CHUNK_BYTES = 1 << 20


# All partition-{locus}.yaml files paired-loci output writes. Subdirs and loci
# are derived from utils.locus_pairs (plus top-level + single-chain/) so a
# future paired-loci pair addition picks up automatically. Originally the
# gate only walked igh+igk/ (or igh+igl/ as fallback) and missed single-chain
# divergence at 50k — see #386.
def _partition_relpaths():
    pairs = utils.locus_pairs['ig']
    loci = sorted({l for lp in pairs for l in lp})
    subs = [''] + ['+'.join(lp) for lp in pairs] + ['single-chain']
    return [(sub, loc, f'{sub}/partition-{loc}.yaml' if sub else f'partition-{loc}.yaml')
            for sub in subs for loc in loci]


PARTITION_RELPATHS = _partition_relpaths()


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(_MD5_CHUNK_BYTES), b''):
            h.update(chunk)
    return h.hexdigest()


def diff_one(label, bp, np_):
    """Compare a single base/new partition-*.yaml pair. Returns 1 on diff, 0 on match, None on skip."""
    if not os.path.exists(bp) and not os.path.exists(np_):
        return None  # neither side has this file — not produced by this run
    if not os.path.exists(bp):
        print(f'{label}: MISSING in base ({bp})')
        return 1
    if not os.path.exists(np_):
        print(f'{label}: MISSING in new ({np_})')
        return 1
    if md5(bp) == md5(np_):
        print(f'{label}: byte-identical')
        return 0
    bpart = best_partition(bp)
    npart = best_partition(np_)
    if bpart != npart:
        print(f'{label}: PARTITION DIFF ({len(bpart)} vs {len(npart)} clusters)')
        return 1
    bd = utils.read_json_yaml(bp)
    nd = utils.read_json_yaml(np_)
    bk = {tuple(sorted(e['unique_ids'])): e for e in bd['events']}
    nk = {tuple(sorted(e['unique_ids'])): e for e in nd['events']}
    if set(bk) != set(nk):
        only_b = set(bk) - set(nk)
        only_n = set(nk) - set(bk)
        print(f'{label}: EVENT-KEY DIFF (only in base: {len(only_b)}, only in new: {len(only_n)})')
        return 1
    field_diffs = 0
    for k in bk:
        for f in FIELDS:
            if bk[k].get(f) != nk[k].get(f):
                field_diffs += 1
                if field_diffs <= 3:
                    print(f'{label} {k[0]}...: {f} diff: {bk[k].get(f)} vs {nk[k].get(f)}')
    if field_diffs:
        print(f'{label}: {field_diffs} FIELD DIFFS across {len(bk)} events')
        return 1
    # Bytes differ but every event keyed by unique_ids has bit-identical content
    # and partitions match — the diff is event-list ordering only (#386).
    cpp_order = [tuple(sorted(e['unique_ids'])) for e in bd['events']]
    zig_order = [tuple(sorted(e['unique_ids'])) for e in nd['events']]
    n_pos_diff = sum(1 for a, b in zip(cpp_order, zig_order) if a != b)
    print(f'{label}: BYTE DIFF (md5 mismatch) but events bit-identical by uid; '
          f'event-list order differs at {n_pos_diff}/{len(bk)} positions')
    return 1


def diff_pair(base, new):
    n_diffs = 0
    n_compared = 0
    for sub, locus, rel in PARTITION_RELPATHS:
        label = f'{sub or "(top)"}/{locus}'
        res = diff_one(label, f'{base}/{rel}', f'{new}/{rel}')
        if res is None:
            continue
        n_compared += 1
        n_diffs += res
    if n_compared == 0:
        print('no partition files found in either dir')
        return 1
    return n_diffs


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    base, new = sys.argv[1], sys.argv[2]
    n = diff_pair(base, new)
    if n:
        print(f'\nFAILED: {n} loci differ')
        sys.exit(1)
    print('\nOK: identical')


if __name__ == '__main__':
    main()
