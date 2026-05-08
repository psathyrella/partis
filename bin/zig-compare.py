#!/usr/bin/env python3
"""Diff two partis paired-loci output dirs for byte-equality on partition + annotations.

Usage: zig-compare.py <baseline_dir> <new_dir>

Per issue #342 validation harness. Exit 0 if identical, 1 otherwise.
"""
import json
import os
import sys

sys.path.insert(0, '/home/matsen/re/partis')
from partis.clusterpath import ClusterPath


def best_partition(path):
    cpath = ClusterPath()
    cpath.readfile(path)
    return sorted(tuple(sorted(c)) for c in cpath.partitions[cpath.i_best])


FIELDS = ('naive_seq', 'v_gene', 'd_gene', 'j_gene',
          'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del',
          'vd_insertion', 'dj_insertion', 'cdr3_length')


def diff_pair(base, new):
    n_diffs = 0
    for locus in ('igh', 'igk', 'igl'):
        bp = f'{base}/igh+igk/partition-{locus}.yaml'
        np_ = f'{new}/igh+igk/partition-{locus}.yaml'
        if not os.path.exists(bp):
            bp = f'{base}/igh+igl/partition-{locus}.yaml'
            np_ = f'{new}/igh+igl/partition-{locus}.yaml'
        if not os.path.exists(bp):
            print(f'{locus}: no baseline file (skipping)')
            continue
        if not os.path.exists(np_):
            print(f'{locus}: MISSING in new ({np_})')
            n_diffs += 1
            continue
        bpart = best_partition(bp)
        npart = best_partition(np_)
        if bpart != npart:
            print(f'{locus}: PARTITION DIFF ({len(bpart)} vs {len(npart)} clusters)')
            n_diffs += 1
            continue
        bd = json.load(open(bp))
        nd = json.load(open(np_))
        bk = {tuple(sorted(e['unique_ids'])): e for e in bd['events']}
        nk = {tuple(sorted(e['unique_ids'])): e for e in nd['events']}
        if set(bk) != set(nk):
            only_b = set(bk) - set(nk)
            only_n = set(nk) - set(bk)
            print(f'{locus}: EVENT-KEY DIFF (only in base: {len(only_b)}, only in new: {len(only_n)})')
            n_diffs += 1
            continue
        field_diffs = 0
        for k in bk:
            for f in FIELDS:
                if bk[k].get(f) != nk[k].get(f):
                    field_diffs += 1
                    if field_diffs <= 3:
                        print(f'{locus} {k[:1]}...: {f} diff: {bk[k].get(f)} vs {nk[k].get(f)}')
        if field_diffs:
            print(f'{locus}: {field_diffs} FIELD DIFFS across {len(bk)} events')
            n_diffs += 1
        else:
            print(f'{locus}: {len(bk)} events identical')
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
