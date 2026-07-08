#!/usr/bin/env python3
"""Standalone CLI for post-partition refinement.

Thin wrapper around partis.refine (the pipeline-agnostic module). Reads a
partition file + SW cache from disk, runs refine_partition(), optionally scores
against ground truth, and writes the refined partition.

  python3 bin/partition-refinement.py \
      --partition-path <partition.yaml> \
      --sw-cache-path <sw-cache.yaml> \
      --simu-path <simu.yaml> \
      --step3-mode relative --ej-margin 0.05 --skip-singleton-merge \
      --min-agreement 0.15 --naive-freq-threshold 10 --rearrangement-guard \
      --cluster-size-csv <cluster_size.csv> --output-path <refined.yaml>
"""
import argparse
import csv
import json
import os
import sys

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from partis import refine  # noqa: E402

# unbuffered output for nohup/redirect progress visibility
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)


def main():
    parser = argparse.ArgumentParser(description='Post-partition refinement (CLI wrapper for partis.refine)')
    parser.add_argument('--partition-path', required=True)
    parser.add_argument('--sw-cache-path', required=True)
    parser.add_argument('--simu-path', default=None)
    parser.add_argument('--jaccard-threshold', type=float, default=None)
    parser.add_argument('--naive-threshold', type=float, default=None)
    parser.add_argument('--min-agreement', type=float, default=0.15,
                        help='min fingerprint agreement for step 2 merge (default 0.15)')
    parser.add_argument('--min-fp-positions', type=int, default=0,
                        help='min strong fingerprint positions for step 2 (0=disabled, try 3)')
    parser.add_argument('--skip-singleton-merge', action='store_true',
                        help='skip merge attempts between two singleton clusters')
    parser.add_argument('--step3-threshold', type=float, default=0.10)
    parser.add_argument('--step3-mode', choices=['fixed', 'gate', 'adaptive', 'relative'],
                        default='fixed',
                        help='Step 3 threshold mode: fixed (original), gate (skip low-signal), '
                             'adaptive (p75 noise floor), relative (median + margin)')
    parser.add_argument('--min-mutations', type=int, default=20,
                        help='For gate mode: skip groups with avg mutations below this')
    parser.add_argument('--ej-margin', type=float, default=0.05,
                        help='For relative mode: threshold = median_EJ + margin')
    parser.add_argument('--rescue-threshold', type=float, default=None,
                        help='EJ threshold for rescuing step 3 singletons (e.g. 0.03)')
    parser.add_argument('--naive-freq-threshold', type=int, default=10)
    parser.add_argument('--max-expansion-ratio', type=float, default=2.0,
                        help='skip splitting if cluster_size/naive_count <= this (expansion guard)')
    parser.add_argument('--output-path', default=None,
                        help='write refined partition to this path (JSON)')
    parser.add_argument('--cluster-size-csv', default=None,
                        help='path to cluster_size.csv for max family size (two-level guard)')
    parser.add_argument('--hmm-annotation-path', default=None,
                        help='HMM annotation file for naives (strips fv_insertion, falls back to SW on length mismatch)')
    parser.add_argument('--min-cluster-size', type=int, default=2)
    parser.add_argument('--rearrangement-guard', action='store_true',
                        help='use V+D+J majority guard for IGH deep-tree protection')
    args = parser.parse_args()

    print('reading partition %s and sw-cache %s' % (args.partition_path, args.sw_cache_path))
    inp = refine.read_refine_inputs(args.partition_path, args.sw_cache_path)
    partition = inp['partition']
    print('  %d clusters, %d seqs, %d per-sequence sw annotations' % (
        len(partition), sum(len(c) for c in partition), len(inp['sw_info'])))

    uid_sw_naives = inp['uid_sw_naives']
    if args.hmm_annotation_path is not None:
        print('loading HMM naives from %s' % args.hmm_annotation_path)
        uid_sw_naives = refine.load_hmm_naives(args.hmm_annotation_path, sw_naives=uid_sw_naives)

    true_partition = None
    if args.simu_path is not None:
        true_partition = refine.load_true_partition(args.simu_path)
        pur, comp = refine.calc_metrics(true_partition, partition)
        print('\nbaseline: purity %.4f, completeness %.4f, %d clusters' % (pur, comp, len(partition)))

    max_family_size = None
    if args.cluster_size_csv is not None:
        with open(args.cluster_size_csv) as f:
            max_family_size = max(int(row['cluster_size']) for row in csv.DictReader(f))
        print('  max family size from %s: %d' % (args.cluster_size_csv, max_family_size))

    final_partition = refine.refine_partition(
        partition, inp['uid_info'], uid_sw_naives,
        uid_rearr_features=(inp['uid_rearr_features'] if args.rearrangement_guard else None),
        max_family_size=max_family_size,
        jaccard_threshold=args.jaccard_threshold, naive_threshold=args.naive_threshold,
        min_agreement=args.min_agreement, min_fp_positions=args.min_fp_positions,
        skip_singleton_merge=args.skip_singleton_merge, step3_threshold=args.step3_threshold,
        step3_mode=args.step3_mode, min_mutations=args.min_mutations, ej_margin=args.ej_margin,
        rescue_threshold=args.rescue_threshold, naive_freq_threshold=args.naive_freq_threshold,
        max_expansion_ratio=args.max_expansion_ratio, min_cluster_size=args.min_cluster_size,
        rearrangement_guard=args.rearrangement_guard)

    if true_partition is not None:
        pur_f, comp_f = refine.calc_metrics(true_partition, final_partition)
        print('\n=== Summary ===')
        print('  %-25s %8.4f %12.4f %8d' % ('baseline', pur, comp, len(partition)))
        print('  %-25s %8.4f %12.4f %8d' % ('refined', pur_f, comp_f, len(final_partition)))

    if args.output_path is not None:
        n = refine.write_full_output(args.output_path, inp['glfo'], final_partition, inp['sw_info'])
        print('\nwrote refined partition (full partis output) to %s (%d clusters)' % (args.output_path, n))


if __name__ == '__main__':
    main()
