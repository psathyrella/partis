#!/usr/bin/env python3
"""Aggregate PERFCOUNTER + PERFCOUNTER_CACHE lines from a perf log file.

Per-query summary: median/p95/max of ksets, fillTrellis calls, addInLogSpace
calls, addInLogSpace+log_exp work calls, active-state distribution stats, and
chunk-cache hit rate. Plus a chunk-cache prefix-list-length histogram across
all queries.

The log file is produced by partis-zig-core when built with
-Dperf-counters=true and run with PARTIS_ZIG_PERF_LOG=<path> in the env. See
packages/zig-core/src/ham/perf_counters.zig for the metric definitions.

These counters pin denominators (active-state distribution, n_ksets,
fillTrellis call counts, addInLogSpace call counts, chunk-cache hit rate)
for the perf-work cost claims tracked in issue #366
(https://github.com/psathyrella/partis/issues/366). Originally landed in
PR #368.

Usage: zig-perf-counters-aggregate.py <log-file>
"""
import re
import sys
from collections import defaultdict


def parse_kv(line):
    return {m.group(1): m.group(2) for m in re.finditer(r"(\w+)=(\S+)", line)}


def stats(xs):
    if not xs:
        return None
    xs = sorted(xs)
    n = len(xs)
    return {
        "n": n,
        "min": xs[0],
        "median": xs[n // 2],
        "p95": xs[int(n * 0.95)] if n > 1 else xs[0],
        "max": xs[-1],
        "mean": sum(xs) / n,
        "total": sum(xs),
    }


def main(path):
    by_alg = defaultdict(lambda: {
        "ksets": [], "fillTrellis": [], "addInLogSpace": [],
        "addInLogSpace_log_exp": [], "chunk_hits": [], "chunk_lookups": [],
        "active_med": [], "active_p95": [], "active_max": [], "active_positions": [],
        "n_seqs": [],
    })
    cache_buckets = defaultdict(lambda: defaultdict(int))
    n_lines = 0
    n_cache_lines = 0

    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith("PERFCOUNTER_CACHE"):
            n_cache_lines += 1
            kv = parse_kv(line[len("PERFCOUNTER_CACHE "):])
            alg = kv.get("alg", "?")
            for k, v in kv.items():
                if k in ("alg", "q"):
                    continue
                cache_buckets[alg][k] += int(v)
            continue
        if not line.startswith("PERFCOUNTER"):
            continue
        n_lines += 1
        kv = parse_kv(line[len("PERFCOUNTER "):])
        alg = kv.get("alg", "?")
        d = by_alg[alg]
        d["n_seqs"].append(int(kv["n_seqs"]))
        d["ksets"].append(int(kv["ksets"]))
        d["fillTrellis"].append(int(kv["fillTrellis"]))
        d["addInLogSpace"].append(int(kv["addInLogSpace"]))
        d["addInLogSpace_log_exp"].append(int(kv["addInLogSpace_log_exp"]))
        h, l = kv["chunk_hits"].split("/")
        d["chunk_hits"].append(int(h))
        d["chunk_lookups"].append(int(l))
        d["active_med"].append(int(kv["active_states_med"]))
        d["active_p95"].append(int(kv["active_states_p95"]))
        d["active_max"].append(int(kv["active_states_max"]))
        d["active_positions"].append(int(kv["active_states_positions"]))

    print(f"# parsed {n_lines} PERFCOUNTER + {n_cache_lines} PERFCOUNTER_CACHE lines from {path}")
    for alg, d in by_alg.items():
        nq = len(d["ksets"])
        print(f"\n## algorithm={alg}, n_queries={nq}")

        def row(name, xs, fmt="{:>12,d}"):
            s = stats(xs)
            if s is None:
                return
            print(f"  {name:<32}  median={fmt.format(s['median'])}  p95={fmt.format(s['p95'])}  max={fmt.format(s['max'])}  mean={s['mean']:>12,.1f}  total={s['total']:>14,d}")

        row("n_seqs (per query)", d["n_seqs"])
        row("ksets (per query)", d["ksets"])
        row("fillTrellis calls (per query)", d["fillTrellis"])
        row("addInLogSpace calls (per query)", d["addInLogSpace"])
        row("addInLogSpace+log_exp (per query)", d["addInLogSpace_log_exp"])

        # active states are already per-query medians/p95/maxes; aggregate across queries
        row("active_states median-of-medians", d["active_med"])
        row("active_states median-of-p95s", d["active_p95"])
        row("active_states median-of-maxes", d["active_max"])
        row("active_positions (per query)", d["active_positions"])

        # chunk-cache rate
        total_hits = sum(d["chunk_hits"])
        total_lookups = sum(d["chunk_lookups"])
        if total_lookups:
            print(f"  {'chunk_cache_rate':<32}  hits/lookups = {total_hits:,} / {total_lookups:,} = {100*total_hits/total_lookups:.1f}%")

    if cache_buckets:
        print("\n## chunk-cache prefix-list-length histogram (across all queries)")
        for alg, buckets in cache_buckets.items():
            print(f"  alg={alg}")
            order = ["0", "1", "2", "3", "4", "5", "6_9", "10_19", "20_49", "50_99", "100_199", "200p"]
            total = sum(buckets.values())
            for k in order:
                v = buckets.get(k, 0)
                pct = 100 * v / total if total else 0
                print(f"    {k:>8}: {v:>10,d}  ({pct:>5.1f}%)")
            print(f"    {'total':>8}: {total:>10,d}")


if __name__ == "__main__":
    main(sys.argv[1])
