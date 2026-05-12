#!/usr/bin/env python3
"""Aggregate PERFCOUNTER / PERFCOUNTER_CACHE / PERFCOUNTER_TIMING lines.

Per-query summary: median/p95/max of ksets, fillTrellis calls, addInLogSpace
calls, addInLogSpace+log_exp work calls, active-state distribution stats, and
chunk-cache hit rate. Plus a chunk-cache prefix-list-length histogram across
all queries.

For Phase-B (added after the issue #366 lever-5 close, PR
https://github.com/psathyrella/partis/pull/383), also aggregates
PERFCOUNTER_TIMING lines: cycle-precise nanosecond accumulators for
fillTrellis sub-phases (chunkScan / init / tmpInit / viterbi / forward /
traceback / put) and the runKSet cache-hit branch (cachedPath). The
timing report is what we should have had before pulling lever 5: it
divides time spent between cache-hit and cache-miss paths, and within
the cache-miss path between trellis init, the DP itself, and traceback.

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


# Fields on PERFCOUNTER_TIMING line that hold per-query ns accumulators.
# fillTrellis_total is the parent (sub-phases sum into it modulo work
# outside the sub-timers, like the model load); cachedPath is a sibling
# (the runKSet branch that bypasses fillTrellis entirely).
TIMING_SUBPHASES = (
    "chunkScan_ns", "init_ns", "tmpInit_ns",
    "viterbi_ns", "forward_ns", "traceback_ns", "put_ns",
)


def main(path):
    by_alg = defaultdict(lambda: {
        "ksets": [], "fillTrellis": [], "addInLogSpace": [],
        "addInLogSpace_log_exp": [], "chunk_hits": [], "chunk_lookups": [],
        "active_med": [], "active_p95": [], "active_max": [], "active_positions": [],
        "n_seqs": [],
    })
    cache_buckets = defaultdict(lambda: defaultdict(int))
    # timing[alg][field] = list of per-query ns values
    timing = defaultdict(lambda: defaultdict(list))
    # glomerator[kind] = list of {glomerator_total_ns, cumul_dpHandler_run_ns}
    glomerator = defaultdict(list)
    n_lines = 0
    n_cache_lines = 0
    n_timing_lines = 0
    n_glom_lines = 0

    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith("PERFCOUNTER_GLOMERATOR"):
            n_glom_lines += 1
            kv = parse_kv(line[len("PERFCOUNTER_GLOMERATOR "):])
            kind = kv.get("kind", "?")
            glomerator[kind].append({
                "glomerator_total_ns": int(kv["glomerator_total_ns"]),
                "cumul_dpHandler_run_ns": int(kv["cumul_dpHandler_run_ns"]),
            })
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
        if line.startswith("PERFCOUNTER_TIMING"):
            n_timing_lines += 1
            kv = parse_kv(line[len("PERFCOUNTER_TIMING "):])
            alg = kv.get("alg", "?")
            t = timing[alg]
            t["fillTrellis_total_ns"].append(int(kv["fillTrellis_total_ns"]))
            t["cachedPath_ns"].append(int(kv["cachedPath_ns"]))
            for k in TIMING_SUBPHASES:
                t[k].append(int(kv[k]))
            # Outer-perimeter accumulators (older logs omit these; default 0).
            t["runKSet_total_ns"].append(int(kv.get("runKSet_total_ns", "0")))
            t["dpHandler_run_total_ns"].append(int(kv.get("dpHandler_run_total_ns", "0")))
            # Phase-B'' deeper-DP timers (older logs omit; default 0).
            for k in ("viterbi_init_ns", "viterbi_ending_ns",
                      "forward_init_ns", "forward_ending_ns"):
                t[k].append(int(kv.get(k, "0")))
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

    print(f"# parsed {n_lines} PERFCOUNTER + {n_cache_lines} PERFCOUNTER_CACHE + {n_timing_lines} PERFCOUNTER_TIMING + {n_glom_lines} PERFCOUNTER_GLOMERATOR lines from {path}")
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

    # Phase-B timing report: where is bcrham wall actually going inside
    # the Zig core? Sub-phases sum into fillTrellis_total; cachedPath is
    # the sibling runKSet branch that bypasses fillTrellis entirely. The
    # denominator we want for "%" is "all time the Zig core spent doing
    # DP-ish work for this query," which is fillTrellis_total + cachedPath.
    # Anything not in that sum is runKSet glue, findPartialCacheMatch,
    # sorted-gene plumbing, etc. — out of the timed perimeter.
    if timing:
        print("\n## sub-phase timing (Phase-B, PERFCOUNTER_TIMING)")
        # Three concentric perimeters, biggest first:
        #   dpHandler_run_total ⊇ Σ runKSet_total ⊇ Σ fillTrellis_total + Σ cachedPath
        # Each "_overhead" row is the residual between a perimeter and the
        # sum of its inner pieces — i.e. work in the outer that isn't in
        # any inner sub-bucket. Denominator is dpHandler_run_total when
        # present (so % adds up over all reported rows); falls back to
        # fillTrellis_total + cachedPath for older logs without outer
        # timers.
        for alg, t in timing.items():
            tot_run = sum(t.get("dpHandler_run_total_ns", []))
            tot_rks = sum(t.get("runKSet_total_ns", []))
            tot_fill = sum(t["fillTrellis_total_ns"])
            tot_cached = sum(t["cachedPath_ns"])
            inner_timed = tot_fill + tot_cached
            denom = tot_run if tot_run > 0 else inner_timed
            print(f"  alg={alg}, n_queries={len(t['fillTrellis_total_ns'])}, "
                  f"dpHandler_run_total={tot_run / 1e9:.3f}s "
                  f"(runKSet={tot_rks / 1e9:.3f}s, "
                  f"fillTrellis+cachedPath={inner_timed / 1e9:.3f}s)")

            def trow(name, ns, denom):
                pct = (100 * ns / denom) if denom else 0
                print(f"    {name:<32}  {ns / 1e9:>9.3f}s  ({pct:>5.1f}% of timed)")

            if tot_run > 0:
                trow("dpHandler_run_total (outer)", tot_run, denom)
                # Pre/post-runKSet work inside run(): Sequences build,
                # only_genes map, rescaleOverallMuteFreqs, fillRecoEvent,
                # Result.finalize, debug summary.
                trow("  (run()-body pre/post-runKSet)", tot_run - tot_rks, denom)
                trow("runKSet_total (Σ ksets)", tot_rks, denom)
            trow("  fillTrellis_total (Σ misses+chunks)", tot_fill, denom)
            for k in TIMING_SUBPHASES:
                trow(f"    {k}", sum(t[k]), denom)
                # Phase-B'' DP-body sub-split (init + ending → position
                # loop derived as residual). Only meaningful for the
                # viterbi/forward rows.
                if k == "viterbi_ns":
                    vi = sum(t.get("viterbi_init_ns", [0]))
                    ve = sum(t.get("viterbi_ending_ns", [0]))
                    vp = sum(t[k]) - vi - ve
                    trow("      viterbi_init_ns", vi, denom)
                    trow("      viterbi_positionLoop_ns (derived)", vp, denom)
                    trow("      viterbi_ending_ns", ve, denom)
                elif k == "forward_ns":
                    fi = sum(t.get("forward_init_ns", [0]))
                    fe = sum(t.get("forward_ending_ns", [0]))
                    fp = sum(t[k]) - fi - fe
                    trow("      forward_init_ns", fi, denom)
                    trow("      forward_positionLoop_ns (derived)", fp, denom)
                    trow("      forward_ending_ns", fe, denom)
            sub_sum = sum(sum(t[k]) for k in TIMING_SUBPHASES)
            unaccounted_fill = tot_fill - sub_sum
            # The model load (`self.hmms.get(gene)`) sits between the
            # chunk scan and the init/tmpInit; it's not timed. So
            # fillTrellis_total − Σ sub-phases is the un-timed remainder,
            # mostly that lookup plus the addWithMinusInfinities math.
            trow("    (fillTrellis untimed)", unaccounted_fill, denom)
            trow("  cachedPath (sibling, Σ hits)", tot_cached, denom)
            if tot_rks > 0:
                # runKSet plumbing: sorted_genes alloc/sort, regions loop,
                # findPartialCacheMatch, regional_total/best updates, the
                # debug-output and fillRecoEvent assembly. Anything inside
                # runKSet body but outside fillTrellis() and cachedPath_ns.
                trow("  (runKSet body, excl. fillTrellis+cachedPath)",
                     tot_rks - inner_timed, denom)

    # Phase-B' outermost perimeter: one PERFCOUNTER_GLOMERATOR per bcrham
    # invocation that ran Glomerator.cluster() or cacheNaiveSeqs(). The
    # `cumul_dpHandler_run_ns` is the sum of every DPHandler.run() body
    # time during that bcrham process; the residual is Glomerator driver
    # overhead — merge decisions, partition I/O, hfrac calculations.
    if glomerator:
        print("\n## Glomerator wrapper (Phase-B', PERFCOUNTER_GLOMERATOR)")
        for kind, rows in glomerator.items():
            tot_glom = sum(r["glomerator_total_ns"] for r in rows)
            tot_dp = sum(r["cumul_dpHandler_run_ns"] for r in rows)
            n_proc = len(rows)
            # `glomerator_total` is summed across bcrham procs, not elapsed
            # wall — without the mean-per-proc column readers can mis-read
            # a 152s sum as a single-process wall claim. Both reported.
            mean_glom = tot_glom / n_proc if n_proc else 0
            print(f"  kind={kind}, n_bcrham_procs={n_proc}, "
                  f"mean_glom_per_proc={mean_glom / 1e9:.3f}s, "
                  f"sum_glomerator_total={tot_glom / 1e9:.3f}s "
                  f"(cumul_dpHandler_run={tot_dp / 1e9:.3f}s, "
                  f"driver_overhead={(tot_glom - tot_dp) / 1e9:.3f}s, "
                  f"driver_overhead_pct={100 * (tot_glom - tot_dp) / tot_glom if tot_glom else 0:.1f}%)")


if __name__ == "__main__":
    main(sys.argv[1])
