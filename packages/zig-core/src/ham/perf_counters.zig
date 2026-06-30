/// ham/perf_counters.zig — Phase-0 diagnostics counters for issue #366.
///
/// The whole module is gated on the `-Dperf-counters` build option (see
/// `build.zig`). When the option is `false` (default), `enabled` is a
/// comptime `false` and every `if (enabled) { ... }` block in the hot
/// path compiles to nothing — bit-equality with the un-instrumented build
/// is preserved.
///
/// Counters cover the metrics requested in #366 phase 0.1 / 0.2:
///
///   - addInLogSpace_calls          — `mathutils.zig:addInLogSpace` entry
///   - addInLogSpace_log_exp_calls  — calls past both NEG_INF early-outs
///                                    (each = 1 log + 1 exp of work)
///   - fillTrellis_calls            — `dp_handler.zig:fillTrellis` entry
///   - n_ksets                      — `dp_handler.zig:run` kset loop
///   - chunk_cache_lookups/hits     — `dp_handler.zig:fillTrellis` cache
///                                    branch entry / prefix-match success
///   - active_state_hist[count]     — `trellis.zig` middleViterbiVals/
///                                    middleForwardVals
///
/// Phase-B sub-phase nanosecond accumulators (added 2026-05 after the
/// lever-5 close — see `/tmp/perf-1.1-results/next-steps-plan.md`). The
/// Phase-0 counts told us the chunk-cache HIT RATE; they did not tell us
/// how time is split between the cache-hit path (the `partial_match`
/// branch in `runKSet`) and the various sub-phases inside `fillTrellis`
/// on the cache-miss path. The lever-5 close was the proximate motivator:
/// a 24.63% perf-record symbol bucket turned out to over-attribute work
/// inside an amortized hot loop, and eliminating ~70% of those calls only
/// shifted bcrham wall by ~2%. The next-steps plan calls for cycle-
/// precise per-sub-phase timing before any further lever pull. These
/// `*_ns` accumulators provide that data:
///
///   - fillTrellis_total_ns         — whole `fillTrellis` body
///   - fillTrellis_chunkScan_ns     — chunk-cache prefix-match for-loop
///   - fillTrellis_init_ns          — cache-miss path: create+clone+initWithSeqs+append
///   - fillTrellis_tmpInit_ns       — cache-hit (trellis) path: tmptrell init
///   - fillTrellis_viterbi_ns       — `trell_ptr.viterbi()` call
///   - fillTrellis_forward_ns       — `trell_ptr.forward()` call
///   - fillTrellis_traceback_ns     — `trell_ptr.traceback(...)` call
///   - fillTrellis_put_ns           — `gc.paths.put` + `gc.scores.put`
///   - cachedPath_ns                — `runKSet`'s `partial_match` cache-hit
///                                    branch (does NOT enter `fillTrellis`)
///
/// All timers use `std.time.Instant` (CLOCK_MONOTONIC). They compile to
/// nothing when `-Dperf-counters=false` (the `Tick` type is `void`, and
/// `tick`/`addElapsed` early-return). The corresponding output line is
/// `PERFCOUNTER_TIMING` (separate from `PERFCOUNTER` so the existing
/// aggregator keeps parsing the original line unchanged).
///
/// The `_calls` vs `_log_exp_calls` split for `addInLogSpace` matters for
/// item-3.2 costing: at function entry the call may still hit a cheap
/// `NEG_INF` short-circuit and do zero transcendental work. The first
/// counter is the call rate; the second is the rate of work that an
/// item-3.2 log/exp approximation would actually replace.
///
/// State is process-global and reset at start of each `DPHandler.run()`.
/// `formatPerfCounters` returns one heap-allocated PERFCOUNTER line per
/// query; caller writes it to the destination it picks (current dump site
/// is the file at `PARTIS_ZIG_PERF_LOG`, opened O_APPEND to keep parallel
/// bcrham processes non-clobbering — see `dp_handler.zig:dumpPerfCounters`).
///
/// SAFETY: the counters below are plain `pub var` globals, NOT atomics.
/// Every caller must run single-threaded within a process. bcrham is
/// single-threaded by design (#342 item 1 investigated and rejected
/// `smp_allocator`); cross-process parallelism via `partis --n-procs N`
/// runs each bcrham in its own address space, so the globals stay
/// per-process. If a future caller introduces threading inside a single
/// bcrham process, this module must grow atomic counters or per-thread
/// shadow state.
const std = @import("std");
const build_options = @import("build_options");

pub const enabled: bool = build_options.perf_counters;

/// Length of the active-state histogram. State indices in this codebase
/// are bounded by `state.zig:STATE_MAX = 500`; one extra slot for 500 itself.
const HIST_LEN: usize = 501;

pub var addInLogSpace_calls: u64 = 0;
pub var addInLogSpace_log_exp_calls: u64 = 0;
pub var fillTrellis_calls: u64 = 0;
pub var n_ksets: u64 = 0;
pub var chunk_cache_lookups: u64 = 0;
pub var chunk_cache_hits: u64 = 0;
pub var active_state_hist: [HIST_LEN]u64 = [_]u64{0} ** HIST_LEN;

// Phase-B sub-phase ns accumulators. Wrap-around (`+%=`) is fine: a u64
// of nanoseconds is ~584 years.
pub var fillTrellis_total_ns: u64 = 0;
pub var fillTrellis_chunkScan_ns: u64 = 0;
pub var fillTrellis_init_ns: u64 = 0;
pub var fillTrellis_tmpInit_ns: u64 = 0;
pub var fillTrellis_viterbi_ns: u64 = 0;
pub var fillTrellis_forward_ns: u64 = 0;
pub var fillTrellis_traceback_ns: u64 = 0;
pub var fillTrellis_put_ns: u64 = 0;
pub var cachedPath_ns: u64 = 0;

// Phase-B outer-perimeter accumulators. These pair with the inner ones
// above to bound the un-timed gap inside one `DPHandler.run()`:
//   runKSet_overhead       = runKSet_total       − (fillTrellis_total + cachedPath)
//   dpHandler_run_overhead = dpHandler_run_total − Σ runKSet_total
// The latter covers the run()-body pre/post work: Sequences build, the
// only_genes map, `hmms.rescaleOverallMuteFreqs`, `fillRecoEvent` for
// viterbi, `Result.finalize`, debug logging, etc. Phase B saw the
// inner perimeter (fillTrellis+cachedPath) is ~95% DP; these outer
// accumulators let us check whether the remaining bcrham wall really
// is Glomerator+driver, or whether there's an addressable bucket
// hiding in run() pre/post or runKSet plumbing.
pub var runKSet_total_ns: u64 = 0;
pub var dpHandler_run_total_ns: u64 = 0;

// Phase-B'' deeper-DP timers: split `viterbi()` and `forward()` (in
// `trellis.zig`) into init / position-loop / ending sub-phases. The
// position-loop body is derived in the aggregator as the residual:
//   viterbi_positionLoop_ns ≈ fillTrellis_viterbi_ns − init − ending
//   forward_positionLoop_ns ≈ fillTrellis_forward_ns − init − ending
// Expected outcome: init + ending sum to <5%, so the position loop
// — i.e. `middleViterbiVals`/`middleForwardVals` + `swapColumnsActive`
// — accounts for ~95%+ of the DP time. That confirms the only
// remaining lever surface is the inner DP body, which is not localized
// (it's the algorithm itself).
pub var viterbi_init_ns: u64 = 0;
pub var viterbi_ending_ns: u64 = 0;
pub var forward_init_ns: u64 = 0;
pub var forward_ending_ns: u64 = 0;

// Process-lifetime accumulator: total wall in `Glomerator.cluster()` OR
// `Glomerator.cacheNaiveSeqs()`. The dumped `PERFCOUNTER_GLOMERATOR`
// line tags which mode via its `kind=` field. Distinct from the
// per-DPHandler.run() counters above — this one survives the per-query
// `reset()` so it accumulates across all the nested DPHandler.run()
// invocations inside one Glomerator entry-point call. Bcrham emits one
// `PERFCOUNTER_GLOMERATOR` line at the end of each call reporting
// this value plus the cumulative DPHandler.run time across all calls
// inside it (the latter via a paired `cumul_dpHandler_run_ns`
// accumulator that ALSO survives reset()).
//
// Together with the `bcrham time:` line printed by partis, the gap
//   bcrham_CPU − cumul_dpHandler_run_ns × n_procs
// tells us whether the remaining bcrham CPU is in Glomerator driver
// code (merge decisions, partition writes) or in bcrham startup /
// model loading. Phase B' shows ~28% of bcrham CPU lives outside
// DPHandler.run(); this counter lets us split it.
pub var glomerator_ns: u64 = 0;
pub var cumul_dpHandler_run_ns: u64 = 0;

pub fn reset() void {
    if (!enabled) return;
    addInLogSpace_calls = 0;
    addInLogSpace_log_exp_calls = 0;
    fillTrellis_calls = 0;
    n_ksets = 0;
    chunk_cache_lookups = 0;
    chunk_cache_hits = 0;
    fillTrellis_total_ns = 0;
    fillTrellis_chunkScan_ns = 0;
    fillTrellis_init_ns = 0;
    fillTrellis_tmpInit_ns = 0;
    fillTrellis_viterbi_ns = 0;
    fillTrellis_forward_ns = 0;
    fillTrellis_traceback_ns = 0;
    fillTrellis_put_ns = 0;
    cachedPath_ns = 0;
    runKSet_total_ns = 0;
    dpHandler_run_total_ns = 0;
    viterbi_init_ns = 0;
    viterbi_ending_ns = 0;
    forward_init_ns = 0;
    forward_ending_ns = 0;
    // NB: explicit loop instead of `@memset(&active_state_hist, 0)` —
    // Zig 0.15.2 hits an x86_64 codegen bug ("no encoding found for:
    // none mov m64 m64 none none") on a Debug @memset over a fixed-size
    // global u64 array. The loop sidesteps the bug and any optimizer
    // will turn it back into a memset under ReleaseFast/ReleaseSafe.
    var i: usize = 0;
    while (i < HIST_LEN) : (i += 1) active_state_hist[i] = 0;
}

/// Cycle-precise sub-phase timer pair, gated on the `-Dperf-counters`
/// build option. When disabled, `Tick` is `void`, both helpers early-
/// return, and the call site reduces to nothing under the optimizer —
/// preserving the default-build bit-equality contract.
///
/// Usage:
///   const t = perf_counters.tick();
///   // ... work to time ...
///   perf_counters.addElapsed(&perf_counters.fillTrellis_chunkScan_ns, t);
pub const Tick = if (enabled) std.time.Instant else void;

pub inline fn tick() Tick {
    if (!enabled) return {};
    // `Instant.now()` only returns `error.Unsupported` on platforms
    // without a monotonic clock. Linux/macOS/Windows all provide one.
    return std.time.Instant.now() catch unreachable;
}

pub inline fn addElapsed(accum: *u64, start: Tick) void {
    if (!enabled) return;
    const now_inst = std.time.Instant.now() catch unreachable;
    accum.* +%= now_inst.since(start);
}

pub inline fn bumpAddInLogSpace() void {
    if (enabled) addInLogSpace_calls += 1;
}

pub inline fn bumpAddInLogSpaceLogExp() void {
    if (enabled) addInLogSpace_log_exp_calls += 1;
}

pub inline fn bumpFillTrellis() void {
    if (enabled) fillTrellis_calls += 1;
}

pub inline fn bumpKSet() void {
    if (enabled) n_ksets += 1;
}

pub inline fn bumpChunkCacheLookup() void {
    if (enabled) chunk_cache_lookups += 1;
}

pub inline fn bumpChunkCacheHit() void {
    if (enabled) chunk_cache_hits += 1;
}

pub inline fn recordActiveStates(n: usize) void {
    if (!enabled) return;
    const idx = if (n >= HIST_LEN) HIST_LEN - 1 else n;
    active_state_hist[idx] += 1;
}

/// Compute (median, p95, max, total_positions) over the active-state histogram.
fn activeStateStats() struct { median: usize, p95: usize, max: usize, total: u64 } {
    var total: u64 = 0;
    for (active_state_hist) |c| total += c;
    if (total == 0) return .{ .median = 0, .p95 = 0, .max = 0, .total = 0 };
    const median_target = (total + 1) / 2;
    const p95_target = (total * 95 + 99) / 100;
    var cum: u64 = 0;
    var median: usize = 0;
    var p95: usize = 0;
    var max_idx: usize = 0;
    var i: usize = 0;
    while (i < HIST_LEN) : (i += 1) {
        if (active_state_hist[i] == 0) continue;
        max_idx = i;
        const next_cum = cum + active_state_hist[i];
        if (cum < median_target and median_target <= next_cum) median = i;
        if (cum < p95_target and p95_target <= next_cum) p95 = i;
        cum = next_cum;
    }
    return .{ .median = median, .p95 = p95, .max = max_idx, .total = total };
}

/// Format the per-query counter summary as a heap-allocated `[]u8` (caller frees).
/// Returns null when counters are disabled (caller skips the write).
/// Format is intentionally trivial-to-grep:
///   PERFCOUNTER alg=... q=... n_seqs=... ksets=... fillTrellis=... \
///       chunk_hits=H/L addInLogSpace=... addInLogSpace_log_exp=... \
///       active_states_med=... active_states_p95=... active_states_max=... \
///       active_states_positions=...
///
/// `chunk_hits=H/L`: H = prefix-match successes, L = lookup attempts (the
/// `if (!no_chunk_cache)` branch entries in fillTrellis). Hit rate is H/L.
pub fn formatPerfCounters(
    allocator: std.mem.Allocator,
    algorithm: []const u8,
    query_name: []const u8,
    n_seqs: usize,
) !?[]u8 {
    if (!enabled) return null;
    const stats = activeStateStats();
    return try std.fmt.allocPrint(
        allocator,
        "PERFCOUNTER alg={s} q={s} n_seqs={d} ksets={d} fillTrellis={d} chunk_hits={d}/{d} addInLogSpace={d} addInLogSpace_log_exp={d} active_states_med={d} active_states_p95={d} active_states_max={d} active_states_positions={d}\n",
        .{
            algorithm,
            query_name,
            n_seqs,
            n_ksets,
            fillTrellis_calls,
            chunk_cache_hits,
            chunk_cache_lookups,
            addInLogSpace_calls,
            addInLogSpace_log_exp_calls,
            stats.median,
            stats.p95,
            stats.max,
            stats.total,
        },
    );
}

/// Format the per-query sub-phase ns accumulators as a separate
/// `PERFCOUNTER_TIMING` line. Kept distinct from `PERFCOUNTER` so the
/// existing aggregator parses unchanged; the timing aggregator reads only
/// these lines.
///
///   PERFCOUNTER_TIMING alg=... q=... n_seqs=... \
///       fillTrellis_total_ns=T chunkScan_ns=... init_ns=... \
///       tmpInit_ns=... viterbi_ns=... forward_ns=... traceback_ns=... \
///       put_ns=... cachedPath_ns=...
pub fn formatPerfCountersTiming(
    allocator: std.mem.Allocator,
    algorithm: []const u8,
    query_name: []const u8,
    n_seqs: usize,
) !?[]u8 {
    if (!enabled) return null;
    return try std.fmt.allocPrint(
        allocator,
        "PERFCOUNTER_TIMING alg={s} q={s} n_seqs={d} dpHandler_run_total_ns={d} runKSet_total_ns={d} fillTrellis_total_ns={d} chunkScan_ns={d} init_ns={d} tmpInit_ns={d} viterbi_ns={d} forward_ns={d} traceback_ns={d} put_ns={d} cachedPath_ns={d} viterbi_init_ns={d} viterbi_ending_ns={d} forward_init_ns={d} forward_ending_ns={d}\n",
        .{
            algorithm,
            query_name,
            n_seqs,
            dpHandler_run_total_ns,
            runKSet_total_ns,
            fillTrellis_total_ns,
            fillTrellis_chunkScan_ns,
            fillTrellis_init_ns,
            fillTrellis_tmpInit_ns,
            fillTrellis_viterbi_ns,
            fillTrellis_forward_ns,
            fillTrellis_traceback_ns,
            fillTrellis_put_ns,
            cachedPath_ns,
            viterbi_init_ns,
            viterbi_ending_ns,
            forward_init_ns,
            forward_ending_ns,
        },
    );
}
