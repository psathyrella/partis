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

pub fn reset() void {
    if (!enabled) return;
    addInLogSpace_calls = 0;
    addInLogSpace_log_exp_calls = 0;
    fillTrellis_calls = 0;
    n_ksets = 0;
    chunk_cache_lookups = 0;
    chunk_cache_hits = 0;
    // NB: explicit loop instead of `@memset(&active_state_hist, 0)` —
    // Zig 0.15.2 hits an x86_64 codegen bug ("no encoding found for:
    // none mov m64 m64 none none") on a Debug @memset over a fixed-size
    // global u64 array. The loop sidesteps the bug and any optimizer
    // will turn it back into a memset under ReleaseFast/ReleaseSafe.
    var i: usize = 0;
    while (i < HIST_LEN) : (i += 1) active_state_hist[i] = 0;
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
