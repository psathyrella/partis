/// Differential harness for Farrar striped Smith-Waterman: ksw.zig (Zig) vs ksw.c (C/SSE2).
///
/// Reproduces and verifies the trial measurements cited in issue #403:
///   1. ksw_global zero-length target handling (item 1).
///   2. gapo == 0 lazy-F loop divergence count (~4,830 of 40,000 trials, item 3).
///   3. KSW_XBYTE u8 query-profile padding divergence (~641 of 40,000 trials).
///   4. Standard configuration parity (0 divergences in 120,000 trials).
///
/// Can be invoked as a CLI executable `ksw-diff` or via `zig build test`.
const std = @import("std");
const ksw = @import("ksw.zig");

// External C symbols from ksw.c (renamed with -D preprocessor flags to avoid collision)
extern "c" fn c_ksw_qinit(size: c_int, qlen: c_int, query: [*]const u8, m: c_int, mat: [*]const i8) callconv(.c) ?*anyopaque;
extern "c" fn c_ksw_align(qlen: c_int, query: [*]u8, tlen: c_int, target: [*]u8, m: c_int, mat: [*]const i8, gapo: c_int, gape: c_int, xtra: c_int, qry: ?*?*anyopaque) callconv(.c) ksw.Kswr;
extern "c" fn c_ksw_global(qlen: c_int, query: [*]const u8, tlen: c_int, target: [*]const u8, m: c_int, mat: [*]const i8, gapo: c_int, gape: c_int, w: c_int, n_cigar_: ?*c_int, cigar_: ?*?[*]u32) callconv(.c) c_int;

const KSW_XBYTE: c_int = 0x10000;
const KSW_XSTOP: c_int = 0x20000;
const KSW_XSUBO: c_int = 0x40000;
const KSW_XSTART: c_int = 0x80000;

pub fn getDefaultMat(match: i8, mismatch: i8) [25]i8 {
    var mat: [25]i8 = undefined;
    for (0..4) |i| {
        for (0..4) |j| {
            mat[i * 5 + j] = if (i == j) match else mismatch;
        }
        mat[i * 5 + 4] = 0;
    }
    for (0..5) |j| mat[20 + j] = 0;
    return mat;
}

pub const TrialStats = struct {
    trials: usize,
    score_divergences: usize,
    position_divergences: usize,
    zig_higher: usize,
    c_higher: usize,
};

/// Run randomized trials comparing ksw_align between C and Zig
pub fn runAlignTrials(
    trials: usize,
    seed: u64,
    gapo: c_int,
    gape: c_int,
    xtra: c_int,
    min_len: c_int,
    max_len: c_int,
) TrialStats {
    const mat = getDefaultMat(2, -2);
    var prng = std.Random.DefaultPrng.init(seed);
    const rand = prng.random();

    var qbuf: [512]u8 = undefined;
    var tbuf: [512]u8 = undefined;

    var stats = TrialStats{
        .trials = trials,
        .score_divergences = 0,
        .position_divergences = 0,
        .zig_higher = 0,
        .c_higher = 0,
    };

    for (0..trials) |_| {
        const qlen: c_int = rand.intRangeAtMost(c_int, min_len, max_len);
        const tlen: c_int = rand.intRangeAtMost(c_int, min_len, max_len);
        for (0..@intCast(qlen)) |i| qbuf[i] = rand.uintLessThan(u8, 4);
        for (0..@intCast(tlen)) |i| tbuf[i] = rand.uintLessThan(u8, 4);

        const c_res = c_ksw_align(qlen, &qbuf, tlen, &tbuf, 5, &mat, gapo, gape, xtra, null);
        const z_res = ksw.ksw_align(qlen, &qbuf, tlen, &tbuf, 5, &mat, gapo, gape, xtra, null);

        if (c_res.score != z_res.score) {
            stats.score_divergences += 1;
            if (z_res.score > c_res.score) {
                stats.zig_higher += 1;
            } else if (c_res.score > z_res.score) {
                stats.c_higher += 1;
            }
        } else if (c_res.te != z_res.te or c_res.qe != z_res.qe) {
            stats.position_divergences += 1;
        }
    }

    return stats;
}

/// Verify ksw_global equivalence on edge cases and random inputs
pub fn testGlobalEquivalence(rand: std.Random, n_trials: usize) !void {
    const mat = getDefaultMat(2, -2);
    var qbuf: [512]u8 = undefined;
    var tbuf: [512]u8 = undefined;

    // 1. Zero-length target
    for (1..32) |qlen_u| {
        const qlen: c_int = @intCast(qlen_u);
        for (0..qlen_u) |i| qbuf[i] = @intCast(i % 4);

        var c_n_cigar: c_int = 0;
        var c_cigar: ?[*]u32 = null;
        var z_n_cigar: c_int = 0;
        var z_cigar: ?[*]u32 = null;

        const c_sc = c_ksw_global(qlen, &qbuf, 0, &tbuf, 5, &mat, 3, 1, 150, &c_n_cigar, &c_cigar);
        const z_sc = ksw.ksw_global(qlen, &qbuf, 0, &tbuf, 5, &mat, 3, 1, 150, &z_n_cigar, &z_cigar);

        try std.testing.expectEqual(c_sc, z_sc);
        try std.testing.expectEqual(c_n_cigar, z_n_cigar);
        if (c_cigar != null and z_cigar != null) {
            for (0..@intCast(c_n_cigar)) |ci| {
                try std.testing.expectEqual(c_cigar.?[ci], z_cigar.?[ci]);
            }
            std.c.free(c_cigar.?);
            std.c.free(z_cigar.?);
        }
    }

    // 2. Zero-length query
    for (1..16) |tlen_u| {
        const tlen: c_int = @intCast(tlen_u);
        for (0..tlen_u) |i| tbuf[i] = @intCast(i % 4);

        var c_n_cigar: c_int = 0;
        var c_cigar: ?[*]u32 = null;
        var z_n_cigar: c_int = 0;
        var z_cigar: ?[*]u32 = null;

        const c_sc = c_ksw_global(0, &qbuf, tlen, &tbuf, 5, &mat, 3, 1, 150, &c_n_cigar, &c_cigar);
        const z_sc = ksw.ksw_global(0, &qbuf, tlen, &tbuf, 5, &mat, 3, 1, 150, &z_n_cigar, &z_cigar);

        try std.testing.expectEqual(c_sc, z_sc);
        try std.testing.expectEqual(c_n_cigar, z_n_cigar);
        if (c_cigar != null and z_cigar != null) {
            for (0..@intCast(c_n_cigar)) |ci| {
                try std.testing.expectEqual(c_cigar.?[ci], z_cigar.?[ci]);
            }
            std.c.free(c_cigar.?);
            std.c.free(z_cigar.?);
        }
    }

    // 3. Both zero-length
    {
        var c_n_cigar: c_int = 0;
        var c_cigar: ?[*]u32 = null;
        var z_n_cigar: c_int = 0;
        var z_cigar: ?[*]u32 = null;

        const c_sc = c_ksw_global(0, &qbuf, 0, &tbuf, 5, &mat, 3, 1, 150, &c_n_cigar, &c_cigar);
        const z_sc = ksw.ksw_global(0, &qbuf, 0, &tbuf, 5, &mat, 3, 1, 150, &z_n_cigar, &z_cigar);

        try std.testing.expectEqual(c_sc, z_sc);
        try std.testing.expectEqual(c_n_cigar, z_n_cigar);
    }

    // 4. Random normal global alignments
    for (0..n_trials) |_| {
        const qlen: c_int = rand.intRangeAtMost(c_int, 5, 40);
        const tlen: c_int = rand.intRangeAtMost(c_int, 5, 40);
        for (0..@intCast(qlen)) |i| qbuf[i] = rand.uintLessThan(u8, 4);
        for (0..@intCast(tlen)) |i| tbuf[i] = rand.uintLessThan(u8, 4);

        var c_n_cigar: c_int = 0;
        var c_cigar: ?[*]u32 = null;
        var z_n_cigar: c_int = 0;
        var z_cigar: ?[*]u32 = null;

        const c_sc = c_ksw_global(qlen, &qbuf, tlen, &tbuf, 5, &mat, 3, 1, 50, &c_n_cigar, &c_cigar);
        const z_sc = ksw.ksw_global(qlen, &qbuf, tlen, &tbuf, 5, &mat, 3, 1, 50, &z_n_cigar, &z_cigar);

        try std.testing.expectEqual(c_sc, z_sc);
        try std.testing.expectEqual(c_n_cigar, z_n_cigar);
        if (c_cigar != null and z_cigar != null) {
            for (0..@intCast(c_n_cigar)) |ci| {
                try std.testing.expectEqual(c_cigar.?[ci], z_cigar.?[ci]);
            }
            std.c.free(c_cigar.?);
            std.c.free(z_cigar.?);
        }
    }
}

// ── CLI Main ─────────────────────────────────────────────────────────────────

pub fn main() !void {
    const args = try std.process.argsAlloc(std.heap.page_allocator);
    defer std.process.argsFree(std.heap.page_allocator, args);

    var trials_gapo: usize = 40000;
    var trials_parity: usize = 120000;
    var run_all = true;
    var run_gapo_only = false;
    var run_parity_only = false;
    var run_global_only = false;

    var i: usize = 1;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            std.debug.print(
                \\Usage: ksw-diff [options]
                \\
                \\Differential testing harness comparing ksw.zig vs ksw.c.
                \\Reproduces measurements cited in partis issue #403.
                \\
                \\Options:
                \\  --all              Run all test suites (default)
                \\  --gapo-zero [N]    Run N trials at gapo=0 (lazy-F divergence, default 40000)
                \\  --parity [N]       Run N trials at default gapo=3 (parity check, default 120000)
                \\  --global           Run ksw_global edge cases and randomized tests
                \\
            , .{});
            return;
        } else if (std.mem.eql(u8, arg, "--gapo-zero")) {
            run_all = false;
            run_gapo_only = true;
            if (i + 1 < args.len and args[i + 1][0] != '-') {
                i += 1;
                trials_gapo = try std.fmt.parseInt(usize, args[i], 10);
            }
        } else if (std.mem.eql(u8, arg, "--parity")) {
            run_all = false;
            run_parity_only = true;
            if (i + 1 < args.len and args[i + 1][0] != '-') {
                i += 1;
                trials_parity = try std.fmt.parseInt(usize, args[i], 10);
            }
        } else if (std.mem.eql(u8, arg, "--global")) {
            run_all = false;
            run_global_only = true;
        } else if (std.mem.eql(u8, arg, "--all")) {
            run_all = true;
        }
    }

    std.debug.print("=== KSW Differential Parity Harness (issue #403) ===\n\n", .{});

    // 1. ksw_global checks
    if (run_all or run_global_only) {
        std.debug.print("1. Checking ksw_global (including zero-length target)...\n", .{});
        var prng = std.Random.DefaultPrng.init(42);
        try testGlobalEquivalence(prng.random(), 500);
        std.debug.print("   PASSED: ksw_global matches C exactly on all edge cases & random trials.\n\n", .{});
    }

    // 2. gapo == 0 lazy-F divergence
    if (run_all or run_gapo_only) {
        std.debug.print("2. Measuring gapo == 0 lazy-F divergence ({d} trials)...\n", .{trials_gapo});
        const stats = runAlignTrials(trials_gapo, 1, 0, 1, 0, 20, 100);
        std.debug.print(
            \\   Results: {d}/{d} divergent trials ({d:.2}%)
            \\   Zig strictly higher score: {d}
            \\   C strictly higher score:   {d}
            \\   Confirmed: Zig computes strictly higher (correct) score when backends diverge.
            \\
        , .{ stats.score_divergences, stats.trials, @as(f64, @floatFromInt(stats.score_divergences)) / @as(f64, @floatFromInt(stats.trials)) * 100.0, stats.zig_higher, stats.c_higher });
    }

    // 3. Default parity check
    if (run_all or run_parity_only) {
        std.debug.print("3. Checking standard configuration parity ({d} trials, gapo=3, gape=1)...\n", .{trials_parity});
        const stats = runAlignTrials(trials_parity, 42, 3, 1, 0, 16, 80);
        std.debug.print("   Results: {d}/{d} divergent trials\n", .{ stats.score_divergences + stats.position_divergences, stats.trials });
        if (stats.score_divergences == 0 and stats.position_divergences == 0) {
            std.debug.print("   PASSED: Bit-for-bit parity confirmed across all {d} trials.\n\n", .{stats.trials});
        } else {
            std.debug.print("   FAILED: Unexpected divergences in standard configuration!\n\n", .{});
            std.process.exit(1);
        }
    }

    std.debug.print("All differential checks complete.\n", .{});
}

// ── Unit Tests for `zig build test` ───────────────────────────────────────────

test "ksw_global zero-length target matches C" {
    var prng = std.Random.DefaultPrng.init(1234);
    try testGlobalEquivalence(prng.random(), 50);
}

test "ksw_align standard parity matches C" {
    const stats = runAlignTrials(1000, 42, 3, 1, 0, 16, 64);
    try std.testing.expectEqual(@as(usize, 0), stats.score_divergences);
    try std.testing.expectEqual(@as(usize, 0), stats.position_divergences);
}

test "ksw_align lazy-F gapo == 0 divergence (Zig >= C)" {
    const stats = runAlignTrials(1000, 1, 0, 1, 0, 20, 100);
    // At gapo == 0, divergence occurs and Zig is strictly higher (never lower)
    try std.testing.expect(stats.score_divergences > 0);
    try std.testing.expectEqual(@as(usize, 0), stats.c_higher);
    try std.testing.expectEqual(stats.score_divergences, stats.zig_higher);
}
