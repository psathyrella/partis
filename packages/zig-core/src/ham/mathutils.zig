/// ham/mathutils.zig — Zig port of ham/include/mathutils.h + ham/src/mathutils.cc
///
/// Math utility functions used throughout ham.
/// All functions in namespace ham (free functions, no class).
///
/// C++ source: packages/ham/src/mathutils.cc, packages/ham/include/mathutils.h
/// C++ author: psathyrella/ham
///
/// Checkpoint strategy: per-call checkpoints are impractical for HMM inner-loop
/// functions (e.g. AddWithMinusInfinities is called ~117k times per sequence).
/// Checkpoint coverage is provided by higher-level callers (DpHandler, Trellis)
/// once those are ported.

const std = @import("std");
const math = std.math;

// Use libc log/exp to match C++ bcrham bit-for-bit across all platforms.
// Zig @log/@exp map to LLVM intrinsics which may differ from glibc on some CPUs.
extern fn log(x: f64) f64;
extern fn exp(x: f64) f64;

/// Negative infinity constant for log-probability computations.
/// Replaces the verbose `-std.math.inf(f64)` throughout the codebase.
pub const NEG_INF: f64 = -math.inf(f64);

// ── POWER table ───────────────────────────────────────────────────────────────
//
// Corresponds to C++ `POWER[b][a-1] = a**b` in mathutils.h.
// Used by LexicalTable for alphabet size lookups.
// In C++ this is a compile-time static const uint64_t[32][128].
// We represent it the same way: comptime-evaluated lookup table.
//
// Rather than replicate the 32×128 literal (which is large), we compute it
// at comptime. C++ defines POWER[b][a-1] = a**b (note: a-1 indexing).
// We expose a function power(b, a) = a**b matching the C++ semantics.
pub fn power(b: usize, a: usize) u64 {
    if (b == 0) return 1;
    var result: u64 = 1;
    var i: usize = 0;
    while (i < b) : (i += 1) {
        result *= @as(u64, a);
    }
    return result;
}

// ── AddWithMinusInfinities ────────────────────────────────────────────────────

/// Add two log-space values, treating -inf as the absorbing zero element.
/// Computes log(a * b) = log(a) + log(b) where -inf represents probability 0.
/// Corresponds to C++ `ham::AddWithMinusInfinities(double first, double second)`.
pub fn addWithMinusInfinities(first: f64, second: f64) f64 {
    if (first == NEG_INF or second == NEG_INF)
        return NEG_INF;
    return first + second;
}

// ── AddInLogSpace ─────────────────────────────────────────────────────────────

/// Log-sum-exp: computes log(a + b) from (log a, log b).
/// Implements the log-space *or* operation (a OR b = a + b in probability space).
/// Corresponds to C++ `ham::AddInLogSpace<T>(T first, T second)`.
pub fn addInLogSpace(first: f64, second: f64) f64 {
    if (first == NEG_INF) return second;
    if (second == NEG_INF) return first;
    if (first > second) {
        return first + log(1.0 + exp(second - first));
    } else {
        return second + log(1.0 + exp(first - second));
    }
}

// ── Vector utilities ──────────────────────────────────────────────────────────
//
// These correspond to the C++ template functions in mathutils.h.
// All operate on slices (the Zig equivalent of vector<T>).

/// Sum all elements of a slice.
/// Corresponds to C++ `ham::sumVector(vector<T>& data)`.
pub fn sumVector(comptime T: type, data: []const T) T {
    var s: T = 0;
    for (data) |v| s += v;
    return s;
}

/// Product of all elements in a slice.
/// Corresponds to C++ `ham::productVector(vector<T>& data)`.
pub fn productVector(comptime T: type, data: []const T) T {
    if (data.len == 0) return 1;
    var p: T = data[0];
    for (data[1..]) |v| p *= v;
    return p;
}

/// Minimum element.
/// Corresponds to C++ `ham::minVector(vector<T>& data)`.
pub fn minVector(comptime T: type, data: []const T) T {
    return std.mem.min(T, data);
}

/// Maximum element.
/// Corresponds to C++ `ham::maxVector(vector<T>& data)`.
pub fn maxVector(comptime T: type, data: []const T) T {
    return std.mem.max(T, data);
}

/// Average of all elements.
/// Corresponds to C++ `ham::avgVector(vector<T>& data)`.
pub fn avgVector(comptime T: type, data: []const T) T {
    return sumVector(T, data) / @as(T, @floatFromInt(data.len));
}

/// Take log of each element in-place.
/// Corresponds to C++ `ham::logVector(vector<T>& data)`.
pub fn logVector(data: []f64) void {
    for (data) |*v| v.* = @log(v.*);
}

/// Return a new slice of exp(x) for each element in data.
/// Caller owns the returned slice.
/// Corresponds to C++ `ham::get_exp_vector(vector<T> data)`.
pub fn getExpVector(allocator: std.mem.Allocator, data: []const f64) ![]f64 {
    const result = try allocator.alloc(f64, data.len);
    for (data, 0..) |v, i| result[i] = @exp(v);
    return result;
}

/// Normalize data in-place to sum to 1 (convert to probabilities).
/// Corresponds to C++ `ham::probVector(vector<T>& data)`.
pub fn probVector(data: []f64) void {
    const s = sumVector(f64, data);
    if (s == 0.0) {
        for (data) |*v| v.* = 0.0;
    } else {
        for (data) |*v| v.* /= s;
    }
}

/// Add rhs element-wise into lhs in-place.
/// Corresponds to C++ `ham::addVector(vector<T>& lhs, vector<T>& rhs)`.
pub fn addVector(comptime T: type, lhs: []T, rhs: []const T) void {
    std.debug.assert(lhs.len == rhs.len);
    for (lhs, rhs) |*l, r| l.* += r;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

test "addWithMinusInfinities: both finite" {
    const result = addWithMinusInfinities(-1.5, -2.3);
    try std.testing.expectApproxEqAbs(-3.8, result, 1e-10);
}

test "addWithMinusInfinities: first -inf" {
    const result = addWithMinusInfinities(NEG_INF, -2.3);
    try std.testing.expect(result == NEG_INF);
}

test "addWithMinusInfinities: second -inf" {
    const result = addWithMinusInfinities(-1.5, NEG_INF);
    try std.testing.expect(result == NEG_INF);
}

test "addInLogSpace: log(e^0 + e^0) = log(2)" {
    const result = addInLogSpace(0.0, 0.0);
    try std.testing.expectApproxEqAbs(@log(2.0), result, 1e-10);
}

test "addInLogSpace: first -inf" {
    const result = addInLogSpace(NEG_INF, -1.5);
    try std.testing.expectApproxEqAbs(-1.5, result, 1e-10);
}

test "addInLogSpace: second -inf" {
    const result = addInLogSpace(-1.5, NEG_INF);
    try std.testing.expectApproxEqAbs(-1.5, result, 1e-10);
}

test "sumVector" {
    const data = [_]f64{ 1.0, 2.0, 3.0 };
    try std.testing.expectApproxEqAbs(6.0, sumVector(f64, &data), 1e-10);
}

test "productVector" {
    const data = [_]f64{ 2.0, 3.0, 4.0 };
    try std.testing.expectApproxEqAbs(24.0, productVector(f64, &data), 1e-10);
}

test "power" {
    try std.testing.expectEqual(@as(u64, 8), power(3, 2)); // 2^3 = 8
    try std.testing.expectEqual(@as(u64, 1), power(0, 5)); // 5^0 = 1
    try std.testing.expectEqual(@as(u64, 4), power(2, 2)); // 2^2 = 4
}
