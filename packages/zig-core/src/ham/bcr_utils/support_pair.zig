/// ham/bcr_utils/support_pair.zig — Zig port of ham::SupportPair
///
/// A (gene-name, log-prob) pair used to rank per-gene support.
/// C++ source: packages/ham/include/bcrutils.h
/// C++ author: psathyrella/ham

const std = @import("std");

/// Corresponds to C++ `ham::SupportPair`.
pub const SupportPair = struct {
    gene: []const u8,
    logprob: f64,

    /// Less-than: returns true when `self` is LESS likely than `rhs`.
    /// Matches C++ `operator<`: return logprob() < rhs.logprob().
    pub fn lessThan(_: void, lhs: SupportPair, rhs: SupportPair) bool {
        return lhs.logprob < rhs.logprob;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "SupportPair: lessThan ordering" {
    const a = SupportPair{ .gene = "IGHV3-15*01", .logprob = -5.0 };
    const b = SupportPair{ .gene = "IGHV3-23*01", .logprob = -3.0 };
    try std.testing.expect(SupportPair.lessThan({}, a, b));
    try std.testing.expect(!SupportPair.lessThan({}, b, a));
    try std.testing.expect(!SupportPair.lessThan({}, a, a));
}
