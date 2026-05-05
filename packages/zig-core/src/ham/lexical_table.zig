/// ham/lexical_table.zig — Zig port of ham/src/lexicaltable.cc + ham/include/lexicaltable.h
///
/// Emission probability lookup table: maps a digitized symbol index → log-probability.
/// Used by the Emission model and ultimately by Trellis (Forward/Viterbi).
///
/// C++ source: packages/ham/src/lexicaltable.cc, packages/ham/include/lexicaltable.h
/// C++ author: psathyrella/ham

const std = @import("std");
const Track = @import("track.zig").Track;
const Sequence = @import("sequences.zig").Sequence;

/// EPS for normalization check (matches C++ build define -DEPS=1e-6).
const EPS: f64 = 1e-6;

/// Emission probability lookup table indexed by digitized symbol.
/// Corresponds to C++ `ham::LexicalTable`.
pub const LexicalTable = struct {
    /// Pointer to the track (not owned; must outlive this LexicalTable).
    track: ?*const Track,
    /// Log-probabilities for each symbol in the track's alphabet.
    /// Owned by this struct, allocated by setLogProbs and reused thereafter.
    log_probs: []f64,
    /// Saved copy of log_probs from the last replaceLogProbs call. Sized
    /// the same as log_probs and allocated together by setLogProbs (item 15
    /// of #342). Treated as a scratch buffer: replace/unReplace @memcpy
    /// rather than realloc, so per-query rescale/un-rescale cycles do not
    /// hit the allocator (~50k calls per 5k-seq run).
    original_log_probs: []f64,

    /// Create an empty LexicalTable with no track.
    /// Corresponds to C++ `LexicalTable()`.
    pub fn init() LexicalTable {
        return LexicalTable{
            .track = null,
            .log_probs = &[_]f64{},
            .original_log_probs = &[_]f64{},
        };
    }

    /// Set the track pointer.
    /// Corresponds to C++ `LexicalTable::Init(Track*)`.
    pub fn setTrack(self: *LexicalTable, trk: *const Track) void {
        self.track = trk;
    }

    pub fn deinit(self: *LexicalTable, allocator: std.mem.Allocator) void {
        if (self.log_probs.len > 0) allocator.free(self.log_probs);
        if (self.original_log_probs.len > 0) allocator.free(self.original_log_probs);
    }

    /// Set (replace) log_probs from a slice, taking a copy. Also (re)allocates
    /// the `original_log_probs` scratch buffer to the same size, so subsequent
    /// replace/unReplace calls can @memcpy without hitting the allocator.
    /// Corresponds to C++ `LexicalTable::SetLogProbs(vector<double>)`.
    pub fn setLogProbs(self: *LexicalTable, allocator: std.mem.Allocator, logprobs: []const f64) !void {
        if (self.log_probs.len > 0) allocator.free(self.log_probs);
        if (self.original_log_probs.len > 0) allocator.free(self.original_log_probs);
        self.log_probs = try allocator.dupe(f64, logprobs);
        self.original_log_probs = try allocator.alloc(f64, logprobs.len);
    }

    /// Replace log_probs with new values, saving the originals for unReplaceLogProbs.
    /// Asserts that sizes match. Checks that exp(new probs) sum to ~1.0 within EPS.
    /// Corresponds to C++ `LexicalTable::ReplaceLogProbs(vector<double>)`.
    /// `allocator` is unused now that the scratch buffer is allocated up front
    /// in setLogProbs (item 15 of #342); kept for API symmetry with setLogProbs.
    pub fn replaceLogProbs(self: *LexicalTable, allocator: std.mem.Allocator, new_log_probs: []const f64) !void {
        _ = allocator;
        std.debug.assert(self.log_probs.len == new_log_probs.len);
        std.debug.assert(self.original_log_probs.len == new_log_probs.len);
        @memcpy(self.original_log_probs, self.log_probs);
        @memcpy(self.log_probs, new_log_probs);
        var total: f64 = 0.0;
        for (self.log_probs) |lp| total += @exp(lp);
        if (@abs(total - 1.0) >= EPS) {
            return error.BadNormalization;
        }
    }

    /// Revert log_probs to the values saved by replaceLogProbs.
    /// Corresponds to C++ `LexicalTable::UnReplaceLogProbs()`.
    /// `allocator` is unused now that the scratch buffer is allocated up front
    /// in setLogProbs (item 15 of #342); kept for API symmetry with setLogProbs.
    pub fn unReplaceLogProbs(self: *LexicalTable, allocator: std.mem.Allocator) void {
        _ = allocator;
        std.debug.assert(self.original_log_probs.len == self.log_probs.len);
        @memcpy(self.log_probs, self.original_log_probs);
    }

    /// Return the log-probability for a digitized symbol index.
    /// Corresponds to C++ `inline double LogProb(uint8_t index)`.
    pub fn logProb(self: *const LexicalTable, index: u8) f64 {
        return self.log_probs[index];
    }

    /// Return the log-probability for the symbol at `pos` in `seq`.
    /// Corresponds to C++ `double LogProb(Sequence *seq, size_t pos)`.
    pub fn logProbSeq(self: *const LexicalTable, seq: *const Sequence, pos: usize) f64 {
        std.debug.assert(pos < seq.size());
        const idx = seq.seqq[pos];
        std.debug.assert(idx < self.log_probs.len);
        return self.log_probs[idx];
    }

    /// Return a copy of the log_probs slice.
    /// Corresponds to C++ `vector<double> log_probs()` (which returns a copy).
    pub fn logProbs(self: *const LexicalTable, allocator: std.mem.Allocator) ![]f64 {
        return allocator.dupe(f64, self.log_probs);
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "LexicalTable: init and setLogProbs" {
    const allocator = std.testing.allocator;
    var lt = LexicalTable.init();
    defer lt.deinit(allocator);

    const lps = [_]f64{ @log(0.25), @log(0.25), @log(0.25), @log(0.25) };
    try lt.setLogProbs(allocator, &lps);

    try std.testing.expectEqual(@as(usize, 4), lt.log_probs.len);
    try std.testing.expectApproxEqAbs(@log(0.25), lt.logProb(0), 1e-9);
    try std.testing.expectApproxEqAbs(@log(0.25), lt.logProb(3), 1e-9);
}

test "LexicalTable: replaceLogProbs and unReplaceLogProbs" {
    const allocator = std.testing.allocator;
    var lt = LexicalTable.init();
    defer lt.deinit(allocator);

    const orig = [_]f64{ @log(0.25), @log(0.25), @log(0.25), @log(0.25) };
    try lt.setLogProbs(allocator, &orig);

    // Replace with new valid distribution
    const new_lps = [_]f64{ @log(0.5), @log(0.3), @log(0.1), @log(0.1) };
    try lt.replaceLogProbs(allocator, &new_lps);
    try std.testing.expectApproxEqAbs(@log(0.5), lt.logProb(0), 1e-9);

    // Revert
    lt.unReplaceLogProbs(allocator);
    try std.testing.expectApproxEqAbs(@log(0.25), lt.logProb(0), 1e-9);
}

test "LexicalTable: replaceLogProbs rejects bad normalization" {
    const allocator = std.testing.allocator;
    var lt = LexicalTable.init();
    defer lt.deinit(allocator);

    const orig = [_]f64{ @log(0.25), @log(0.25), @log(0.25), @log(0.25) };
    try lt.setLogProbs(allocator, &orig);

    // Distribution that sums to 0.9, not 1.0
    const bad = [_]f64{ @log(0.3), @log(0.3), @log(0.2), @log(0.1) };
    const result = lt.replaceLogProbs(allocator, &bad);
    try std.testing.expectError(error.BadNormalization, result);
}

test "LexicalTable: logProbSeq" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try @import("track.zig").Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var lt = LexicalTable.init();
    defer lt.deinit(allocator);
    lt.setTrack(&track);

    const lps = [_]f64{ @log(0.1), @log(0.2), @log(0.3), @log(0.4) };
    try lt.setLogProbs(allocator, &lps);

    var seq = try @import("sequences.zig").Sequence.initFromString(allocator, &track, "s1", "ACGT");
    defer seq.deinit(allocator);

    // A → index 0 → log(0.1)
    try std.testing.expectApproxEqAbs(@log(0.1), lt.logProbSeq(&seq, 0), 1e-9);
    // T → index 3 → log(0.4)
    try std.testing.expectApproxEqAbs(@log(0.4), lt.logProbSeq(&seq, 3), 1e-9);
}
