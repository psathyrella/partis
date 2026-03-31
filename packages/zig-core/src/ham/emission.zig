/// ham/emission.zig — Zig port of ham/src/emission.cc + ham/include/emission.h
///
/// Emission probability model for a single HMM state.
/// Wraps a LexicalTable and provides Parse (from YAML-extracted data), score, and rescaling.
///
/// C++ source: packages/ham/src/emission.cc, packages/ham/include/emission.h
/// C++ author: psathyrella/ham

const std = @import("std");
const Track = @import("track.zig").Track;
const Sequence = @import("sequences.zig").Sequence;
const LexicalTable = @import("lexical_table.zig").LexicalTable;

extern fn log(x: f64) f64;

/// EPS for normalization check (matches C++ build define -DEPS=1e-6).
const EPS: f64 = 1e-6;

/// Emission probability model for a single HMM state.
/// Corresponds to C++ `ham::Emission`.
pub const Emission = struct {
    /// Sum of raw probabilities (checked during parse; kept for diagnostics).
    total: f64,
    /// Inner lookup table (indexed by digitized symbol).
    scores: LexicalTable,
    /// Pointer to the track (not owned; must outlive this Emission).
    track: ?*const Track,

    /// Create an empty Emission with no data.
    /// Corresponds to C++ `Emission()`.
    pub fn init() Emission {
        return Emission{
            .total = 0.0,
            .scores = LexicalTable.init(),
            .track = null,
        };
    }

    pub fn deinit(self: *Emission, allocator: std.mem.Allocator) void {
        self.scores.deinit(allocator);
    }

    /// Parse emission probabilities from an already-extracted prob map.
    ///
    /// `probs_map` is a `StringHashMap(f64)` mapping each symbol name (e.g. "A") to its
    /// raw probability.  This corresponds to the `probs` sub-node extracted from the YAML
    /// `emissions` node by the caller (ham/Model).
    ///
    /// Corresponds to C++ `Emission::Parse(YAML::Node config, Track *track)`.
    pub fn parse(
        self: *Emission,
        allocator: std.mem.Allocator,
        probs_map: std.StringHashMap(f64),
        trk: *const Track,
    ) !void {
        self.track = trk;
        self.scores.setTrack(trk);

        const alphabet_size = trk.alphabetSize();
        if (probs_map.count() != alphabet_size) {
            return error.EmissionProbsSizeMismatch;
        }

        var log_probs = try allocator.alloc(f64, alphabet_size);
        defer allocator.free(log_probs);

        self.total = 0.0;
        for (0..alphabet_size) |ip| {
            const sym = trk.symbol(ip);
            const prob = probs_map.get(sym) orelse return error.MissingEmissionSymbol;
            log_probs[ip] = log(prob);
            self.total += prob;
        }
        if (@abs(self.total - 1.0) >= EPS) {
            return error.BadNormalization;
        }
        try self.scores.setLogProbs(allocator, log_probs);
    }

    /// Replace log probs (for mute-freq rescaling).
    /// Corresponds to C++ `void ReplaceLogProbs(vector<double>)`.
    pub fn replaceLogProbs(self: *Emission, allocator: std.mem.Allocator, new_log_probs: []const f64) !void {
        try self.scores.replaceLogProbs(allocator, new_log_probs);
    }

    /// Revert to original log probs.
    /// Corresponds to C++ `void UnReplaceLogProbs()`.
    pub fn unReplaceLogProbs(self: *Emission, allocator: std.mem.Allocator) void {
        self.scores.unReplaceLogProbs(allocator);
    }

    /// Score (log-probability) for a digitized symbol index.
    /// Corresponds to C++ `inline double score(uint8_t index)`.
    pub inline fn scoreIndex(self: *const Emission, index: u8) f64 {
        return self.scores.logProb(index);
    }

    /// Score (log-probability) for the symbol at `pos` in `seq`.
    /// Corresponds to C++ `double score(Sequence *seq, size_t pos)`.
    pub fn scoreSeq(self: *const Emission, seq: *const Sequence, pos: usize) f64 {
        return self.scores.logProbSeq(seq, pos);
    }

    /// Return a copy of the log_probs slice.
    /// Corresponds to C++ `vector<double> log_probs()`.
    pub fn logProbs(self: *const Emission, allocator: std.mem.Allocator) ![]f64 {
        return self.scores.logProbs(allocator);
    }

    /// Return the track pointer.
    pub fn getTrack(self: *const Emission) ?*const Track {
        return self.track;
    }

    /// Print emission probabilities to stdout.
    /// Corresponds to C++ `Emission::Print()`.
    pub fn print(self: *const Emission) void {
        const trk = self.track orelse return;
        const n = trk.alphabetSize();
        std.debug.print("    {s}     (normed to within at least {e})\n", .{ trk.name, EPS });
        for (0..n) |i| {
            const prob = @exp(self.scores.logProb(@intCast(i)));
            if (prob < 0.01) {
                std.debug.print("{e:>22}", .{prob});
            } else {
                std.debug.print("{d:>22.9}", .{prob});
            }
        }
        std.debug.print("\n", .{});
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Emission: parse uniform distribution" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var probs = std.StringHashMap(f64).init(allocator);
    defer probs.deinit();
    try probs.put("A", 0.25);
    try probs.put("C", 0.25);
    try probs.put("G", 0.25);
    try probs.put("T", 0.25);

    var em = Emission.init();
    defer em.deinit(allocator);

    try em.parse(allocator, probs, &track);

    try std.testing.expectApproxEqAbs(1.0, em.total, EPS);
    try std.testing.expectApproxEqAbs(@log(0.25), em.scoreIndex(0), 1e-9);
    try std.testing.expectApproxEqAbs(@log(0.25), em.scoreIndex(3), 1e-9);
}

test "Emission: parse non-uniform distribution" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var probs = std.StringHashMap(f64).init(allocator);
    defer probs.deinit();
    try probs.put("A", 0.7);
    try probs.put("C", 0.1);
    try probs.put("G", 0.1);
    try probs.put("T", 0.1);

    var em = Emission.init();
    defer em.deinit(allocator);

    try em.parse(allocator, probs, &track);
    try std.testing.expectApproxEqAbs(@log(0.7), em.scoreIndex(0), 1e-9);
    try std.testing.expectApproxEqAbs(@log(0.1), em.scoreIndex(1), 1e-9);
}

test "Emission: parse rejects bad normalization" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var probs = std.StringHashMap(f64).init(allocator);
    defer probs.deinit();
    try probs.put("A", 0.3);
    try probs.put("C", 0.3);
    try probs.put("G", 0.2);
    try probs.put("T", 0.1); // sums to 0.9

    var em = Emission.init();
    defer em.deinit(allocator);

    const result = em.parse(allocator, probs, &track);
    try std.testing.expectError(error.BadNormalization, result);
}

test "Emission: scoreSeq" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var probs = std.StringHashMap(f64).init(allocator);
    defer probs.deinit();
    try probs.put("A", 0.1);
    try probs.put("C", 0.2);
    try probs.put("G", 0.3);
    try probs.put("T", 0.4);

    var em = Emission.init();
    defer em.deinit(allocator);
    try em.parse(allocator, probs, &track);

    var seq = try @import("sequences.zig").Sequence.initFromString(allocator, &track, "s1", "ACGT");
    defer seq.deinit(allocator);

    // A at pos 0 → log(0.1)
    try std.testing.expectApproxEqAbs(@log(0.1), em.scoreSeq(&seq, 0), 1e-9);
    // T at pos 3 → log(0.4)
    try std.testing.expectApproxEqAbs(@log(0.4), em.scoreSeq(&seq, 3), 1e-9);
}
