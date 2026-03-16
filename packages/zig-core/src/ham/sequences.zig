/// ham/sequences.zig — Zig port of ham/src/sequences.cc + ham/include/sequences.h
///
/// Sequence container: stores a nucleotide string with its digitized (u8 index) form.
///
/// C++ source: packages/ham/src/sequences.cc, packages/ham/include/sequences.h
/// C++ author: psathyrella/ham

const std = @import("std");
const ham_text = @import("text.zig");
const Track = @import("track.zig").Track;

/// A single sequence with its name, undigitized string, and digitized u8 slice.
/// Corresponds to C++ `ham::Sequence`.
pub const Sequence = struct {
    name: []u8,
    /// Optional header line (e.g. FASTA header).
    header: []u8,
    /// The original nucleotide string (with any whitespace stripped).
    undigitized: []u8,
    /// Digitized form: each character → u8 index via the track's alphabet.
    /// Owned by this struct.
    seqq: []u8,
    /// Pointer to the track (not owned; must outlive this Sequence).
    track: ?*const Track,

    /// Initialize an empty Sequence (not generally useful; exists for container compatibility).
    /// Corresponds to C++ `Sequence()`.
    pub fn init() Sequence {
        return Sequence{
            .name = &[_]u8{},
            .header = &[_]u8{},
            .undigitized = &[_]u8{},
            .seqq = &[_]u8{},
            .track = null,
        };
    }

    /// Create a Sequence from a track, name, and undigitized string.
    /// Strips `\n` from the string (matching C++ ClearWhitespace call).
    /// Digitizes the sequence via the track's symbol_index.
    /// Corresponds to C++ `Sequence(Track*, string name, string& undigitized)`.
    pub fn initFromString(
        allocator: std.mem.Allocator,
        trk: *const Track,
        name: []const u8,
        undigitized_in: []const u8,
    ) !Sequence {
        var seq = Sequence{
            .name = try allocator.dupe(u8, name),
            .header = try allocator.dupe(u8, ""),
            .undigitized = try allocator.dupe(u8, undigitized_in),
            .seqq = try allocator.dupe(u8, ""),
            .track = trk,
        };
        // ClearWhitespace("\n", &undigitized_)
        var ud_list = std.ArrayList(u8).fromOwnedSlice(seq.undigitized);
        ham_text.clearWhitespace(allocator, "\n", &ud_list);
        seq.undigitized = try ud_list.toOwnedSlice(allocator);
        // Digitize (free the placeholder empty allocation)
        allocator.free(seq.seqq);
        seq.seqq = try digitize(allocator, trk, seq.undigitized);
        return seq;
    }

    /// Create a subsequence from position `pos` of length `len`.
    /// Corresponds to C++ `Sequence(Track*, string name, string& undigitized, size_t pos, size_t len)`.
    pub fn initSubstring(
        allocator: std.mem.Allocator,
        trk: *const Track,
        name: []const u8,
        undigitized_in: []const u8,
        pos: usize,
        len: usize,
    ) !Sequence {
        try checkPosLen(name, undigitized_in, pos, len);
        const sub = undigitized_in[pos .. pos + len];
        return initFromString(allocator, trk, name, sub);
    }

    /// Create a subsequence of an existing Sequence (slice without re-digitizing from string).
    /// Corresponds to C++ `Sequence(Sequence& rhs, size_t pos, size_t len)`.
    pub fn initSlice(
        allocator: std.mem.Allocator,
        rhs: *const Sequence,
        pos: usize,
        len: usize,
    ) !Sequence {
        try checkPosLen(rhs.name, rhs.undigitized, pos, len);
        return Sequence{
            .name = try allocator.dupe(u8, rhs.name),
            .header = try allocator.dupe(u8, rhs.header),
            .undigitized = try allocator.dupe(u8, rhs.undigitized[pos .. pos + len]),
            .seqq = try allocator.dupe(u8, rhs.seqq[pos .. pos + len]),
            .track = rhs.track,
        };
    }

    /// Clone a Sequence (deep copy).
    /// Corresponds to C++ `Sequence(const Sequence& rhs)`.
    pub fn clone(self: *const Sequence, allocator: std.mem.Allocator) !Sequence {
        return Sequence{
            .name = try allocator.dupe(u8, self.name),
            .header = try allocator.dupe(u8, self.header),
            .undigitized = try allocator.dupe(u8, self.undigitized),
            .seqq = try allocator.dupe(u8, self.seqq),
            .track = self.track,
        };
    }

    pub fn deinit(self: *Sequence, allocator: std.mem.Allocator) void {
        allocator.free(self.name);
        allocator.free(self.header);
        allocator.free(self.undigitized);
        allocator.free(self.seqq);
    }

    pub fn size(self: *const Sequence) usize {
        return self.seqq.len;
    }

    pub fn setName(self: *Sequence, allocator: std.mem.Allocator, nm: []const u8) !void {
        allocator.free(self.name);
        self.name = try allocator.dupe(u8, nm);
    }
};

/// Validate that `[pos, pos+len)` is in bounds for `undigitized`.
fn checkPosLen(name: []const u8, undigitized: []const u8, pos: usize, len: usize) !void {
    if (pos >= undigitized.len or pos + len > undigitized.len) {
        std.debug.print("ERROR: len {d} too large for '{s}' in '{s}'\n", .{ len, undigitized, name });
        return error.InvalidPosLen;
    }
}

/// Digitize `undigitized` using the track's symbol_index.
fn digitize(allocator: std.mem.Allocator, trk: *const Track, undigitized: []const u8) ![]u8 {
    const result = try allocator.alloc(u8, undigitized.len);
    for (undigitized, 0..) |_, i| {
        const sym = undigitized[i .. i + 1];
        result[i] = trk.symbolIndex(sym) catch |err| {
            allocator.free(result);
            return err;
        };
    }
    return result;
}

/// A collection of equal-length Sequences.
/// Corresponds to C++ `ham::Sequences`.
pub const Sequences = struct {
    seqs: std.ArrayListUnmanaged(Sequence),
    /// Length of each sequence (all must be equal).
    sequence_length: usize,

    pub fn init() Sequences {
        return Sequences{
            .seqs = .{},
            .sequence_length = 0,
        };
    }

    pub fn deinit(self: *Sequences, allocator: std.mem.Allocator) void {
        for (self.seqs.items) |*s| s.deinit(allocator);
        self.seqs.deinit(allocator);
    }

    pub fn nSeqs(self: *const Sequences) usize {
        return self.seqs.items.len;
    }

    pub fn get(self: *const Sequences, index: usize) *const Sequence {
        return &self.seqs.items[index];
    }

    pub fn getPtr(self: *Sequences, index: usize) *Sequence {
        return &self.seqs.items[index];
    }

    /// Deep-copy this Sequences (matches C++ copy constructor semantics).
    pub fn clone(self: *const Sequences, allocator: std.mem.Allocator) !Sequences {
        var result = Sequences.init();
        errdefer result.deinit(allocator);
        for (self.seqs.items) |*sq| {
            try result.addSeq(allocator, try sq.clone(allocator));
        }
        return result;
    }

    /// Add a sequence, taking ownership.
    /// All sequences must have the same length.
    /// Corresponds to C++ `Sequences::AddSeq(Sequence sq)`.
    pub fn addSeq(self: *Sequences, allocator: std.mem.Allocator, sq: Sequence) !void {
        if (self.seqs.items.len == 0) {
            self.sequence_length = sq.size();
        } else if (sq.size() != self.sequence_length) {
            return error.SequenceLengthMismatch;
        }
        try self.seqs.append(allocator, sq);
    }

    /// Return a new Sequences with the union of self and otherseqs.
    /// Detects duplicate names.
    pub fn unionWith(
        self: *const Sequences,
        allocator: std.mem.Allocator,
        other: *const Sequences,
    ) !Sequences {
        var result = Sequences.init();
        var seen = std.StringHashMap(void).init(allocator);
        defer seen.deinit();

        for (self.seqs.items) |*s| {
            if (seen.contains(s.name)) return error.DuplicateSequenceName;
            try seen.put(s.name, {});
            try result.addSeq(allocator, try s.clone(allocator));
        }
        for (other.seqs.items) |*s| {
            if (seen.contains(s.name)) return error.DuplicateSequenceName;
            try seen.put(s.name, {});
            try result.addSeq(allocator, try s.clone(allocator));
        }
        return result;
    }

    /// Return a comma-separated string of sequence names.
    /// Corresponds to C++ `Sequences::name_str(string delimiter)`.
    pub fn nameStr(self: *const Sequences, allocator: std.mem.Allocator, delimiter: []const u8) ![]u8 {
        var names = std.ArrayList([]const u8).init(allocator);
        defer names.deinit();
        for (self.seqs.items) |*s| try names.append(s.name);
        return ham_text.joinStrings(allocator, names.items, delimiter);
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Sequence: initFromString digitizes correctly" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var seq = try Sequence.initFromString(allocator, &track, "test_seq", "ACGT");
    defer seq.deinit(allocator);

    try std.testing.expectEqualStrings("test_seq", seq.name);
    try std.testing.expectEqual(@as(usize, 4), seq.size());
    try std.testing.expectEqual(@as(u8, 0), seq.seqq[0]); // A
    try std.testing.expectEqual(@as(u8, 1), seq.seqq[1]); // C
    try std.testing.expectEqual(@as(u8, 2), seq.seqq[2]); // G
    try std.testing.expectEqual(@as(u8, 3), seq.seqq[3]); // T
}

test "Sequences: addSeq and length check" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var seqs = Sequences.init();
    defer seqs.deinit(allocator);

    const s1 = try Sequence.initFromString(allocator, &track, "s1", "ACGT");
    try seqs.addSeq(allocator, s1);
    try std.testing.expectEqual(@as(usize, 4), seqs.sequence_length);

    // Adding a sequence of different length should error
    const s2 = try Sequence.initFromString(allocator, &track, "s2", "AC");
    var bad_seq = s2;
    const result = seqs.addSeq(allocator, bad_seq);
    _ = result catch {};
    bad_seq.deinit(allocator);
}
