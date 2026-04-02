/// ham/track.zig — Zig port of ham/src/track.cc + ham/include/track.h
///
/// Alphabet/symbol-index mapping for HMM tracks (nucleotide sequences).
///
/// C++ source: packages/ham/src/track.cc, packages/ham/include/track.h
/// C++ author: psathyrella/ham

const std = @import("std");

/// Maximum number of symbols in an alphabet (matches C++ `max_alphabet_size_`).
pub const MAX_ALPHABET_SIZE: usize = 255;

/// The reserved index for the ambiguous character (matches C++ `ambiguous_index_`).
pub const AMBIGUOUS_INDEX: u8 = MAX_ALPHABET_SIZE - 1;

/// Symbol alphabet and index map for HMM tracks.
/// Corresponds to C++ `ham::Track`.
pub const Track = struct {
    name: []u8,
    /// Ordered list of symbols in the alphabet (owned slices).
    alphabet: std.ArrayListUnmanaged([]u8),
    /// Map from symbol string to u8 index (keys are slices owned by `alphabet`).
    symbol_indices: std.StringHashMapUnmanaged(u8),
    /// The string representing an ambiguous/wildcard character (e.g. "N").
    /// Empty string means no ambiguous character set.
    ambiguous_char: []u8,

    pub fn init(allocator: std.mem.Allocator, name: []const u8) !Track {
        return Track{
            .name = try allocator.dupe(u8, name),
            .alphabet = .{},
            .symbol_indices = .{},
            .ambiguous_char = try allocator.dupe(u8, ""),
        };
    }

    /// Create a Track with a name, initial symbol list, and optional ambiguous char.
    /// Corresponds to C++ `Track(string name, vector<string> symbols, string ambiguous_char)`.
    pub fn initWithSymbols(
        allocator: std.mem.Allocator,
        name: []const u8,
        symbols: []const []const u8,
        ambiguous_char: []const u8,
    ) !Track {
        var t = try Track.init(allocator, name);
        allocator.free(t.ambiguous_char);
        t.ambiguous_char = try allocator.dupe(u8, ambiguous_char);
        try t.addSymbols(allocator, symbols);
        return t;
    }

    pub fn deinit(self: *Track, allocator: std.mem.Allocator) void {
        allocator.free(self.name);
        allocator.free(self.ambiguous_char);
        // alphabet owns all the symbol strings
        for (self.alphabet.items) |sym| allocator.free(sym);
        self.alphabet.deinit(allocator);
        self.symbol_indices.deinit(allocator);
    }

    pub fn setName(self: *Track, allocator: std.mem.Allocator, nm: []const u8) !void {
        allocator.free(self.name);
        self.name = try allocator.dupe(u8, nm);
    }

    /// Add a single symbol to the alphabet.
    /// Corresponds to C++ `Track::AddSymbol(string symbol)`.
    pub fn addSymbol(self: *Track, allocator: std.mem.Allocator, sym: []const u8) !void {
        std.debug.assert(self.alphabet.items.len < MAX_ALPHABET_SIZE);
        const owned = try allocator.dupe(u8, sym);
        try self.alphabet.append(allocator, owned);
        const index: u8 = @intCast(self.alphabet.items.len - 1);
        try self.symbol_indices.put(allocator, owned, index);
    }

    /// Add multiple symbols.
    /// Corresponds to C++ `Track::AddSymbols(vector<string>& symbols)`.
    pub fn addSymbols(self: *Track, allocator: std.mem.Allocator, symbols: []const []const u8) !void {
        for (symbols) |s| try self.addSymbol(allocator, s);
    }

    pub fn setAmbiguous(self: *Track, allocator: std.mem.Allocator, amb: []const u8) !void {
        allocator.free(self.ambiguous_char);
        self.ambiguous_char = try allocator.dupe(u8, amb);
    }

    pub fn alphabetSize(self: *const Track) usize {
        return self.alphabet.items.len;
    }

    /// Get the symbol string at position `iter`.
    /// Corresponds to C++ `Track::symbol(size_t iter)`.
    pub fn symbol(self: *const Track, iter: usize) []const u8 {
        return self.alphabet.items[iter];
    }

    /// Return the u8 index for a symbol string.
    /// Returns AMBIGUOUS_INDEX for the ambiguous character.
    /// Corresponds to C++ `Track::symbol_index(const string &symbol)`.
    pub fn symbolIndex(self: *const Track, s: []const u8) !u8 {
        if (self.ambiguous_char.len > 0 and std.mem.eql(u8, s, self.ambiguous_char)) {
            return AMBIGUOUS_INDEX;
        }
        return self.symbol_indices.get(s) orelse
            return error.SymbolNotFound;
    }

    /// Return a human-readable string listing all symbols.
    /// Caller owns the returned string.
    /// Corresponds to C++ `Track::Stringify()`.
    pub fn stringify(self: *const Track, allocator: std.mem.Allocator) ![]u8 {
        var buf: std.ArrayListUnmanaged(u8) = .{};
        errdefer buf.deinit(allocator);
        for (self.alphabet.items) |sym| {
            try buf.append(allocator, ' ');
            try buf.appendSlice(allocator, sym);
        }
        if (self.ambiguous_char.len > 0) {
            try buf.appendSlice(allocator, "   ambiguous: ");
            try buf.appendSlice(allocator, self.ambiguous_char);
        }
        return buf.toOwnedSlice(allocator);
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Track: init and addSymbol" {
    const allocator = std.testing.allocator;
    var track = try Track.init(allocator, "nucs");
    defer track.deinit(allocator);

    try track.addSymbol(allocator, "A");
    try track.addSymbol(allocator, "C");
    try track.addSymbol(allocator, "G");
    try track.addSymbol(allocator, "T");

    try std.testing.expectEqual(@as(usize, 4), track.alphabetSize());
    try std.testing.expectEqualStrings("G", track.symbol(2));
    try std.testing.expectEqual(@as(u8, 0), try track.symbolIndex("A"));
    try std.testing.expectEqual(@as(u8, 3), try track.symbolIndex("T"));
}

test "Track: ambiguous character" {
    const allocator = std.testing.allocator;
    const syms = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &syms, "N");
    defer track.deinit(allocator);

    try std.testing.expectEqual(AMBIGUOUS_INDEX, try track.symbolIndex("N"));
}

test "Track: symbolIndex error on unknown" {
    const allocator = std.testing.allocator;
    var track = try Track.init(allocator, "nucs");
    defer track.deinit(allocator);
    try track.addSymbol(allocator, "A");

    try std.testing.expectError(error.SymbolNotFound, track.symbolIndex("Z"));
}

test "Track: stringify" {
    const allocator = std.testing.allocator;
    const syms = [_][]const u8{ "A", "C" };
    var track = try Track.initWithSymbols(allocator, "nucs", &syms, "");
    defer track.deinit(allocator);
    const s = try track.stringify(allocator);
    defer allocator.free(s);
    try std.testing.expectEqualStrings(" A C", s);
}
