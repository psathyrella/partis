/// ham/bcr_utils/germ_lines.zig — Zig port of ham::GermLines
///
/// Loads germline FASTA sequences + cyst/tryp positions for a locus.
/// C++ source: packages/ham/src/bcrutils.cc

const std = @import("std");
const ham_text = @import("../text.zig");

pub const regions = [_][]const u8{ "v", "d", "j" };

/// Corresponds to C++ `ham::GermLines`.
pub const GermLines = struct {
    locus: []u8,
    /// Region gene name lists (v/d/j keys).
    names: std.StringHashMapUnmanaged(std.ArrayListUnmanaged([]u8)),
    /// Gene-name → sequence mapping.
    seqs: std.StringHashMapUnmanaged([]u8),
    /// Gene-name → cyst position.
    cyst_positions: std.StringHashMapUnmanaged(i32),
    /// Gene-name → tryp position.
    tryp_positions: std.StringHashMapUnmanaged(i32),
    /// For loci without D genes, the synthetic dummy D gene name.
    dummy_d_gene: ?[]u8,

    allocator: std.mem.Allocator,

    /// Load germline data from `gldir` for `locus`.
    /// Corresponds to C++ `GermLines::GermLines(gldir, locus)`.
    pub fn init(allocator: std.mem.Allocator, gldir: []const u8, locus: []const u8) !GermLines {
        var self = GermLines{
            .locus = try allocator.dupe(u8, locus),
            .names = .{},
            .seqs = .{},
            .cyst_positions = .{},
            .tryp_positions = .{},
            .dummy_d_gene = null,
            .allocator = allocator,
        };
        errdefer self.deinit();

        // Initialise empty name lists for each region
        for (regions) |r| {
            try self.names.put(allocator, try allocator.dupe(u8, r), .{});
        }

        const has_d = hasDGene(locus);

        if (!has_d) {
            // Build dummy D gene name e.g. "IGKDx-x*x"
            const upper = try allocator.dupe(u8, locus);
            defer allocator.free(upper);
            for (upper) |*c| c.* = std.ascii.toUpper(c.*);
            const dummy = try std.fmt.allocPrint(allocator, "{s}Dx-x*x", .{upper});
            self.dummy_d_gene = dummy;

            if (self.names.getPtr("d")) |d_list| {
                try d_list.append(allocator, try allocator.dupe(u8, dummy));
            }
            try self.seqs.put(allocator, try allocator.dupe(u8, dummy), try allocator.dupe(u8, "A"));
        }

        // Read FASTA files for each region
        for (regions) |r| {
            if (!has_d and std.mem.eql(u8, r, "d")) continue;

            const infname = try std.fmt.allocPrint(allocator, "{s}/{s}/{s}{s}.fasta", .{ gldir, locus, locus, r });
            defer allocator.free(infname);

            const content = std.fs.cwd().readFileAlloc(allocator, infname, 64 * 1024 * 1024) catch |err| {
                return err; // file must exist
            };
            defer allocator.free(content);

            var name_buf: ?[]u8 = null;
            var seq_buf: std.ArrayListUnmanaged(u8) = .{};
            defer seq_buf.deinit(allocator);

            var line_iter = std.mem.splitScalar(u8, content, '\n');
            while (line_iter.next()) |raw| {
                const line = std.mem.trimRight(u8, raw, "\r");
                if (line.len == 0) continue;

                if (line[0] == '>') {
                    // flush previous sequence
                    if (name_buf) |nm| {
                        const seq_owned = try allocator.dupe(u8, seq_buf.items);
                        try self.seqs.put(allocator, nm, seq_owned);
                        seq_buf.clearRetainingCapacity();
                        name_buf = null;
                    }
                    // parse new header: take up to first space after '>'
                    const rest = line[1..];
                    const sp = std.mem.indexOfScalar(u8, rest, ' ') orelse rest.len;
                    const gene_name = try allocator.dupe(u8, rest[0..sp]);
                    name_buf = gene_name;
                    if (self.names.getPtr(r)) |lst| {
                        try lst.append(allocator, try allocator.dupe(u8, gene_name));
                    }
                } else {
                    try seq_buf.appendSlice(allocator, line);
                }
            }
            // flush last sequence
            if (name_buf) |nm| {
                const seq_owned = try allocator.dupe(u8, seq_buf.items);
                try self.seqs.put(allocator, nm, seq_owned);
            }
        }

        // Read extras.csv for cyst/tryp positions
        const extras_fname = try std.fmt.allocPrint(allocator, "{s}/{s}/extras.csv", .{ gldir, locus });
        defer allocator.free(extras_fname);

        const extras = std.fs.cwd().readFileAlloc(allocator, extras_fname, 4 * 1024 * 1024) catch null;
        if (extras) |content| {
            defer allocator.free(content);
            var line_iter = std.mem.splitScalar(u8, content, '\n');

            // Read and checkpoint the header line (matches C++ SplitString call).
            if (line_iter.next()) |header_raw| {
                const header_line = std.mem.trim(u8, header_raw, " \t\r");
                if (header_line.len > 0) {
                    var header_parts = try ham_text.split_string(allocator, header_line, ",");
                    defer {
                        for (header_parts.items) |s| allocator.free(s);
                        header_parts.deinit(allocator);
                    }
                }
            }

            while (line_iter.next()) |raw| {
                const line = std.mem.trim(u8, raw, " \t\r");
                if (line.len == 0) continue;

                // Use split_string to emit checkpoint (matches C++ instrumentation).
                var parts = try ham_text.split_string(allocator, line, ",");
                defer {
                    for (parts.items) |s| allocator.free(s);
                    parts.deinit(allocator);
                }

                if (parts.items.len < 2) continue;
                const gene_name = parts.items[0];
                const cyst_str = parts.items[1];
                const tryp_str = if (parts.items.len > 2) parts.items[2] else "";
                const phen_str = if (parts.items.len > 3) parts.items[3] else "";

                if (cyst_str.len > 0) {
                    const val = std.fmt.parseInt(i32, cyst_str, 10) catch continue;
                    try self.cyst_positions.put(allocator, try allocator.dupe(u8, gene_name), val);
                } else if (tryp_str.len > 0) {
                    const val = std.fmt.parseInt(i32, tryp_str, 10) catch continue;
                    try self.tryp_positions.put(allocator, try allocator.dupe(u8, gene_name), val);
                } else if (phen_str.len > 0) {
                    const val = std.fmt.parseInt(i32, phen_str, 10) catch continue;
                    try self.tryp_positions.put(allocator, try allocator.dupe(u8, gene_name), val);
                }
            }
        }

        return self;
    }

    pub fn deinit(self: *GermLines) void {
        const allocator = self.allocator;

        allocator.free(self.locus);
        if (self.dummy_d_gene) |d| allocator.free(d);

        // free names map
        var names_it = self.names.iterator();
        while (names_it.next()) |entry| {
            for (entry.value_ptr.items) |s| allocator.free(s);
            entry.value_ptr.deinit(allocator);
            allocator.free(entry.key_ptr.*);
        }
        self.names.deinit(allocator);

        // free seqs map
        var seqs_it = self.seqs.iterator();
        while (seqs_it.next()) |entry| {
            allocator.free(entry.key_ptr.*);
            allocator.free(entry.value_ptr.*);
        }
        self.seqs.deinit(allocator);

        // free cyst positions
        var cyst_it = self.cyst_positions.iterator();
        while (cyst_it.next()) |entry| allocator.free(entry.key_ptr.*);
        self.cyst_positions.deinit(allocator);

        // free tryp positions
        var tryp_it = self.tryp_positions.iterator();
        while (tryp_it.next()) |entry| allocator.free(entry.key_ptr.*);
        self.tryp_positions.deinit(allocator);
    }

    /// Sanitize gene name: replace '*' with '_star_' and '/' with '_slash_'.
    /// Corresponds to C++ `GermLines::SanitizeName`.
    pub fn sanitizeName(allocator: std.mem.Allocator, gene_name: []const u8) ![]u8 {
        // Replace * with _star_, then / with _slash_ (chaining ownership)
        const step1 = try std.mem.replaceOwned(u8, allocator, gene_name, "*", "_star_");
        errdefer allocator.free(step1);
        const step2 = try std.mem.replaceOwned(u8, allocator, step1, "/", "_slash_");
        allocator.free(step1);
        return step2;
    }

    /// Get region letter from gene name.
    /// Corresponds to C++ `GermLines::GetRegion`.
    pub fn getRegion(gene: []const u8) !u8 {
        if (!std.mem.startsWith(u8, gene, "IG") and !std.mem.startsWith(u8, gene, "TR"))
            return error.BadGeneName;
        if (gene.len < 4) return error.BadGeneName;
        const ch = std.ascii.toLower(gene[3]);
        if (ch != 'v' and ch != 'd' and ch != 'j') return error.BadRegion;
        return ch;
    }
};

/// Whether the locus has a D gene segment.
/// Corresponds to C++ `ham::HasDGene`.
pub fn hasDGene(locus: []const u8) bool {
    return std.mem.eql(u8, locus, "igh") or
        std.mem.eql(u8, locus, "trb") or
        std.mem.eql(u8, locus, "trd");
}

// ── Tests ─────────────────────────────────────────────────────────────────────

test "hasDGene" {
    try std.testing.expect(hasDGene("igh"));
    try std.testing.expect(hasDGene("trb"));
    try std.testing.expect(!hasDGene("igk"));
    try std.testing.expect(!hasDGene("igl"));
}

test "GermLines.sanitizeName" {
    const allocator = std.testing.allocator;
    const s = try GermLines.sanitizeName(allocator, "IGHV3-15*01");
    defer allocator.free(s);
    try std.testing.expectEqualStrings("IGHV3-15_star_01", s);

    const s2 = try GermLines.sanitizeName(allocator, "IGHV3-15*01/OR1");
    defer allocator.free(s2);
    try std.testing.expect(std.mem.indexOf(u8, s2, "_slash_") != null);
}
