/// ham/bcr_utils/term_colors.zig — Zig port of ham::TermColors
///
/// ANSI terminal color helpers used for debug printing in bcrham.
/// C++ source: packages/ham/src/bcrutils.cc, packages/ham/include/bcrutils.h

const std = @import("std");
const Region = @import("germ_lines.zig").Region;
const regionStr = @import("germ_lines.zig").regionStr;

/// ANSI escape code table.  Corresponds to C++ `TermColors::codes_`.
const color_codes = std.StaticStringMap([]const u8).initComptime(.{
    .{ "bold", "\x1b[1m" },
    .{ "reverse", "\x1b[7m" },
    .{ "purple", "\x1b[95m" },
    .{ "blue", "\x1b[94m" },
    .{ "light_blue", "\x1b[1;34m" },
    .{ "green", "\x1b[92m" },
    .{ "yellow", "\x1b[93m" },
    .{ "red", "\x1b[91m" },
    .{ "end", "\x1b[0m" },
});

/// Corresponds to C++ `ham::TermColors`.
pub const TermColors = struct {
    /// Wrap `seq` with ANSI color escape codes.
    /// Caller owns the returned slice (allocated with `allocator`).
    /// Corresponds to C++ `TermColors::Color(col, seq)`.
    pub fn color(allocator: std.mem.Allocator, col: []const u8, seq: []const u8) ![]u8 {
        const start = color_codes.get(col) orelse return error.UnknownColor;
        const end_code = color_codes.get("end").?;
        return std.fmt.allocPrint(allocator, "{s}{s}{s}", .{ start, seq, end_code });
    }

    /// Return `nuc` wrapped in red if it differs from `germline_nuc`.
    /// Corresponds to C++ `TermColors::RedifyIfMuted`.
    pub fn redifyIfMuted(allocator: std.mem.Allocator, germline_nuc: u8, nuc: u8) ![]u8 {
        if (nuc == germline_nuc) {
            const s = try allocator.alloc(u8, 1);
            s[0] = nuc;
            return s;
        }
        const single = [_]u8{nuc};
        return color(allocator, "red", &single);
    }

    /// Extract region from a gene name (e.g. "IGHD1-1*01" → .d).
    /// Corresponds to C++ `TermColors::GetRegion(gene)`.
    pub fn getRegion(gene: []const u8) !Region {
        if (!std.mem.startsWith(u8, gene, "IG") and !std.mem.startsWith(u8, gene, "TR"))
            return error.BadGeneName;
        if (gene.len < 4) return error.BadGeneName;
        const ch = std.ascii.toLower(gene[3]);
        return switch (ch) {
            'v' => .v,
            'd' => .d,
            'j' => .j,
            else => error.BadRegion,
        };
    }

    /// Return `seq` with every occurrence of `ch_to_color` colored.
    /// Corresponds to C++ `TermColors::ColorChars`.
    pub fn colorChars(allocator: std.mem.Allocator, ch_to_color: u8, col: []const u8, seq: []const u8) ![]u8 {
        var buf: std.ArrayListUnmanaged(u8) = .{};
        defer buf.deinit(allocator);
        for (seq) |c| {
            if (c == ch_to_color) {
                const colored = try color(allocator, col, &[_]u8{c});
                defer allocator.free(colored);
                try buf.appendSlice(allocator, colored);
            } else {
                try buf.append(allocator, c);
            }
        }
        return buf.toOwnedSlice(allocator);
    }

    /// Colorize a gene name (e.g. "IGHV3-15*01") with region, version, allele parts.
    /// Caller owns the returned slice (allocated with `allocator`).
    /// Corresponds to C++ `TermColors::ColorGene(gene)`.
    pub fn colorGene(allocator: std.mem.Allocator, gene: []const u8) ![]u8 {
        const region = try getRegion(gene);
        const rs = regionStr(region);
        const region_colored = try color(allocator, "red", rs);
        defer allocator.free(region_colored);

        var buf: std.ArrayListUnmanaged(u8) = .{};
        defer buf.deinit(allocator);
        try buf.appendSlice(allocator, region_colored);

        // Find positions of '-' and '*' and '_' in the gene name
        const dash_pos = std.mem.indexOfScalar(u8, gene, '-');
        const star_pos = std.mem.indexOfScalar(u8, gene, '*');
        const underscore_pos = std.mem.indexOfScalar(u8, gene, '_');

        if (region == .j) {
            // For j genes: n_version = gene[4..star_pos], no subversion
            const version_end = star_pos orelse gene.len;
            const n_version = gene[4..version_end];
            const version_colored = try color(allocator, "purple", n_version);
            defer allocator.free(version_colored);
            try buf.appendSlice(allocator, version_colored);
        } else {
            // For v and d genes: n_version = gene[4..dash_pos], n_subversion = gene[dash_pos+1..star_pos]
            const d_pos = dash_pos orelse gene.len;
            const n_version = gene[4..d_pos];
            const version_colored = try color(allocator, "purple", n_version);
            defer allocator.free(version_colored);
            try buf.appendSlice(allocator, version_colored);
            try buf.append(allocator, '-');
            const s_pos = star_pos orelse gene.len;
            if (d_pos + 1 <= s_pos) {
                const n_subversion = gene[d_pos + 1 .. s_pos];
                const subversion_colored = try color(allocator, "purple", n_subversion);
                defer allocator.free(subversion_colored);
                try buf.appendSlice(allocator, subversion_colored);
            }
        }

        // Allele: gene[star_pos+1..underscore_pos]
        if (star_pos) |sp| {
            const allele_end = underscore_pos orelse gene.len;
            const allele = gene[sp + 1 .. allele_end];
            const allele_colored = try color(allocator, "yellow", allele);
            defer allocator.free(allele_colored);
            try buf.appendSlice(allocator, allele_colored);
        }

        // Suffix after '_' (if present)
        if (underscore_pos) |up| {
            try buf.appendSlice(allocator, gene[up..]);
        }

        return buf.toOwnedSlice(allocator);
    }

    /// Return `seq` with mutant positions highlighted.
    /// Corresponds to C++ `TermColors::ColorMutants`.
    pub fn colorMutants(
        allocator: std.mem.Allocator,
        col: []const u8,
        seq: []const u8,
        ref1: []const u8,
        other_refs: []const []const u8,
        ambiguous_char: []const u8,
    ) ![]u8 {
        // Build full ref list: other_refs + ref1 (if non-empty)
        var refs: std.ArrayListUnmanaged([]const u8) = .{};
        defer refs.deinit(allocator);
        try refs.appendSlice(allocator, other_refs);
        if (ref1.len > 0) try refs.append(allocator, ref1);

        var buf: std.ArrayListUnmanaged(u8) = .{};
        defer buf.deinit(allocator);

        for (seq, 0..) |nuc, inuc| {
            if (nuc == 'i') {
                const colored = try color(allocator, "yellow", &[_]u8{nuc});
                defer allocator.free(colored);
                try buf.appendSlice(allocator, colored);
            } else {
                var ndiff: usize = 0;
                for (refs.items) |ref| {
                    if (inuc >= ref.len) continue;
                    const ambig = ambiguous_char.len > 0 and (ref[inuc] == ambiguous_char[0] or nuc == ambiguous_char[0]);
                    if (!ambig and nuc != ref[inuc]) ndiff += 1;
                }
                if (ndiff == 0) {
                    try buf.append(allocator, nuc);
                } else if (ndiff == 1) {
                    const colored = try color(allocator, col, &[_]u8{nuc});
                    defer allocator.free(colored);
                    try buf.appendSlice(allocator, colored);
                } else {
                    const inner = try color(allocator, col, &[_]u8{nuc});
                    defer allocator.free(inner);
                    const outer = try color(allocator, "reverse", inner);
                    defer allocator.free(outer);
                    try buf.appendSlice(allocator, outer);
                }
            }
        }
        return buf.toOwnedSlice(allocator);
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "TermColors: getRegion" {
    try std.testing.expectEqual(Region.v, try TermColors.getRegion("IGHV3-15*01"));
    try std.testing.expectEqual(Region.d, try TermColors.getRegion("IGHD1-1*01"));
    try std.testing.expectEqual(Region.j, try TermColors.getRegion("IGHJ4*02"));
    try std.testing.expectError(error.BadGeneName, TermColors.getRegion("bad"));
}

test "TermColors: color wraps with ANSI codes" {
    const allocator = std.testing.allocator;
    const s = try TermColors.color(allocator, "red", "ABC");
    defer allocator.free(s);
    try std.testing.expect(std.mem.startsWith(u8, s, "\x1b[91m"));
    try std.testing.expect(std.mem.endsWith(u8, s, "\x1b[0m"));
    try std.testing.expect(std.mem.indexOf(u8, s, "ABC") != null);
}
