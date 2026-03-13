/// ham/bcr_utils/insertions.zig — Zig port of ham::Insertions
///
/// Maps each VDJ region to its associated insertion names.
/// C++ source: packages/ham/include/bcrutils.h

const std = @import("std");

/// Corresponds to C++ `ham::Insertions`.
pub const Insertions = struct {
    /// v region: fv insertion
    v: [1][]const u8 = .{"fv"},
    /// d region: vd insertion
    d: [1][]const u8 = .{"vd"},
    /// j region: dj and jf insertions
    j: [2][]const u8 = .{ "dj", "jf" },

    /// Return a slice of insertion names for a given region ("v", "d", or "j").
    /// Corresponds to C++ `operator[](string region)`.
    pub fn forRegion(self: *const Insertions, region: []const u8) []const []const u8 {
        if (std.mem.eql(u8, region, "v")) return &self.v;
        if (std.mem.eql(u8, region, "d")) return &self.d;
        if (std.mem.eql(u8, region, "j")) return &self.j;
        return &.{};
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Insertions: forRegion" {
    const ins = Insertions{};
    const v_ins = ins.forRegion("v");
    try std.testing.expectEqual(@as(usize, 1), v_ins.len);
    try std.testing.expectEqualStrings("fv", v_ins[0]);

    const j_ins = ins.forRegion("j");
    try std.testing.expectEqual(@as(usize, 2), j_ins.len);
    try std.testing.expectEqualStrings("dj", j_ins[0]);
    try std.testing.expectEqualStrings("jf", j_ins[1]);

    const unknown = ins.forRegion("x");
    try std.testing.expectEqual(@as(usize, 0), unknown.len);
}
