/// ham/bcr_utils/k_bounds.zig — Zig port of ham::KSet and ham::KBounds
///
/// k_v / k_d cut-point pairs and their bounding box.
/// C++ source: packages/ham/include/bcrutils.h

const std = @import("std");

/// A pair of (k_v, k_d) values specifying how to chop a query sequence.
/// Corresponds to C++ `ham::KSet`.
pub const KSet = struct {
    v: usize,
    d: usize,

    pub fn isNull(self: KSet) bool {
        return self.v == 0 and self.d == 0;
    }

    pub fn eql(self: KSet, rhs: KSet) bool {
        return self.v == rhs.v and self.d == rhs.d;
    }

    /// Less-than for use as map key (v first, then d).
    pub fn lessThan(_: void, lhs: KSet, rhs: KSet) bool {
        if (lhs.v != rhs.v) return lhs.v < rhs.v;
        return lhs.d < rhs.d;
    }
};

/// Bounding box over KSet space.
/// Corresponds to C++ `ham::KBounds`.
pub const KBounds = struct {
    vmin: usize = 0,
    dmin: usize = 0,
    vmax: usize = 0,
    dmax: usize = 0,

    pub fn initFromSets(kmin: KSet, kmax: KSet) KBounds {
        return KBounds{
            .vmin = kmin.v,
            .dmin = kmin.d,
            .vmax = kmax.v,
            .dmax = kmax.d,
        };
    }

    pub fn eql(self: KBounds, rhs: KBounds) bool {
        return self.vmin == rhs.vmin and self.vmax == rhs.vmax and
            self.dmin == rhs.dmin and self.dmax == rhs.dmax;
    }

    /// Return the union of `self` and `rhs`.
    /// Corresponds to C++ `KBounds::LogicalOr`.
    pub fn logicalOr(self: KBounds, rhs: KBounds) KBounds {
        return KBounds{
            .vmin = @min(self.vmin, rhs.vmin),
            .dmin = @min(self.dmin, rhs.dmin),
            .vmax = @max(self.vmax, rhs.vmax),
            .dmax = @max(self.dmax, rhs.dmax),
        };
    }

    pub fn stringify(self: KBounds, allocator: std.mem.Allocator) ![]u8 {
        return std.fmt.allocPrint(allocator, "{d}-{d}, {d}-{d}", .{
            self.vmin, self.vmax, self.dmin, self.dmax,
        });
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "KSet: isNull and eql" {
    const null_ks = KSet{ .v = 0, .d = 0 };
    const real_ks = KSet{ .v = 3, .d = 2 };
    try std.testing.expect(null_ks.isNull());
    try std.testing.expect(!real_ks.isNull());
    try std.testing.expect(real_ks.eql(KSet{ .v = 3, .d = 2 }));
    try std.testing.expect(!real_ks.eql(null_ks));
}

test "KBounds: logicalOr" {
    const a = KBounds{ .vmin = 1, .vmax = 5, .dmin = 1, .dmax = 3 };
    const b = KBounds{ .vmin = 2, .vmax = 7, .dmin = 0, .dmax = 4 };
    const u = a.logicalOr(b);
    try std.testing.expectEqual(@as(usize, 1), u.vmin);
    try std.testing.expectEqual(@as(usize, 7), u.vmax);
    try std.testing.expectEqual(@as(usize, 0), u.dmin);
    try std.testing.expectEqual(@as(usize, 4), u.dmax);
}
