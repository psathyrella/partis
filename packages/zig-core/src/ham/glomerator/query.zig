/// ham/glomerator/query.zig — Zig port of ham::Query
///
/// A named cluster with associated sequences, kbounds, and gene lists.
/// C++ source: packages/ham/include/glomerator.h

const std = @import("std");
const Sequence = @import("../sequences.zig").Sequence;
const KBounds = @import("../bcr_utils/k_bounds.zig").KBounds;

/// Corresponds to C++ `ham::Query`.
///
/// Sequence ownership (issue #342, item 6): `seqs` holds **borrowed** Sequence
/// values — shallow copies whose slice fields (`name`, `header`, `undigitized`,
/// `seqq`) point at buffers owned by `Glomerator.single_seqs`. The Query must
/// not outlive the Glomerator that owns those buffers, and `Query.deinit` must
/// not call `Sequence.deinit` on its entries.
pub const Query = struct {
    name: []u8,
    /// Borrowed shallow copies of Sequences. Buffers are owned by
    /// `Glomerator.single_seqs`; do not free per-item in deinit.
    seqs: std.ArrayListUnmanaged(Sequence),
    seed_missing: bool,
    only_genes: std.ArrayListUnmanaged([]u8),
    kbounds: KBounds,
    mute_freq: f32,
    cdr3_length: usize,
    /// Names of the two parent queries (may be empty when not set).
    parents: ?[2][]u8,

    pub fn init(allocator: std.mem.Allocator) !Query {
        _ = allocator;
        return Query{
            .name = &.{},
            .seqs = .{},
            .seed_missing = false,
            .only_genes = .{},
            .kbounds = .{},
            .mute_freq = 0.0,
            .cdr3_length = 0,
            .parents = null,
        };
    }

    /// Create a fully-specified Query. `seqs` is borrowed: each Sequence is
    /// stored as a shallow copy (slice descriptors), and the buffers must
    /// outlive this Query. See struct doc comment.
    pub fn create(
        allocator: std.mem.Allocator,
        name: []const u8,
        seqs: []const Sequence,
        seed_missing: bool,
        only_genes: []const []const u8,
        kbounds: KBounds,
        mute_freq: f32,
        cdr3_length: usize,
        parent1: ?[]const u8,
        parent2: ?[]const u8,
    ) !Query {
        var q = Query{
            .name = try allocator.dupe(u8, name),
            .seqs = .{},
            .seed_missing = seed_missing,
            .only_genes = .{},
            .kbounds = kbounds,
            .mute_freq = mute_freq,
            .cdr3_length = cdr3_length,
            .parents = null,
        };
        errdefer q.deinit(allocator);

        try q.seqs.ensureUnusedCapacity(allocator, seqs.len);
        for (seqs) |sq| {
            q.seqs.appendAssumeCapacity(sq);
        }
        for (only_genes) |g| {
            try q.only_genes.append(allocator, try allocator.dupe(u8, g));
        }
        if (parent1 != null and parent2 != null) {
            q.parents = .{
                try allocator.dupe(u8, parent1.?),
                try allocator.dupe(u8, parent2.?),
            };
        }
        return q;
    }

    pub fn deinit(self: *Query, allocator: std.mem.Allocator) void {
        if (self.name.len > 0) allocator.free(self.name);
        // seqs entries are borrowed shallow copies (see struct doc); do not
        // call Sequence.deinit here. Only the ArrayList itself is owned.
        self.seqs.deinit(allocator);
        for (self.only_genes.items) |g| allocator.free(g);
        self.only_genes.deinit(allocator);
        if (self.parents) |*ps| {
            allocator.free(ps[0]);
            allocator.free(ps[1]);
        }
    }

    /// Return number of sequences.
    pub fn nSeqs(self: *const Query) usize {
        return self.seqs.items.len;
    }
};
