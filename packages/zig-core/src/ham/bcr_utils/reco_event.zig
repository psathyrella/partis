/// ham/bcr_utils/reco_event.zig — Zig port of ham::RecoEvent
///
/// A single VDJ recombination event annotation.
/// C++ source: packages/ham/src/bcrutils.cc

const std = @import("std");
const SupportPair = @import("support_pair.zig").SupportPair;
const GermLines = @import("germ_lines.zig").GermLines;

/// Corresponds to C++ `ham::RecoEvent`.
pub const RecoEvent = struct {
    /// Region gene assignments: "v", "d", "j" → gene name.
    genes: std.StringHashMapUnmanaged([]u8),
    /// Deletion lengths: "v_5p", "v_3p", "d_5p", "d_3p", "j_5p", "j_3p".
    deletions: std.StringHashMapUnmanaged(usize),
    /// Insertion sequences: "fv", "vd", "dj", "jf".
    insertions: std.StringHashMapUnmanaged([]u8),
    naive_seq: []u8,
    score: f32,
    cyst_position: i32,
    tryp_position: i32,
    cdr3_length: i32,
    /// Per-region sorted (gene, logprob) support.
    per_gene_support: std.StringHashMapUnmanaged(std.ArrayListUnmanaged(SupportPair)),

    /// Create an empty RecoEvent (score = 999, as in C++).
    pub fn init(allocator: std.mem.Allocator) !RecoEvent {
        _ = allocator;
        return RecoEvent{
            .genes = .{},
            .deletions = .{},
            .insertions = .{},
            .naive_seq = &.{},
            .score = 999,
            .cyst_position = 0,
            .tryp_position = 0,
            .cdr3_length = 0,
            .per_gene_support = .{},
        };
    }

    pub fn deinit(self: *RecoEvent, allocator: std.mem.Allocator) void {
        var git = self.genes.iterator();
        while (git.next()) |e| {
            allocator.free(e.key_ptr.*);
            allocator.free(e.value_ptr.*);
        }
        self.genes.deinit(allocator);

        var dit = self.deletions.iterator();
        while (dit.next()) |e| allocator.free(e.key_ptr.*);
        self.deletions.deinit(allocator);

        var iit = self.insertions.iterator();
        while (iit.next()) |e| {
            allocator.free(e.key_ptr.*);
            allocator.free(e.value_ptr.*);
        }
        self.insertions.deinit(allocator);

        if (self.naive_seq.len > 0) allocator.free(self.naive_seq);

        var pit = self.per_gene_support.iterator();
        while (pit.next()) |e| {
            allocator.free(e.key_ptr.*);
            e.value_ptr.deinit(allocator);
        }
        self.per_gene_support.deinit(allocator);
    }

    pub fn setGenes(self: *RecoEvent, allocator: std.mem.Allocator, vgene: []const u8, dgene: []const u8, jgene: []const u8) !void {
        try self.setGene(allocator, "v", vgene);
        try self.setGene(allocator, "d", dgene);
        try self.setGene(allocator, "j", jgene);
    }

    pub fn setGene(self: *RecoEvent, allocator: std.mem.Allocator, region: []const u8, gene: []const u8) !void {
        const key = try allocator.dupe(u8, region);
        errdefer allocator.free(key);
        const val = try allocator.dupe(u8, gene);
        errdefer allocator.free(val);
        // Free old entry if present
        if (self.genes.fetchRemove(key)) |old| {
            allocator.free(old.key);
            allocator.free(old.value);
        }
        try self.genes.put(allocator, key, val);
    }

    pub fn setDeletion(self: *RecoEvent, allocator: std.mem.Allocator, name: []const u8, len: usize) !void {
        const key = try allocator.dupe(u8, name);
        errdefer allocator.free(key);
        if (self.deletions.fetchRemove(key)) |old| allocator.free(old.key);
        try self.deletions.put(allocator, key, len);
    }

    pub fn setInsertion(self: *RecoEvent, allocator: std.mem.Allocator, name: []const u8, seq: []const u8) !void {
        const key = try allocator.dupe(u8, name);
        errdefer allocator.free(key);
        const val = try allocator.dupe(u8, seq);
        errdefer allocator.free(val);
        if (self.insertions.fetchRemove(key)) |old| {
            allocator.free(old.key);
            allocator.free(old.value);
        }
        try self.insertions.put(allocator, key, val);
    }

    pub fn setScore(self: *RecoEvent, s: f64) void {
        self.score = @floatCast(s);
    }

    pub fn clear(self: *RecoEvent, allocator: std.mem.Allocator) void {
        var git = self.genes.iterator();
        while (git.next()) |e| {
            allocator.free(e.key_ptr.*);
            allocator.free(e.value_ptr.*);
        }
        self.genes.clearRetainingCapacity();

        var dit = self.deletions.iterator();
        while (dit.next()) |e| allocator.free(e.key_ptr.*);
        self.deletions.clearRetainingCapacity();

        var iit = self.insertions.iterator();
        while (iit.next()) |e| {
            allocator.free(e.key_ptr.*);
            allocator.free(e.value_ptr.*);
        }
        self.insertions.clearRetainingCapacity();
    }

    /// Less-than by score (for sorting).
    pub fn lessThan(_: void, lhs: RecoEvent, rhs: RecoEvent) bool {
        return lhs.score < rhs.score;
    }

    /// Build naive_seq from germline sequences and indel lengths.
    /// Corresponds to C++ `RecoEvent::SetNaiveSeq`.
    pub fn setNaiveSeq(self: *RecoEvent, allocator: std.mem.Allocator, gl: *const GermLines) !void {
        const vgene = self.genes.get("v") orelse return error.MissingVGene;
        const dgene = self.genes.get("d") orelse return error.MissingDGene;
        const jgene = self.genes.get("j") orelse return error.MissingJGene;

        const v_seq = gl.seqs.get(vgene) orelse return error.MissingSeq;
        const d_seq = gl.seqs.get(dgene) orelse return error.MissingSeq;
        const j_seq = gl.seqs.get(jgene) orelse return error.MissingSeq;

        const v_5p: usize = self.deletions.get("v_5p") orelse 0;
        const v_3p: usize = self.deletions.get("v_3p") orelse 0;
        const d_5p: usize = self.deletions.get("d_5p") orelse 0;
        const d_3p: usize = self.deletions.get("d_3p") orelse 0;
        const j_5p: usize = self.deletions.get("j_5p") orelse 0;
        const j_3p: usize = self.deletions.get("j_3p") orelse 0;

        const fv = self.insertions.get("fv") orelse "";
        const vd = self.insertions.get("vd") orelse "";
        const dj = self.insertions.get("dj") orelse "";
        const jf = self.insertions.get("jf") orelse "";

        const eroded_v = if (v_5p + v_3p <= v_seq.len) v_seq[v_5p .. v_seq.len - v_3p] else "";
        const eroded_d = if (d_5p + d_3p <= d_seq.len) d_seq[d_5p .. d_seq.len - d_3p] else "";
        const eroded_j = if (j_5p + j_3p <= j_seq.len) j_seq[j_5p .. j_seq.len - j_3p] else "";

        if (self.naive_seq.len > 0) allocator.free(self.naive_seq);
        self.naive_seq = try std.fmt.allocPrint(allocator, "{s}{s}{s}{s}{s}{s}{s}", .{
            fv, eroded_v, vd, eroded_d, dj, eroded_j, jf,
        });

        // Set CDR3 positions
        const v_cyst = gl.cyst_positions.get(vgene) orelse 0;
        const j_tryp = gl.tryp_positions.get(jgene) orelse 0;
        const eroded_cpos: i32 = v_cyst - @as(i32, @intCast(v_5p)) + @as(i32, @intCast(fv.len));
        const tpos_in_joined: i32 = j_tryp - @as(i32, @intCast(j_5p)) +
            @as(i32, @intCast(fv.len + eroded_v.len + vd.len + eroded_d.len + dj.len));
        self.cyst_position = eroded_cpos;
        self.tryp_position = tpos_in_joined;
        self.cdr3_length = tpos_in_joined - eroded_cpos + 3;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "RecoEvent: init, setGene, setDeletion, clear" {
    const allocator = std.testing.allocator;
    var ev = try RecoEvent.init(allocator);
    defer ev.deinit(allocator);

    try ev.setGene(allocator, "v", "IGHV3-15*01");
    try ev.setDeletion(allocator, "v_3p", 2);
    try std.testing.expectEqualStrings("IGHV3-15*01", ev.genes.get("v").?);
    try std.testing.expectEqual(@as(usize, 2), ev.deletions.get("v_3p").?);

    ev.clear(allocator);
    try std.testing.expectEqual(@as(usize, 0), ev.genes.count());
}
