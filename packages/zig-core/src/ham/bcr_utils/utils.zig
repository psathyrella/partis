/// ham/bcr_utils/utils.zig — Free functions from ham::bcrutils
///
/// CSV/output streaming, sequence utilities, memory query helpers.
/// C++ source: packages/ham/src/bcrutils.cc

const std = @import("std");
const c = @cImport(@cInclude("stdio.h"));
const SupportPair = @import("support_pair.zig").SupportPair;
const RecoEvent = @import("reco_event.zig").RecoEvent;
const Sequence = @import("../sequences.zig").Sequence;

/// Whether the locus has a D gene segment (re-exported for use outside this module).
pub const hasDGene = @import("germ_lines.zig").hasDGene;

// ── C++-compatible float formatting ─────────────────────────────────────────

/// Format a float like C++ `ostream << double` (6 significant digits, `%.6g`).
fn fmtG6(buf: *[32]u8, val: f64) []const u8 {
    const n: usize = @intCast(c.snprintf(buf, 32, "%.6g", val));
    return buf[0..n];
}

/// Format a float like C++ `std::to_string(double)` (6 decimal places, `%.6f`).
fn fmtF6(buf: *[64]u8, val: f64) []const u8 {
    const n: usize = @intCast(c.snprintf(buf, 64, "%.6f", val));
    return buf[0..n];
}

// ── Output streaming ──────────────────────────────────────────────────────────

/// Write CSV header line for the given algorithm.
/// Corresponds to C++ `ham::StreamHeader`.
pub fn streamHeader(writer: anytype, algorithm: []const u8) !void {
    if (std.mem.eql(u8, algorithm, "viterbi")) {
        try writer.writeAll(
            "unique_ids,v_gene,d_gene,j_gene,fv_insertion,vd_insertion,dj_insertion,jf_insertion," ++
                "v_5p_del,v_3p_del,d_5p_del,d_3p_del,j_5p_del,j_3p_del,logprob,seqs," ++
                "v_per_gene_support,d_per_gene_support,j_per_gene_support,errors\n",
        );
    } else if (std.mem.eql(u8, algorithm, "forward")) {
        try writer.writeAll("unique_ids,logprob,errors\n");
    } else {
        return error.BadAlgorithm;
    }
}

/// Build per-gene support string: "gene:logprob;gene:logprob;..."
/// Corresponds to C++ `ham::PerGeneSupportString`.
pub fn perGeneSupportString(
    allocator: std.mem.Allocator,
    support: []const SupportPair,
) ![]u8 {
    var buf: std.ArrayListUnmanaged(u8) = .{};
    defer buf.deinit(allocator);
    for (support, 0..) |sp, i| {
        if (i > 0) try buf.append(allocator, ';');
        var fbuf: [64]u8 = undefined;
        const fstr = fmtF6(&fbuf, sp.logprob);
        try buf.appendSlice(allocator, sp.gene);
        try buf.append(allocator, ':');
        try buf.appendSlice(allocator, fstr);
    }
    return buf.toOwnedSlice(allocator);
}

/// Write error output row to the CSV file.
/// Corresponds to C++ `ham::StreamErrorput`.
pub fn streamErrorput(
    writer: anytype,
    algorithm: []const u8,
    allocator: std.mem.Allocator,
    seqs: []const Sequence,
    errors: []const u8,
) !void {
    const names = try seqNameStr(allocator, seqs, ":");
    defer allocator.free(names);
    const seq_str = try seqStr(allocator, seqs, ":");
    defer allocator.free(seq_str);

    if (std.mem.eql(u8, algorithm, "viterbi")) {
        try writer.print("{s},,,,,,,,,,,,,,,{s},,,,{s}\n", .{ names, seq_str, errors });
    } else {
        try writer.print("{s},,{s}\n", .{ names, errors });
    }
}

/// Write Viterbi output row to the CSV file.
/// Corresponds to C++ `ham::StreamViterbiOutput`.
pub fn streamViterbiOutput(
    writer: anytype,
    allocator: std.mem.Allocator,
    event: *const RecoEvent,
    seqs: []const Sequence,
    errors: []const u8,
) !void {
    const names = try seqNameStr(allocator, seqs, ":");
    defer allocator.free(names);
    const seq_str = try seqStr(allocator, seqs, ":");
    defer allocator.free(seq_str);

    const vgene = event.genes.get("v") orelse "";
    const dgene = event.genes.get("d") orelse "";
    const jgene = event.genes.get("j") orelse "";
    const fv_ins = event.insertions.get("fv") orelse "";
    const vd_ins = event.insertions.get("vd") orelse "";
    const dj_ins = event.insertions.get("dj") orelse "";
    const jf_ins = event.insertions.get("jf") orelse "";
    const v5p = event.deletions.get("v_5p") orelse 0;
    const v3p = event.deletions.get("v_3p") orelse 0;
    const d5p = event.deletions.get("d_5p") orelse 0;
    const d3p = event.deletions.get("d_3p") orelse 0;
    const j5p = event.deletions.get("j_5p") orelse 0;
    const j3p = event.deletions.get("j_3p") orelse 0;

    const empty_support: []const SupportPair = &.{};
    const vsupp = if (event.per_gene_support.get("v")) |s| s.items else empty_support;
    const dsupp = if (event.per_gene_support.get("d")) |s| s.items else empty_support;
    const jsupp = if (event.per_gene_support.get("j")) |s| s.items else empty_support;

    const vs = try perGeneSupportString(allocator, vsupp);
    defer allocator.free(vs);
    const ds = try perGeneSupportString(allocator, dsupp);
    defer allocator.free(ds);
    const js = try perGeneSupportString(allocator, jsupp);
    defer allocator.free(js);

    var score_buf: [32]u8 = undefined;
    const score_str = fmtG6(&score_buf, event.score);
    try writer.print("{s},{s},{s},{s},{s},{s},{s},{s},{d},{d},{d},{d},{d},{d},{s},{s},{s},{s},{s},{s}\n", .{
        names, vgene, dgene, jgene,
        fv_ins, vd_ins, dj_ins, jf_ins,
        v5p, v3p, d5p, d3p, j5p, j3p,
        score_str, seq_str, vs, ds, js, errors,
    });
}

/// Write Forward output row.
/// Corresponds to C++ `ham::StreamForwardOutput`.
pub fn streamForwardOutput(
    writer: anytype,
    allocator: std.mem.Allocator,
    seqs: []const Sequence,
    total_score: f64,
    errors: []const u8,
) !void {
    const names = try seqNameStr(allocator, seqs, ":");
    defer allocator.free(names);
    var score_buf: [32]u8 = undefined;
    const score_str = fmtG6(&score_buf, total_score);
    try writer.print("{s},{s},{s}\n", .{ names, score_str, errors });
}

// ── Sequence string helpers ────────────────────────────────────────────────────

/// Join undigitized sequences with a delimiter.
/// Corresponds to C++ `ham::SeqStr`.
pub fn seqStr(
    allocator: std.mem.Allocator,
    seqs: []const Sequence,
    delimiter: []const u8,
) ![]u8 {
    var buf: std.ArrayListUnmanaged(u8) = .{};
    defer buf.deinit(allocator);
    for (seqs, 0..) |*seq, i| {
        if (i > 0) try buf.appendSlice(allocator, delimiter);
        try buf.appendSlice(allocator, seq.undigitized);
    }
    return buf.toOwnedSlice(allocator);
}

/// Join sequence names with a delimiter.
/// Corresponds to C++ `ham::SeqNameStr`.
pub fn seqNameStr(
    allocator: std.mem.Allocator,
    seqs: []const Sequence,
    delimiter: []const u8,
) ![]u8 {
    var buf: std.ArrayListUnmanaged(u8) = .{};
    defer buf.deinit(allocator);
    for (seqs, 0..) |*seq, i| {
        if (i > 0) try buf.appendSlice(allocator, delimiter);
        try buf.appendSlice(allocator, seq.name);
    }
    return buf.toOwnedSlice(allocator);
}

// ── Multi-sequence annotation quality check ────────────────────────────────────

/// Return true if the annotation looks unreliable for multi-seq input.
/// Corresponds to C++ `ham::FishyMultiSeqAnnotation`.
pub fn fishyMultiSeqAnnotation(n_seqs: usize, event: *const RecoEvent) bool {
    if (n_seqs < 3) return false;
    const real_deletions = [_][]const u8{ "v_3p", "d_5p", "d_3p", "j_5p" };
    for (real_deletions) |delname| {
        if ((event.deletions.get(delname) orelse 0) > 2) return true;
    }
    return false;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

test "streamHeader: viterbi and forward" {
    const allocator = std.testing.allocator;
    var buf: std.ArrayListUnmanaged(u8) = .{};
    defer buf.deinit(allocator);
    const writer = buf.writer(allocator);
    try streamHeader(writer, "viterbi");
    try std.testing.expect(std.mem.startsWith(u8, buf.items, "unique_ids,v_gene"));

    buf.clearRetainingCapacity();
    try streamHeader(writer, "forward");
    try std.testing.expect(std.mem.startsWith(u8, buf.items, "unique_ids,logprob"));
}

test "perGeneSupportString" {
    const allocator = std.testing.allocator;
    const support = [_]SupportPair{
        .{ .gene = "IGHV3-15*01", .logprob = -5.0 },
        .{ .gene = "IGHV3-23*01", .logprob = -3.0 },
    };
    const s = try perGeneSupportString(allocator, &support);
    defer allocator.free(s);
    try std.testing.expect(std.mem.indexOf(u8, s, "IGHV3-15*01") != null);
    try std.testing.expect(std.mem.indexOf(u8, s, ";") != null);
}
