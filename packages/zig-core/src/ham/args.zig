/// ham/args.zig — Zig port of ham/src/args.cc + ham/include/args.h
///
/// Command-line argument holder and CSV input file reader for bcrham.
/// In the C++ code this uses tclap for argument parsing; in Zig we use
/// std.process.args directly (matching the existing main.zig pattern).
///
/// C++ source: packages/ham/src/args.cc, packages/ham/include/args.h
/// C++ author: psathyrella/ham

const std = @import("std");
const ham_text = @import("text.zig");

/// Valid algorithm values (matches C++ `algo_vals_`).
pub const Algorithm = enum { viterbi, forward };

/// Valid loci (matches C++ loci vector).
const valid_loci = [_][]const u8{ "igh", "igk", "igl", "tra", "trb", "trg", "trd" };

/// Column headers from the input CSV file (matches C++ *_headers_ sets).
pub const str_list_headers = [_][]const u8{ "names", "seqs", "only_genes" };
pub const int_headers = [_][]const u8{ "k_v_min", "k_v_max", "k_d_min", "k_d_max", "cdr3_length" };
pub const float_headers = [_][]const u8{"mut_freq"};

/// One row of per-sequence query data (derived from one row of the input CSV).
pub const QueryRow = struct {
    /// Sequence names (one per sequence in this query).
    names: std.ArrayListUnmanaged([]u8),
    /// Sequences (one per sequence in this query).
    seqs: std.ArrayListUnmanaged([]u8),
    /// Gene names to restrict search to (may be empty = no restriction).
    only_genes: std.ArrayListUnmanaged([]u8),
    /// VDJ k-boundaries.
    k_v_min: i32,
    k_v_max: i32,
    k_d_min: i32,
    k_d_max: i32,
    cdr3_length: i32,
    /// Mutation frequency.
    mut_freq: f64,

    pub fn deinit(self: *QueryRow, allocator: std.mem.Allocator) void {
        for (self.names.items) |s| allocator.free(s);
        self.names.deinit(allocator);
        for (self.seqs.items) |s| allocator.free(s);
        self.seqs.deinit(allocator);
        for (self.only_genes.items) |s| allocator.free(s);
        self.only_genes.deinit(allocator);
    }
};

/// Parsed bcrham arguments and CSV data.
/// Corresponds to C++ `ham::Args`.
pub const Args = struct {
    // ── CLI args ──────────────────────────────────────────────────────────────
    hmmdir: []u8,
    datadir: []u8,
    infile: []u8,
    outfile: []u8,
    annotationfile: []u8,
    input_cachefname: []u8,
    output_cachefname: []u8,
    locus: []u8,
    algorithm: []u8,
    ambig_base: []u8,
    seed_unique_id: []u8,

    hamming_fraction_bound_lo: f32,
    hamming_fraction_bound_hi: f32,
    logprob_ratio_threshold: f32,
    max_logprob_drop: f32,

    debug: i32,
    naive_hamming_cluster: i32,
    biggest_naive_seq_cluster_to_calculate: i32,
    biggest_logprob_cluster_to_calculate: i32,
    n_partitions_to_write: i32,

    n_final_clusters: u32,
    min_largest_cluster_size: u32,
    max_cluster_size: u32,
    random_seed: u32,

    no_chunk_cache: bool,
    partition: bool,
    dont_rescale_emissions: bool,
    cache_naive_seqs: bool,
    cache_naive_hfracs: bool,
    only_cache_new_vals: bool,
    write_logprob_for_each_partition: bool,

    // ── CSV data ──────────────────────────────────────────────────────────────
    /// One entry per query row in the input CSV.
    queries: std.ArrayListUnmanaged(QueryRow),

    /// Create an empty Args with default values.
    pub fn initDefaults(allocator: std.mem.Allocator) !Args {
        return Args{
            .hmmdir = try allocator.dupe(u8, ""),
            .datadir = try allocator.dupe(u8, ""),
            .infile = try allocator.dupe(u8, ""),
            .outfile = try allocator.dupe(u8, ""),
            .annotationfile = try allocator.dupe(u8, ""),
            .input_cachefname = try allocator.dupe(u8, ""),
            .output_cachefname = try allocator.dupe(u8, ""),
            .locus = try allocator.dupe(u8, ""),
            .algorithm = try allocator.dupe(u8, ""),
            .ambig_base = try allocator.dupe(u8, ""),
            .seed_unique_id = try allocator.dupe(u8, ""),
            .hamming_fraction_bound_lo = 0.0,
            .hamming_fraction_bound_hi = 1.0,
            .logprob_ratio_threshold = -std.math.inf(f32),
            .max_logprob_drop = -1.0,
            .debug = 0,
            .naive_hamming_cluster = 0,
            .biggest_naive_seq_cluster_to_calculate = 99999,
            .biggest_logprob_cluster_to_calculate = 99999,
            .n_partitions_to_write = 99999,
            .n_final_clusters = 0,
            .min_largest_cluster_size = 0,
            .max_cluster_size = 0,
            .random_seed = @intCast(std.time.timestamp() & 0xFFFFFFFF),
            .no_chunk_cache = false,
            .partition = false,
            .dont_rescale_emissions = false,
            .cache_naive_seqs = false,
            .cache_naive_hfracs = false,
            .only_cache_new_vals = false,
            .write_logprob_for_each_partition = false,
            .queries = .{},
        };
    }

    pub fn deinit(self: *Args, allocator: std.mem.Allocator) void {
        allocator.free(self.hmmdir);
        allocator.free(self.datadir);
        allocator.free(self.infile);
        allocator.free(self.outfile);
        allocator.free(self.annotationfile);
        allocator.free(self.input_cachefname);
        allocator.free(self.output_cachefname);
        allocator.free(self.locus);
        allocator.free(self.algorithm);
        allocator.free(self.ambig_base);
        allocator.free(self.seed_unique_id);
        for (self.queries.items) |*q| q.deinit(allocator);
        self.queries.deinit(allocator);
    }

    /// Parse the CSV input file (infile must already be set).
    /// Corresponds to the file-reading portion of C++ `Args::Args(argc, argv)`.
    pub fn readInfile(self: *Args, allocator: std.mem.Allocator) !void {
        const content = try std.fs.cwd().readFileAlloc(allocator, self.infile, 64 * 1024 * 1024);
        defer allocator.free(content);

        var line_iter = std.mem.splitScalar(u8, content, '\n');

        // Parse header line
        const header_line = line_iter.next() orelse return error.EmptyInputFile;
        var headers: std.ArrayListUnmanaged([]const u8) = .{};
        defer headers.deinit(allocator);
        {
            var tok = std.mem.splitScalar(u8, std.mem.trim(u8, header_line, " \t\r"), ' ');
            while (tok.next()) |h| {
                if (h.len > 0) try headers.append(allocator, h);
            }
        }

        // Parse data rows
        while (line_iter.next()) |raw_line| {
            const line = std.mem.trim(u8, raw_line, " \t\r");
            if (line.len < 10) continue; // skip blank/short lines

            var q = QueryRow{
                .names = .{},
                .seqs = .{},
                .only_genes = .{},
                .k_v_min = 0,
                .k_v_max = 0,
                .k_d_min = 0,
                .k_d_max = 0,
                .cdr3_length = 0,
                .mut_freq = 0.0,
            };
            errdefer q.deinit(allocator);

            var tok = std.mem.splitScalar(u8, line, ' ');
            for (headers.items) |head| {
                const field = tok.next() orelse break;
                if (isStrListHeader(head)) {
                    // Split on ':' — matches C++ SplitString(tmpstr, ":")
                    var parts = try ham_text.splitString(allocator, field, ":");
                    defer {
                        for (parts.items) |p| allocator.free(p);
                        parts.deinit(allocator);
                    }
                    if (std.mem.eql(u8, head, "names")) {
                        for (parts.items) |p| try q.names.append(allocator, try allocator.dupe(u8, p));
                    } else if (std.mem.eql(u8, head, "seqs")) {
                        for (parts.items) |p| {
                            // Strip newlines from each sequence
                            const owned = try allocator.dupe(u8, p);
                            try q.seqs.append(allocator, owned);
                        }
                    } else if (std.mem.eql(u8, head, "only_genes")) {
                        for (parts.items) |p| try q.only_genes.append(allocator, try allocator.dupe(u8, p));
                    }
                } else if (isIntHeader(head)) {
                    const val = try std.fmt.parseInt(i32, field, 10);
                    if (std.mem.eql(u8, head, "k_v_min")) q.k_v_min = val
                    else if (std.mem.eql(u8, head, "k_v_max")) q.k_v_max = val
                    else if (std.mem.eql(u8, head, "k_d_min")) q.k_d_min = val
                    else if (std.mem.eql(u8, head, "k_d_max")) q.k_d_max = val
                    else if (std.mem.eql(u8, head, "cdr3_length")) q.cdr3_length = val;
                } else if (isFloatHeader(head)) {
                    const val = try std.fmt.parseFloat(f64, field);
                    if (std.mem.eql(u8, head, "mut_freq")) q.mut_freq = val;
                } else {
                    return error.UnexpectedHeader;
                }
            }
            try self.queries.append(allocator, q);
        }
    }

    fn isStrListHeader(head: []const u8) bool {
        for (str_list_headers) |h| if (std.mem.eql(u8, h, head)) return true;
        return false;
    }

    fn isIntHeader(head: []const u8) bool {
        for (int_headers) |h| if (std.mem.eql(u8, h, head)) return true;
        return false;
    }

    fn isFloatHeader(head: []const u8) bool {
        for (float_headers) |h| if (std.mem.eql(u8, h, head)) return true;
        return false;
    }

    /// Validate that locus is one of the known valid values.
    pub fn validateLocus(self: *const Args) !void {
        for (valid_loci) |l| if (std.mem.eql(u8, l, self.locus)) return;
        return error.InvalidLocus;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Args: initDefaults" {
    const allocator = std.testing.allocator;
    var args = try Args.initDefaults(allocator);
    defer args.deinit(allocator);

    try std.testing.expectEqualStrings("", args.hmmdir);
    try std.testing.expectEqual(@as(i32, 0), args.debug);
    try std.testing.expectEqual(@as(i32, 99999), args.biggest_naive_seq_cluster_to_calculate);
    try std.testing.expect(!args.partition);
}

test "Args: readInfile with minimal CSV" {
    const allocator = std.testing.allocator;

    // Write a minimal input CSV to a temp file
    const csv = "names seqs only_genes k_v_min k_v_max k_d_min k_d_max cdr3_length mut_freq\n" ++
        "seq1 ACGT gene1:gene2 1 10 1 5 30 0.05\n";

    const tmp = "/tmp/test_args_input.csv";
    {
        const f = try std.fs.createFileAbsolute(tmp, .{});
        defer f.close();
        try f.writeAll(csv);
    }
    defer std.fs.deleteFileAbsolute(tmp) catch {};

    var args = try Args.initDefaults(allocator);
    defer args.deinit(allocator);
    allocator.free(args.infile);
    args.infile = try allocator.dupe(u8, tmp);

    try args.readInfile(allocator);

    try std.testing.expectEqual(@as(usize, 1), args.queries.items.len);
    const q = &args.queries.items[0];
    try std.testing.expectEqual(@as(usize, 1), q.names.items.len);
    try std.testing.expectEqualStrings("seq1", q.names.items[0]);
    try std.testing.expectEqualStrings("ACGT", q.seqs.items[0]);
    try std.testing.expectEqual(@as(usize, 2), q.only_genes.items.len);
    try std.testing.expectEqual(@as(i32, 1), q.k_v_min);
    try std.testing.expectApproxEqAbs(0.05, q.mut_freq, 1e-9);
}
