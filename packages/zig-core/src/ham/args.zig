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

    hamming_fraction_bound_lo: f64,
    hamming_fraction_bound_hi: f64,
    logprob_ratio_threshold: f64,
    max_logprob_drop: f64,

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

    max_infile_bytes: usize,

    // ── CSV data ──────────────────────────────────────────────────────────────
    /// One entry per query row in the input CSV.
    queries: std.ArrayListUnmanaged(QueryRow),

    pub const default_max_infile_bytes: usize = 4 * 1024 * 1024 * 1024; // 4 GiB

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
            .logprob_ratio_threshold = -std.math.inf(f64),
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
            .max_infile_bytes = default_max_infile_bytes,
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

    pub var initial_read_buffer_size: usize = 2 * 1024 * 1024; // 2 MiB initial streaming buffer (dynamically grows on StreamTooLong)

    /// Read the next line from `file_reader`, dynamically growing `buffer` if the line exceeds
    /// the current buffer capacity, up to `max_line_bytes`.
    fn readNextLine(file_reader: *std.fs.File.Reader, buffer: *[]u8, allocator: std.mem.Allocator, max_line_bytes: usize, filename: []const u8) !?[]const u8 {
        while (true) {
            if (file_reader.interface.takeDelimiter('\n')) |maybe_line| {
                return maybe_line;
            } else |err| switch (err) {
                error.StreamTooLong => {
                    if (buffer.*.len >= max_line_bytes) {
                        std.debug.print("error: line in '{s}' ({d} bytes) exceeds maximum supported line length ({d} bytes)\n", .{
                            filename,
                            buffer.*.len,
                            max_line_bytes,
                        });
                        return error.StreamTooLong;
                    }
                    var new_size = buffer.*.len * 2;
                    if (new_size > max_line_bytes) new_size = max_line_bytes;

                    const unconsumed_len = file_reader.interface.end - file_reader.interface.seek;
                    const new_buf = try allocator.alloc(u8, new_size);
                    @memcpy(new_buf[0..unconsumed_len], file_reader.interface.buffer[file_reader.interface.seek..file_reader.interface.end]);
                    allocator.free(buffer.*);
                    buffer.* = new_buf;

                    file_reader.interface.buffer = buffer.*;
                    file_reader.interface.seek = 0;
                    file_reader.interface.end = unconsumed_len;
                },
                else => return err,
            }
        }
    }

    /// Parse the CSV input file (infile must already be set).
    /// Corresponds to the file-reading portion of C++ `Args::Args(argc, argv)`.
    pub fn readInfile(self: *Args, allocator: std.mem.Allocator) !void {
        const file = std.fs.cwd().openFile(self.infile, .{}) catch |err| {
            std.debug.print("error: unable to open input file '{s}': {}\n", .{ self.infile, err });
            return err;
        };
        defer file.close();

        // Check file size against max_infile_bytes up front
        const file_stat = file.stat() catch |err| {
            std.debug.print("error: unable to stat input file '{s}': {}\n", .{ self.infile, err });
            return err;
        };
        if (file_stat.size > self.max_infile_bytes) {
            ham_text.printFileTooBig("input file", self.infile, file_stat.size, self.max_infile_bytes);
            return error.FileTooBig;
        }

        var read_buf = try allocator.alloc(u8, initial_read_buffer_size);
        defer allocator.free(read_buf);

        var file_reader = file.readerStreaming(read_buf);

        // Parse header line
        const header_raw = (try readNextLine(&file_reader, &read_buf, allocator, self.max_infile_bytes, self.infile)) orelse return error.EmptyInputFile;

        var headers: std.ArrayListUnmanaged([]const u8) = .{};
        defer {
            for (headers.items) |h| allocator.free(h);
            headers.deinit(allocator);
        }
        {
            var tok = std.mem.splitScalar(u8, std.mem.trim(u8, header_raw, " \t\r"), ' ');
            while (tok.next()) |h| {
                if (h.len > 0) {
                    const owned = try allocator.dupe(u8, h);
                    errdefer allocator.free(owned);
                    try headers.append(allocator, owned);
                }
            }
        }

        errdefer {
            for (self.queries.items) |*q| q.deinit(allocator);
            self.queries.clearRetainingCapacity();
        }

        // Parse data rows
        while (try readNextLine(&file_reader, &read_buf, allocator, self.max_infile_bytes, self.infile)) |raw_line| {
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
                    // Split on ':' — matches C++ SplitString(tmpstr, ":") without intermediate vector allocations
                    var parts_iter = std.mem.splitScalar(u8, field, ':');
                    if (std.mem.eql(u8, head, "names")) {
                        while (parts_iter.next()) |p| {
                            const owned = try allocator.dupe(u8, p);
                            errdefer allocator.free(owned);
                            try q.names.append(allocator, owned);
                        }
                    } else if (std.mem.eql(u8, head, "seqs")) {
                        while (parts_iter.next()) |p| {
                            // Strip newlines from each sequence
                            const owned = try allocator.dupe(u8, p);
                            errdefer allocator.free(owned);
                            try q.seqs.append(allocator, owned);
                        }
                    } else if (std.mem.eql(u8, head, "only_genes")) {
                        while (parts_iter.next()) |p| {
                            const owned = try allocator.dupe(u8, p);
                            errdefer allocator.free(owned);
                            try q.only_genes.append(allocator, owned);
                        }
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

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const abs_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/test_args_input.csv", .{abs_path});
    defer allocator.free(full_path);

    const csv = "names seqs only_genes k_v_min k_v_max k_d_min k_d_max cdr3_length mut_freq\n" ++
        "seq1 ACGT gene1:gene2 1 10 1 5 30 0.05\n";
    try tmp_dir.dir.writeFile(.{ .sub_path = "test_args_input.csv", .data = csv });

    var args = try Args.initDefaults(allocator);
    defer args.deinit(allocator);
    allocator.free(args.infile);
    args.infile = try allocator.dupe(u8, full_path);

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

test "Args: readInfile FileTooBig diagnostics" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const abs_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/test_args_input_large.csv", .{abs_path});
    defer allocator.free(full_path);

    const csv = "names seqs only_genes k_v_min k_v_max k_d_min k_d_max cdr3_length mut_freq\n" ++
        "seq1 ACGT gene1:gene2 1 10 1 5 30 0.05\n";
    try tmp_dir.dir.writeFile(.{ .sub_path = "test_args_input_large.csv", .data = csv });

    var args = try Args.initDefaults(allocator);
    defer args.deinit(allocator);
    allocator.free(args.infile);
    args.infile = try allocator.dupe(u8, full_path);
    args.max_infile_bytes = 20;

    try std.testing.expectError(error.FileTooBig, args.readInfile(allocator));
}

test "Args: readInfile line longer than initial buffer (dynamic buffer growth)" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const abs_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/test_args_long_line.csv", .{abs_path});
    defer allocator.free(full_path);

    // Create a 200+ byte row with colon-separated names
    const csv = "names seqs only_genes k_v_min k_v_max k_d_min k_d_max cdr3_length mut_freq\n" ++
        "n1:n2:n3:n4:n5:n6:n7:n8:n9:n10:n11:n12:n13:n14:n15:n16:n17:n18:n19:n20 " ++
        "s1:s2:s3:s4:s5:s6:s7:s8:s9:s10:s11:s12:s13:s14:s15:s16:s17:s18:s19:s20 " ++
        "gene1:gene2 1 10 1 5 30 0.05\n";
    try tmp_dir.dir.writeFile(.{ .sub_path = "test_args_long_line.csv", .data = csv });

    const orig_buf_size = Args.initial_read_buffer_size;
    Args.initial_read_buffer_size = 64; // Force small 64-byte initial buffer
    defer Args.initial_read_buffer_size = orig_buf_size;

    var args = try Args.initDefaults(allocator);
    defer args.deinit(allocator);
    allocator.free(args.infile);
    args.infile = try allocator.dupe(u8, full_path);

    try args.readInfile(allocator);

    try std.testing.expectEqual(@as(usize, 1), args.queries.items.len);
    const q = &args.queries.items[0];
    try std.testing.expectEqual(@as(usize, 20), q.names.items.len);
    try std.testing.expectEqualStrings("n1", q.names.items[0]);
    try std.testing.expectEqualStrings("n20", q.names.items[19]);
    try std.testing.expectEqual(@as(usize, 20), q.seqs.items.len);
    try std.testing.expectEqualStrings("s1", q.seqs.items[0]);
    try std.testing.expectEqualStrings("s20", q.seqs.items[19]);
}

test "Args: readInfile input without trailing newline" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const abs_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/test_args_no_newline.csv", .{abs_path});
    defer allocator.free(full_path);

    // Two rows, no trailing newline after row 2
    const csv = "names seqs only_genes k_v_min k_v_max k_d_min k_d_max cdr3_length mut_freq\n" ++
        "seq1 ACGT gene1 1 10 1 5 30 0.05\n" ++
        "seq2 TGCA gene2 2 12 2 6 32 0.08";
    try tmp_dir.dir.writeFile(.{ .sub_path = "test_args_no_newline.csv", .data = csv });

    var args = try Args.initDefaults(allocator);
    defer args.deinit(allocator);
    allocator.free(args.infile);
    args.infile = try allocator.dupe(u8, full_path);

    try args.readInfile(allocator);

    try std.testing.expectEqual(@as(usize, 2), args.queries.items.len);
    try std.testing.expectEqualStrings("seq1", args.queries.items[0].names.items[0]);
    try std.testing.expectEqualStrings("seq2", args.queries.items[1].names.items[0]);
    try std.testing.expectApproxEqAbs(0.08, args.queries.items[1].mut_freq, 1e-9);
}

test "Args: readInfile row spanning buffer refill boundary" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const abs_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/test_args_refill.csv", .{abs_path});
    defer allocator.free(full_path);

    // Each row is ~45 bytes, with 64-byte buffer, row 2 crosses boundary
    const csv = "names seqs only_genes k_v_min k_v_max k_d_min k_d_max cdr3_length mut_freq\n" ++
        "seq1 ACGT gene1 1 10 1 5 30 0.05\n" ++
        "seq2 TGCA gene2 2 12 2 6 32 0.08\n" ++
        "seq3 GGGG gene3 3 15 3 7 35 0.10\n";
    try tmp_dir.dir.writeFile(.{ .sub_path = "test_args_refill.csv", .data = csv });

    const orig_buf_size = Args.initial_read_buffer_size;
    Args.initial_read_buffer_size = 64;
    defer Args.initial_read_buffer_size = orig_buf_size;

    var args = try Args.initDefaults(allocator);
    defer args.deinit(allocator);
    allocator.free(args.infile);
    args.infile = try allocator.dupe(u8, full_path);

    try args.readInfile(allocator);

    try std.testing.expectEqual(@as(usize, 3), args.queries.items.len);
    try std.testing.expectEqualStrings("seq1", args.queries.items[0].names.items[0]);
    try std.testing.expectEqualStrings("seq2", args.queries.items[1].names.items[0]);
    try std.testing.expectEqualStrings("seq3", args.queries.items[2].names.items[0]);
}
