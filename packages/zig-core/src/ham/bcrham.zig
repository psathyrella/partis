/// ham/bcrham.zig — Zig port of ham/src/bcrham.cc
///
/// bcrham entry point: wires together Args, Track, GermLines, HMMHolder,
/// Sequences, DPHandler, and Glomerator.  The C++ `main()` is replaced here
/// by `run()`, which is called from cabi.zig's bcrham_run / bcrham_forward /
/// bcrham_viterbi stubs once they are implemented.
///
/// C++ source: packages/ham/src/bcrham.cc
/// C++ author: psathyrella/ham

const std = @import("std");
const Args = @import("args.zig").Args;
const Track = @import("track.zig").Track;
const Sequence = @import("sequences.zig").Sequence;
const bcr = @import("bcr_utils/root.zig");
const GermLines = bcr.GermLines;
const HMMHolder = bcr.HMMHolder;
const KSet = bcr.KSet;
const KBounds = bcr.KBounds;
const DPHandler = @import("dp_handler.zig").DPHandler;
const glomerator_mod = @import("glomerator/root.zig");
const Glomerator = glomerator_mod.Glomerator;
const bcr_utils = @import("bcr_utils/root.zig");

// ─── helpers ─────────────────────────────────────────────────────────────────

/// Build the flat list-of-lists of Sequence from parsed Args.
/// Corresponds to C++ `GetSeqs(args, trk)`.
pub fn getSeqs(
    allocator: std.mem.Allocator,
    args: *const Args,
    trk: *const Track,
) !std.ArrayListUnmanaged(std.ArrayListUnmanaged(Sequence)) {
    var all_seqs: std.ArrayListUnmanaged(std.ArrayListUnmanaged(Sequence)) = .{};
    errdefer {
        for (all_seqs.items) |*sv| {
            for (sv.items) |*sq| sq.deinit(allocator);
            sv.deinit(allocator);
        }
        all_seqs.deinit(allocator);
    }
    for (args.queries.items) |*qrow| {
        var seqs: std.ArrayListUnmanaged(Sequence) = .{};
        errdefer {
            for (seqs.items) |*sq| sq.deinit(allocator);
            seqs.deinit(allocator);
        }
        const n = qrow.names.items.len;
        for (0..n) |iseq| {
            const sq = try Sequence.initFromString(
                allocator,
                trk,
                qrow.names.items[iseq],
                qrow.seqs.items[iseq],
            );
            try seqs.append(allocator, sq);
        }
        try all_seqs.append(allocator, seqs);
    }
    return all_seqs;
}

// ─── run_algorithm ────────────────────────────────────────────────────────────

/// Run the Forward or Viterbi algorithm for each query and stream output.
/// Corresponds to C++ `run_algorithm(hmms, gl, qry_seq_list, args)`.
pub fn runAlgorithm(
    allocator: std.mem.Allocator,
    hmms: *HMMHolder,
    gl: *GermLines,
    qry_seq_list: *std.ArrayListUnmanaged(std.ArrayListUnmanaged(Sequence)),
    args: *Args,
    writer: anytype,
) !void {
    const bcr_text = @import("bcr_utils/utils.zig");
    try bcr_text.streamHeader(writer, args.algorithm);

    var n_vtb: i32 = 0;
    var n_fwd: i32 = 0;

    for (qry_seq_list.items, 0..) |*qseqs, iqry| {
        const qrow = &args.queries.items[iqry];
        const kmin = KSet{ .v = @intCast(qrow.k_v_min), .d = @intCast(qrow.k_d_min) };
        const kmax = KSet{ .v = @intCast(qrow.k_v_max), .d = @intCast(qrow.k_d_max) };
        const kbounds = KBounds.initFromSets(kmin, kmax);

        // Build only_genes slice ([]const []const u8) from qrow.only_genes
        const only_genes = try allocator.alloc([]const u8, qrow.only_genes.items.len);
        defer allocator.free(only_genes);
        for (qrow.only_genes.items, only_genes) |g, *out| out.* = g;

        // Build pointer view into qseqs for the DP and stream calls (item 16).
        const qseq_ptrs = try allocator.alloc(*const Sequence, qseqs.items.len);
        defer allocator.free(qseq_ptrs);
        for (qseqs.items, qseq_ptrs) |*sq, *out| out.* = sq;

        var dph = try DPHandler.init(allocator, args.algorithm, args, gl, hmms);
        defer dph.deinit();

        var result = try dph.run(qseq_ptrs, kbounds, only_genes, @floatCast(qrow.mut_freq), true);
        defer result.deinit();

        if (result.no_path) {
            try bcr_text.streamErrorput(writer, args.algorithm, allocator, qseq_ptrs, "no_path");
        } else if (std.mem.eql(u8, args.algorithm, "viterbi")) {
            try bcr_text.streamViterbiOutput(writer, allocator, &result.best_event.?, qseq_ptrs, "");
        } else if (std.mem.eql(u8, args.algorithm, "forward")) {
            try bcr_text.streamForwardOutput(writer, allocator, qseq_ptrs, result.total_score, "");
        } else {
            return error.UnknownAlgorithm;
        }

        if (std.mem.eql(u8, args.algorithm, "viterbi"))
            n_vtb += 1
        else if (std.mem.eql(u8, args.algorithm, "forward"))
            n_fwd += 1;
    }
    {
        var buf: [128]u8 = undefined;
        const s = try std.fmt.bufPrint(&buf, "        calcd:   vtb {: <4}  fwd {: <4}\n", .{ @as(u32, @intCast(n_vtb)), @as(u32, @intCast(n_fwd)) });
        try std.fs.File.stdout().writeAll(s);
    }
}

// ─── top-level run ────────────────────────────────────────────────────────────

/// Top-level bcrham run.  Corresponds to C++ `main()`.
/// `argv` is a slice of null-terminated argument strings (excluding argv[0]).
pub fn run(allocator: std.mem.Allocator, argv: []const [*:0]const u8) !void {
    const start = std.time.milliTimestamp();

    // Parse args from argv: build defaults then populate each --flag value
    var args = try Args.initDefaults(allocator);
    defer args.deinit(allocator);
    {
        // Table-driven arg parsing: each tuple maps "--flag-name" → "field_name".
        const string_flags = .{
            .{ "--hmmdir", "hmmdir" },
            .{ "--datadir", "datadir" },
            .{ "--infile", "infile" },
            .{ "--outfile", "outfile" },
            .{ "--annotationfile", "annotationfile" },
            .{ "--input-cachefname", "input_cachefname" },
            .{ "--output-cachefname", "output_cachefname" },
            .{ "--locus", "locus" },
            .{ "--algorithm", "algorithm" },
            .{ "--ambig-base", "ambig_base" },
            .{ "--seed-unique-id", "seed_unique_id" },
        };
        const i32_flags = .{
            .{ "--debug", "debug" },
            .{ "--n-partitions-to-write", "n_partitions_to_write" },
            .{ "--biggest-naive-seq-cluster-to-calculate", "biggest_naive_seq_cluster_to_calculate" },
            .{ "--biggest-logprob-cluster-to-calculate", "biggest_logprob_cluster_to_calculate" },
        };
        const u32_flags = .{
            .{ "--n-final-clusters", "n_final_clusters" },
            .{ "--max-cluster-size", "max_cluster_size" },
            .{ "--random-seed", "random_seed" },
            .{ "--min-largest-cluster-size", "min_largest_cluster_size" },
        };
        const f32_flags = .{
            .{ "--hamming-fraction-bound-lo", "hamming_fraction_bound_lo" },
            .{ "--hamming-fraction-bound-hi", "hamming_fraction_bound_hi" },
            .{ "--logprob-ratio-threshold", "logprob_ratio_threshold" },
            .{ "--max-logprob-drop", "max_logprob_drop" },
        };
        const bool_flags = .{
            .{ "--no-chunk-cache", "no_chunk_cache" },
            .{ "--partition", "partition" },
            .{ "--dont-rescale-emissions", "dont_rescale_emissions" },
            .{ "--cache-naive-seqs", "cache_naive_seqs" },
            .{ "--cache-naive-hfracs", "cache_naive_hfracs" },
            .{ "--only-cache-new-vals", "only_cache_new_vals" },
            .{ "--write-logprob-for-each-partition", "write_logprob_for_each_partition" },
        };

        var i: usize = 0;
        while (i < argv.len) : (i += 1) {
            const flag = std.mem.span(argv[i]);
            var matched = false;

            // Boolean flags: no value argument consumed
            inline for (bool_flags) |entry| {
                if (std.mem.eql(u8, flag, entry[0])) {
                    @field(args, entry[1]) = true;
                    matched = true;
                }
            }
            if (matched) continue;

            // All remaining flags consume a value argument
            if (i + 1 >= argv.len) break;
            const val = std.mem.span(argv[i + 1]);

            inline for (string_flags) |entry| {
                if (std.mem.eql(u8, flag, entry[0])) {
                    allocator.free(@field(args, entry[1]));
                    @field(args, entry[1]) = try allocator.dupe(u8, val);
                    matched = true;
                }
            }
            inline for (i32_flags) |entry| {
                if (std.mem.eql(u8, flag, entry[0])) {
                    @field(args, entry[1]) = try std.fmt.parseInt(i32, val, 10);
                    matched = true;
                }
            }
            inline for (u32_flags) |entry| {
                if (std.mem.eql(u8, flag, entry[0])) {
                    @field(args, entry[1]) = try std.fmt.parseInt(u32, val, 10);
                    matched = true;
                }
            }
            inline for (f32_flags) |entry| {
                if (std.mem.eql(u8, flag, entry[0])) {
                    @field(args, entry[1]) = try std.fmt.parseFloat(f32, val);
                    matched = true;
                }
            }
            if (matched) i += 1; // consume the value argument
        }
    }
    // Validate required arguments
    var missing = false;
    if (args.hmmdir.len == 0) { std.debug.print("ERROR: missing required argument --hmmdir\n", .{}); missing = true; }
    if (args.datadir.len == 0) { std.debug.print("ERROR: missing required argument --datadir\n", .{}); missing = true; }
    if (args.infile.len == 0) { std.debug.print("ERROR: missing required argument --infile\n", .{}); missing = true; }
    if (args.outfile.len == 0) { std.debug.print("ERROR: missing required argument --outfile\n", .{}); missing = true; }
    if (args.locus.len == 0) { std.debug.print("ERROR: missing required argument --locus\n", .{}); missing = true; }
    if (args.algorithm.len == 0) { std.debug.print("ERROR: missing required argument --algorithm\n", .{}); missing = true; }
    if (missing) {
        std.debug.print(
            \\
            \\Usage: partis-zig-core --hmmdir <dir> --datadir <dir> --infile <file>
            \\                       --outfile <file> --locus <igh|igk|igl|tra|trb|trg|trd>
            \\                       --algorithm <viterbi|forward> [--ambig-base <base>]
            \\                       [--debug <0|1|2>] [--random-seed <n>] [--partition]
            \\
        , .{});
        std.process.exit(1);
    }

    if (args.infile.len > 0) try args.readInfile(allocator);

    // Build track
    const characters = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.init(allocator, "NUKES");
    defer track.deinit(allocator);
    for (characters) |c| {
        try track.addSymbol(allocator, c);
    }
    if (args.ambig_base.len > 0) try track.setAmbiguous(allocator, args.ambig_base);

    // Load germlines and HMMs
    var gl = try GermLines.init(allocator, args.datadir, args.locus);
    defer gl.deinit();

    var hmms = try HMMHolder.init(allocator, args.hmmdir, &gl);
    defer hmms.deinit();

    // Build query sequences
    var qry_seq_list = try getSeqs(allocator, &args, &track);
    defer {
        for (qry_seq_list.items) |*sv| {
            for (sv.items) |*sq| sq.deinit(allocator);
            sv.deinit(allocator);
        }
        qry_seq_list.deinit(allocator);
    }

    // Open outfile
    const outfile = try std.fs.createFileAbsolute(args.outfile, .{});
    defer outfile.close();
    var bw_buf: [65536]u8 = undefined;
    var bw = outfile.writer(&bw_buf);

    if (args.cache_naive_seqs or args.partition) {
        // Convert ArrayList(ArrayList(Sequence)) → [][]Sequence for Glomerator
        const seq_slices = try allocator.alloc([]Sequence, qry_seq_list.items.len);
        defer allocator.free(seq_slices);
        for (qry_seq_list.items, 0..) |sv, i| seq_slices[i] = sv.items;
        var glom = try Glomerator.init(allocator, &hmms, &gl, seq_slices, &args, &track);
        defer glom.deinit();
        if (args.cache_naive_seqs) {
            try glom.cacheNaiveSeqs();
        } else {
            try glom.cluster();
        }
    } else {
        try runAlgorithm(allocator, &hmms, &gl, &qry_seq_list, &args, &bw.interface);
    }

    try bw.interface.flush();

    const elapsed_s = @as(f64, @floatFromInt(std.time.milliTimestamp() - start)) / 1000.0;
    {
        var buf: [64]u8 = undefined;
        const s = try std.fmt.bufPrint(&buf, "        time: bcrham {d:.1}\n", .{elapsed_s});
        try std.fs.File.stdout().writeAll(s);
    }
}
