/// ham/glomerator/glomerator.zig — Zig port of ham/src/glomerator.cc
///
/// Hierarchical clustering (glomerating) of BCR sequences.
///
/// C++ source: packages/ham/src/glomerator.cc, packages/ham/include/glomerator.h
/// C++ author: psathyrella/ham

const std = @import("std");
const Args = @import("../args.zig").Args;
const DPHandler = @import("../dp_handler.zig").DPHandler;
const track_mod = @import("../track.zig");
const Track = track_mod.Track;
const ClusterPath = @import("../cluster_path.zig").ClusterPath;
const Partition = @import("../cluster_path.zig").Partition;
const sortPartition = @import("../cluster_path.zig").sortPartition;
const Sequence = @import("../sequences.zig").Sequence;
const mathutils = @import("../mathutils.zig");
const ham_text = @import("../text.zig");
const bcr = @import("../bcr_utils/root.zig");
const GermLines = bcr.GermLines;
const HMMHolder = bcr.HMMHolder;
const KSet = bcr.KSet;
const KBounds = bcr.KBounds;
const RecoEvent = bcr.RecoEvent;
const Result = bcr.Result;
const Query = @import("query.zig").Query;

/// Corresponds to C++ `ham::Glomerator`.
pub const Glomerator = struct {
    track: *Track,
    args: *Args,
    gl: *GermLines,
    hmms: *HMMHolder,
    allocator: std.mem.Allocator,

    /// The starting partition (one cluster per input query).
    initial_partition: Partition,

    /// All actual sequences, keyed by sequence name.
    single_seqs: std.StringHashMapUnmanaged(Sequence),
    /// Single-sequence cache entries.
    single_seq_cachefo: std.StringHashMapUnmanaged(Query),
    /// Permanent cache entries for clusters we've actually merged.
    cachefo: std.StringHashMapUnmanaged(Query),
    /// Temporary cache entries for clusters we're considering merging.
    tmp_cachefo: std.StringHashMapUnmanaged(Query),

    /// Cached log-probabilities.
    log_probs: std.StringHashMapUnmanaged(f64),
    /// Cached naive hamming fractions (joint key).
    naive_hfracs: std.StringHashMapUnmanaged(f64),
    /// Cached log-probability ratios.
    lratios: std.StringHashMapUnmanaged(f64),
    /// Cached naive sequences.
    naive_seqs: std.StringHashMapUnmanaged([]u8),
    /// Error strings.
    errors: std.StringHashMapUnmanaged([]u8),

    /// Queries that failed with no valid path.
    failed_queries: std.StringHashMapUnmanaged(void),

    /// Keys present in the initial cache file.
    initial_log_probs: std.StringHashMapUnmanaged(void),
    initial_naive_hfracs: std.StringHashMapUnmanaged(void),
    initial_naive_seqs: std.StringHashMapUnmanaged(void),

    /// Name translation maps.
    naive_seq_name_translations: std.StringHashMapUnmanaged([]u8),
    logprob_name_translations: std.StringHashMapUnmanaged([2][]u8),
    logprob_asymetric_translations: std.StringHashMapUnmanaged([]u8),
    name_subsets: std.StringHashMapUnmanaged([]u8),

    /// Counters.
    n_fwd_calculated: i32,
    n_vtb_calculated: i32,
    n_hfrac_calculated: i32,
    n_hfrac_merges: i32,
    n_lratio_merges: i32,

    asym_factor: f64,
    force_merge: bool,

    /// Asymmetry scaling factor (C++ asym_factor_ = 4.0).
    const ASYM_FACTOR: f64 = 4.0;

    /// C++: Glomerator::Glomerator() — glomerator.cc:6
    pub fn init(
        allocator: std.mem.Allocator,
        hmms: *HMMHolder,
        gl: *GermLines,
        qry_seq_lists: []const []Sequence,
        args: *Args,
        track: *Track,
    ) !Glomerator {
        var self = Glomerator{
            .track = track,
            .args = args,
            .gl = gl,
            .hmms = hmms,
            .allocator = allocator,
            .initial_partition = .{},
            .single_seqs = .{},
            .single_seq_cachefo = .{},
            .cachefo = .{},
            .tmp_cachefo = .{},
            .log_probs = .{},
            .naive_hfracs = .{},
            .lratios = .{},
            .naive_seqs = .{},
            .errors = .{},
            .failed_queries = .{},
            .initial_log_probs = .{},
            .initial_naive_hfracs = .{},
            .initial_naive_seqs = .{},
            .naive_seq_name_translations = .{},
            .logprob_name_translations = .{},
            .logprob_asymetric_translations = .{},
            .name_subsets = .{},
            .n_fwd_calculated = 0,
            .n_vtb_calculated = 0,
            .n_hfrac_calculated = 0,
            .n_hfrac_merges = 0,
            .n_lratio_merges = 0,
            .asym_factor = ASYM_FACTOR,
            .force_merge = false,
        };
        errdefer self.deinit();

        try self.readCacheFile();

        for (qry_seq_lists, 0..) |seq_list, iqry| {
            // Build the colon-separated key for this query group
            var names: std.ArrayListUnmanaged([]const u8) = .{};
            defer names.deinit(allocator);
            for (seq_list) |*sq| {
                try names.append(allocator, sq.name);
            }
            const key = try ham_text.joinStrings(allocator, names.items, ":");
            defer allocator.free(key);

            // Stash all sequences in single_seqs
            for (seq_list) |*sq| {
                const name_key = try allocator.dupe(u8, sq.name);
                errdefer allocator.free(name_key);
                if (!self.single_seqs.contains(name_key)) {
                    try self.single_seqs.put(allocator, name_key, try sq.clone(allocator));
                } else {
                    allocator.free(name_key);
                }
            }

            // Build bounds from args queries
            const qrow = &args.queries.items[iqry];
            const kmin = KSet{ .v = @intCast(@max(0, qrow.k_v_min)), .d = @intCast(@max(0, qrow.k_d_min)) };
            const kmax = KSet{ .v = @intCast(@max(1, qrow.k_v_max)), .d = @intCast(@max(1, qrow.k_d_max)) };
            const kbounds = KBounds.initFromSets(kmin, kmax);

            // Add to initial partition
            const key_owned = try allocator.dupe(u8, key);
            errdefer allocator.free(key_owned);
            try self.initial_partition.append(allocator, key_owned);

            // Build only_genes list
            var only_genes: std.ArrayListUnmanaged([]const u8) = .{};
            defer only_genes.deinit(allocator);
            for (qrow.only_genes.items) |g| {
                try only_genes.append(allocator, g);
            }

            // Build single_seq_cachefo entries
            const names_in_key = try splitName(allocator, key);
            defer {
                for (names_in_key.items) |n| allocator.free(n);
                var nl = names_in_key;
                nl.deinit(allocator);
            }
            for (names_in_key.items) |uid| {
                if (!self.single_seq_cachefo.contains(uid)) {
                    const uid_seqs = try self.getSeqs(uid);
                    defer {
                        for (uid_seqs.items) |*sq| sq.deinit(allocator);
                        var s = uid_seqs;
                        s.deinit(allocator);
                    }
                    const is_seed_missing = !ham_text.inString(args.seed_unique_id, uid, ":");
                    var uid_q = try Query.create(
                        allocator,
                        uid,
                        uid_seqs.items,
                        is_seed_missing,
                        only_genes.items,
                        kbounds,
                        @floatCast(qrow.mut_freq),
                        @intCast(@max(0, qrow.cdr3_length)),
                        null, null,
                    );
                    errdefer uid_q.deinit(allocator);
                    const uid_key = try allocator.dupe(u8, uid);
                    errdefer allocator.free(uid_key);
                    try self.single_seq_cachefo.put(allocator, uid_key, uid_q);
                }
            }

            // Build cachefo entry for the full query group
            const key_seqs = try self.getSeqs(key);
            defer {
                for (key_seqs.items) |*ks| ks.deinit(allocator);
                var ksl = key_seqs;
                ksl.deinit(allocator);
            }
            const is_seed_missing = !ham_text.inString(args.seed_unique_id, key, ":");
            var q = try Query.create(
                allocator,
                key,
                key_seqs.items,
                is_seed_missing,
                only_genes.items,
                kbounds,
                @floatCast(qrow.mut_freq),
                @intCast(@max(0, qrow.cdr3_length)),
                null, null,
            );
            errdefer q.deinit(allocator);
            const cf_key = try allocator.dupe(u8, key);
            errdefer allocator.free(cf_key);
            try self.cachefo.put(allocator, cf_key, q);
        }

        return self;
    }

    /// C++: ~Glomerator() — glomerator.cc:58 (cache write + stats moved to cluster())
    pub fn deinit(self: *Glomerator) void {
        const allocator = self.allocator;

        // Free initial_partition
        for (self.initial_partition.items) |s| allocator.free(s);
        self.initial_partition.deinit(allocator);

        // Free single_seqs
        var ssit = self.single_seqs.iterator();
        while (ssit.next()) |e| {
            allocator.free(e.key_ptr.*);
            e.value_ptr.deinit(allocator);
        }
        self.single_seqs.deinit(allocator);

        freeQueryMap(allocator, &self.single_seq_cachefo);
        freeQueryMap(allocator, &self.cachefo);
        freeQueryMap(allocator, &self.tmp_cachefo);

        freeStringF64Map(allocator, &self.log_probs);
        freeStringF64Map(allocator, &self.naive_hfracs);
        freeStringF64Map(allocator, &self.lratios);

        var nsit = self.naive_seqs.iterator();
        while (nsit.next()) |e| {
            allocator.free(e.key_ptr.*);
            allocator.free(e.value_ptr.*);
        }
        self.naive_seqs.deinit(allocator);

        var eit = self.errors.iterator();
        while (eit.next()) |e| {
            allocator.free(e.key_ptr.*);
            allocator.free(e.value_ptr.*);
        }
        self.errors.deinit(allocator);

        freeStringVoidMap(allocator, &self.failed_queries);
        freeStringVoidMap(allocator, &self.initial_log_probs);
        freeStringVoidMap(allocator, &self.initial_naive_hfracs);
        freeStringVoidMap(allocator, &self.initial_naive_seqs);

        freeStringStringMap(allocator, &self.naive_seq_name_translations);
        freeStringStringMap(allocator, &self.logprob_asymetric_translations);
        freeStringStringMap(allocator, &self.name_subsets);

        var lntit = self.logprob_name_translations.iterator();
        while (lntit.next()) |e| {
            allocator.free(e.key_ptr.*);
            allocator.free(e.value_ptr[0]);
            allocator.free(e.value_ptr[1]);
        }
        self.logprob_name_translations.deinit(allocator);
    }

    // ── Helper free functions ─────────────────────────────────────────────────

    fn freeQueryMap(allocator: std.mem.Allocator, map: *std.StringHashMapUnmanaged(Query)) void {
        var it = map.iterator();
        while (it.next()) |e| {
            allocator.free(e.key_ptr.*);
            e.value_ptr.deinit(allocator);
        }
        map.deinit(allocator);
    }

    fn freeStringF64Map(allocator: std.mem.Allocator, map: *std.StringHashMapUnmanaged(f64)) void {
        var it = map.iterator();
        while (it.next()) |e| allocator.free(e.key_ptr.*);
        map.deinit(allocator);
    }

    fn freeStringVoidMap(allocator: std.mem.Allocator, map: *std.StringHashMapUnmanaged(void)) void {
        var it = map.iterator();
        while (it.next()) |e| allocator.free(e.key_ptr.*);
        map.deinit(allocator);
    }

    fn freeStringStringMap(allocator: std.mem.Allocator, map: *std.StringHashMapUnmanaged([]u8)) void {
        var it = map.iterator();
        while (it.next()) |e| {
            allocator.free(e.key_ptr.*);
            allocator.free(e.value_ptr.*);
        }
        map.deinit(allocator);
    }

    // ── Main clustering entry points ──────────────────────────────────────────

    /// Run full hierarchical clustering and write output.
    /// Corresponds to C++ `Glomerator::Cluster`.
    pub fn cluster(self: *Glomerator) !void {
        if (self.args.logprob_ratio_threshold == -std.math.inf(f32)) {
            return error.LogprobRatioThresholdNotSet;
        }

        if (self.args.debug > 0) {
            var n_seeded: usize = 0;
            if (self.args.seed_unique_id.len > 0) {
                for (self.initial_partition.items) |cl| {
                    if (ham_text.inString(self.args.seed_unique_id, cl, ":")) n_seeded += 1;
                }
            }
            const seeded_str = if (self.args.seed_unique_id.len > 0)
                try std.fmt.allocPrint(self.allocator, "  ({d} seeded)", .{n_seeded})
            else
                try std.fmt.allocPrint(self.allocator, "", .{});
            defer self.allocator.free(seeded_str);
            const s = try std.fmt.allocPrint(self.allocator, "   hieragloming {d} clusters{s}\n", .{ self.initial_partition.items.len, seeded_str });
            defer self.allocator.free(s);
            try std.fs.File.stdout().writeAll(s);
        }

        // Clone partition entries and sort to match C++ set<string> ordering
        var init_part: Partition = .{};
        for (self.initial_partition.items) |name| {
            try init_part.append(self.allocator, try self.allocator.dupe(u8, name));
        }
        sortPartition(&init_part);
        var cp = try ClusterPath.initWithPartition(self.allocator, init_part, mathutils.NEG_INF);
        defer cp.deinit(self.allocator);

        while (!cp.finished) {
            try self.merge(&cp);
        }

        try self.writePartitions(&cp);
        try self.writeCacheFile();
        if (self.args.annotationfile.len > 0) {
            // WriteAnnotations is deprecated in C++; skip
        }
        {
            var buf: [256]u8 = undefined;
            const s = try std.fmt.bufPrint(&buf,
                "        calcd:   vtb {: <4}  fwd {: <4}  hfrac {: <8}\n        merged:  hfrac {: <4} lratio {: <4}\n",
                .{ @as(u32, @intCast(self.n_vtb_calculated)), @as(u32, @intCast(self.n_fwd_calculated)), @as(u32, @intCast(self.n_hfrac_calculated)), @as(u32, @intCast(self.n_hfrac_merges)), @as(u32, @intCast(self.n_lratio_merges)) },
            );
            try std.fs.File.stdout().writeAll(s);
        }
    }

    /// Cache all naive sequences.
    /// Corresponds to C++ `Glomerator::CacheNaiveSeqs`.
    pub fn cacheNaiveSeqs(self: *Glomerator) !void {
        try std.fs.File.stdout().writeAll("      caching all naive sequences\n");
        // Sort keys to match C++ map<string,...> iteration order
        const sorted_keys = try self.allocator.alloc([]const u8, self.cachefo.count());
        defer self.allocator.free(sorted_keys);
        {
            var it = self.cachefo.iterator();
            var idx: usize = 0;
            while (it.next()) |e| : (idx += 1) sorted_keys[idx] = e.key_ptr.*;
        }
        std.mem.sort([]const u8, sorted_keys, {}, struct {
            fn lessThan(_: void, a: []const u8, b: []const u8) bool {
                return std.mem.order(u8, a, b) == .lt;
            }
        }.lessThan);
        for (sorted_keys) |key| {
            _ = try self.getNaiveSeq(key, null);
        }
        try self.writeCacheFile();
        // Signal completion by opening and closing the outfile
        const f = try std.fs.cwd().createFile(self.args.outfile, .{});
        f.close();
        {
            var buf: [256]u8 = undefined;
            const s = try std.fmt.bufPrint(&buf,
                "        calcd:   vtb {: <4}  fwd {: <4}  hfrac {: <8}\n        merged:  hfrac {: <4} lratio {: <4}\n",
                .{ @as(u32, @intCast(self.n_vtb_calculated)), @as(u32, @intCast(self.n_fwd_calculated)), @as(u32, 0), @as(u32, 0), @as(u32, 0) },
            );
            try std.fs.File.stdout().writeAll(s);
        }
    }

    // ── Cache file I/O ────────────────────────────────────────────────────────

    /// C++: Glomerator::ReadCacheFile() — glomerator.cc:114
    fn readCacheFile(self: *Glomerator) !void {
        const allocator = self.allocator;
        if (self.args.input_cachefname.len == 0) {
            try std.fs.File.stdout().writeAll("        read-cache:  logprobs 0   naive-seqs 0\n");
            return;
        }

        const content = try std.fs.cwd().readFileAlloc(allocator, self.args.input_cachefname, 256 * 1024 * 1024);
        defer allocator.free(content);

        var line_iter = std.mem.splitScalar(u8, content, '\n');
        const header = line_iter.next() orelse {
            try std.fs.File.stdout().writeAll("        empty cachefile\n");
            return;
        };
        _ = header; // just skip it

        while (line_iter.next()) |raw| {
            const line = std.mem.trim(u8, raw, " \t\r");
            if (line.len == 0) continue;

            var fields = std.mem.splitScalar(u8, line, ',');
            const query = fields.next() orelse continue;
            const logprob_str = fields.next() orelse "";
            const naive_seq = fields.next() orelse "";
            const naive_hfrac_str = fields.next() orelse "";
            const errors_str = fields.next() orelse "";

            if (std.mem.indexOf(u8, errors_str, "no_path") != null) {
                const fq_key = try allocator.dupe(u8, query);
                try self.failed_queries.put(allocator, fq_key, {});
                continue;
            }

            if (logprob_str.len > 0) {
                // C++ uses stof() which parses to float32 then widens to double.
                // We must match that precision loss.
                const val32: f32 = std.fmt.parseFloat(f32, logprob_str) catch continue;
                const val: f64 = @as(f64, val32);
                const gop_lp = try self.log_probs.getOrPut(allocator, query);
                if (!gop_lp.found_existing) gop_lp.key_ptr.* = try allocator.dupe(u8, query);
                gop_lp.value_ptr.* = val;
                const gop_ilp = try self.initial_log_probs.getOrPut(allocator, query);
                if (!gop_ilp.found_existing) gop_ilp.key_ptr.* = try allocator.dupe(u8, query);
            }

            if (naive_hfrac_str.len > 0) {
                // C++ uses stof() which parses to float32 then widens to double.
                const val32: f32 = std.fmt.parseFloat(f32, naive_hfrac_str) catch continue;
                const val: f64 = @as(f64, val32);
                const gop_nh = try self.naive_hfracs.getOrPut(allocator, query);
                if (!gop_nh.found_existing) gop_nh.key_ptr.* = try allocator.dupe(u8, query);
                gop_nh.value_ptr.* = val;
                const gop_inh = try self.initial_naive_hfracs.getOrPut(allocator, query);
                if (!gop_inh.found_existing) gop_inh.key_ptr.* = try allocator.dupe(u8, query);
            }

            if (naive_seq.len > 0) {
                const gop_ns = try self.naive_seqs.getOrPut(allocator, query);
                if (gop_ns.found_existing) {
                    allocator.free(gop_ns.value_ptr.*);
                } else {
                    gop_ns.key_ptr.* = try allocator.dupe(u8, query);
                }
                gop_ns.value_ptr.* = try allocator.dupe(u8, naive_seq);
                const gop_ins = try self.initial_naive_seqs.getOrPut(allocator, query);
                if (!gop_ins.found_existing) gop_ins.key_ptr.* = try allocator.dupe(u8, query);
            }
        }

        {
            var buf: [128]u8 = undefined;
            const s = try std.fmt.bufPrint(&buf, "        read-cache:  logprobs {d}   naive-seqs {d}\n", .{ self.log_probs.count(), self.naive_seqs.count() });
            try std.fs.File.stdout().writeAll(s);
        }
    }

    /// C++: Glomerator::WriteCacheFile() — glomerator.cc:190
    fn writeCacheFile(self: *Glomerator) !void {
        if (self.args.output_cachefname.len == 0) return;

        const f = try std.fs.cwd().createFile(self.args.output_cachefname, .{});
        defer f.close();
        var wr_buf: [65536]u8 = undefined;
        var wr = f.writer(&wr_buf);
        const writer = &wr.interface;

        try writer.writeAll("unique_ids,logprob,naive_seq,naive_hfrac,errors\n");

        var keys_to_cache: std.StringHashMapUnmanaged(void) = .{};
        defer {
            var it = keys_to_cache.iterator();
            while (it.next()) |e| self.allocator.free(e.key_ptr.*);
            keys_to_cache.deinit(self.allocator);
        }

        {
            var it = self.log_probs.iterator();
            while (it.next()) |e| {
                if (self.args.only_cache_new_vals and self.initial_log_probs.contains(e.key_ptr.*)) continue;
                if (!keys_to_cache.contains(e.key_ptr.*))
                    try keys_to_cache.put(self.allocator, try self.allocator.dupe(u8, e.key_ptr.*), {});
            }
        }
        {
            var it = self.naive_seqs.iterator();
            while (it.next()) |e| {
                if (self.args.only_cache_new_vals and self.initial_naive_seqs.contains(e.key_ptr.*)) continue;
                if (!keys_to_cache.contains(e.key_ptr.*))
                    try keys_to_cache.put(self.allocator, try self.allocator.dupe(u8, e.key_ptr.*), {});
            }
        }
        if (self.args.cache_naive_hfracs) {
            var it = self.naive_hfracs.iterator();
            while (it.next()) |e| {
                if (self.args.only_cache_new_vals and self.initial_naive_hfracs.contains(e.key_ptr.*)) continue;
                if (!keys_to_cache.contains(e.key_ptr.*))
                    try keys_to_cache.put(self.allocator, try self.allocator.dupe(u8, e.key_ptr.*), {});
            }
        }

        // Sort cache keys to match C++ map<string,...> iteration order
        const sorted_cache_keys = try self.allocator.alloc([]const u8, keys_to_cache.count());
        defer self.allocator.free(sorted_cache_keys);
        {
            var kit = keys_to_cache.iterator();
            var kidx: usize = 0;
            while (kit.next()) |e| : (kidx += 1) sorted_cache_keys[kidx] = e.key_ptr.*;
        }
        std.mem.sort([]const u8, sorted_cache_keys, {}, struct {
            fn lessThan(_: void, a: []const u8, b: []const u8) bool {
                return std.mem.order(u8, a, b) == .lt;
            }
        }.lessThan);
        for (sorted_cache_keys) |key| {
            try writer.writeAll(key);
            try writer.writeByte(',');
            if (self.log_probs.get(key)) |lp| try writer.print("{d}", .{lp});
            try writer.writeByte(',');
            if (self.naive_seqs.get(key)) |ns| try writer.writeAll(ns);
            try writer.writeByte(',');
            if (self.args.cache_naive_hfracs) {
                if (self.naive_hfracs.get(key)) |nh| try writer.print("{d}", .{nh});
            }
            try writer.writeByte(',');
            if (self.errors.get(key)) |err| try writer.writeAll(err);
            try writer.writeByte('\n');
        }
        try wr.interface.flush();
    }

    // ── Partition I/O ─────────────────────────────────────────────────────────

    /// C++: Glomerator::WritePartitions() — glomerator.cc:227
    fn writePartitions(self: *Glomerator, cp: *ClusterPath) !void {
        const run_start = std.time.milliTimestamp();
        if (self.args.debug > 0) {
            try std.fs.File.stdout().writeAll("        writing partitions\n");
        }
        const f = try std.fs.cwd().createFile(self.args.outfile, .{});
        defer f.close();
        var wr_buf: [65536]u8 = undefined;
        var wr = f.writer(&wr_buf);
        const writer = &wr.interface;

        try writer.writeAll("partition,logprob\n");

        const n_parts = cp.partitions.items.len;
        const istart: usize = if (n_parts > @as(usize, @intCast(self.args.n_partitions_to_write)))
            n_parts - @as(usize, @intCast(self.args.n_partitions_to_write))
        else
            0;

        for (istart..n_parts) |ipart| {
            if (self.args.write_logprob_for_each_partition) {
                const lp = try self.logProbOfPartition(&cp.partitions.items[ipart]);
                cp.setLogprob(ipart, lp);
            }

            var first_cluster = true;
            for (cp.partitions.items[ipart].items) |cl| {
                if (!first_cluster) try writer.writeByte(';');
                try writer.writeAll(cl);
                first_cluster = false;
            }
            try writer.print(",{d}\n", .{cp.logprobs.items[ipart]});
        }

        try writer.flush();

        if (self.args.write_logprob_for_each_partition) {
            const elapsed = @as(f64, @floatFromInt(std.time.milliTimestamp() - run_start)) / 1000.0;
            const s = try std.fmt.allocPrint(self.allocator, "        partition writing time (probably includes calculating a bunch of new logprobs) {d:.1}\n", .{elapsed});
            defer self.allocator.free(s);
            try std.fs.File.stdout().writeAll(s);
        }
    }

    // ── Merge step ────────────────────────────────────────────────────────────

    /// Perform one merge step.
    /// Corresponds to C++ `Glomerator::Merge`.
    pub fn merge(self: *Glomerator, path: *ClusterPath) !void {
        var qpair = try self.findHfracMerge(path);
        defer if (qpair) |*qp| qp.deinit(self.allocator);

        // if no hfrac merge, try lratio
        if (qpair == null) {
            qpair = try self.findLRatioMerge(path);
        }

        // max_cluster_size check
        if (self.args.max_cluster_size > 0) {
            const cur_part = path.currentPartition();
            for (cur_part.items) |cl| {
                if (countMembers(cl) > self.args.max_cluster_size) {
                    const s = try std.fmt.allocPrint(self.allocator, "    --max-cluster-size: stopping with a cluster of size {d} (> {d})\n", .{ @as(u32, @intCast(countMembers(cl))), self.args.max_cluster_size });
                    defer self.allocator.free(s);
                    try std.fs.File.stdout().writeAll(s);
                    path.finished = true;
                }
            }
            if (path.finished) return;
        }

        if (qpair == null) {
            // No valid merge found
            if (self.args.n_final_clusters == 0 and self.args.min_largest_cluster_size == 0) {
                path.finished = true;
            } else {
                const cur_size = path.currentPartition().items.len;
                const largest = largestClusterSize(path.currentPartition());
                const need_more = (self.args.n_final_clusters > 0 and cur_size > self.args.n_final_clusters) or
                    (self.args.min_largest_cluster_size > 0 and largest < self.args.min_largest_cluster_size);
                if (need_more) {
                    if (self.force_merge) {
                        path.finished = true;
                        if (self.args.debug > 0) {
                            if (self.args.n_final_clusters > 0) {
                                const s = try std.fmt.allocPrint(self.allocator, "    couldn't merge beyond {d} clusters despite setting force merge\n", .{cur_size});
                                defer self.allocator.free(s);
                                try std.fs.File.stdout().writeAll(s);
                            }
                            if (self.args.min_largest_cluster_size > 0) {
                                const s = try std.fmt.allocPrint(self.allocator, "    couldn't merge past a biggest cluster of {d} despite setting force merge\n", .{largest});
                                defer self.allocator.free(s);
                                try std.fs.File.stdout().writeAll(s);
                            }
                        }
                    } else {
                        self.force_merge = true;
                        if (self.args.debug > 0) {
                            if (self.args.n_final_clusters > 0) {
                                const s = try std.fmt.allocPrint(self.allocator, "    setting force merge (currently have {d} clusters, requested {d})\n", .{ cur_size, self.args.n_final_clusters });
                                defer self.allocator.free(s);
                                try std.fs.File.stdout().writeAll(s);
                            }
                            if (self.args.min_largest_cluster_size > 0) {
                                const s = try std.fmt.allocPrint(self.allocator, "    setting force merge (current biggest cluster {d}, requested {d})\n", .{ largest, self.args.min_largest_cluster_size });
                                defer self.allocator.free(s);
                                try std.fs.File.stdout().writeAll(s);
                            }
                        }
                    }
                } else {
                    path.finished = true;
                }
            }
            return;
        }

        const chosen = qpair.?;
        const parents = chosen.parents orelse return;
        const p1 = parents[0];
        const p2 = parents[1];

        // Move to permanent cache
        {
            const key = try self.allocator.dupe(u8, chosen.name);
            errdefer self.allocator.free(key);
            var q_copy = try cloneQuery(self.allocator, &chosen);
            errdefer q_copy.deinit(self.allocator);
            try self.cachefo.put(self.allocator, key, q_copy);
        }

        _ = try self.getNaiveSeq(chosen.name, if (chosen.parents != null) &[2][]const u8{ p1, p2 } else null);
        try self.updateLogProbTranslationsForAsymetrics(&chosen);
        try self.moveSubsetsFromTmpCache(chosen.name);

        // Build new partition
        var new_partition: Partition = .{};
        errdefer {
            for (new_partition.items) |s| self.allocator.free(s);
            new_partition.deinit(self.allocator);
        }
        const cur_part = path.currentPartition();
        for (cur_part.items) |cl| {
            if (!std.mem.eql(u8, cl, p1) and !std.mem.eql(u8, cl, p2)) {
                try new_partition.append(self.allocator, try self.allocator.dupe(u8, cl));
            }
        }
        try new_partition.append(self.allocator, try self.allocator.dupe(u8, chosen.name));
        sortPartition(&new_partition);
        try path.addPartition(self.allocator, new_partition, mathutils.NEG_INF, @intCast(self.args.n_partitions_to_write));

        if (self.args.debug > 0) {
            const s = try std.fmt.allocPrint(self.allocator, "       merged   {s}  {s}\n          removing {d} entries from tmp cache\n", .{ p1, p2, self.tmp_cachefo.count() });
            defer self.allocator.free(s);
            try std.fs.File.stdout().writeAll(s);
        }

        // Clear tmp cache
        var tcit = self.tmp_cachefo.iterator();
        while (tcit.next()) |e| {
            self.allocator.free(e.key_ptr.*);
            e.value_ptr.deinit(self.allocator);
        }
        self.tmp_cachefo.clearRetainingCapacity();

        // Check if we're done
        const new_size = path.currentPartition().items.len;
        const new_largest = largestClusterSize(path.currentPartition());
        if ((self.args.n_final_clusters > 0 and new_size <= self.args.n_final_clusters) or
            (self.args.min_largest_cluster_size > 0 and new_largest >= self.args.min_largest_cluster_size))
        {
            path.finished = true;
            if (self.args.n_final_clusters > 0) {
                const s = try std.fmt.allocPrint(self.allocator, "    finished glomerating to {d} clusters (requested {d}), force merge status {d}\n", .{ new_size, self.args.n_final_clusters, @as(u8, if (self.force_merge) 1 else 0) });
                defer self.allocator.free(s);
                try std.fs.File.stdout().writeAll(s);
            }
            if (self.args.min_largest_cluster_size > 0) {
                const s = try std.fmt.allocPrint(self.allocator, "    finished glomerating to a biggest cluster of {d} (requested {d}), force merge status {d}\n", .{ new_largest, self.args.min_largest_cluster_size, @as(u8, if (self.force_merge) 1 else 0) });
                defer self.allocator.free(s);
                try std.fs.File.stdout().writeAll(s);
            }
        }
    }

    // ── Log-probability calculation ───────────────────────────────────────────

    /// C++: Glomerator::LogProbOfPartition() — glomerator.cc:282
    fn logProbOfPartition(self: *Glomerator, partition: *const Partition) !f64 {
        var total: f64 = 0.0;
        for (partition.items) |key| {
            const lp = try self.getLogProb(key);
            total = mathutils.addWithMinusInfinities(total, lp);
        }
        return total;
    }

    /// C++: Glomerator::GetLogProb() — glomerator.cc:681
    fn getLogProb(self: *Glomerator, queries: []const u8) !f64 {
        if (self.log_probs.get(queries)) |lp| return lp;
        const lp = try self.calculateLogProb(queries);
        const key = try self.allocator.dupe(u8, queries);
        try self.log_probs.put(self.allocator, key, lp);
        return lp;
    }

    /// C++: Glomerator::CalculateLogProb() — glomerator.cc:754
    fn calculateLogProb(self: *Glomerator, queries: []const u8) !f64 {
        std.debug.assert(!self.log_probs.contains(queries));
        self.n_fwd_calculated += 1;

        var dph = try DPHandler.init(self.allocator, "forward", self.args, self.gl, self.hmms);
        defer dph.deinit();

        const cacheref = try self.getCachefo(queries);
        const only_genes_slice = try self.geneListToSlice(cacheref.only_genes.items);
        defer self.allocator.free(only_genes_slice);

        var result = try dph.run(cacheref.seqs.items, cacheref.kbounds, only_genes_slice, cacheref.mute_freq, true);
        defer result.deinit();

        if (result.no_path) {
            try self.addFailedQuery(queries, "no_path");
            return mathutils.NEG_INF;
        }
        return result.total_score;
    }

    /// C++: Glomerator::GetLogProbRatio() — glomerator.cc:692
    fn getLogProbRatio(self: *Glomerator, key_a: []const u8, key_b: []const u8) !f64 {
        const joint_name = try self.joinNames(key_a, key_b);
        defer self.allocator.free(joint_name);

        if (self.lratios.get(joint_name)) |lr| return lr;

        var full_qmerged = try self.getMergedQuery(key_a, key_b);
        defer full_qmerged.deinit(self.allocator);
        const parents_to_calc = try self.getLogProbPairOfNamesToCalculate(joint_name, if (full_qmerged.parents != null) &full_qmerged.parents.? else null);
        defer {
            self.allocator.free(parents_to_calc[0]);
            self.allocator.free(parents_to_calc[1]);
        }

        const lp_a = try self.getLogProb(parents_to_calc[0]);
        const lp_b = try self.getLogProb(parents_to_calc[1]);

        var qmerged_calc = try self.getMergedQuery(parents_to_calc[0], parents_to_calc[1]);
        defer qmerged_calc.deinit(self.allocator);
        const lp_ab = try self.getLogProb(qmerged_calc.name);

        const lratio = lp_ab - lp_a - lp_b;

        if (self.args.debug > 0) {
            const dbg = try std.fmt.allocPrint(self.allocator, "             {d: >8.3} = {s} - {s} - {s}", .{ lratio, joint_name, key_a, key_b });
            defer self.allocator.free(dbg);
            if (!std.mem.eql(u8, qmerged_calc.name, joint_name) or !std.mem.eql(u8, parents_to_calc[0], key_a) or !std.mem.eql(u8, parents_to_calc[1], key_b)) {
                const calcd = try std.fmt.allocPrint(self.allocator, " (calcd  {s} - {s} - {s})", .{ qmerged_calc.name, parents_to_calc[0], parents_to_calc[1] });
                defer self.allocator.free(calcd);
                const combined = try std.fmt.allocPrint(self.allocator, "{s}{s}\n", .{ dbg, calcd });
                defer self.allocator.free(combined);
                try std.fs.File.stdout().writeAll(combined);
            } else {
                const with_nl = try std.fmt.allocPrint(self.allocator, "{s}\n", .{dbg});
                defer self.allocator.free(with_nl);
                try std.fs.File.stdout().writeAll(with_nl);
            }
        }

        const lr_key = try self.allocator.dupe(u8, joint_name);
        try self.lratios.put(self.allocator, lr_key, lratio);
        return lratio;
    }

    // ── Naive sequence calculation ────────────────────────────────────────────

    /// C++: Glomerator::GetNaiveSeq() — glomerator.cc:641
    fn getNaiveSeq(self: *Glomerator, queries: []const u8, parents: ?*const [2][]const u8) anyerror![]const u8 {
        if (self.naive_seqs.get(queries)) |ns| return ns;

        // maybe use parent naive seq directly
        if (parents != null) {
            const replace = try self.findNaiveSeqNameReplace(parents.?);
            if (replace) |rep| {
                defer self.allocator.free(rep);
                const parent_ns = try self.getNaiveSeq(rep, null);
                const key = try self.allocator.dupe(u8, queries);
                const val = try self.allocator.dupe(u8, parent_ns);
                try self.naive_seqs.put(self.allocator, key, val);
                return self.naive_seqs.get(queries).?;
            }
        }

        const queries_to_calc = try self.getNaiveSeqNameToCalculate(queries);
        defer self.allocator.free(queries_to_calc);

        if (!self.naive_seqs.contains(queries_to_calc)) {
            const tmp_ns = try self.calculateNaiveSeq(queries_to_calc, null);
            if (tmp_ns.len == 0) {
                self.allocator.free(tmp_ns);
                {
                    const warn = try std.fmt.allocPrint(self.allocator, "  warning: zero length naive sequence for {s}\n", .{queries_to_calc});
                    defer self.allocator.free(warn);
                    try std.fs.File.stdout().writeAll(warn);
                }
                return "";
            }
            const key = try self.allocator.dupe(u8, queries_to_calc);
            errdefer self.allocator.free(key);
            // tmp_ns is already an owned allocation from calculateNaiveSeq; store directly
            errdefer self.allocator.free(tmp_ns);
            try self.naive_seqs.put(self.allocator, key, tmp_ns);
        }

        if (!std.mem.eql(u8, queries_to_calc, queries)) {
            const translated_ns = self.naive_seqs.get(queries_to_calc).?;
            const key = try self.allocator.dupe(u8, queries);
            const val = try self.allocator.dupe(u8, translated_ns);
            try self.naive_seqs.put(self.allocator, key, val);
        }

        return self.naive_seqs.get(queries).?;
    }

    /// C++: Glomerator::CalculateNaiveSeq() — glomerator.cc:727
    fn calculateNaiveSeq(self: *Glomerator, queries: []const u8, event_out: ?*RecoEvent) ![]u8 {
        std.debug.assert(!self.naive_seqs.contains(queries));
        self.n_vtb_calculated += 1;

        var dph = try DPHandler.init(self.allocator, "viterbi", self.args, self.gl, self.hmms);
        defer dph.deinit();

        const cacheref = try self.getCachefo(queries);
        const only_genes_slice = try self.geneListToSlice(cacheref.only_genes.items);
        defer self.allocator.free(only_genes_slice);

        var result = try dph.run(cacheref.seqs.items, cacheref.kbounds, only_genes_slice, cacheref.mute_freq, true);
        defer result.deinit();

        if (result.no_path) {
            try self.addFailedQuery(queries, "no_path");
            return self.allocator.dupe(u8, "");
        }

        if (event_out != null) {
            // NOTE: transferring event ownership requires nulling result.best_event
            // so that result.deinit() doesn't free the slices the caller now owns.
            if (result.best_event) |ev| {
                event_out.?.* = ev;
                result.best_event = null;
                return self.allocator.dupe(u8, ev.naive_seq);
            }
        }

        // Must dupe before defer result.deinit() frees the memory
        return self.allocator.dupe(u8, result.best_event.?.naive_seq);
    }

    // ── Hamming fraction ──────────────────────────────────────────────────────

    /// C++: Glomerator::CalculateHfrac() — glomerator.cc:455
    fn calculateHfrac(self: *Glomerator, seq_a: []const u8, seq_b: []const u8) !f64 {
        self.n_hfrac_calculated += 1;
        if (seq_a.len != seq_b.len) return error.SequenceLengthMismatch;

        var distance: usize = 0;
        var len_excl_ambig: usize = 0;
        const ambig_idx = track_mod.AMBIGUOUS_INDEX;

        for (seq_a, seq_b) |ca, cb| {
            const ia = self.track.symbolIndex(&[_]u8{ca}) catch ambig_idx;
            const ib = self.track.symbolIndex(&[_]u8{cb}) catch ambig_idx;
            if (ia == ambig_idx or ib == ambig_idx) continue;
            len_excl_ambig += 1;
            if (ia != ib) distance += 1;
        }

        if (len_excl_ambig == 0) return std.math.inf(f64);
        return @as(f64, @floatFromInt(distance)) / @as(f64, @floatFromInt(len_excl_ambig));
    }

    /// C++: Glomerator::NaiveHfrac() — glomerator.cc:474
    fn naiveHfrac(self: *Glomerator, key_a: []const u8, key_b: []const u8) !f64 {
        const joint_key = try self.joinNames(key_a, key_b);
        defer self.allocator.free(joint_key);

        if (self.naive_hfracs.get(joint_key)) |hf| return hf;

        const ns_a = try self.getNaiveSeq(key_a, null);
        const ns_b = try self.getNaiveSeq(key_b, null);

        if (ns_a.len == 0 or ns_b.len == 0 or
            self.failed_queries.contains(key_a) or self.failed_queries.contains(key_b))
            return std.math.inf(f64);

        const hf = try self.calculateHfrac(ns_a, ns_b);
        const jk = try self.allocator.dupe(u8, joint_key);
        try self.naive_hfracs.put(self.allocator, jk, hf);
        return hf;
    }

    // ── Merge finders ────────────────────────────────────────────────────────

    /// C++: Glomerator::FindHfracMerge() — glomerator.cc:1000
    fn findHfracMerge(self: *Glomerator, path: *ClusterPath) !?Query {
        var min_hf: f64 = std.math.inf(f64);
        var found = false;
        var best_query: ?Query = null;

        const cur_part = path.currentPartition();
        const has_seed = self.args.seed_unique_id.len > 0;

        var ia: usize = 0;
        while (ia < cur_part.items.len) : (ia += 1) {
            const key_a = cur_part.items[ia];
            // C++: when seed is set, outer loop only covers seeded clusters (GetSeededClusters)
            if (has_seed) {
                const ca = try self.getCachefo(key_a);
                if (ca.seed_missing) continue;
            }

            // C++: without seed, inner starts at ia+1; with seed, inner starts at 0 (all clusters)
            var ib: usize = if (has_seed) 0 else ia + 1;
            while (ib < cur_part.items.len) : (ib += 1) {
                const key_b = cur_part.items[ib];

                if (std.mem.eql(u8, key_a, key_b)) continue; // skip self-pairs (needed with seed)
                if (self.failed_queries.contains(key_a) or self.failed_queries.contains(key_b)) continue;

                _ = try self.getCachefo(key_a);
                _ = try self.getCachefo(key_b);
                const ca = self.cachefo.getPtr(key_a) orelse self.tmp_cachefo.getPtr(key_a) orelse return error.MissingSeqCache;
                const cb = self.cachefo.getPtr(key_b) orelse self.tmp_cachefo.getPtr(key_b) orelse return error.MissingSeqCache;
                if (ca.cdr3_length != cb.cdr3_length) continue;

                const hf = try self.naiveHfrac(key_a, key_b);
                if (hf > self.args.hamming_fraction_bound_hi) continue;
                if (self.args.hamming_fraction_bound_lo <= 0.0 or hf >= self.args.hamming_fraction_bound_lo) continue;

                if (hf < min_hf) {
                    min_hf = hf;
                    if (best_query) |*bq| bq.deinit(self.allocator);
                    var merged = try self.getMergedQuery(key_a, key_b);
                    defer merged.deinit(self.allocator);
                    best_query = try cloneQuery(self.allocator, &merged);
                    found = true;
                }
            }
        }

        if (found) {
            self.n_hfrac_merges += 1;
            if (self.args.debug > 0) {
                if (best_query) |bq| {
                    if (bq.parents) |parents| {
                        const s = try std.fmt.allocPrint(self.allocator, "          hfrac merge  {d:.3}   {s}  {s}\n", .{ min_hf, printStr(parents[0]), printStr(parents[1]) });
                        defer self.allocator.free(s);
                        try std.fs.File.stdout().writeAll(s);
                    }
                }
            }
            return best_query;
        }
        return null;
    }

    /// C++: Glomerator::FindLRatioMerge() — glomerator.cc:1063
    fn findLRatioMerge(self: *Glomerator, path: *ClusterPath) !?Query {
        var max_lratio: f64 = mathutils.NEG_INF;
        var best_query: ?Query = null;

        const cur_part = path.currentPartition();
        const has_seed = self.args.seed_unique_id.len > 0;

        var ia: usize = 0;
        while (ia < cur_part.items.len) : (ia += 1) {
            const key_a = cur_part.items[ia];
            // C++: when seed is set, outer loop only covers seeded clusters (GetSeededClusters)
            if (has_seed) {
                const ca = try self.getCachefo(key_a);
                if (ca.seed_missing) continue;
            }

            // C++: without seed, inner starts at ia+1; with seed, inner starts at 0 (all clusters)
            var ib: usize = if (has_seed) 0 else ia + 1;
            while (ib < cur_part.items.len) : (ib += 1) {
                const key_b = cur_part.items[ib];

                if (std.mem.eql(u8, key_a, key_b)) continue; // skip self-pairs (needed with seed)
                if (self.failed_queries.contains(key_a) or self.failed_queries.contains(key_b)) continue;

                _ = try self.getCachefo(key_a);
                _ = try self.getCachefo(key_b);
                const ca = self.cachefo.getPtr(key_a) orelse self.tmp_cachefo.getPtr(key_a) orelse return error.MissingSeqCache;
                const cb = self.cachefo.getPtr(key_b) orelse self.tmp_cachefo.getPtr(key_b) orelse return error.MissingSeqCache;
                if (ca.cdr3_length != cb.cdr3_length) continue;

                const hf = try self.naiveHfrac(key_a, key_b);
                if (hf > self.args.hamming_fraction_bound_hi) continue;

                const lratio = try self.getLogProbRatio(key_a, key_b);
                const cluster_size: i32 = @intCast(countMembers(key_a) + countMembers(key_b));
                if (!self.force_merge and self.likelihoodRatioTooSmall(lratio, cluster_size)) continue;

                if (lratio > max_lratio) {
                    max_lratio = lratio;
                    if (best_query) |*bq| bq.deinit(self.allocator);
                    var merged = try self.getMergedQuery(key_a, key_b);
                    defer merged.deinit(self.allocator);
                    best_query = try cloneQuery(self.allocator, &merged);
                }
            }
        }

        if (max_lratio != mathutils.NEG_INF) {
            self.n_lratio_merges += 1;
            if (self.args.debug > 0) {
                if (best_query) |bq| {
                    if (bq.parents) |parents| {
                        const s = try std.fmt.allocPrint(self.allocator, "          lratio merge  {d:.3}   {s}  {s}\n", .{ max_lratio, printStr(parents[0]), printStr(parents[1]) });
                        defer self.allocator.free(s);
                        try std.fs.File.stdout().writeAll(s);
                    }
                }
            }
            return best_query;
        }
        return null;
    }

    /// Check whether the log-likelihood ratio is below a cluster-size-dependent threshold.
    /// Ported from C++ `Glomerator::LikelihoodRatioTooSmall()`.
    /// The threshold loosens by 1.0 for each additional sequence beyond 2, capped at ccs>=6.
    /// `ccs` is the combined cluster size (number of sequences in the merged cluster).
    fn likelihoodRatioTooSmall(self: *Glomerator, lratio: f64, combined_cluster_size: i32) bool {
        const threshold = self.args.logprob_ratio_threshold;
        if (combined_cluster_size == 2) return lratio < threshold;
        if (combined_cluster_size == 3) return lratio < threshold - 2.0;
        if (combined_cluster_size == 4) return lratio < threshold - 3.0;
        if (combined_cluster_size == 5) return lratio < threshold - 4.0;
        return lratio < threshold - 5.0;
    }

    // ── Cache utilities ───────────────────────────────────────────────────────

    /// C++: Glomerator::cachefo() — glomerator.cc:782
    fn getCachefo(self: *Glomerator, queries: []const u8) !*Query {
        if (self.cachefo.getPtr(queries)) |q| return q;
        if (self.tmp_cachefo.getPtr(queries)) |q| return q;

        // Build on the fly from single_seq_cachefo
        var only_gene_set: std.StringHashMapUnmanaged(void) = .{};
        defer {
            var it = only_gene_set.iterator();
            while (it.next()) |e| self.allocator.free(e.key_ptr.*);
            only_gene_set.deinit(self.allocator);
        }

        var kbounds = KBounds{};
        var mute_freq_total: f64 = 0.0;
        var cdr3_length: usize = 0;

        const names = try splitName(self.allocator, queries);
        defer {
            for (names.items) |n| self.allocator.free(n);
            var nl = names;
            nl.deinit(self.allocator);
        }

        for (names.items, 0..) |uid, is| {
            const scache = self.single_seq_cachefo.getPtr(uid) orelse return error.MissingSeqCache;
            for (scache.only_genes.items) |g| {
                if (!only_gene_set.contains(g)) {
                    try only_gene_set.put(self.allocator, try self.allocator.dupe(u8, g), {});
                }
            }
            if (is == 0) {
                kbounds = scache.kbounds;
                cdr3_length = scache.cdr3_length;
            } else {
                kbounds = kbounds.logicalOr(scache.kbounds);
                if (cdr3_length != scache.cdr3_length) return error.Cdr3LengthMismatch;
            }
            mute_freq_total += scache.mute_freq;
        }

        var only_genes_list: std.ArrayListUnmanaged([]const u8) = .{};
        defer only_genes_list.deinit(self.allocator);
        {
            var ogit = only_gene_set.iterator();
            while (ogit.next()) |e| {
                try only_genes_list.append(self.allocator, e.key_ptr.*);
            }
        }
        // Sort to match C++ set<string> iteration order
        std.mem.sort([]const u8, only_genes_list.items, {}, struct {
            fn lessThan(_: void, a: []const u8, b: []const u8) bool {
                return std.mem.order(u8, a, b) == .lt;
            }
        }.lessThan);

        const key_seqs = try self.getSeqs(queries);
        defer {
            var ks = key_seqs;
            ks.deinit(self.allocator);
        }

        const is_seed_missing = !ham_text.inString(self.args.seed_unique_id, queries, ":");
        var new_q = try Query.create(
            self.allocator,
            queries,
            key_seqs.items,
            is_seed_missing,
            only_genes_list.items,
            kbounds,
            @floatCast(mute_freq_total / @as(f64, @floatFromInt(names.items.len))),
            cdr3_length,
            null, null,
        );
        errdefer new_q.deinit(self.allocator);

        const cf_key = try self.allocator.dupe(u8, queries);
        errdefer self.allocator.free(cf_key);
        try self.tmp_cachefo.put(self.allocator, cf_key, new_q);
        return self.tmp_cachefo.getPtr(queries).?;
    }

    /// C++: Glomerator::GetMergedQuery() — glomerator.cc:908
    fn getMergedQuery(self: *Glomerator, name_a: []const u8, name_b: []const u8) !Query {
        const joint_name = try self.joinNames(name_a, name_b);
        defer self.allocator.free(joint_name);

        if (self.cachefo.getPtr(joint_name)) |q| return try cloneQuery(self.allocator, q);
        if (self.tmp_cachefo.getPtr(joint_name)) |q| return try cloneQuery(self.allocator, q);

        // NOTE: getCachefo returns a pointer into tmp_cachefo; the second call
        // may insert and resize the map, invalidating the first pointer.
        // So we must call both first, then re-lookup to get stable pointers.
        _ = try self.getCachefo(name_a);
        _ = try self.getCachefo(name_b);
        const ref_a = self.cachefo.getPtr(name_a) orelse self.tmp_cachefo.getPtr(name_a) orelse return error.MissingSeqCache;
        const ref_b = self.cachefo.getPtr(name_b) orelse self.tmp_cachefo.getPtr(name_b) orelse return error.MissingSeqCache;

        if (ref_a.cdr3_length != ref_b.cdr3_length)
            return error.Cdr3LengthMismatch;

        // Union of only_genes
        var joint_genes: std.StringHashMapUnmanaged(void) = .{};
        defer {
            var it = joint_genes.iterator();
            while (it.next()) |e| self.allocator.free(e.key_ptr.*);
            joint_genes.deinit(self.allocator);
        }
        for (ref_a.only_genes.items) |g| {
            if (!joint_genes.contains(g))
                try joint_genes.put(self.allocator, try self.allocator.dupe(u8, g), {});
        }
        for (ref_b.only_genes.items) |g| {
            if (!joint_genes.contains(g))
                try joint_genes.put(self.allocator, try self.allocator.dupe(u8, g), {});
        }
        var genes_list: std.ArrayListUnmanaged([]const u8) = .{};
        defer genes_list.deinit(self.allocator);
        {
            var git = joint_genes.iterator();
            while (git.next()) |e| {
                try genes_list.append(self.allocator, e.key_ptr.*);
            }
        }
        // Sort to match C++ set<string> iteration order
        std.mem.sort([]const u8, genes_list.items, {}, struct {
            fn lessThan(_: void, a: []const u8, b: []const u8) bool {
                return std.mem.order(u8, a, b) == .lt;
            }
        }.lessThan);

        const joint_kbounds = ref_a.kbounds.logicalOr(ref_b.kbounds);
        const n_a: f64 = @floatFromInt(ref_a.nSeqs());
        const n_b: f64 = @floatFromInt(ref_b.nSeqs());
        const joint_mute_freq: f32 = @floatCast((n_a * ref_a.mute_freq + n_b * ref_b.mute_freq) / (n_a + n_b));
        const is_seed_missing = !ham_text.inString(self.args.seed_unique_id, joint_name, ":");

        const key_seqs = try self.getSeqs(joint_name);
        defer {
            for (key_seqs.items) |*ks| ks.deinit(self.allocator);
            var ksl = key_seqs;
            ksl.deinit(self.allocator);
        }

        var q = try Query.create(
            self.allocator,
            joint_name,
            key_seqs.items,
            is_seed_missing,
            genes_list.items,
            joint_kbounds,
            joint_mute_freq,
            ref_a.cdr3_length,
            name_a, name_b,
        );
        errdefer q.deinit(self.allocator);

        const tmp_key = try self.allocator.dupe(u8, joint_name);
        errdefer self.allocator.free(tmp_key);
        const q_copy = try cloneQuery(self.allocator, &q);
        try self.tmp_cachefo.put(self.allocator, tmp_key, q);

        return q_copy;
    }

    // ── Name translation helpers ─────────────────────────────────────────────

    /// C++: Glomerator::GetNaiveSeqNameToCalculate() — glomerator.cc:549
    fn getNaiveSeqNameToCalculate(self: *Glomerator, actual_queries: []const u8) ![]u8 {
        if (self.naive_seq_name_translations.get(actual_queries)) |t| return try self.allocator.dupe(u8, t);

        const n = countMembersWithExcludeSeeds(actual_queries, self.args.seed_unique_id);
        const threshold = @as(i32, @intCast(self.args.biggest_naive_seq_cluster_to_calculate));
        // C++: CountMembers(...) < 1.5 * biggest_naive_seq_cluster_to_calculate()
        // Must use float comparison to match C++ behavior exactly.
        if (@as(f64, @floatFromInt(n)) < 1.5 * @as(f64, @floatFromInt(threshold))) {
            return try self.allocator.dupe(u8, actual_queries);
        }

        const sub = try self.chooseSubsetOfNames(actual_queries, threshold);
        if (self.args.debug > 0) {
            const dbg = try std.fmt.allocPrint(self.allocator, "                translate for naive seq  {s}  -->  {s}\n", .{ actual_queries, sub });
            defer self.allocator.free(dbg);
            try std.fs.File.stdout().writeAll(dbg);
        }
        const key = try self.allocator.dupe(u8, actual_queries);
        errdefer self.allocator.free(key);
        const val = try self.allocator.dupe(u8, sub);
        try self.naive_seq_name_translations.put(self.allocator, key, val);
        return try self.allocator.dupe(u8, sub);
    }

    /// C++: Glomerator::GetLogProbPairOfNamesToCalculate() — glomerator.cc:584
    fn getLogProbPairOfNamesToCalculate(self: *Glomerator, actual_queries: []const u8, actual_parents: ?*const [2][]u8) ![2][]u8 {
        if (self.logprob_name_translations.get(actual_queries)) |t| {
            return .{
                try self.allocator.dupe(u8, t[0]),
                try self.allocator.dupe(u8, t[1]),
            };
        }

        const n_max: i32 = @intCast(self.args.biggest_logprob_cluster_to_calculate);
        if (countMembers(actual_queries) <= 2 * n_max) {
            if (actual_parents) |parents| {
                return .{
                    try self.allocator.dupe(u8, parents[0]),
                    try self.allocator.dupe(u8, parents[1]),
                };
            }
            return .{
                try self.allocator.dupe(u8, actual_queries),
                try self.allocator.dupe(u8, actual_queries),
            };
        }

        const ap1 = if (actual_parents) |p| p[0] else actual_queries;
        const ap2 = if (actual_parents) |p| p[1] else actual_queries;
        const p1 = if (actual_parents) |p| try self.getLogProbNameToCalculate(p[0], n_max) else try self.allocator.dupe(u8, actual_queries);
        const p2 = if (actual_parents) |p| try self.getLogProbNameToCalculate(p[1], n_max) else try self.allocator.dupe(u8, actual_queries);
        const pair = [2][]u8{ p1, p2 };

        if (self.args.debug > 0) {
            const dbg = try std.fmt.allocPrint(self.allocator, "                translate for lratio ({s})   {s}  {s}  -->  {s}  {s}\n", .{ actual_queries, ap1, ap2, p1, p2 });
            defer self.allocator.free(dbg);
            try std.fs.File.stdout().writeAll(dbg);
        }

        const key = try self.allocator.dupe(u8, actual_queries);
        errdefer self.allocator.free(key);
        try self.logprob_name_translations.put(self.allocator, key, .{
            try self.allocator.dupe(u8, p1),
            try self.allocator.dupe(u8, p2),
        });

        return pair;
    }

    /// C++: Glomerator::GetLogProbNameToCalculate() — glomerator.cc:568
    fn getLogProbNameToCalculate(self: *Glomerator, queries: []const u8, n_max: i32) ![]u8 {
        var q = queries;
        if (self.logprob_asymetric_translations.get(queries)) |t| {
            if (self.args.debug > 0) {
                const dbg = try std.fmt.allocPrint(self.allocator, "             using asymetric translation  {s}  -->  {s}\n", .{ queries, t });
                defer self.allocator.free(dbg);
                try std.fs.File.stdout().writeAll(dbg);
            }
            q = t;
        }
        if (countMembersWithExcludeSeeds(q, self.args.seed_unique_id) > n_max) {
            return self.chooseSubsetOfNames(q, n_max);
        }
        return self.allocator.dupe(u8, q);
    }

    /// C++: Glomerator::ChooseSubsetOfNames() — glomerator.cc:490
    fn chooseSubsetOfNames(self: *Glomerator, queries: []const u8, n_max_in: i32) ![]u8 {
        if (self.name_subsets.get(queries)) |ns| return try self.allocator.dupe(u8, ns);

        var namevec = try splitName(self.allocator, queries);
        defer {
            for (namevec.items) |n| self.allocator.free(n);
            namevec.deinit(self.allocator);
        }

        // C++: set<string> nameset(namevector.begin(), namevector.end());
        // Deduplicate names — with seed clustering you can get many duplicate seqs
        var nameset = std.StringHashMapUnmanaged(void){};
        defer nameset.deinit(self.allocator);
        for (namevec.items) |name| try nameset.put(self.allocator, name, {});
        var n_max: usize = @intCast(n_max_in);
        if (nameset.count() < n_max) // C++: if(nameset.size() < unsigned(n_max))
            n_max = nameset.count();

        // C++: srand(hash<string>{}(queries)) — GCC's std::hash<string> is a 64-bit
        // MurmurHash variant, and srand() truncates to unsigned int (32-bit)
        const hash_val = gccStringHash(queries);
        const seed_val: u32 = @truncate(hash_val);
        var rand = GlibcRand.init(seed_val);

        // C++ uses set<int> ichosen + set<string> chosen_strs for deduplication,
        // plus vector<int> ichosen_vec for ordering
        var ichosen = std.AutoHashMapUnmanaged(usize, void){};
        defer ichosen.deinit(self.allocator);
        var ichosen_vec = std.ArrayListUnmanaged(usize){};
        defer ichosen_vec.deinit(self.allocator);
        var chosen_strs = std.StringHashMapUnmanaged(void){};
        defer chosen_strs.deinit(self.allocator);

        for (0..n_max) |_| {
            var n_tries: usize = 0;
            while (true) {
                // C++: ich = rand() % namevector.size()
                const ich: usize = rand.next() % namevec.items.len;
                n_tries += 1;
                if (ichosen.contains(ich) or chosen_strs.contains(namevec.items[ich])) {
                    if (n_tries > 1_000_000) return error.TooManyTriesInChooseSubset;
                    continue;
                }
                try ichosen.put(self.allocator, ich, {});
                try chosen_strs.put(self.allocator, namevec.items[ich], {});
                try ichosen_vec.append(self.allocator, ich);
                break;
            }
        }

        // C++: sort(ichosen_vec.begin(), ichosen_vec.end())
        std.mem.sort(usize, ichosen_vec.items, {}, std.sort.asc(usize));

        var sub_names: std.ArrayListUnmanaged([]const u8) = .{};
        defer sub_names.deinit(self.allocator);
        for (ichosen_vec.items) |i| try sub_names.append(self.allocator, namevec.items[i]);

        const subqueries = try ham_text.joinStrings(self.allocator, sub_names.items, ":");
        errdefer self.allocator.free(subqueries);

        // C++ caches the subset in tmp_cachefo_ using the parent cluster's attributes
        // (only_genes, kbounds, mute_freq, cdr3_length). This is important because
        // later getCachefo() calls for the subset will find this entry instead of
        // rebuilding from singles, and the parent's mute_freq (which carries f32
        // merge cascade truncation) must be used for consistency with C++.
        const cacheref = try self.getCachefo(queries);
        {
            const sub_seqs = try self.getSeqs(subqueries);
            defer {
                var ss = sub_seqs;
                ss.deinit(self.allocator);
            }
            const is_seed_missing = !ham_text.inString(self.args.seed_unique_id, subqueries, ":");
            var sub_q = try Query.create(
                self.allocator,
                subqueries,
                sub_seqs.items,
                is_seed_missing,
                cacheref.only_genes.items,
                cacheref.kbounds,
                cacheref.mute_freq,
                cacheref.cdr3_length,
                null, null,
            );
            errdefer sub_q.deinit(self.allocator);
            const cf_key = try self.allocator.dupe(u8, subqueries);
            errdefer self.allocator.free(cf_key);
            try self.tmp_cachefo.put(self.allocator, cf_key, sub_q);
        }

        if (self.args.debug > 0) {
            const dbg = try std.fmt.allocPrint(self.allocator, "                chose subset  {s}  -->  {s}\n", .{ queries, subqueries });
            defer self.allocator.free(dbg);
            try std.fs.File.stdout().writeAll(dbg);
        }

        const key = try self.allocator.dupe(u8, queries);
        errdefer self.allocator.free(key);
        try self.name_subsets.put(self.allocator, key, subqueries);

        return try self.allocator.dupe(u8, subqueries);
    }

    /// C++: Glomerator::FindNaiveSeqNameReplace() — glomerator.cc:623
    fn findNaiveSeqNameReplace(self: *Glomerator, parents: *const [2][]const u8) !?[]u8 {
        const ns_a = try self.getNaiveSeq(parents[0], null);
        const ns_b = try self.getNaiveSeq(parents[1], null);
        if (std.mem.eql(u8, ns_a, ns_b)) return try self.allocator.dupe(u8, parents[0]);

        const nmax: i32 = @intCast(@as(u64, @intCast(self.args.biggest_naive_seq_cluster_to_calculate)) * 3 / 2);
        if (try self.firstParentMuchBigger(parents[0], parents[1], nmax)) return try self.allocator.dupe(u8, parents[0]);
        if (try self.firstParentMuchBigger(parents[1], parents[0], nmax)) return try self.allocator.dupe(u8, parents[1]);
        return null;
    }

    /// C++: Glomerator::FirstParentMuchBigger() — glomerator.cc:608
    fn firstParentMuchBigger(self: *Glomerator, queries: []const u8, queries_other: []const u8, nmax: i32) !bool {
        const nseq: i32 = countMembers(queries);
        const nseq_other: i32 = countMembers(queries_other);
        const result = nseq > nmax and @as(f64, @floatFromInt(nseq)) / @as(f64, @floatFromInt(nseq_other)) > self.asym_factor;
        if (result and self.args.debug > 0) {
            const joined = try self.joinNames(queries, queries_other);
            defer self.allocator.free(joined);
            const dbg = try std.fmt.allocPrint(self.allocator, "                asymetric  {d} {d}  use {s}  instead of {s}\n", .{ nseq, nseq_other, queries, joined });
            defer self.allocator.free(dbg);
            try std.fs.File.stdout().writeAll(dbg);
            if (self.naive_seq_name_translations.get(queries)) |t| {
                const dbg2 = try std.fmt.allocPrint(self.allocator, "                    naive seq translates to {s}\n", .{t});
                defer self.allocator.free(dbg2);
                try std.fs.File.stdout().writeAll(dbg2);
            }
        }
        return result;
    }

    /// C++: Glomerator::UpdateLogProbTranslationsForAsymetrics() — glomerator.cc:1121
    fn updateLogProbTranslationsForAsymetrics(self: *Glomerator, qmerge: *const Query) !void {
        const parents = qmerge.parents orelse return;
        const nmax: i32 = @intCast(@as(u64, @intCast(self.args.biggest_logprob_cluster_to_calculate)) * 3 / 2);

        var big_parent: ?[]const u8 = null;
        if (try self.firstParentMuchBigger(parents[0], parents[1], nmax))
            big_parent = parents[0]
        else if (try self.firstParentMuchBigger(parents[1], parents[0], nmax))
            big_parent = parents[1];

        if (big_parent == null) return;
        const bp = big_parent.?;

        // Follow chain of translations
        var subqueries = bp;
        while (self.logprob_asymetric_translations.get(subqueries)) |t| {
            if (self.args.debug > 0) {
                const dbg = try std.fmt.allocPrint(self.allocator, "                  turtles {s}  -->  {s}\n", .{ subqueries, t });
                defer self.allocator.free(dbg);
                try std.fs.File.stdout().writeAll(dbg);
            }
            subqueries = t;
        }

        const ratio = @as(f64, @floatFromInt(countMembers(bp))) /
            @as(f64, @floatFromInt(countMembers(subqueries)));
        if (ratio < 2.0) {
            if (self.args.debug > 0) {
                const dbg = try std.fmt.allocPrint(self.allocator, "                logprob asymetric translation  {s}  -->  {s}\n", .{ qmerge.name, subqueries });
                defer self.allocator.free(dbg);
                try std.fs.File.stdout().writeAll(dbg);
            }
            const key = try self.allocator.dupe(u8, qmerge.name);
            errdefer self.allocator.free(key);
            const val = try self.allocator.dupe(u8, subqueries);
            errdefer self.allocator.free(val);
            try self.logprob_asymetric_translations.put(self.allocator, key, val);
        } else {
            if (self.args.debug > 0) {
                const dbg = try std.fmt.allocPrint(self.allocator, "                  ratio too big for asymetric {d} {d}\n", .{ countMembers(bp), countMembers(subqueries) });
                defer self.allocator.free(dbg);
                try std.fs.File.stdout().writeAll(dbg);
            }
        }
    }

    /// C++: Glomerator::MoveSubsetsFromTmpCache() — glomerator.cc:865
    fn moveSubsetsFromTmpCache(self: *Glomerator, query: []const u8) !void {
        if (self.naive_seq_name_translations.get(query)) |tq| {
            try self.copyToPermanentCache(tq, query);
        }
        if (self.logprob_name_translations.get(query)) |tpair| {
            try self.copyToPermanentCache(tpair[0], query);
            try self.copyToPermanentCache(tpair[1], query);
        }
        if (self.logprob_asymetric_translations.get(query)) |tq| {
            try self.copyToPermanentCache(tq, query);
        }
    }

    /// C++: Glomerator::CopyToPermanentCache() — glomerator.cc:889
    fn copyToPermanentCache(self: *Glomerator, translated_query: []const u8, superquery: []const u8) !void {
        if (self.tmp_cachefo.getPtr(translated_query)) |tq| {
            const key = try self.allocator.dupe(u8, translated_query);
            errdefer self.allocator.free(key);
            const val = try cloneQuery(self.allocator, tq);
            try self.cachefo.put(self.allocator, key, val);
        } else {
            const supercache = try self.getCachefo(superquery);
            const key_seqs = try self.getSeqs(translated_query);
            defer {
                var ks = key_seqs;
                ks.deinit(self.allocator);
            }
            const is_seed_missing = !ham_text.inString(self.args.seed_unique_id, translated_query, ":");
            var q = try Query.create(
                self.allocator,
                translated_query,
                key_seqs.items,
                is_seed_missing,
                supercache.only_genes.items,
                supercache.kbounds,
                supercache.mute_freq,
                supercache.cdr3_length,
                null, null,
            );
            errdefer q.deinit(self.allocator);
            const cf_key = try self.allocator.dupe(u8, translated_query);
            errdefer self.allocator.free(cf_key);
            try self.cachefo.put(self.allocator, cf_key, q);
        }
    }

    // ── Utility helpers ───────────────────────────────────────────────────────

    /// C++: Glomerator::GetSeqs() — glomerator.cc:850
    fn getSeqs(self: *Glomerator, query: []const u8) !std.ArrayListUnmanaged(Sequence) {
        const names = try splitName(self.allocator, query);
        defer {
            for (names.items) |n| self.allocator.free(n);
            var nl = names;
            nl.deinit(self.allocator);
        }

        var result: std.ArrayListUnmanaged(Sequence) = .{};
        errdefer {
            for (result.items) |*sq| sq.deinit(self.allocator);
            result.deinit(self.allocator);
        }
        for (names.items) |uid| {
            const sq = self.single_seqs.getPtr(uid) orelse return error.SequenceNotFound;
            try result.append(self.allocator, try sq.clone(self.allocator));
        }
        return result;
    }

    /// C++: Glomerator::AddFailedQuery() — glomerator.cc:776
    fn addFailedQuery(self: *Glomerator, queries: []const u8, error_str: []const u8) !void {
        // Append error string
        const old_err = self.errors.get(queries) orelse "";
        const new_err = try std.fmt.allocPrint(self.allocator, "{s}:{s}", .{ old_err, error_str });
        errdefer self.allocator.free(new_err);

        if (self.errors.fetchRemove(queries)) |old| {
            self.allocator.free(old.key);
            self.allocator.free(old.value);
        }
        const err_key = try self.allocator.dupe(u8, queries);
        errdefer self.allocator.free(err_key);
        try self.errors.put(self.allocator, err_key, new_err);

        if (!self.failed_queries.contains(queries)) {
            const fq_key = try self.allocator.dupe(u8, queries);
            try self.failed_queries.put(self.allocator, fq_key, {});
        }
    }

    /// C++: Glomerator::JoinNames() — glomerator.cc:404
    fn joinNames(self: *Glomerator, name1: []const u8, name2: []const u8) ![]u8 {
        var names = [_][]const u8{ name1, name2 };
        std.mem.sort([]const u8, &names, {}, struct {
            fn lessThan(_: void, a: []const u8, b: []const u8) bool {
                return std.mem.lessThan(u8, a, b);
            }
        }.lessThan);
        return std.fmt.allocPrint(self.allocator, "{s}:{s}", .{ names[0], names[1] });
    }

    fn geneListToSlice(self: *Glomerator, genes: []const []u8) ![]const []const u8 {
        const result = try self.allocator.alloc([]const u8, genes.len);
        for (genes, result) |g, *r| r.* = g;
        return result;
    }
};

// ── Module-level helpers ──────────────────────────────────────────────────────

fn splitName(allocator: std.mem.Allocator, query: []const u8) !std.ArrayListUnmanaged([]u8) {
    var result: std.ArrayListUnmanaged([]u8) = .{};
    var it = std.mem.splitScalar(u8, query, ':');
    while (it.next()) |part| {
        if (part.len > 0) try result.append(allocator, try allocator.dupe(u8, part));
    }
    return result;
}

/// Return "len(N)" if the cluster has >= 10 members, otherwise the full string.
/// Corresponds to C++ `Glomerator::PrintStr`.
/// Uses alternating static buffers since this may be called twice per print.
fn printStr(queries: []const u8) []const u8 {
    if (countMembers(queries) < 10) return queries;
    const S = struct {
        var bufs: [2][32]u8 = .{ undefined, undefined };
        var idx: u1 = 0;
    };
    const i = S.idx;
    S.idx +%= 1;
    const s = std.fmt.bufPrint(&S.bufs[i], "len({d})", .{countMembers(queries)}) catch return queries;
    return s;
}

/// C++: Glomerator::CountMembers() — glomerator.cc:360
fn countMembers(namestr: []const u8) i32 {
    var n: i32 = 1;
    for (namestr) |c| {
        if (c == ':') n += 1;
    }
    return n;
}

fn countMembersWithExcludeSeeds(namestr: []const u8, seed_uid: []const u8) i32 {
    const n_members = countMembers(namestr);
    if (seed_uid.len == 0) return n_members;
    // Count how many times seed_uid appears as a member, then subtract extras
    // (keep one copy of the seed in the count, exclude duplicates)
    var n_seeds: i32 = 0;
    var it = std.mem.splitScalar(u8, namestr, ':');
    while (it.next()) |part| {
        if (std.mem.eql(u8, part, seed_uid)) n_seeds += 1;
    }
    if (n_seeds > 0) return n_members - (n_seeds - 1);
    return n_members;
}

/// C++: Glomerator::LargestClusterSize() — glomerator.cc:374
fn largestClusterSize(partition: *const Partition) u32 {
    var max: u32 = 0;
    for (partition.items) |cluster| {
        const sz: u32 = @intCast(countMembers(cluster));
        if (sz > max) max = sz;
    }
    return max;
}

/// GCC libstdc++ `std::hash<std::string>` on 64-bit Linux.
/// Based on MurmurHashUnaligned2 variant from gcc/libstdc++-v3/libsupc++/hash_bytes.cc,
/// `_Hash_bytes` function and `load_bytes` helper. Seed is 0xc70f6907 (from functional_hash.h).
/// Only valid on little-endian architectures (x86-64).
fn gccStringHash(s: []const u8) u64 {
    comptime if (@import("builtin").cpu.arch.endian() != .little)
        @compileError("gccStringHash assumes little-endian (x86-64)");
    const mul: u64 = 0xc6a4a7935bd1e995;
    const seed: u64 = 0xc70f6907;
    const len = s.len;

    var hash: u64 = seed ^ (len *% mul);

    // Process 8 bytes at a time
    var i: usize = 0;
    while (i + 8 <= len) : (i += 8) {
        const data = shiftMix(std.mem.readInt(u64, s[i..][0..8], .little) *% mul) *% mul;
        hash ^= data;
        hash *%= mul;
    }

    // Tail: remaining bytes
    const remaining = len & 7;
    if (remaining != 0) {
        // Matches libstdc++ load_bytes (hash_bytes.cc): packs remaining
        // bytes into u64 with last byte most-significant.
        var data: u64 = 0;
        var j: usize = remaining;
        while (j > 0) {
            j -= 1;
            data = (data << 8) + @as(u64, s[i + j]);
        }
        hash ^= data;
        hash *%= mul;
    }

    // Finalization
    hash = shiftMix(hash) *% mul;
    hash = shiftMix(hash);
    return hash;
}

fn shiftMix(v: u64) u64 {
    return v ^ (v >> 47);
}

/// glibc TYPE_3 PRNG — matches srand()/rand() on Linux.
/// State: 31 words, separation 3, trinomial x^31 + x^3 + 1.
const GlibcRand = struct {
    state: [31]i32,
    fptr: usize, // index into state (initially 3)
    rptr: usize, // index into state (initially 0)

    fn init(seed_val: u32) GlibcRand {
        var self = GlibcRand{ .state = undefined, .fptr = 3, .rptr = 0 };
        // glibc treats seed=0 as seed=1. @bitCast reinterprets u32 as i32 (may
        // be negative for seeds > 0x7FFFFFFF), matching glibc's cast from
        // unsigned int to int32_t in __srandom_r.
        const seed: i32 = if (seed_val == 0) 1 else @bitCast(seed_val);
        self.state[0] = seed;

        // Fill state[1..30] with Park-Miller LCG: x = 16807*x mod 2^31-1
        var word: i32 = seed;
        for (1..31) |i| {
            const lo: i64 = @as(i64, @rem(word, 127773)) * 16807;
            const hi: i64 = @as(i64, @divTrunc(word, 127773)) * 2836;
            var val: i64 = lo - hi;
            if (val < 0) val += 2147483647;
            word = @intCast(val);
            self.state[i] = word;
        }

        // Warm up: discard 310 values (DEG_3 * 10)
        for (0..310) |_| _ = self.next();

        return self;
    }

    /// Returns a 31-bit non-negative value (0..2^31-1), matching glibc rand().
    fn next(self: *GlibcRand) u32 {
        // Additive feedback: state[fptr] += state[rptr]
        const val: u32 = @as(u32, @bitCast(self.state[self.fptr])) +% @as(u32, @bitCast(self.state[self.rptr]));
        self.state[self.fptr] = @bitCast(val);
        const result: u32 = val >> 1; // drop LSB, return 31-bit value

        // Advance pointers cyclically
        self.fptr += 1;
        if (self.fptr >= 31) {
            self.fptr = 0;
            self.rptr += 1;
        } else {
            self.rptr += 1;
            if (self.rptr >= 31) self.rptr = 0;
        }

        return result;
    }
};

fn cloneQuery(allocator: std.mem.Allocator, q: *const Query) !Query {
    return Query.create(
        allocator,
        q.name,
        q.seqs.items,
        q.seed_missing,
        q.only_genes.items,
        q.kbounds,
        q.mute_freq,
        q.cdr3_length,
        if (q.parents) |p| p[0] else null,
        if (q.parents) |p| p[1] else null,
    );
}
