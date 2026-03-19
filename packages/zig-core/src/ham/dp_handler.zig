/// ham/dp_handler.zig — Zig port of ham/src/dphandler.cc + ham/include/dphandler.h
///
/// Dynamic programming handler: runs Forward + Viterbi over a kset grid.
/// Core bcrham algorithm.
///
/// C++ source: packages/ham/src/dphandler.cc, packages/ham/include/dphandler.h
/// C++ author: psathyrella/ham

const std = @import("std");
const Trellis = @import("trellis.zig").Trellis;
const TracebackPath = @import("traceback_path.zig").TracebackPath;
const Model = @import("model.zig").Model;
const Sequences = @import("sequences.zig").Sequences;
const Sequence = @import("sequences.zig").Sequence;
const mathutils = @import("mathutils.zig");
const Args = @import("args.zig").Args;
const bcr = @import("bcr_utils/root.zig");
const GermLines = bcr.GermLines;
const HMMHolder = bcr.HMMHolder;
const KSet = bcr.KSet;
const KBounds = bcr.KBounds;
const RecoEvent = bcr.RecoEvent;
const Result = bcr.Result;
const Insertions = bcr.Insertions;
const Region = bcr.Region;
const regionStr = bcr.regionStr;
const TermColors = bcr.TermColors;

/// Key for the per-gene score/path caches.
const GeneKSetKey = struct {
    gene: []const u8,
    kset: KSet,
};

/// Cached trellis entry: query strings → Trellis
const TrellisEntry = struct {
    query_strs: std.ArrayListUnmanaged([]u8),
    trellis: Trellis,
};

/// Corresponds to C++ `ham::DPHandler`.
pub const DPHandler = struct {
    algorithm: []const u8,
    args: *Args,
    gl: *GermLines,
    hmms: *HMMHolder,
    allocator: std.mem.Allocator,

    /// scratch_cachefo_: gene → list of (query-string-vec, Trellis) entries
    scratch_cachefo: std.StringHashMapUnmanaged(std.ArrayListUnmanaged(TrellisEntry)),
    /// paths_: gene → (KSet → TracebackPath)
    paths: std.StringHashMapUnmanaged(std.AutoHashMapUnmanaged(KSet, TracebackPath)),
    /// scores_: gene → (KSet → f64)
    scores: std.StringHashMapUnmanaged(std.AutoHashMapUnmanaged(KSet, f64)),
    /// per_gene_support_: gene → best full-annotation log-prob
    per_gene_support: std.StringHashMapUnmanaged(f64),

    pub fn init(
        allocator: std.mem.Allocator,
        algorithm: []const u8,
        args: *Args,
        gl: *GermLines,
        hmms: *HMMHolder,
    ) !DPHandler {
        return DPHandler{
            .algorithm = algorithm,
            .args = args,
            .gl = gl,
            .hmms = hmms,
            .allocator = allocator,
            .scratch_cachefo = .{},
            .paths = .{},
            .scores = .{},
            .per_gene_support = .{},
        };
    }

    pub fn deinit(self: *DPHandler) void {
        self.clear();
    }

    /// Clear all cached trellis data.
    /// Corresponds to C++ `DPHandler::Clear`.
    pub fn clear(self: *DPHandler) void {
        const allocator = self.allocator;

        // Free scratch_cachefo
        var cit = self.scratch_cachefo.iterator();
        while (cit.next()) |entry| {
            allocator.free(entry.key_ptr.*);
            for (entry.value_ptr.items) |*te| {
                for (te.query_strs.items) |s| allocator.free(s);
                te.query_strs.deinit(allocator);
                te.trellis.deinit();
            }
            entry.value_ptr.deinit(allocator);
        }
        self.scratch_cachefo.deinit(allocator);
        self.scratch_cachefo = .{};

        // Free paths
        var pit = self.paths.iterator();
        while (pit.next()) |entry| {
            allocator.free(entry.key_ptr.*);
            var inner_it = entry.value_ptr.iterator();
            while (inner_it.next()) |kv| kv.value_ptr.deinit(allocator);
            entry.value_ptr.deinit(allocator);
        }
        self.paths.deinit(allocator);
        self.paths = .{};

        // Free scores
        var sit = self.scores.iterator();
        while (sit.next()) |entry| {
            allocator.free(entry.key_ptr.*);
            entry.value_ptr.deinit(allocator);
        }
        self.scores.deinit(allocator);
        self.scores = .{};

        // Free per_gene_support
        var pgit = self.per_gene_support.iterator();
        while (pgit.next()) |entry| allocator.free(entry.key_ptr.*);
        self.per_gene_support.deinit(allocator);
        self.per_gene_support = .{};
    }

    /// Run the DP for a single Sequence.
    /// Corresponds to C++ `DPHandler::Run(Sequence, ...)`.
    pub fn runSeq(
        self: *DPHandler,
        seq: Sequence,
        kbounds: KBounds,
        only_gene_list: []const []const u8,
        overall_mute_freq: f64,
        clear_cache: bool,
    ) !Result {
        // Clone into a Sequences container (single sequence)
        var seqvec: std.ArrayListUnmanaged(Sequence) = .{};
        defer seqvec.deinit(self.allocator);
        try seqvec.append(self.allocator, try seq.clone(self.allocator));
        return self.run(seqvec.items, kbounds, only_gene_list, overall_mute_freq, clear_cache);
    }

    /// Run the DP over a vector of sequences.
    /// Corresponds to C++ `DPHandler::Run(vector<Sequence>, ...)`.
    pub fn run(
        self: *DPHandler,
        seqvector: []const Sequence,
        kbounds: KBounds,
        only_gene_list: []const []const u8,
        overall_mute_freq: f64,
        clear_cache: bool,
    ) !Result {
        const allocator = self.allocator;
        const run_start = std.time.milliTimestamp();

        var seqs = Sequences.init();
        defer seqs.deinit(allocator);
        for (seqvector) |*sq| {
            try seqs.addSeq(allocator, try sq.clone(allocator));
        }

        // Build only_genes map: region → set of gene names
        var only_genes: std.AutoHashMapUnmanaged(Region, std.StringHashMapUnmanaged(void)) = .{};
        defer {
            var oit = only_genes.iterator();
            while (oit.next()) |entry| {
                var inner_it = entry.value_ptr.iterator();
                while (inner_it.next()) |kv| allocator.free(kv.key_ptr.*);
                entry.value_ptr.deinit(allocator);
            }
            only_genes.deinit(allocator);
        }

        if (only_gene_list.len > 0) {
            for (bcr.germ_lines.regions) |r| {
                try only_genes.put(allocator, r, .{});
            }
            for (only_gene_list) |gene| {
                const r = try GermLines.getRegion(gene);
                if (only_genes.getPtr(r)) |set| {
                    const gkey = try allocator.dupe(u8, gene);
                    errdefer allocator.free(gkey);
                    try set.put(allocator, gkey, {});
                }
            }
            // Validate each region has at least one gene
            for (bcr.germ_lines.regions) |r| {
                const set = only_genes.get(r) orelse return error.NoGenesForRegion;
                if (set.count() == 0) return error.NoGenesForRegion;
            }
        } else {
            // Populate from GermLines
            for (bcr.germ_lines.regions) |r| {
                var set: std.StringHashMapUnmanaged(void) = .{};
                if (self.gl.names.get(regionStr(r))) |name_list| {
                    for (name_list.items) |gene| {
                        try set.put(allocator, try allocator.dupe(u8, gene), {});
                    }
                }
                try only_genes.put(allocator, r, set);
            }
        }

        if (kbounds.vmin == 0 or kbounds.dmin == 0 or
            kbounds.vmax <= kbounds.vmin or kbounds.dmax <= kbounds.dmin)
            return error.TrivialKBounds;

        if (clear_cache) self.clear();

        var best_scores: std.AutoHashMapUnmanaged(KSet, f64) = .{};
        defer best_scores.deinit(allocator);
        var total_scores: std.AutoHashMapUnmanaged(KSet, f64) = .{};
        defer total_scores.deinit(allocator);
        var best_genes: std.AutoHashMapUnmanaged(KSet, std.AutoHashMapUnmanaged(Region, []u8)) = .{};
        defer {
            var bgit = best_genes.iterator();
            while (bgit.next()) |entry| {
                var inner = entry.value_ptr.iterator();
                while (inner.next()) |kv| {
                    allocator.free(kv.value_ptr.*);
                }
                entry.value_ptr.deinit(allocator);
            }
            best_genes.deinit(allocator);
        }

        if (!self.args.dont_rescale_emissions) {
            try self.hmms.rescaleOverallMuteFreqs(&only_genes, overall_mute_freq);
        }

        var result = try Result.init(allocator, kbounds, self.args.locus);

        var best_score: f64 = mathutils.NEG_INF;
        var best_kset = KSet{ .v = 0, .d = 0 };
        var n_too_long: usize = 0;
        var n_run: usize = 0;
        var n_total: usize = 0;

        // Loop over k-space in reverse order (for chunk caching)
        var k_v: usize = kbounds.vmax - 1;
        while (true) {
            var k_d: usize = kbounds.dmax - 1;
            while (true) {
                n_total += 1;
                if (k_v + k_d >= seqs.sequence_length) {
                    n_too_long += 1;
                } else {
                    const kset = KSet{ .v = k_v, .d = k_d };
                    try self.runKSet(&seqs, kset, &only_genes, &best_scores, &total_scores, &best_genes);
                    n_run += 1;
                    const total_kset = total_scores.get(kset) orelse mathutils.NEG_INF;
                    result.total_score = mathutils.addInLogSpace(total_kset, result.total_score);
                    if (self.args.debug == 2 and std.mem.eql(u8, self.algorithm, "forward")) {
                        const fwd_msg = try std.fmt.allocPrint(allocator, "            {d: >9.2} ({e:.1})  tot: {d: >7.2}\n", .{ total_kset, @as(f64, @exp(total_kset)), result.total_score });
                        defer allocator.free(fwd_msg);
                        try std.fs.File.stdout().writeAll(fwd_msg);
                    }

                    const best_kset_score = best_scores.get(kset) orelse mathutils.NEG_INF;
                    if (best_kset_score > best_score) {
                        best_score = best_kset_score;
                        best_kset = kset;
                    }

                    if (std.mem.eql(u8, self.algorithm, "viterbi") and best_kset_score != mathutils.NEG_INF) {
                        const bg = best_genes.get(kset) orelse continue;
                        const event = try self.fillRecoEvent(&seqs, kset, &bg, best_kset_score);
                        try result.pushBackRecoEvent(event);
                    }
                }
                if (k_d == kbounds.dmin) break;
                k_d -= 1;
            }
            if (k_v == kbounds.vmin) break;
            k_v -= 1;
        }
        if (self.args.debug > 0 and n_too_long > 0) {
            const skip_msg = try std.fmt.allocPrint(allocator, "      skipped {d} (of {d}) k sets 'cause they were longer than the sequence (ran {d})\n", .{ n_too_long, n_total, n_run });
            defer allocator.free(skip_msg);
            try std.fs.File.stdout().writeAll(skip_msg);
        }

        // No valid path
        if (best_kset.v == 0 and best_kset.d == 0) {
            const name_str = try seqs.nameStr(allocator, ":");
            defer allocator.free(name_str);
            const no_path_msg = try std.fmt.allocPrint(allocator, "    no valid paths for query {s}\n", .{name_str});
            defer allocator.free(no_path_msg);
            try std.fs.File.stdout().writeAll(no_path_msg);
            result.no_path = true;
            if (!self.args.dont_rescale_emissions) {
                try self.hmms.unRescaleOverallMuteFreqs(&only_genes);
            }
            return result;
        }

        if (std.mem.eql(u8, self.algorithm, "viterbi")) {
            try result.finalize(self.gl, &self.per_gene_support, best_kset, kbounds);
        }

        // Debug: summary line (debug >= 1)
        if (self.args.debug > 0) {
            const prob: f64 = if (std.mem.eql(u8, self.algorithm, "viterbi")) best_score else result.total_score;
            const alg_str: []const u8 = if (std.mem.eql(u8, self.algorithm, "viterbi")) "vtb" else "fwd";
            var kstr_buf: [300]u8 = undefined;
            const kstr = if (std.mem.eql(u8, self.algorithm, "viterbi"))
                try std.fmt.bufPrint(&kstr_buf, "{d} [{d}-{d})  {d} [{d}-{d})", .{ best_kset.v, kbounds.vmin, kbounds.vmax, best_kset.d, kbounds.dmin, kbounds.dmax })
            else
                try std.fmt.bufPrint(&kstr_buf, "    [{d}-{d})     [{d}-{d})", .{ kbounds.vmin, kbounds.vmax, kbounds.dmin, kbounds.dmax });
            const n_v = if (only_genes.get(.v)) |s| s.count() else 0;
            const n_d = if (only_genes.get(.d)) |s| s.count() else 0;
            const n_j = if (only_genes.get(.j)) |s| s.count() else 0;
            const name_str = try seqs.nameStr(allocator, ":");
            defer allocator.free(name_str);
            const elapsed_ms = std.time.milliTimestamp() - run_start;
            const cpu_seconds = @as(f64, @floatFromInt(elapsed_ms)) / 1000.0;
            const summary = try std.fmt.allocPrint(allocator, "           {s} {d: >12.3}   {s: <25}  {d: >2}v {d: >2}d {d: >2}j  {d: >5.2}s   {d: >4}  {s}\n", .{ alg_str, prob, kstr, n_v, n_d, n_j, cpu_seconds, seqs.nSeqs(), name_str });
            defer allocator.free(summary);
            try std.fs.File.stdout().writeAll(summary);
        }

        if (!self.args.dont_rescale_emissions) {
            try self.hmms.unRescaleOverallMuteFreqs(&only_genes);
        }

        return result;
    }

    /// Get subsequences for one region.
    /// Corresponds to C++ `DPHandler::GetSubSeqs(seqs, kset, region)`.
    fn getSubSeqs(self: *DPHandler, seqs: *const Sequences, kset: KSet, region: Region) !Sequences {
        const allocator = self.allocator;
        var result = Sequences.init();
        errdefer result.deinit(allocator);

        const k_v = kset.v;
        const k_d = kset.d;
        const start: usize = switch (region) {
            .v => 0,
            .d => k_v,
            .j => k_v + k_d,
        };
        const length: usize = switch (region) {
            .v => k_v,
            .d => k_d,
            .j => seqs.sequence_length - k_v - k_d,
        };

        for (seqs.seqs.items) |*sq| {
            var sub = try Sequence.initSlice(allocator, sq, start, length);
            errdefer sub.deinit(allocator);
            try result.addSeq(allocator, sub);
        }
        return result;
    }

    /// Collect undigitized strings for the given region's subsequences.
    fn getQueryStrs(self: *DPHandler, seqs: *const Sequences, kset: KSet, region: Region) !std.ArrayListUnmanaged([]u8) {
        const allocator = self.allocator;
        var query_seqs = try self.getSubSeqs(seqs, kset, region);
        defer query_seqs.deinit(allocator);

        var strs: std.ArrayListUnmanaged([]u8) = .{};
        errdefer {
            for (strs.items) |s| allocator.free(s);
            strs.deinit(allocator);
        }
        for (query_seqs.seqs.items) |*sq| {
            try strs.append(allocator, try allocator.dupe(u8, sq.undigitized));
        }
        return strs;
    }

    /// Ensure per-gene caches are initialised for `gene`.
    /// Three separate key copies are needed because each HashMap owns its key independently.
    fn initCache(self: *DPHandler, gene: []const u8) !void {
        if (self.scores.contains(gene)) return;

        const key = try self.allocator.dupe(u8, gene);
        errdefer self.allocator.free(key);
        try self.scratch_cachefo.put(self.allocator, key, .{});

        const key2 = try self.allocator.dupe(u8, gene);
        errdefer self.allocator.free(key2);
        try self.paths.put(self.allocator, key2, .{});

        const key3 = try self.allocator.dupe(u8, gene);
        errdefer self.allocator.free(key3);
        try self.scores.put(self.allocator, key3, .{});
    }

    /// Look for a cached kset whose region sequences are a prefix of `query_strs`.
    /// Returns null-kset if none found.
    fn findPartialCacheMatch(self: *DPHandler, region: Region, gene: []const u8, kset: KSet) KSet {
        const gene_scores = self.scores.get(gene) orelse return KSet{ .v = 0, .d = 0 };
        if (gene_scores.get(kset) != null) return kset;

        switch (region) {
            .v => {
                var it = gene_scores.iterator();
                while (it.next()) |kv| {
                    if (kv.key_ptr.v == kset.v) return kv.key_ptr.*;
                }
            },
            .j => {
                var it = gene_scores.iterator();
                while (it.next()) |kv| {
                    if (kv.key_ptr.v + kv.key_ptr.d == kset.v + kset.d) return kv.key_ptr.*;
                }
            },
            .d => {},
        }
        return KSet{ .v = 0, .d = 0 };
    }

    /// Fill the trellis for the given (gene, kset, query_seqs).
    /// Stores result in scores_[gene][kset] and (for viterbi) paths_[gene][kset].
    fn fillTrellis(
        self: *DPHandler,
        kset: KSet,
        query_seqs: Sequences,
        query_strs: []const []const u8,
        gene: []const u8,
    ) ![]const u8 {
        const allocator = self.allocator;

        // Look for a chunk-cached trellis
        var cached_trellis: ?*const Trellis = null;
        if (!self.args.no_chunk_cache) {
            if (self.scratch_cachefo.getPtr(gene)) |cache_list| {
                for (cache_list.items) |*te| {
                    if (te.query_strs.items.len != query_strs.len) continue;
                    var all_match = true;
                    for (te.query_strs.items, query_strs) |cached, current| {
                        // cached must START WITH current (prefix match)
                        if (!std.mem.startsWith(u8, cached, current)) {
                            all_match = false;
                            break;
                        }
                    }
                    if (all_match) {
                        cached_trellis = &te.trellis;
                        break;
                    }
                }
            }
        }

        const model = try self.hmms.get(gene);

        // Match C++ DPHandler::FillTrellis logic exactly:
        //   - When no cached trellis: create scratch, store it, run algorithm ON THE SCRATCH directly.
        //   - When cached trellis exists: create a temp trellis that borrows from the cache,
        //     run algorithm on the temp trellis.
        // This is critical: when no cache, do NOT create a second trellis that borrows from
        // the (empty) scratch — the scratch itself is the working trellis.

        // Build the TrellisEntry query strings (needed if storing a new scratch).
        var qstrs: std.ArrayListUnmanaged([]u8) = .{};
        defer {
            for (qstrs.items) |s| allocator.free(s);
            qstrs.deinit(allocator);
        }
        for (query_strs) |s| {
            try qstrs.append(allocator, try allocator.dupe(u8, s));
        }

        // trell_ptr points to the trellis on which the algorithm is run.
        // tmptrell is only used when borrowing from a cached trellis.
        var tmptrell: ?Trellis = null;
        defer if (tmptrell) |*t| t.deinit();

        var origin: []const u8 = "scratch";
        const trell_ptr: *Trellis = if (cached_trellis == null) blk: {
            // No cache: create scratch WITHOUT cached_trellis, store it, use it directly.
            var fresh = try Trellis.initWithSeqs(allocator, model, query_seqs, null);
            errdefer fresh.deinit();
            const entry = TrellisEntry{ .query_strs = qstrs, .trellis = fresh };
            qstrs = .{}; // ownership transferred to entry
            origin = "scratch";
            if (self.scratch_cachefo.getPtr(gene)) |cache_list| {
                try cache_list.append(allocator, entry);
                break :blk &cache_list.items[cache_list.items.len - 1].trellis;
            } else {
                // gene not in scratch_cachefo — should never happen after initCache
                unreachable;
            }
        } else blk: {
            // Have cache: create temp trellis that borrows from the cached scratch.
            tmptrell = try Trellis.initWithSeqs(allocator, model, query_seqs, cached_trellis);
            origin = "chunk";
            break :blk &tmptrell.?;
        };

        // Run the algorithm on trell_ptr
        const uncorrected_score: f64 = if (std.mem.eql(u8, self.algorithm, "viterbi")) blk: {
            try trell_ptr.viterbi();
            const sc = trell_ptr.ending_viterbi_log_prob;
            // Store traceback path
            var path = TracebackPath.initWithModel(model);
            errdefer path.deinit(allocator);
            if (sc != mathutils.NEG_INF) {
                try trell_ptr.traceback(&path);
            }
            if (self.paths.getPtr(gene)) |gene_paths| {
                try gene_paths.put(allocator, kset, path);
            } else {
                path.deinit(allocator);
            }
            break :blk sc;
        } else blk: {
            try trell_ptr.forward();
            break :blk trell_ptr.ending_forward_log_prob;
        };

        const gene_choice_score = @log(model.overall_prob);
        const final_score = mathutils.addWithMinusInfinities(uncorrected_score, gene_choice_score);

        if (self.scores.getPtr(gene)) |gene_scores| {
            try gene_scores.put(allocator, kset, final_score);
        }

        return origin;
    }

    /// Print colored germline alignment for a gene path (debug==2 viterbi output).
    /// Corresponds to C++ `DPHandler::PrintPath`.
    fn printPath(self: *DPHandler, kset: KSet, query_strs: []const []const u8, gene: []const u8, score: f64, extra_str: []const u8) !void {
        const allocator = self.allocator;
        if (score == mathutils.NEG_INF) return;

        const path_names_list = try self.getPathNames(gene, kset, allocator);
        defer {
            for (path_names_list.items) |s| allocator.free(s);
            var pl = path_names_list;
            pl.deinit(allocator);
        }
        if (path_names_list.items.len == 0) {
            if (self.args.debug > 0) {
                const msg = try std.fmt.allocPrint(allocator, "                     {s} has no valid path\n", .{gene});
                defer allocator.free(msg);
                try std.fs.File.stdout().writeAll(msg);
            }
            return;
        }

        const path_names = path_names_list.items;
        var buf_l: [200]u8 = undefined;
        const left_insert = getInsertion("left", path_names, &buf_l);
        var buf_r: [200]u8 = undefined;
        const right_insert = getInsertion("right", path_names, &buf_r);
        const left_erosion_length = self.getErosionLength("left", path_names, gene);
        const right_erosion_length = self.getErosionLength("right", path_names, gene);

        // Build modified germline: left_insert + trimmed_germline + right_insert
        const germline = self.gl.seqs.get(gene) orelse return;
        const mid_end = if (germline.len >= right_erosion_length) germline.len - right_erosion_length else 0;
        const mid = if (left_erosion_length <= mid_end) germline[left_erosion_length..mid_end] else "";

        var mg_buf: std.ArrayListUnmanaged(u8) = .{};
        defer mg_buf.deinit(allocator);
        try mg_buf.appendSlice(allocator, left_insert);
        try mg_buf.appendSlice(allocator, mid);
        try mg_buf.appendSlice(allocator, right_insert);
        const modified_germline = mg_buf.items;

        // Get ambiguous char from args (matches C++ hmms_.track()->ambiguous_char())
        const ambig = self.args.ambig_base;

        // Color the match: mutants in red relative to query strings, ambiguous in light_blue
        const match_colored = try TermColors.colorMutants(allocator, "red", modified_germline, "", query_strs, ambig);
        defer allocator.free(match_colored);
        const ambig_ch: u8 = if (ambig.len > 0) ambig[0] else 0;
        const match_str_inner = if (ambig_ch != 0)
            try TermColors.colorChars(allocator, ambig_ch, "light_blue", match_colored)
        else
            try allocator.dupe(u8, match_colored);
        defer allocator.free(match_str_inner);

        // Build output with erosion decorations
        var out_buf: std.ArrayListUnmanaged(u8) = .{};
        defer out_buf.deinit(allocator);

        if (left_erosion_length > 0) {
            var le_buf: [16]u8 = undefined;
            const le_str = try std.fmt.bufPrint(&le_buf, "{d}", .{left_erosion_length});
            const pad = if (le_str.len < 2) @as(usize, 2) - le_str.len else @as(usize, 0);
            try out_buf.appendNTimes(allocator, ' ', pad);
            try out_buf.append(allocator, '.');
            try out_buf.appendSlice(allocator, le_str);
            try out_buf.append(allocator, '.');
        } else {
            try out_buf.appendSlice(allocator, "    ");
        }
        try out_buf.appendSlice(allocator, match_str_inner);
        if (right_erosion_length > 0) {
            var re_buf: [16]u8 = undefined;
            const re_str = try std.fmt.bufPrint(&re_buf, "{d}", .{right_erosion_length});
            const pad = if (re_str.len < 2) @as(usize, 2) - re_str.len else @as(usize, 0);
            try out_buf.appendNTimes(allocator, ' ', pad);
            try out_buf.append(allocator, '.');
            try out_buf.appendSlice(allocator, re_str);
            try out_buf.append(allocator, '.');
        } else {
            try out_buf.appendSlice(allocator, "    ");
        }

        // "                    " + match_str + "  " + extra_str + score(w12) + gene(w25)
        const score_str = try fmtSigFigs(allocator, score, 6);
        defer allocator.free(score_str);
        // Pad score to 12 chars right-justified
        const score_pad = if (score_str.len < 12) @as(usize, 12) - score_str.len else @as(usize, 0);
        // Pad gene to 25 chars right-justified
        const gene_pad = if (gene.len < 25) @as(usize, 25) - gene.len else @as(usize, 0);

        var line_buf: std.ArrayListUnmanaged(u8) = .{};
        defer line_buf.deinit(allocator);
        try line_buf.appendSlice(allocator, "                    ");
        try line_buf.appendSlice(allocator, out_buf.items);
        try line_buf.appendSlice(allocator, "  ");
        try line_buf.appendSlice(allocator, extra_str);
        try line_buf.appendNTimes(allocator, ' ', score_pad);
        try line_buf.appendSlice(allocator, score_str);
        try line_buf.appendNTimes(allocator, ' ', gene_pad);
        try line_buf.appendSlice(allocator, gene);
        try line_buf.append(allocator, '\n');
        try std.fs.File.stdout().writeAll(line_buf.items);

        // Print insertion indicator line if there were insertions
        if (left_insert.len + right_insert.len > 0) {
            var ins_buf: std.ArrayListUnmanaged(u8) = .{};
            defer ins_buf.deinit(allocator);
            try ins_buf.appendNTimes(allocator, 'i', left_insert.len);
            try ins_buf.appendNTimes(allocator, ' ', mid.len);
            try ins_buf.appendNTimes(allocator, 'i', right_insert.len);
            const ins_colored = try TermColors.color(allocator, "yellow", ins_buf.items);
            defer allocator.free(ins_colored);
            const ins_line = try std.fmt.allocPrint(allocator, "                        {s}  \n", .{ins_colored});
            defer allocator.free(ins_line);
            try std.fs.File.stdout().writeAll(ins_line);
        }
    }

    /// Run all genes for one kset.
    fn runKSet(
        self: *DPHandler,
        seqs: *Sequences,
        kset: KSet,
        only_genes: *std.AutoHashMapUnmanaged(Region, std.StringHashMapUnmanaged(void)),
        best_scores: *std.AutoHashMapUnmanaged(KSet, f64),
        total_scores: *std.AutoHashMapUnmanaged(KSet, f64),
        best_genes: *std.AutoHashMapUnmanaged(KSet, std.AutoHashMapUnmanaged(Region, []u8)),
    ) !void {
        const allocator = self.allocator;

        try best_scores.put(allocator, kset, mathutils.NEG_INF);
        try total_scores.put(allocator, kset, mathutils.NEG_INF);
        const bg_entry: std.AutoHashMapUnmanaged(Region, []u8) = .{};
        try best_genes.put(allocator, kset, bg_entry);

        // Debug: kset header
        if (self.args.debug == 2) {
            if (std.mem.eql(u8, self.algorithm, "forward")) {
                var hdr_buf: [256]u8 = undefined;
                const s = try std.fmt.bufPrint(&hdr_buf, "         {d: >3}{d: >3} {s: >6} {s: >9}  {s: >7}  {s: >7} {s}\n", .{ kset.v, kset.d, "prob", "logprob", "total", "origin", "---------------" });
                try std.fs.File.stdout().writeAll(s);
            } else {
                var hdr_buf: [256]u8 = undefined;
                const s = try std.fmt.bufPrint(&hdr_buf, "         {d: >3}{d: >3} {s}\n", .{ kset.v, kset.d, "---------------" });
                try std.fs.File.stdout().writeAll(s);
            }
        }

        var regional_best: std.AutoHashMapUnmanaged(Region, f64) = .{};
        defer regional_best.deinit(allocator);
        var regional_total: std.AutoHashMapUnmanaged(Region, f64) = .{};
        defer regional_total.deinit(allocator);
        var per_gene_support_this_kset: std.StringHashMapUnmanaged(f64) = .{};
        defer {
            var it = per_gene_support_this_kset.iterator();
            while (it.next()) |e| allocator.free(e.key_ptr.*);
            per_gene_support_this_kset.deinit(allocator);
        }

        for (bcr.germ_lines.regions) |region| {
            try regional_best.put(allocator, region, mathutils.NEG_INF);
            try regional_total.put(allocator, region, mathutils.NEG_INF);

            var query_strs = try self.getQueryStrs(seqs, kset, region);
            defer {
                for (query_strs.items) |s| allocator.free(s);
                query_strs.deinit(allocator);
            }

            var query_seqs = try self.getSubSeqs(seqs, kset, region);
            defer query_seqs.deinit(allocator);

            // Debug: region query display
            if (self.args.debug == 2) {
                const rs = regionStr(region);
                if (std.mem.eql(u8, self.algorithm, "viterbi")) {
                    // Get ambiguous char from args (matches C++ hmms_.track()->ambiguous_char())
                    const ambig = self.args.ambig_base;
                    const ambig_ch: u8 = if (ambig.len > 0) ambig[0] else 0;
                    if (query_strs.items.len > 0) {
                        const colored = if (ambig_ch != 0)
                            try TermColors.colorChars(allocator, ambig_ch, "light_blue", query_strs.items[0])
                        else
                            try allocator.dupe(u8, query_strs.items[0]);
                        defer allocator.free(colored);
                        const line = try std.fmt.allocPrint(allocator, "                {s} query {s}\n", .{ rs, colored });
                        defer allocator.free(line);
                        try std.fs.File.stdout().writeAll(line);
                    }
                    // Additional query strings colored as mutants relative to first
                    if (query_strs.items.len > 1) {
                        // Build const slice for colorMutants
                        const const_strs = try allocator.alloc([]const u8, query_strs.items.len);
                        defer allocator.free(const_strs);
                        for (query_strs.items, 0..) |s, i| const_strs[i] = s;
                        for (query_strs.items[1..]) |qs| {
                            const mutant = try TermColors.colorMutants(allocator, "purple", qs, "", const_strs, ambig);
                            defer allocator.free(mutant);
                            const colored2 = if (ambig_ch != 0)
                                try TermColors.colorChars(allocator, ambig_ch, "light_blue", mutant)
                            else
                                try allocator.dupe(u8, mutant);
                            defer allocator.free(colored2);
                            const line = try std.fmt.allocPrint(allocator, "                {s} query {s}\n", .{ rs, colored2 });
                            defer allocator.free(line);
                            try std.fs.File.stdout().writeAll(line);
                        }
                    }
                } else {
                    const line = try std.fmt.allocPrint(allocator, "              {s}\n", .{rs});
                    defer allocator.free(line);
                    try std.fs.File.stdout().writeAll(line);
                }
            }

            const gene_set = only_genes.get(region) orelse continue;
            // Collect and sort gene names to match C++ set<string> iteration order
            const sorted_genes = try allocator.alloc([]const u8, gene_set.count());
            defer allocator.free(sorted_genes);
            {
                var gi = gene_set.iterator();
                var idx: usize = 0;
                while (gi.next()) |e| : (idx += 1) sorted_genes[idx] = e.key_ptr.*;
            }
            std.mem.sort([]const u8, sorted_genes, {}, struct {
                fn lessThan(_: void, a: []const u8, b: []const u8) bool {
                    return std.mem.order(u8, a, b) == .lt;
                }
            }.lessThan);
            for (sorted_genes) |gene| {
                try self.initCache(gene);

                var origin: []const u8 = "";
                const partial_match = self.findPartialCacheMatch(region, gene, kset);
                if (!partial_match.isNull()) {
                    // Copy path and score from cached kset
                    if (self.paths.getPtr(gene)) |gp| {
                        if (gp.get(partial_match)) |cached_path| {
                            var path_copy = TracebackPath.init();
                            errdefer path_copy.deinit(allocator);
                            try path_copy.path.appendSlice(allocator, cached_path.path.items);
                            path_copy.setScore(cached_path.score);
                            path_copy.setModel(cached_path.hmm.?);
                            try gp.put(allocator, kset, path_copy);
                        }
                    }
                    if (self.scores.getPtr(gene)) |gs| {
                        const cached_score = gs.get(partial_match) orelse mathutils.NEG_INF;
                        try gs.put(allocator, kset, cached_score);
                    }
                    origin = "cached";
                } else {
                    origin = try self.fillTrellis(kset, query_seqs, query_strs.items, gene);
                }

                const gene_score = if (self.scores.getPtr(gene)) |gs| gs.get(kset) orelse mathutils.NEG_INF else mathutils.NEG_INF;

                // Debug: per-gene viterbi output
                if (self.args.debug == 2 and std.mem.eql(u8, self.algorithm, "viterbi")) {
                    // Build const slice for printPath
                    const const_strs = try allocator.alloc([]const u8, query_strs.items.len);
                    defer allocator.free(const_strs);
                    for (query_strs.items, 0..) |s, i| const_strs[i] = s;
                    try self.printPath(kset, const_strs, gene, gene_score, origin);
                }

                // Update regional totals
                if (regional_total.getPtr(region)) |rt| {
                    rt.* = mathutils.addInLogSpace(gene_score, rt.*);
                }

                // Debug: forward per-gene line
                if (self.args.debug == 2 and std.mem.eql(u8, self.algorithm, "forward")) {
                    const rt_val = regional_total.get(region) orelse mathutils.NEG_INF;
                    const colored_gene = try TermColors.colorGene(allocator, gene);
                    defer allocator.free(colored_gene);
                    const fwd_line = try std.fmt.allocPrint(allocator, "                {e: >6.0} {d: >9.2}  {d: >7.2}  {s}  {s}\n", .{ @as(f64, @exp(gene_score)), gene_score, rt_val, origin, colored_gene });
                    defer allocator.free(fwd_line);
                    try std.fs.File.stdout().writeAll(fwd_line);
                }

                // Update regional best + best_genes for this kset
                const reg_best = regional_best.get(region) orelse mathutils.NEG_INF;
                if (gene_score > reg_best) {
                    if (regional_best.getPtr(region)) |rb| rb.* = gene_score;
                    if (best_genes.getPtr(kset)) |bg_map| {
                        const gene_val = try allocator.dupe(u8, gene);
                        errdefer allocator.free(gene_val);
                        if (bg_map.fetchRemove(region)) |old| {
                            allocator.free(old.value);
                        }
                        try bg_map.put(allocator, region, gene_val);
                    }
                }

                // per_gene_support for this kset
                const pg_key = try allocator.dupe(u8, gene);
                errdefer allocator.free(pg_key);
                try per_gene_support_this_kset.put(allocator, pg_key, gene_score);
            }

            // If no gene found for this region, return early
            if (best_genes.getPtr(kset)) |bg_map| {
                if (!bg_map.contains(region)) {
                    if (self.args.debug == 2) {
                        const msg = try std.fmt.allocPrint(allocator, "                  found no gene for {s} so skip\n", .{regionStr(region)});
                        defer allocator.free(msg);
                        try std.fs.File.stdout().writeAll(msg);
                    }
                    return;
                }
            }
        }

        // Store combined best and total scores
        const rb_v = regional_best.get(.v) orelse mathutils.NEG_INF;
        const rb_d = regional_best.get(.d) orelse mathutils.NEG_INF;
        const rb_j = regional_best.get(.j) orelse mathutils.NEG_INF;
        const rt_v = regional_total.get(.v) orelse mathutils.NEG_INF;
        const rt_d = regional_total.get(.d) orelse mathutils.NEG_INF;
        const rt_j = regional_total.get(.j) orelse mathutils.NEG_INF;

        try best_scores.put(allocator, kset, mathutils.addWithMinusInfinities(rb_v, mathutils.addWithMinusInfinities(rb_d, rb_j)));
        try total_scores.put(allocator, kset, mathutils.addWithMinusInfinities(rt_v, mathutils.addWithMinusInfinities(rt_d, rt_j)));

        // Compute per_gene_support across ksets
        for (bcr.germ_lines.regions) |region| {
            const gene_set = only_genes.get(region) orelse continue;
            // Sort gene names to match C++ set<string> iteration order
            const sorted_support_genes = try allocator.alloc([]const u8, gene_set.count());
            defer allocator.free(sorted_support_genes);
            {
                var si = gene_set.iterator();
                var sidx: usize = 0;
                while (si.next()) |e| : (sidx += 1) sorted_support_genes[sidx] = e.key_ptr.*;
            }
            std.mem.sort([]const u8, sorted_support_genes, {}, struct {
                fn lessThan(_: void, a: []const u8, b: []const u8) bool {
                    return std.mem.order(u8, a, b) == .lt;
                }
            }.lessThan);
            for (sorted_support_genes) |gene| {
                var score_this_kset: f64 = 0.0;
                for (bcr.germ_lines.regions) |tmpreg| {
                    if (tmpreg == region) {
                        score_this_kset = mathutils.addWithMinusInfinities(score_this_kset, per_gene_support_this_kset.get(gene) orelse mathutils.NEG_INF);
                    } else {
                        score_this_kset = mathutils.addWithMinusInfinities(score_this_kset, regional_best.get(tmpreg) orelse mathutils.NEG_INF);
                    }
                }

                const existing = self.per_gene_support.get(gene) orelse mathutils.NEG_INF;
                if (score_this_kset > existing) {
                    if (self.per_gene_support.contains(gene)) {
                        self.per_gene_support.getPtr(gene).?.* = score_this_kset;
                    } else {
                        const pg_key = try self.allocator.dupe(u8, gene);
                        errdefer self.allocator.free(pg_key);
                        try self.per_gene_support.put(self.allocator, pg_key, score_this_kset);
                    }
                }
            }
        }
    }

    /// Build a RecoEvent for the given kset and best gene assignments.
    fn fillRecoEvent(
        self: *DPHandler,
        seqs: *Sequences,
        kset: KSet,
        best_genes_for_kset: *const std.AutoHashMapUnmanaged(Region, []u8),
        score: f64,
    ) !RecoEvent {
        const allocator = self.allocator;
        var event = try RecoEvent.init(allocator);
        errdefer event.deinit(allocator);

        for (bcr.germ_lines.regions) |region| {
            const rs = regionStr(region);
            const gene = best_genes_for_kset.get(region) orelse {
                event.setScore(mathutils.NEG_INF);
                return event;
            };

            var query_strs = try self.getQueryStrs(seqs, kset, region);
            defer {
                for (query_strs.items) |s| allocator.free(s);
                query_strs.deinit(allocator);
            }
            if (query_strs.items.len == 0) {
                event.setScore(mathutils.NEG_INF);
                return event;
            }

            // Get traceback path names
            const path_names_list = try self.getPathNames(gene, kset, allocator);
            defer {
                for (path_names_list.items) |s| allocator.free(s);
                var pl = path_names_list;
                pl.deinit(allocator);
            }

            if (path_names_list.items.len == 0) {
                event.setScore(mathutils.NEG_INF);
                return event;
            }

            try event.setGene(allocator, rs, gene);
            const del_3p = try std.fmt.allocPrint(allocator, "{s}_3p", .{rs});
            defer allocator.free(del_3p);
            const del_5p = try std.fmt.allocPrint(allocator, "{s}_5p", .{rs});
            defer allocator.free(del_5p);
            try event.setDeletion(allocator, del_3p, self.getErosionLength("right", path_names_list.items, gene));
            try event.setDeletion(allocator, del_5p, self.getErosionLength("left", path_names_list.items, gene));
            try self.setInsertions(region, path_names_list.items, &event);
        }

        event.setScore(score);
        try event.setNaiveSeq(allocator, self.gl);
        return event;
    }

    /// Get the stored path name vector for a gene+kset combination.
    fn getPathNames(
        self: *DPHandler,
        gene: []const u8,
        kset: KSet,
        allocator: std.mem.Allocator,
    ) !std.ArrayListUnmanaged([]u8) {
        const gene_paths = self.paths.getPtr(gene) orelse return .{};
        const path = gene_paths.getPtr(kset) orelse return .{};
        return path.nameVector(allocator);
    }

    /// Set insertions on `event` based on path state names for `region`.
    /// Corresponds to C++ `DPHandler::SetInsertions`.
    fn setInsertions(self: *DPHandler, region: Region, path_names: []const []const u8, event: *RecoEvent) !void {
        const allocator = self.allocator;
        const ins = Insertions{};
        var buf: [200]u8 = undefined;
        for (ins.forRegion(region)) |insertion| {
            const side: []const u8 = if (std.mem.eql(u8, insertion, "jf")) "right" else "left";
            const inserted = getInsertion(side, path_names, &buf);
            try event.setInsertion(allocator, insertion, inserted);
        }
    }

    /// Extract inserted bases from path state names for the given side.
    /// Writes into caller-provided `buf` and returns a slice into it.
    /// Corresponds to C++ `DPHandler::GetInsertion`.
    fn getInsertion(side: []const u8, names: []const []const u8, buf: *[200]u8) []const u8 {
        if (std.mem.eql(u8, side, "left")) {
            // Scan left-to-right: collect last char of each "insert..." state
            const max = @min(names.len, 200);
            var len: usize = 0;
            for (names[0..max]) |name| {
                if (!std.mem.startsWith(u8, name, "insert")) break;
                buf[len] = name[name.len - 1];
                len += 1;
            }
            return buf[0..len];
        } else {
            // right: scan right-to-left
            var len: usize = 0;
            var i: usize = names.len;
            while (i > 0) {
                i -= 1;
                if (!std.mem.startsWith(u8, names[i], "insert")) break;
                if (len < 200) {
                    // prepend
                    std.mem.copyBackwards(u8, buf[1 .. len + 1], buf[0..len]);
                    buf[0] = names[i][names[i].len - 1];
                    len += 1;
                }
            }
            return buf[0..len];
        }
    }

    /// Get the number of eroded bases on the given side.
    /// Corresponds to C++ `DPHandler::GetErosionLength`.
    fn getErosionLength(self: *DPHandler, side: []const u8, names: []const []const u8, gene_name: []const u8) usize {
        const germline = self.gl.seqs.get(gene_name) orelse return 0;

        // Check if all states are inserts
        var all_inserts = true;
        for (names) |name| {
            if (!std.mem.startsWith(u8, name, "insert")) {
                all_inserts = false;
                break;
            }
        }
        if (all_inserts) {
            if (std.mem.eql(u8, side, "left"))
                return germline.len / 2
            else
                return (germline.len + 1) / 2;
        }

        // Find the relevant boundary state
        if (std.mem.eql(u8, side, "left")) {
            // leftmost non-insert state
            var istate: usize = 0;
            for (names, 0..) |name, il| {
                if (!std.mem.startsWith(u8, name, "insert")) {
                    istate = il;
                    break;
                }
            }
            const state_idx = parseStateIndex(names[istate]) orelse return 0;
            return state_idx;
        } else {
            // rightmost non-insert state
            var istate: usize = 0;
            var found = false;
            var il: usize = names.len;
            while (il > 0) {
                il -= 1;
                if (!std.mem.startsWith(u8, names[il], "insert")) {
                    istate = il;
                    found = true;
                    break;
                }
            }
            if (!found) return 0;
            const state_idx = parseStateIndex(names[istate]) orelse return 0;
            if (germline.len == 0) return 0;
            return germline.len - state_idx - 1;
        }
    }

    /// Handle fishy multi-seq annotations by re-running with naive seq.
    /// Corresponds to C++ `DPHandler::HandleFishyAnnotations`.
    pub fn handleFishyAnnotations(
        self: *DPHandler,
        multi_seq_result: *Result,
        qry_seqs: []const Sequence,
        kbounds: KBounds,
        only_gene_list: []const []const u8,
        overall_mute_freq: f64,
    ) !void {
        const allocator = self.allocator;
        const best_ev = multi_seq_result.best_event orelse return;

        // Create naive sequence (use first query's track)
        const first_track = qry_seqs[0].track orelse return;
        var naive_seq = try Sequence.initFromString(allocator, first_track, "naive-seq", best_ev.naive_seq);
        defer naive_seq.deinit(allocator);

        const naive_seqs = [_]Sequence{naive_seq};
        var naive_result = try self.run(&naive_seqs, kbounds, only_gene_list, overall_mute_freq, true);
        defer naive_result.deinit();

        const naive_ev = naive_result.best_event orelse return;
        if (multi_seq_result.best_event) |*me| {
            for (bcr.germ_lines.regions) |region| {
                const rs = regionStr(region);
                const ng = naive_ev.genes.get(rs) orelse continue;
                try me.setGene(allocator, rs, ng);
            }
            const all_dels = [_][]const u8{ "v_5p", "v_3p", "d_5p", "d_3p", "j_5p", "j_3p" };
            for (all_dels) |delname| {
                const nd = naive_ev.deletions.get(delname) orelse 0;
                try me.setDeletion(allocator, delname, nd);
            }
            const all_ins = [_][]const u8{ "fv", "vd", "dj", "jf" };
            for (all_ins) |ins_name| {
                const ni = naive_ev.insertions.get(ins_name) orelse "";
                try me.setInsertion(allocator, ins_name, ni);
            }
        }
    }
};

/// Parse the state index from a state name like "IGHV3-15*01_42".
/// Returns the integer after the last underscore, or null.
fn parseStateIndex(name: []const u8) ?usize {
    const last_under = std.mem.lastIndexOfScalar(u8, name, '_') orelse return null;
    const idx_str = name[last_under + 1 ..];
    return std.fmt.parseInt(usize, idx_str, 10) catch null;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

/// Format a float with N significant digits, matching C++ `cout << setprecision(N) << value`
/// (i.e. %g format with rounding).  Caller owns the returned slice.
fn fmtSigFigs(allocator: std.mem.Allocator, value: f64, sig_figs: u8) ![]u8 {
    const abs_val = @abs(value);
    // Determine number of integer digits (digits before decimal point)
    const int_digits: u8 = if (abs_val < 1.0) 0 else blk: {
        var d: u8 = 0;
        var v = abs_val;
        while (v >= 1.0) : (d += 1) v /= 10.0;
        break :blk d;
    };
    // Number of decimal places needed for sig_figs significant digits
    const dec_places: u8 = if (sig_figs > int_digits) sig_figs - int_digits else 0;
    // Round to the right number of decimal places
    const factor = std.math.pow(f64, 10.0, @as(f64, @floatFromInt(dec_places)));
    const rounded = @round(value * factor) / factor;
    // Format with exact decimal places
    const raw = try std.fmt.allocPrint(allocator, "{d}", .{rounded});
    defer allocator.free(raw);
    // Find decimal point
    const dot_pos = std.mem.indexOfScalar(u8, raw, '.') orelse return try allocator.dupe(u8, raw);
    // Compute how many chars to keep after decimal point
    const end = @min(dot_pos + 1 + @as(usize, dec_places), raw.len);
    // Trim trailing zeros after decimal point (matching %g behavior)
    var trim_end = end;
    while (trim_end > dot_pos + 1 and raw[trim_end - 1] == '0') trim_end -= 1;
    if (trim_end > 0 and raw[trim_end - 1] == '.') trim_end -= 1;
    return try allocator.dupe(u8, raw[0..trim_end]);
}

test "fmtSigFigs matches C++ cout default" {
    const a = std.testing.allocator;
    // -103.3918 rounded to 6 sig figs → -103.392
    const s1 = try fmtSigFigs(a, -103.3918, 6);
    defer a.free(s1);
    try std.testing.expectEqualStrings("-103.392", s1);
    // -72.9176 = 6 sig figs (exact)
    const s2 = try fmtSigFigs(a, -72.9176, 6);
    defer a.free(s2);
    try std.testing.expectEqualStrings("-72.9176", s2);
    // -10.4 = trailing zeros trimmed
    const s3 = try fmtSigFigs(a, -10.4, 6);
    defer a.free(s3);
    try std.testing.expectEqualStrings("-10.4", s3);
}

test "DPHandler: parseStateIndex" {
    try std.testing.expectEqual(@as(?usize, 42), parseStateIndex("IGHV3-15_star_01_42"));
    try std.testing.expectEqual(@as(?usize, 0), parseStateIndex("IGHV3-15_star_01_0"));
    try std.testing.expectEqual(@as(?usize, null), parseStateIndex("insert_left_A"));
}
