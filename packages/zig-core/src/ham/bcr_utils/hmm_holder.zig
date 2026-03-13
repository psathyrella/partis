/// ham/bcr_utils/hmm_holder.zig — Zig port of ham::HMMHolder
///
/// Lazy-loading cache of Model objects keyed by gene name.
/// C++ source: packages/ham/src/bcrutils.cc

const std = @import("std");
const Model = @import("../model.zig").Model;
const germ_lines_mod = @import("germ_lines.zig");
const GermLines = germ_lines_mod.GermLines;
const regions = germ_lines_mod.regions;

/// Corresponds to C++ `ham::HMMHolder`.
pub const HMMHolder = struct {
    hmm_dir: []u8,
    gl: *GermLines,
    /// gene name → owned Model pointer
    hmms: std.StringHashMapUnmanaged(*Model),
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator, hmm_dir: []const u8, gl: *GermLines) !HMMHolder {
        return HMMHolder{
            .hmm_dir = try allocator.dupe(u8, hmm_dir),
            .gl = gl,
            .hmms = .{},
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *HMMHolder) void {
        const allocator = self.allocator;
        allocator.free(self.hmm_dir);
        var it = self.hmms.iterator();
        while (it.next()) |entry| {
            allocator.free(entry.key_ptr.*);
            entry.value_ptr.*.deinit(allocator);
            allocator.destroy(entry.value_ptr.*);
        }
        self.hmms.deinit(allocator);
    }

    /// Get (or lazy-load) a Model for `gene`.
    /// Corresponds to C++ `HMMHolder::Get`.
    pub fn get(self: *HMMHolder, gene: []const u8) !*Model {
        if (self.hmms.get(gene)) |m| return m;

        // Not yet cached — load from disk
        const sanitized = try GermLines.sanitizeName(self.allocator, gene);
        defer self.allocator.free(sanitized);
        const infname = try std.fmt.allocPrint(self.allocator, "{s}/{s}.yaml", .{ self.hmm_dir, sanitized });
        defer self.allocator.free(infname);

        const m = try self.allocator.create(Model);
        errdefer self.allocator.destroy(m);
        m.* = try Model.init(self.allocator);
        errdefer m.deinit(self.allocator);
        try m.parse(self.allocator, infname);

        const key = try self.allocator.dupe(u8, gene);
        errdefer self.allocator.free(key);
        try self.hmms.put(self.allocator, key, m);
        return m;
    }

    /// Read all available HMMs into memory.
    /// Corresponds to C++ `HMMHolder::CacheAll`.
    pub fn cacheAll(self: *HMMHolder) !void {
        for (regions) |region| {
            const names = self.gl.names.get(region) orelse continue;
            for (names.items) |gene| {
                const sanitized = try GermLines.sanitizeName(self.allocator, gene);
                defer self.allocator.free(sanitized);
                const infname = try std.fmt.allocPrint(self.allocator, "{s}/{s}.yaml", .{ self.hmm_dir, sanitized });
                defer self.allocator.free(infname);
                // Only cache if file exists
                std.fs.cwd().access(infname, .{}) catch continue;
                _ = try self.get(gene);
            }
        }
    }

    /// Rescale overall_mute_freq for all genes in `only_genes`.
    pub fn rescaleOverallMuteFreqs(
        self: *HMMHolder,
        only_genes: *std.StringHashMapUnmanaged(std.StringHashMapUnmanaged(void)),
        overall_mute_freq: f64,
    ) !void {
        for (regions) |region| {
            const gene_set = only_genes.get(region) orelse continue;
            var git = gene_set.iterator();
            while (git.next()) |entry| {
                const m = try self.get(entry.key_ptr.*);
                try m.rescaleOverallMuteFreq(self.allocator, overall_mute_freq);
            }
        }
    }

    /// Undo rescaling.
    pub fn unRescaleOverallMuteFreqs(
        self: *HMMHolder,
        only_genes: *std.StringHashMapUnmanaged(std.StringHashMapUnmanaged(void)),
    ) !void {
        for (regions) |region| {
            const gene_set = only_genes.get(region) orelse continue;
            var git = gene_set.iterator();
            while (git.next()) |entry| {
                const m = try self.get(entry.key_ptr.*);
                m.unRescaleOverallMuteFreq(self.allocator);
            }
        }
    }
};
