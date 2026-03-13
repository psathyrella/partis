/// ham/bcr_utils/result.zig — Zig port of ham::Result
///
/// Holds the collection of RecoEvents for all ksets, finalizes to best.
/// C++ source: packages/ham/src/bcrutils.cc

const std = @import("std");
const SupportPair = @import("support_pair.zig").SupportPair;
const RecoEvent = @import("reco_event.zig").RecoEvent;
const GermLines = @import("germ_lines.zig").GermLines;
const k_bounds_mod = @import("k_bounds.zig");
const KSet = k_bounds_mod.KSet;
const KBounds = k_bounds_mod.KBounds;
const germ_lines_mod = @import("germ_lines.zig");
const hasDGene = germ_lines_mod.hasDGene;
const regions = germ_lines_mod.regions;

/// Corresponds to C++ `ham::Result`.
pub const Result = struct {
    total_score: f64,
    no_path: bool,
    locus: []u8,
    better_kbounds: KBounds,
    boundary_error: bool,
    could_not_expand: bool,
    finalized: bool,
    events: std.ArrayListUnmanaged(RecoEvent),
    best_event: ?RecoEvent,

    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator, kbounds: KBounds, locus: []const u8) !Result {
        return Result{
            .total_score = -std.math.inf(f64),
            .no_path = false,
            .locus = try allocator.dupe(u8, locus),
            .better_kbounds = kbounds,
            .boundary_error = false,
            .could_not_expand = false,
            .finalized = false,
            .events = .{},
            .best_event = null,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Result) void {
        const allocator = self.allocator;
        allocator.free(self.locus);
        for (self.events.items) |*ev| ev.deinit(allocator);
        self.events.deinit(allocator);
        if (self.best_event) |*ev| ev.deinit(allocator);
    }

    pub fn pushBackRecoEvent(self: *Result, event: RecoEvent) !void {
        try self.events.append(self.allocator, event);
    }

    /// Sort events by score, set best_event, set per_gene_support on best.
    /// Corresponds to C++ `Result::Finalize`.
    pub fn finalize(
        self: *Result,
        _gl: *const GermLines,
        unsorted_per_gene_support: *const std.StringHashMapUnmanaged(f64),
        best_kset: KSet,
        kbounds: KBounds,
    ) !void {
        _ = _gl;
        std.debug.assert(!self.finalized);

        // Sort events descending by score
        std.mem.sort(RecoEvent, self.events.items, {}, struct {
            fn desc(_: void, a: RecoEvent, b: RecoEvent) bool {
                return a.score > b.score;
            }
        }.desc);

        // Take ownership of the best event by removing it from the events list.
        // This avoids a double-free: deinit frees both best_event and every item
        // in events, so the best event must appear in exactly one of the two.
        self.best_event = self.events.orderedRemove(0);

        // Build per_gene_support for best event
        for (regions) |region| {
            var support: std.ArrayListUnmanaged(SupportPair) = .{};
            var it = unsorted_per_gene_support.iterator();
            while (it.next()) |entry| {
                const gene = entry.key_ptr.*;
                const logprob = entry.value_ptr.*;
                const r = GermLines.getRegion(gene) catch continue;
                if (r != region[0]) continue;
                try support.append(self.allocator, SupportPair{ .gene = gene, .logprob = logprob });
            }
            // Sort descending by logprob
            std.mem.sort(SupportPair, support.items, {}, struct {
                fn desc(_: void, a: SupportPair, b: SupportPair) bool {
                    return a.logprob > b.logprob;
                }
            }.desc);
            const key = try self.allocator.dupe(u8, region);
            errdefer self.allocator.free(key);
            if (self.best_event) |*be| {
                try be.per_gene_support.put(self.allocator, key, support);
            } else {
                support.deinit(self.allocator);
                self.allocator.free(key);
            }
        }

        self.checkBoundaries(best_kset, kbounds);
        self.finalized = true;
    }

    fn checkBoundaries(self: *Result, best: KSet, kbounds: KBounds) void {
        const delta: usize = 2;
        if (best.v == kbounds.vmin) {
            self.boundary_error = true;
            self.better_kbounds.vmin = if (kbounds.vmin > delta) kbounds.vmin - delta else 1;
        }
        if (best.v == kbounds.vmax -| 1) {
            self.boundary_error = true;
            self.better_kbounds.vmax = kbounds.vmax + delta;
        }
        if (hasDGene(self.locus) and best.d == kbounds.dmin) {
            self.boundary_error = true;
            self.better_kbounds.dmin = if (kbounds.dmin > delta) kbounds.dmin - delta else 1;
        }
        if (hasDGene(self.locus) and best.d == kbounds.dmax -| 1) {
            self.boundary_error = true;
            self.better_kbounds.dmax = kbounds.dmax + delta;
        }
        if (self.boundary_error and self.better_kbounds.eql(kbounds))
            self.could_not_expand = true;
    }
};
