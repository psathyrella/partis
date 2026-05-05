/// ham/trellis.zig — Zig port of ham/src/trellis.cc + ham/include/trellis.h
///
/// Forward/Viterbi DP trellis for a single HMM run over a Sequences object.
///
/// C++ source: packages/ham/src/trellis.cc, packages/ham/include/trellis.h
/// C++ author: psathyrella/ham

const std = @import("std");
const Model = @import("model.zig").Model;
const State = @import("state.zig").State;
const state_mod = @import("state.zig");
const Sequences = @import("sequences.zig").Sequences;
const Sequence = @import("sequences.zig").Sequence;
const TracebackPath = @import("traceback_path.zig").TracebackPath;
const mathutils = @import("mathutils.zig");

/// DP trellis for Forward/Viterbi algorithms.
/// Corresponds to C++ `ham::Trellis`.
pub const Trellis = struct {
    /// HMM model (not owned).
    hmm: *const Model,
    /// Sequences to run DP over (BORROWED; not owned).
    /// The pointee must outlive this Trellis. For trellises stored in
    /// `DPHandler.scratch_cachefo`, the owning Sequences lives in the
    /// enclosing TrellisEntry and is freed alongside the trellis.
    /// For temporary chunk-borrowing trellises, the pointee is the caller's
    /// stack-local `query_seqs` and lives only for the duration of fillTrellis.
    seqs: *const Sequences,
    /// Cached at init from `seqs.sequence_length`. Read by `traceback`,
    /// `viterbi`, and `forward`; valid even if the borrowed seqs is no longer
    /// being read post-DP.
    seq_len: usize,

    /// Flat traceback table indexed as [position * traceback_n_states + state].
    /// Each entry is the previous state index (-1 if no valid predecessor).
    /// Corresponds to C++ `typedef vector<vector<int16_t>> int_2D`, but stored
    /// as a single contiguous allocation for locality and one alloc/free per
    /// Trellis instead of seq_len + 1.
    /// Empty (`&.{}`) until `viterbi()` allocates it; `forward()` never touches
    /// this field. Read via `pickTracebackView` (returns a `TracebackView` with
    /// the slice and stride); written via `tracebackSet`.
    traceback_table: []i16,
    /// State width of `traceback_table`. Captured at allocation time and equal
    /// to `hmm.nStates()` for the trellis's lifetime. Stored explicitly so a
    /// chunk-borrowing temp trellis can verify its `cached_trellis` was built
    /// against the same model (see `initWithSeqs` assert).
    traceback_n_states: usize,

    /// Pointer to another trellis with pre-filled DP tables (not owned).
    cached_trellis: ?*const Trellis,

    ending_viterbi_pointer: i16,
    ending_viterbi_log_prob: f64,
    ending_forward_log_prob: f64,

    /// Per-position Viterbi log-probs (including ending transition).
    viterbi_log_probs: std.ArrayListUnmanaged(f64),
    /// Per-position Forward log-probs.
    forward_log_probs: std.ArrayListUnmanaged(f64),
    /// Per-position best-state index (for Viterbi).
    viterbi_indices: std.ArrayListUnmanaged(i32),

    /// Current and previous scoring columns (size = n_states each).
    scoring_current: []f64,
    scoring_previous: []f64,

    /// Scratch bitset for deduplicating next_states additions.
    next_seen: state_mod.StateBitset,

    /// Captured at `initWithSeqs` time. For trellises stored in
    /// `DPHandler.scratch_cachefo`, this is `dph.scratchAllocator()`,
    /// whose `Allocator.ptr` aliases `&dph.query_arena`. The DPHandler
    /// must therefore not be moved after a Trellis is constructed against
    /// it, or every stored `allocator.ptr` would dangle. All current
    /// callers bind `dph` to a stack-local `var` and never reassign,
    /// which preserves the invariant. (Item 4, #342.)
    allocator: std.mem.Allocator,

    /// Create a Trellis borrowing the given Sequences.
    /// The caller retains ownership; `seqs` must outlive this Trellis.
    /// Corresponds to C++ `Trellis(Model*, Sequences, Trellis*)`, but without
    /// the C++ copy-by-value clone — the borrow is sufficient for DP read
    /// access, and for cached/stored trellises only the result arrays are
    /// re-read.
    pub fn initWithSeqs(allocator: std.mem.Allocator, hmm: *const Model, seqs: *const Sequences, cached: ?*const Trellis) !Trellis {
        const n = hmm.nStates();
        // Stale-cache guard: a chunk-borrowing trellis must be built against
        // the same model as its cached source, otherwise its scalar reads from
        // `cached.traceback_table` (via `pickTracebackView`) would index with
        // the wrong stride. The cached trellis can have an empty traceback
        // table (forward-only), so guard only when it has one.
        if (cached) |ct| {
            if (ct.traceback_table.len > 0) {
                std.debug.assert(ct.traceback_n_states == n);
            }
        }

        const scoring_current = try allocator.alloc(f64, n);
        errdefer allocator.free(scoring_current);
        const scoring_previous = try allocator.alloc(f64, n);
        errdefer allocator.free(scoring_previous);
        @memset(scoring_current, mathutils.NEG_INF);
        @memset(scoring_previous, mathutils.NEG_INF);

        return Trellis{
            .hmm = hmm,
            .seqs = seqs,
            .seq_len = seqs.sequence_length,
            .traceback_table = &.{},
            .traceback_n_states = 0,
            .cached_trellis = cached,
            .ending_viterbi_pointer = -1,
            .ending_viterbi_log_prob = mathutils.NEG_INF,
            .ending_forward_log_prob = mathutils.NEG_INF,
            .viterbi_log_probs = .{},
            .forward_log_probs = .{},
            .viterbi_indices = .{},
            .scoring_current = scoring_current,
            .scoring_previous = scoring_previous,
            .next_seen = state_mod.StateBitset.initEmpty(),
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Trellis) void {
        const allocator = self.allocator;
        // seqs is borrowed; not freed here.
        if (self.traceback_table.len > 0) allocator.free(self.traceback_table);
        self.viterbi_log_probs.deinit(allocator);
        self.forward_log_probs.deinit(allocator);
        self.viterbi_indices.deinit(allocator);
        allocator.free(self.scoring_current);
        allocator.free(self.scoring_previous);
    }

    /// Run the Viterbi algorithm.
    /// After calling this, `ending_viterbi_log_prob` and `traceback_table` are set.
    /// Corresponds to C++ `Trellis::Viterbi()`.
    pub fn viterbi(self: *Trellis) !void {
        const allocator = self.allocator;
        const seq_len = self.seq_len;
        const n_states = self.hmm.nStates();
        // Cache the borrowed seqs pointer in a local so the inner DP loop reads
        // it from a register rather than re-loading `self.seqs` (a heap-pointer
        // member after #357) on every emission lookup.
        const seqs = self.seqs;

        if (self.cached_trellis) |ct| {
            // Poach scalar values from cached trellis; array data accessed via accessors.
            self.ending_viterbi_pointer = @intCast(ct.viterbiIndicesAt(seq_len));
            self.ending_viterbi_log_prob = ct.endingViterbiLogProbAt(seq_len);
            return;
        }

        // Initialize storage
        try self.viterbi_log_probs.resize(allocator, seq_len);
        try self.viterbi_indices.resize(allocator, seq_len);
        @memset(self.viterbi_log_probs.items, mathutils.NEG_INF);
        @memset(self.viterbi_indices.items, -1);

        // Initialize flat traceback table [seq_len * n_states] = -1.
        // Single allocation replaces seq_len per-row allocations (item 3 of #342).
        // Note the explicit `traceback_table = &.{}` between free and alloc:
        // if the alloc fails, deinit reads `traceback_table` and would otherwise
        // see the just-freed slice with len > 0 and double-free it.
        std.debug.assert(seq_len <= std.math.maxInt(usize) / @max(n_states, 1));
        const total = seq_len * n_states;
        if (self.traceback_table.len != total) {
            if (self.traceback_table.len > 0) allocator.free(self.traceback_table);
            self.traceback_table = &.{};
            self.traceback_table = try allocator.alloc(i16, total);
        }
        @memset(self.traceback_table, -1);
        self.traceback_n_states = n_states;

        // Reset scoring columns
        @memset(self.scoring_current, mathutils.NEG_INF);
        @memset(self.scoring_previous, mathutils.NEG_INF);

        var current_states = std.ArrayListUnmanaged(usize){};
        defer current_states.deinit(allocator);
        var next_states = std.ArrayListUnmanaged(usize){};
        defer next_states.deinit(allocator);

        // Position 0: transitions from init state
        const init_st = self.hmm.initial orelse return error.NoInitState;
        // Reset next_seen before populating position 0
        self.next_seen = state_mod.StateBitset.initEmpty();
        for (init_st.to_state_indices.items) |i_st| {
            const emission_val = self.hmm.stateByIndex(i_st).emissionLogprobSeqs(seqs, 0);
            const dpval = emission_val + init_st.transitionLogprob(i_st);
            if (std.math.isNegativeInf(dpval)) continue;
            self.scoring_current[i_st] = dpval;
            self.cacheViterbiVals(0, dpval, i_st);
            // Mark outbound transitions for next position (deduplicated)
            for (self.hmm.stateByIndex(i_st).to_state_indices.items) |j| {
                if (!self.next_seen.isSet(j)) {
                    self.next_seen.set(j);
                    try next_states.append(allocator, j);
                }
            }
        }

        // Positions 1..seq_len-1
        var position: usize = 1;
        while (position < seq_len) : (position += 1) {
            try self.swapColumnsActive(allocator, &current_states, &next_states);
            try self.middleViterbiVals(allocator, current_states, &next_states, position);
        }
        try self.swapColumnsActive(allocator, &current_states, &next_states);

        // Compute ending probability
        self.ending_viterbi_pointer = -1;
        self.ending_viterbi_log_prob = mathutils.NEG_INF;
        for (0..n_states) |st_prev| {
            if (std.math.isNegativeInf(self.scoring_previous[st_prev])) continue;
            const dpval = self.scoring_previous[st_prev] + self.hmm.stateByIndex(st_prev).endTransitionLogprob();
            if (dpval > self.ending_viterbi_log_prob) {
                self.ending_viterbi_log_prob = dpval;
                self.ending_viterbi_pointer = @intCast(st_prev);
            }
        }
    }

    /// Run the Forward algorithm.
    /// After calling this, `ending_forward_log_prob` is set.
    /// Corresponds to C++ `Trellis::Forward()`.
    pub fn forward(self: *Trellis) !void {
        const allocator = self.allocator;
        const seq_len = self.seq_len;
        const n_states = self.hmm.nStates();
        // See viterbi(): cache borrowed seqs pointer in a local for the inner DP loop.
        const seqs = self.seqs;

        if (self.cached_trellis) |ct| {
            self.ending_forward_log_prob = ct.endingForwardLogProbAt(seq_len);
            return;
        }

        // Initialize storage
        try self.forward_log_probs.resize(allocator, seq_len);
        @memset(self.forward_log_probs.items, mathutils.NEG_INF);

        // Reset scoring columns
        @memset(self.scoring_current, mathutils.NEG_INF);
        @memset(self.scoring_previous, mathutils.NEG_INF);

        var current_states = std.ArrayListUnmanaged(usize){};
        defer current_states.deinit(allocator);
        var next_states = std.ArrayListUnmanaged(usize){};
        defer next_states.deinit(allocator);

        // Position 0: transitions from init state
        const init_st = self.hmm.initial orelse return error.NoInitState;
        // Reset next_seen before populating position 0
        self.next_seen = state_mod.StateBitset.initEmpty();
        for (init_st.to_state_indices.items) |i_st| {
            const emission_val = self.hmm.stateByIndex(i_st).emissionLogprobSeqs(seqs, 0);
            const dpval = emission_val + init_st.transitionLogprob(i_st);
            if (std.math.isNegativeInf(dpval)) continue;
            self.scoring_current[i_st] = dpval;
            // Mark outbound transitions for next position (deduplicated)
            for (self.hmm.stateByIndex(i_st).to_state_indices.items) |j| {
                if (!self.next_seen.isSet(j)) {
                    self.next_seen.set(j);
                    try next_states.append(allocator, j);
                }
            }
            self.cacheForwardVals(0, dpval, i_st);
        }

        // Positions 1..seq_len-1
        var position: usize = 1;
        while (position < seq_len) : (position += 1) {
            try self.swapColumnsActive(allocator, &current_states, &next_states);
            try self.middleForwardVals(allocator, current_states, &next_states, position);
        }
        try self.swapColumnsActive(allocator, &current_states, &next_states);

        // Compute ending probability
        self.ending_forward_log_prob = mathutils.NEG_INF;
        for (0..n_states) |st_prev| {
            if (std.math.isNegativeInf(self.scoring_previous[st_prev])) continue;
            const dpval = self.scoring_previous[st_prev] + self.hmm.stateByIndex(st_prev).endTransitionLogprob();
            if (std.math.isNegativeInf(dpval)) continue;
            self.ending_forward_log_prob = mathutils.addInLogSpace(self.ending_forward_log_prob, dpval);
        }
    }

    /// Traceback to produce a path.
    /// Corresponds to C++ `Trellis::Traceback(TracebackPath&)`.
    ///
    /// `path_alloc` backs the path's growing item list, which may have a
    /// different lifetime than the trellis itself: in the chunk-cache-hit
    /// path of `DPHandler.fillTrellis`, the temp trellis lives on the
    /// per-kset arena but the path is stored in `DPHandler.paths` (per
    /// query). Passing the allocator explicitly avoids the use-after-free
    /// that would otherwise come from `self.allocator` capturing a
    /// shorter-lived arena. (Item 4, #342.)
    pub fn traceback(self: *const Trellis, path_alloc: std.mem.Allocator, path: *TracebackPath) !void {
        std.debug.assert(self.seq_len != 0);
        path.setModel(self.hmm);
        if (std.math.isNegativeInf(self.ending_viterbi_log_prob)) return;
        path.setScore(self.ending_viterbi_log_prob);
        try path.pushBack(path_alloc, self.ending_viterbi_pointer);

        // Resolve the traceback source once: cached trellis if it has a
        // populated table, else self. If neither does, there's nothing to
        // walk back through.
        const view = self.pickTracebackView() orelse return;
        var pointer: i16 = self.ending_viterbi_pointer;
        var position: usize = self.seq_len - 1;
        while (position > 0) : (position -= 1) {
            pointer = view.tbl[position * view.stride + @as(usize, @intCast(pointer))];
            if (pointer == -1) {
                std.debug.print("No valid path at position {d}\n", .{position});
                return;
            }
            try path.pushBack(path_alloc, pointer);
        }
        std.debug.assert(path.size() > 0);
    }

    // ── Private helpers ───────────────────────────────────────────────────────

    fn swapColumnsActive(
        self: *Trellis,
        allocator: std.mem.Allocator,
        current_states: *std.ArrayListUnmanaged(usize),
        next_states: *std.ArrayListUnmanaged(usize),
    ) !void {
        // Swap scoring_current ↔ scoring_previous
        const tmp = self.scoring_previous;
        self.scoring_previous = self.scoring_current;
        self.scoring_current = tmp;
        @memset(self.scoring_current, mathutils.NEG_INF);

        // Clear next_seen for all indices that were in next_states
        // (they will become current_states; next_seen tracks what's queued for the NEW next)
        for (next_states.items) |j| self.next_seen.unset(j);

        // current_states ← next_states; next_states ← empty
        // Swap the backing allocations to avoid reallocation
        const old_current_items = current_states.items;
        const old_current_cap = current_states.capacity;
        current_states.items = next_states.items;
        current_states.capacity = next_states.capacity;
        next_states.items = old_current_items;
        next_states.capacity = old_current_cap;
        next_states.clearRetainingCapacity();
        // Sort current_states in ascending order to match C++ bitset iteration
        // (0..n_states). This is critical: cacheForwardVals accumulates
        // forward_log_probs[position] across all (state, prev_state) pairs,
        // and addInLogSpace is not associative, so different accumulation orders
        // produce different results. The cached forward_log_probs values are
        // then read by chunk-cached trellises via endingForwardLogProbAt().
        std.mem.sort(usize, current_states.items, {}, std.sort.asc(usize));
        _ = allocator; // backing storage reuse avoids new allocations
    }

    fn middleViterbiVals(
        self: *Trellis,
        allocator: std.mem.Allocator,
        current_states: std.ArrayListUnmanaged(usize),
        next_states: *std.ArrayListUnmanaged(usize),
        position: usize,
    ) !void {
        // See viterbi(): cache borrowed seqs pointer in a local for the inner loop.
        const seqs = self.seqs;
        for (current_states.items) |i_st_cur| {
            const st_cur = self.hmm.stateByIndex(i_st_cur);
            const emission_val = st_cur.emissionLogprobSeqs(seqs, position);
            if (std.math.isNegativeInf(emission_val)) continue;
            for (st_cur.from_state_indices.items) |i_st_prev| {
                const prev_val = self.scoring_previous[i_st_prev];
                if (std.math.isNegativeInf(prev_val)) continue;
                const dpval = prev_val + emission_val + self.hmm.stateByIndex(i_st_prev).transitionLogprob(i_st_cur);
                if (dpval > self.scoring_current[i_st_cur]) {
                    self.scoring_current[i_st_cur] = dpval;
                    self.tracebackSet(position, i_st_cur, @intCast(i_st_prev));
                }
                self.cacheViterbiVals(position, dpval, i_st_cur);
                // Mark outbound transitions inside the i_st_prev loop (as in C++) so we only
                // include states that are really needed, i.e. that have at least one valid predecessor
                for (st_cur.to_state_indices.items) |j| {
                    if (!self.next_seen.isSet(j)) {
                        self.next_seen.set(j);
                        try next_states.append(allocator, j);
                    }
                }
            }
        }
    }

    fn middleForwardVals(
        self: *Trellis,
        allocator: std.mem.Allocator,
        current_states: std.ArrayListUnmanaged(usize),
        next_states: *std.ArrayListUnmanaged(usize),
        position: usize,
    ) !void {
        // See viterbi(): cache borrowed seqs pointer in a local for the inner loop.
        const seqs = self.seqs;
        for (current_states.items) |i_st_cur| {
            const st_cur = self.hmm.stateByIndex(i_st_cur);
            const emission_val = st_cur.emissionLogprobSeqs(seqs, position);
            if (std.math.isNegativeInf(emission_val)) continue;
            for (st_cur.from_state_indices.items) |i_st_prev| {
                const prev_val = self.scoring_previous[i_st_prev];
                if (std.math.isNegativeInf(prev_val)) continue;
                const dpval = prev_val + emission_val + self.hmm.stateByIndex(i_st_prev).transitionLogprob(i_st_cur);
                self.scoring_current[i_st_cur] = mathutils.addInLogSpace(dpval, self.scoring_current[i_st_cur]);
                self.cacheForwardVals(position, dpval, i_st_cur);
                // Mark outbound transitions inside the i_st_prev loop (as in C++) so we only
                // include states that are really needed, i.e. that have at least one valid predecessor
                for (st_cur.to_state_indices.items) |j| {
                    if (!self.next_seen.isSet(j)) {
                        self.next_seen.set(j);
                        try next_states.append(allocator, j);
                    }
                }
            }
        }
    }

    inline fn cacheViterbiVals(self: *Trellis, position: usize, dpval: f64, i_st_cur: usize) void {
        const end_trans = self.hmm.stateByIndex(i_st_cur).endTransitionLogprob();
        const logprob = dpval + end_trans;
        if (logprob > self.viterbi_log_probs.items[position]) {
            self.viterbi_log_probs.items[position] = logprob;
            self.viterbi_indices.items[position] = @intCast(i_st_cur);
        }
    }

    inline fn cacheForwardVals(self: *Trellis, position: usize, dpval: f64, i_st_cur: usize) void {
        const end_trans = self.hmm.stateByIndex(i_st_cur).endTransitionLogprob();
        const logprob = dpval + end_trans;
        self.forward_log_probs.items[position] = mathutils.addInLogSpace(logprob, self.forward_log_probs.items[position]);
    }

    // ── Accessors: dispatch to cached trellis data or own data ──────────────

    /// Resolved view of the traceback source: the cached trellis's table if it
    /// has one populated, else self's, else null. The stride is captured at
    /// the same time as the slice so callers don't have to re-do the dispatch.
    const TracebackView = struct { tbl: []const i16, stride: usize };

    inline fn pickTracebackView(self: *const Trellis) ?TracebackView {
        if (self.cached_trellis) |ct| {
            if (ct.traceback_table.len > 0) {
                return .{ .tbl = ct.traceback_table, .stride = ct.traceback_n_states };
            }
        }
        if (self.traceback_table.len > 0) {
            return .{ .tbl = self.traceback_table, .stride = self.traceback_n_states };
        }
        return null;
    }

    /// Scalar write of *self*'s traceback table at (position, state). Used by
    /// `middleViterbiVals` while filling the table fresh.
    inline fn tracebackSet(self: *Trellis, position: usize, state: usize, v: i16) void {
        self.traceback_table[position * self.traceback_n_states + state] = v;
    }

    /// Ending Viterbi log-prob for a specific sequence length (used by cached trellis).
    pub fn endingViterbiLogProbAt(self: *const Trellis, length: usize) f64 {
        if (length == 0) return mathutils.NEG_INF;
        const items = if (self.cached_trellis) |ct| ct.viterbi_log_probs.items else self.viterbi_log_probs.items;
        if (length - 1 < items.len) return items[length - 1];
        return mathutils.NEG_INF;
    }

    /// Ending Forward log-prob for a specific sequence length.
    pub fn endingForwardLogProbAt(self: *const Trellis, length: usize) f64 {
        if (length == 0) return mathutils.NEG_INF;
        const items = if (self.cached_trellis) |ct| ct.forward_log_probs.items else self.forward_log_probs.items;
        if (length - 1 < items.len) return items[length - 1];
        return mathutils.NEG_INF;
    }

    /// Best state index at a specific sequence length (Viterbi).
    pub fn viterbiIndicesAt(self: *const Trellis, length: usize) i32 {
        if (length == 0) return -1;
        const items = if (self.cached_trellis) |ct| ct.viterbi_indices.items else self.viterbi_indices.items;
        if (length - 1 < items.len) return items[length - 1];
        return -1;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

// Integration test: parse a real model and run Forward on a short sequence.
test "Trellis: forward on real IGHV3-15 model" {
    const allocator = std.testing.allocator;

    const yaml_path = "/fh/fast/matsen_e/shared/partis-zig/partis/test/ref-results-slow/test/parameters/simu/sw/hmms/IGHV3-15_star_07.yaml";
    var model = try @import("model.zig").Model.init(allocator);
    defer model.deinit(allocator);
    model.parse(allocator, yaml_path) catch |err| {
        // Skip test if YAML file not available (e.g. CI without data)
        std.debug.print("skipping trellis test: {}\n", .{err});
        return;
    };

    // Build a 4-symbol sequence using the model's track
    const trk = model.track orelse return;
    var seq = try Sequence.initFromString(allocator, trk, "test_seq", "ACGT");
    defer seq.deinit(allocator);

    // Trellis borrows the Sequences; the caller retains ownership.
    var seqs = Sequences.init();
    defer seqs.deinit(allocator);
    try seqs.addSeq(allocator, try seq.clone(allocator));

    var trellis = try Trellis.initWithSeqs(allocator, &model, &seqs, null);
    defer trellis.deinit();

    try trellis.forward();

    // A 4-symbol query against a ~300nt V-gene HMM may legitimately return
    // -inf (no valid path); a finite log-prob is the success signal but not
    // a hard requirement. The structural success is that forward() ran without
    // crashing through the trellis allocation/teardown path.
    std.debug.print("forward log-prob: {d}\n", .{trellis.ending_forward_log_prob});
}

// Layout test: writing through `tracebackSet` and reading back via
// `pickTracebackView` must round-trip the value at the same (position, state).
// Catches row/col transposition in the indexing helpers without needing a full
// Model — the byte-equality gate at 5k/30k can't catch a transposition when
// seq_len ≈ n_states (square trellis), so this covers the gap.
test "Trellis: flat traceback layout round-trips position/state" {
    const allocator = std.testing.allocator;

    const seq_len: usize = 4;
    const n_states: usize = 5;

    // Hand-build just the indexing fields. The non-indexing fields aren't
    // touched by tracebackSet / pickTracebackView.
    var trellis: Trellis = undefined;
    trellis.allocator = allocator;
    trellis.cached_trellis = null;
    trellis.traceback_n_states = n_states;
    trellis.traceback_table = try allocator.alloc(i16, seq_len * n_states);
    defer allocator.free(trellis.traceback_table);
    @memset(trellis.traceback_table, -1);

    // Write a unique value at every (position, state).
    for (0..seq_len) |pos| {
        for (0..n_states) |st| {
            trellis.tracebackSet(pos, st, @intCast(pos * 10 + st));
        }
    }
    // Read back via the resolved view. A row/col swap in the indexing math
    // would surface as the wrong value (or out-of-range index) at every cell.
    const view = trellis.pickTracebackView() orelse return error.TestUnexpectedNull;
    for (0..seq_len) |pos| {
        for (0..n_states) |st| {
            const expected: i16 = @intCast(pos * 10 + st);
            try std.testing.expectEqual(expected, view.tbl[pos * view.stride + st]);
        }
    }

    // Empty-table fallback: pickTracebackView returns null when neither self
    // nor cached has a table, mirroring the original `?[][]i16` accessor's
    // null path. `traceback()` short-circuits on this.
    var empty_trellis: Trellis = undefined;
    empty_trellis.allocator = allocator;
    empty_trellis.cached_trellis = null;
    empty_trellis.traceback_n_states = 0;
    empty_trellis.traceback_table = &.{};
    try std.testing.expectEqual(@as(?Trellis.TracebackView, null), empty_trellis.pickTracebackView());

    // Cached-fallthrough: when self has no table but cached does, the view
    // should come from cached, with cached's stride (not self's, which is 0).
    var borrow_trellis: Trellis = undefined;
    borrow_trellis.allocator = allocator;
    borrow_trellis.cached_trellis = &trellis;
    borrow_trellis.traceback_n_states = 0;
    borrow_trellis.traceback_table = &.{};
    const borrow_view = borrow_trellis.pickTracebackView() orelse return error.TestUnexpectedNull;
    try std.testing.expectEqual(@as(usize, n_states), borrow_view.stride);
    for (0..seq_len) |pos| {
        for (0..n_states) |st| {
            const expected: i16 = @intCast(pos * 10 + st);
            try std.testing.expectEqual(expected, borrow_view.tbl[pos * borrow_view.stride + st]);
        }
    }
}

// Viterbi → traceback parity: run viterbi on a real model, walk the traceback
// path through the flat table, and assert the path is consistent with the
// per-position best-state record. This is end-to-end coverage of the indexing
// helpers under the actual DP loop, complementing the layout test above.
test "Trellis: viterbi traceback walks through flat table consistently" {
    const allocator = std.testing.allocator;

    const yaml_path = "/fh/fast/matsen_e/shared/partis-zig/partis/test/ref-results-slow/test/parameters/simu/sw/hmms/IGHV3-15_star_07.yaml";
    var model = try @import("model.zig").Model.init(allocator);
    defer model.deinit(allocator);
    model.parse(allocator, yaml_path) catch |err| {
        std.debug.print("skipping viterbi traceback test: {}\n", .{err});
        return;
    };

    const trk = model.track orelse return;
    var seq = try Sequence.initFromString(allocator, trk, "test_seq", "ACGTACGT");
    defer seq.deinit(allocator);

    var seqs = Sequences.init();
    defer seqs.deinit(allocator);
    try seqs.addSeq(allocator, try seq.clone(allocator));

    var trellis = try Trellis.initWithSeqs(allocator, &model, &seqs, null);
    defer trellis.deinit();

    try trellis.viterbi();

    // Traceback table should be populated under the new flat layout.
    try std.testing.expectEqual(seqs.sequence_length * model.nStates(), trellis.traceback_table.len);
    try std.testing.expectEqual(model.nStates(), trellis.traceback_n_states);

    // Skip the path comparison if viterbi found no valid path through the
    // model (depends on YAML data); the structural checks above are enough
    // to fail on a misshapen allocation.
    if (std.math.isNegativeInf(trellis.ending_viterbi_log_prob)) {
        std.debug.print("viterbi found no valid path; skipping traceback walk\n", .{});
        return;
    }

    var path = TracebackPath.initWithModel(&model);
    defer path.deinit(allocator);
    try trellis.traceback(allocator, &path);

    // Path length should equal seq_len: traceback pushes ending pointer, then
    // one entry per position from seq_len-1 down to 1.
    try std.testing.expectEqual(seqs.sequence_length, path.size());

    // Every traceback step must read a valid (>=0) state from the flat table.
    // This is what would fail catastrophically under a row/col transpose: the
    // wrong-strided read would land on -1 entries (or out of bounds — the
    // overflow assert above guards bounds).
    const view = trellis.pickTracebackView() orelse return error.TestUnexpectedNull;
    var pointer: i16 = trellis.ending_viterbi_pointer;
    var position: usize = trellis.seq_len - 1;
    while (position > 0) : (position -= 1) {
        const next = view.tbl[position * view.stride + @as(usize, @intCast(pointer))];
        try std.testing.expect(next >= 0);
        pointer = next;
    }
}
