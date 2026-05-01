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

/// 2D traceback table: [seq_len][n_states], each entry is the previous state index (or -1).
/// Corresponds to C++ `typedef vector<vector<int16_t>> int_2D`.
pub const Int2D = std.ArrayListUnmanaged([]i16);

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

    /// Traceback table [position][state] → previous state (-1 if none).
    traceback_table: Int2D,

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

    allocator: std.mem.Allocator,

    /// Create a Trellis borrowing the given Sequences.
    /// The caller retains ownership; `seqs` must outlive this Trellis.
    /// Corresponds to C++ `Trellis(Model*, Sequences, Trellis*)`, but without
    /// the C++ copy-by-value clone — the borrow is sufficient for DP read
    /// access, and for cached/stored trellises only the result arrays are
    /// re-read.
    pub fn initWithSeqs(allocator: std.mem.Allocator, hmm: *const Model, seqs: *const Sequences, cached: ?*const Trellis) !Trellis {
        const n = hmm.nStates();
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
            .traceback_table = .{},
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
        // Free traceback table rows
        for (self.traceback_table.items) |row| allocator.free(row);
        self.traceback_table.deinit(allocator);
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

        // Initialize traceback table [seq_len][n_states] = -1
        for (self.traceback_table.items) |row| allocator.free(row);
        self.traceback_table.clearRetainingCapacity();
        for (0..seq_len) |_| {
            const row = try allocator.alloc(i16, n_states);
            @memset(row, -1);
            try self.traceback_table.append(allocator, row);
        }

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
            const emission_val = self.hmm.stateByIndex(i_st).emissionLogprobSeqs(self.seqs, 0);
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
            const emission_val = self.hmm.stateByIndex(i_st).emissionLogprobSeqs(self.seqs, 0);
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
    pub fn traceback(self: *const Trellis, path: *TracebackPath) !void {
        std.debug.assert(self.seq_len != 0);
        path.setModel(self.hmm);
        if (std.math.isNegativeInf(self.ending_viterbi_log_prob)) return;
        path.setScore(self.ending_viterbi_log_prob);
        try path.pushBack(self.allocator, self.ending_viterbi_pointer);

        var pointer: i16 = self.ending_viterbi_pointer;
        var position: usize = self.seq_len - 1;
        const tbl_items = self.tracebackTableItems() orelse return;
        while (position > 0) : (position -= 1) {
            pointer = tbl_items[position][@intCast(pointer)];
            if (pointer == -1) {
                std.debug.print("No valid path at position {d}\n", .{position});
                return;
            }
            try path.pushBack(self.allocator, pointer);
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
        for (current_states.items) |i_st_cur| {
            const st_cur = self.hmm.stateByIndex(i_st_cur);
            const emission_val = st_cur.emissionLogprobSeqs(self.seqs, position);
            if (std.math.isNegativeInf(emission_val)) continue;
            for (st_cur.from_state_indices.items) |i_st_prev| {
                const prev_val = self.scoring_previous[i_st_prev];
                if (std.math.isNegativeInf(prev_val)) continue;
                const dpval = prev_val + emission_val + self.hmm.stateByIndex(i_st_prev).transitionLogprob(i_st_cur);
                if (dpval > self.scoring_current[i_st_cur]) {
                    self.scoring_current[i_st_cur] = dpval;
                    self.traceback_table.items[position][i_st_cur] = @intCast(i_st_prev);
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
        for (current_states.items) |i_st_cur| {
            const st_cur = self.hmm.stateByIndex(i_st_cur);
            const emission_val = st_cur.emissionLogprobSeqs(self.seqs, position);
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

    fn tracebackTableItems(self: *const Trellis) ?[][]i16 {
        if (self.cached_trellis) |ct| {
            if (ct.traceback_table.items.len > 0) return ct.traceback_table.items;
        }
        if (self.traceback_table.items.len > 0) return self.traceback_table.items;
        return null;
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

    // The forward probability should be finite (not -inf) for a valid sequence
    try std.testing.expect(!std.math.isNegativeInf(trellis.ending_forward_log_prob));
    // It should be a plausible log-probability (very negative is OK, but not -inf)
    std.debug.print("forward log-prob: {d}\n", .{trellis.ending_forward_log_prob});
}
