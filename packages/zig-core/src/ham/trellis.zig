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
    /// Sequences to run DP over (owned copy).
    seqs: Sequences,

    /// Traceback table [position][state] → previous state (-1 if none).
    traceback_table: Int2D,
    /// Pointer to the active traceback table (self's or cached trellis's).
    traceback_table_ptr: ?*Int2D,

    /// Pointer to another trellis with pre-filled DP tables (not owned).
    cached_trellis: ?*const Trellis,

    ending_viterbi_pointer: i16,
    ending_viterbi_log_prob: f64,
    ending_forward_log_prob: f64,

    /// Per-position Viterbi log-probs (including ending transition).
    viterbi_log_probs: std.ArrayListUnmanaged(f64),
    viterbi_log_probs_ptr: ?*std.ArrayListUnmanaged(f64),
    /// Per-position Forward log-probs.
    forward_log_probs: std.ArrayListUnmanaged(f64),
    forward_log_probs_ptr: ?*std.ArrayListUnmanaged(f64),
    /// Per-position best-state index (for Viterbi).
    viterbi_indices: std.ArrayListUnmanaged(i32),
    viterbi_indices_ptr: ?*std.ArrayListUnmanaged(i32),

    /// Current and previous scoring columns (size = n_states each).
    scoring_current: []f64,
    scoring_previous: []f64,

    /// Scratch boolean array for deduplicating next_states additions.
    next_seen: [state_mod.STATE_MAX]bool,

    allocator: std.mem.Allocator,

    /// Create a Trellis for a single Sequence.
    /// The sequence is copied into an internal Sequences object.
    /// Corresponds to C++ `Trellis(Model*, Sequence, Trellis*)`.
    pub fn initWithSeq(allocator: std.mem.Allocator, hmm: *const Model, seq: Sequence, cached: ?*const Trellis) !Trellis {
        var seqs = Sequences.init();
        try seqs.addSeq(allocator, try seq.clone(allocator));
        return initImpl(allocator, hmm, seqs, cached);
    }

    /// Create a Trellis for a Sequences object.
    /// Clones `seqs` so the caller retains ownership of the original
    /// (matches C++ Trellis(Model*, Sequences, Trellis*) copy-by-value semantics).
    pub fn initWithSeqs(allocator: std.mem.Allocator, hmm: *const Model, seqs: Sequences, cached: ?*const Trellis) !Trellis {
        const seqs_copy = try seqs.clone(allocator);
        return initImpl(allocator, hmm, seqs_copy, cached);
    }

    fn initImpl(allocator: std.mem.Allocator, hmm: *const Model, seqs: Sequences, cached: ?*const Trellis) !Trellis {
        const n = hmm.nStates();
        const scoring_current = try allocator.alloc(f64, n);
        const scoring_previous = try allocator.alloc(f64, n);
        @memset(scoring_current, -std.math.inf(f64));
        @memset(scoring_previous, -std.math.inf(f64));

        const t = Trellis{
            .hmm = hmm,
            .seqs = seqs,
            .traceback_table = .{},
            .traceback_table_ptr = null,
            .cached_trellis = cached,
            .ending_viterbi_pointer = -1,
            .ending_viterbi_log_prob = -std.math.inf(f64),
            .ending_forward_log_prob = -std.math.inf(f64),
            .viterbi_log_probs = .{},
            .viterbi_log_probs_ptr = null,
            .forward_log_probs = .{},
            .forward_log_probs_ptr = null,
            .viterbi_indices = .{},
            .viterbi_indices_ptr = null,
            .scoring_current = scoring_current,
            .scoring_previous = scoring_previous,
            .next_seen = [_]bool{false} ** state_mod.STATE_MAX,
            .allocator = allocator,
        };
        return t;
    }

    /// Fix up self-referential _ptr fields after the trellis has been copied/moved.
    /// Must be called after any operation that may have moved this Trellis in memory.
    pub fn fixupPtrs(self: *Trellis) void {
        if (self.viterbi_log_probs_ptr != null) self.viterbi_log_probs_ptr = &self.viterbi_log_probs;
        if (self.viterbi_indices_ptr != null) self.viterbi_indices_ptr = &self.viterbi_indices;
        if (self.forward_log_probs_ptr != null) self.forward_log_probs_ptr = &self.forward_log_probs;
        if (self.traceback_table_ptr != null) self.traceback_table_ptr = &self.traceback_table;
    }

    pub fn deinit(self: *Trellis) void {
        const allocator = self.allocator;
        self.seqs.deinit(allocator);
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
        const seq_len = self.seqs.sequence_length;
        const n_states = self.hmm.nStates();

        if (self.cached_trellis) |ct| {
            // Poach values from cached trellis
            self.traceback_table_ptr = @constCast(&ct.traceback_table);
            self.ending_viterbi_pointer = @intCast(ct.viterbiIndicesAt(seq_len));
            self.ending_viterbi_log_prob = ct.endingViterbiLogProbAt(seq_len);
            self.viterbi_log_probs_ptr = @constCast(&ct.viterbi_log_probs);
            self.viterbi_indices_ptr = @constCast(&ct.viterbi_indices);
            return;
        }

        // Initialize storage
        try self.viterbi_log_probs.resize(allocator, seq_len);
        try self.viterbi_indices.resize(allocator, seq_len);
        @memset(self.viterbi_log_probs.items, -std.math.inf(f64));
        @memset(self.viterbi_indices.items, -1);
        self.viterbi_log_probs_ptr = &self.viterbi_log_probs;
        self.viterbi_indices_ptr = &self.viterbi_indices;

        // Initialize traceback table [seq_len][n_states] = -1
        for (self.traceback_table.items) |row| allocator.free(row);
        self.traceback_table.clearRetainingCapacity();
        for (0..seq_len) |_| {
            const row = try allocator.alloc(i16, n_states);
            @memset(row, -1);
            try self.traceback_table.append(allocator, row);
        }
        self.traceback_table_ptr = &self.traceback_table;

        // Reset scoring columns
        @memset(self.scoring_current, -std.math.inf(f64));
        @memset(self.scoring_previous, -std.math.inf(f64));

        var current_states = std.ArrayListUnmanaged(usize){};
        defer current_states.deinit(allocator);
        var next_states = std.ArrayListUnmanaged(usize){};
        defer next_states.deinit(allocator);

        // Position 0: transitions from init state
        const init_st = self.hmm.initial orelse return error.NoInitState;
        // Reset next_seen before populating position 0
        @memset(&self.next_seen, false);
        for (init_st.to_state_indices.items) |i_st| {
            const emission_val = self.hmm.stateByIndex(i_st).emissionLogprobSeqs(&self.seqs, 0);
            const dpval = emission_val + init_st.transitionLogprob(i_st);
            if (std.math.isNegativeInf(dpval)) continue;
            self.scoring_current[i_st] = dpval;
            self.cacheViterbiVals(0, dpval, i_st);
            // Mark outbound transitions for next position (deduplicated)
            for (self.hmm.stateByIndex(i_st).to_state_indices.items) |j| {
                if (!self.next_seen[j]) {
                    self.next_seen[j] = true;
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
        self.ending_viterbi_log_prob = -std.math.inf(f64);
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
        const seq_len = self.seqs.sequence_length;
        const n_states = self.hmm.nStates();

        if (self.cached_trellis) |ct| {
            self.ending_forward_log_prob = ct.endingForwardLogProbAt(seq_len);
            self.forward_log_probs_ptr = @constCast(&ct.forward_log_probs);
            return;
        }

        // Initialize storage
        try self.forward_log_probs.resize(allocator, seq_len);
        @memset(self.forward_log_probs.items, -std.math.inf(f64));
        self.forward_log_probs_ptr = &self.forward_log_probs;

        // Reset scoring columns
        @memset(self.scoring_current, -std.math.inf(f64));
        @memset(self.scoring_previous, -std.math.inf(f64));

        var current_states = std.ArrayListUnmanaged(usize){};
        defer current_states.deinit(allocator);
        var next_states = std.ArrayListUnmanaged(usize){};
        defer next_states.deinit(allocator);

        // Position 0: transitions from init state
        const init_st = self.hmm.initial orelse return error.NoInitState;
        // Reset next_seen before populating position 0
        @memset(&self.next_seen, false);
        for (init_st.to_state_indices.items) |i_st| {
            const emission_val = self.hmm.stateByIndex(i_st).emissionLogprobSeqs(&self.seqs, 0);
            const dpval = emission_val + init_st.transitionLogprob(i_st);
            if (std.math.isNegativeInf(dpval)) continue;
            self.scoring_current[i_st] = dpval;
            // Mark outbound transitions for next position (deduplicated)
            for (self.hmm.stateByIndex(i_st).to_state_indices.items) |j| {
                if (!self.next_seen[j]) {
                    self.next_seen[j] = true;
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
        self.ending_forward_log_prob = -std.math.inf(f64);
        for (0..n_states) |st_prev| {
            if (std.math.isNegativeInf(self.scoring_previous[st_prev])) continue;
            const dpval = self.scoring_previous[st_prev] + self.hmm.stateByIndex(st_prev).endTransitionLogprob();
            if (std.math.isNegativeInf(dpval)) continue;
            self.ending_forward_log_prob = mathutils.add_in_log_space(self.ending_forward_log_prob, dpval);
        }
    }

    /// Traceback to produce a path.
    /// Corresponds to C++ `Trellis::Traceback(TracebackPath&)`.
    pub fn traceback(self: *const Trellis, path: *TracebackPath) !void {
        std.debug.assert(self.seqs.sequence_length != 0);
        path.setModel(self.hmm);
        if (std.math.isNegativeInf(self.ending_viterbi_log_prob)) return;
        path.setScore(self.ending_viterbi_log_prob);
        try path.pushBack(self.allocator, self.ending_viterbi_pointer);

        var pointer: i16 = self.ending_viterbi_pointer;
        var position: usize = self.seqs.sequence_length - 1;
        while (position > 0) : (position -= 1) {
            const tbl = self.traceback_table_ptr orelse return;
            pointer = tbl.items[position][@intCast(pointer)];
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
        @memset(self.scoring_current, -std.math.inf(f64));

        // Clear next_seen for all indices that were in next_states
        // (they will become current_states; next_seen tracks what's queued for the NEW next)
        for (next_states.items) |j| self.next_seen[j] = false;

        // current_states ← next_states; next_states ← empty
        // Swap the backing allocations to avoid reallocation
        const old_current_items = current_states.items;
        const old_current_cap = current_states.capacity;
        current_states.items = next_states.items;
        current_states.capacity = next_states.capacity;
        next_states.items = old_current_items;
        next_states.capacity = old_current_cap;
        next_states.clearRetainingCapacity();
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
            const emission_val = st_cur.emissionLogprobSeqs(&self.seqs, position);
            if (std.math.isNegativeInf(emission_val)) continue;
            for (st_cur.from_state_indices.items) |i_st_prev| {
                const prev_val = self.scoring_previous[i_st_prev];
                if (std.math.isNegativeInf(prev_val)) continue;
                const dpval = prev_val + emission_val + self.hmm.stateByIndex(i_st_prev).transitionLogprob(i_st_cur);
                if (dpval > self.scoring_current[i_st_cur]) {
                    self.scoring_current[i_st_cur] = dpval;
                    if (self.traceback_table_ptr) |tbl| {
                        tbl.items[position][i_st_cur] = @intCast(i_st_prev);
                    }
                }
                self.cacheViterbiVals(position, dpval, i_st_cur);
            }
            // Mark outbound transitions (once per current state, deduplicated via next_seen)
            for (st_cur.to_state_indices.items) |j| {
                if (!self.next_seen[j]) {
                    self.next_seen[j] = true;
                    try next_states.append(allocator, j);
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
            const emission_val = st_cur.emissionLogprobSeqs(&self.seqs, position);
            if (std.math.isNegativeInf(emission_val)) continue;
            for (st_cur.from_state_indices.items) |i_st_prev| {
                const prev_val = self.scoring_previous[i_st_prev];
                if (std.math.isNegativeInf(prev_val)) continue;
                const dpval = prev_val + emission_val + self.hmm.stateByIndex(i_st_prev).transitionLogprob(i_st_cur);
                self.scoring_current[i_st_cur] = mathutils.add_in_log_space(dpval, self.scoring_current[i_st_cur]);
                self.cacheForwardVals(position, dpval, i_st_cur);
            }
            // Mark outbound transitions (once per current state, deduplicated via next_seen)
            for (st_cur.to_state_indices.items) |j| {
                if (!self.next_seen[j]) {
                    self.next_seen[j] = true;
                    try next_states.append(allocator, j);
                }
            }
        }
    }

    fn cacheViterbiVals(self: *Trellis, position: usize, dpval: f64, i_st_cur: usize) void {
        const end_trans = self.hmm.stateByIndex(i_st_cur).endTransitionLogprob();
        const logprob = dpval + end_trans;
        if (logprob > self.viterbi_log_probs.items[position]) {
            self.viterbi_log_probs.items[position] = logprob;
            self.viterbi_indices.items[position] = @intCast(i_st_cur);
        }
    }

    fn cacheForwardVals(self: *Trellis, position: usize, dpval: f64, i_st_cur: usize) void {
        const end_trans = self.hmm.stateByIndex(i_st_cur).endTransitionLogprob();
        const logprob = dpval + end_trans;
        self.forward_log_probs.items[position] = mathutils.add_in_log_space(logprob, self.forward_log_probs.items[position]);
    }

    /// Ending Viterbi log-prob for a specific sequence length (used by cached trellis).
    pub fn endingViterbiLogProbAt(self: *const Trellis, length: usize) f64 {
        if (length == 0) return -std.math.inf(f64);
        if (self.viterbi_log_probs_ptr) |ptr| {
            if (length - 1 < ptr.items.len) return ptr.items[length - 1];
        }
        return -std.math.inf(f64);
    }

    /// Ending Forward log-prob for a specific sequence length.
    pub fn endingForwardLogProbAt(self: *const Trellis, length: usize) f64 {
        if (length == 0) return -std.math.inf(f64);
        if (self.forward_log_probs_ptr) |ptr| {
            if (length - 1 < ptr.items.len) return ptr.items[length - 1];
        }
        return -std.math.inf(f64);
    }

    /// Best state index at a specific sequence length (Viterbi).
    pub fn viterbiIndicesAt(self: *const Trellis, length: usize) i32 {
        if (length == 0) return -1;
        if (self.viterbi_indices_ptr) |ptr| {
            if (length - 1 < ptr.items.len) return ptr.items[length - 1];
        }
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

    var trellis = try Trellis.initWithSeq(allocator, &model, seq, null);
    defer trellis.deinit();

    try trellis.forward();

    // The forward probability should be finite (not -inf) for a valid sequence
    try std.testing.expect(!std.math.isNegativeInf(trellis.ending_forward_log_prob));
    // It should be a plausible log-probability (very negative is OK, but not -inf)
    std.debug.print("forward log-prob: {d}\n", .{trellis.ending_forward_log_prob});
}
