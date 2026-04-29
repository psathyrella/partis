/// ham/state.zig — Zig port of ham/src/state.cc + ham/include/state.h
///
/// HMM state: name, transitions, emission probabilities, and graph connectivity.
///
/// C++ source: packages/ham/src/state.cc, packages/ham/include/state.h
/// C++ author: psathyrella/ham

const std = @import("std");
const Track = @import("track.zig").Track;
const track_mod = @import("track.zig");
const Emission = @import("emission.zig").Emission;
const Transition = @import("transitions.zig").Transition;
const Sequences = @import("sequences.zig").Sequences;
const mathutils = @import("mathutils.zig");

const log = mathutils.log;
const exp = mathutils.exp;

/// Matches C++ build define -DSTATE_MAX=500.
pub const STATE_MAX: usize = 500;

/// EPS for normalization check (matches C++ -DEPS=1e-6).
const EPS: f64 = 1e-6;

/// Bitset over STATE_MAX bits. Uses std.StaticBitSet for compact representation.
/// Corresponds to C++ `std::bitset<STATE_MAX>`.
pub const StateBitset = std.StaticBitSet(STATE_MAX);

/// A single HMM state.
/// Corresponds to C++ `ham::State`.
pub const State = struct {
    name: []u8,
    /// The germline nucleotide for this position (empty if none).
    germline_nuc: []u8,
    /// Log-probability for the ambiguous character emission (-inf if not set).
    ambiguous_emission_logprob: f64,
    /// The ambiguous character string (empty if not set).
    ambiguous_char: []u8,

    /// Transitions to other states (not the "end" state), in parse order.
    /// After reorderTransitions(), reindexed to model order (with null sentinels for absent).
    /// Each pointer is owned by this State.
    transitions: std.ArrayListUnmanaged(?*Transition),
    /// Transition to the "end" state (null if none).
    /// Owned by this State.
    trans_to_end: ?*Transition,

    /// Emission model for this state.
    emission: Emission,

    /// Index of this state in the model's state vector (set by model finalize).
    index: usize,

    /// Bitset of states this state can transition *to* (set by model finalize).
    to_states: StateBitset,
    /// Bitset of states that can transition *to* this state (set by model finalize).
    from_states: StateBitset,
    /// Indices of states in `from_states` (faster iteration than bitset scan).
    from_state_indices: std.ArrayListUnmanaged(usize),
    /// Indices of states in `to_states` (faster iteration than bitset scan).
    to_state_indices: std.ArrayListUnmanaged(usize),

    /// Create an empty State.
    /// Corresponds to C++ `State()`.
    pub fn init() State {
        return State{
            .name = &[_]u8{},
            .germline_nuc = &[_]u8{},
            .ambiguous_emission_logprob = mathutils.NEG_INF,
            .ambiguous_char = &[_]u8{},
            .transitions = .{},
            .trans_to_end = null,
            .emission = Emission.init(),
            .index = std.math.maxInt(usize),
            .to_states = StateBitset.initEmpty(),
            .from_states = StateBitset.initEmpty(),
            .from_state_indices = .{},
            .to_state_indices = .{},
        };
    }

    pub fn deinit(self: *State, allocator: std.mem.Allocator) void {
        if (self.name.len > 0) allocator.free(self.name);
        if (self.germline_nuc.len > 0) allocator.free(self.germline_nuc);
        if (self.ambiguous_char.len > 0) allocator.free(self.ambiguous_char);
        for (self.transitions.items) |maybe_t| {
            if (maybe_t) |t| {
                t.deinit(allocator);
                allocator.destroy(t);
            }
        }
        self.transitions.deinit(allocator);
        if (self.trans_to_end) |t| {
            t.deinit(allocator);
            allocator.destroy(t);
        }
        self.emission.deinit(allocator);
        self.from_state_indices.deinit(allocator);
        self.to_state_indices.deinit(allocator);
    }

    /// Parse state data from already-extracted fields.
    ///
    /// `transitions_map` maps to_state_name → raw probability.
    /// `emission_probs` maps symbol_name → raw emission probability
    ///   (pass null for the "init" state which has no emissions).
    /// `state_names` is the list of valid non-"end" state names (for validation).
    ///
    /// Corresponds to C++ `State::Parse(YAML::Node, vector<string>, Track*)`.
    pub fn parse(
        self: *State,
        allocator: std.mem.Allocator,
        name: []const u8,
        germline_nuc: ?[]const u8,
        ambiguous_emission_prob: ?f64,
        ambiguous_char: ?[]const u8,
        transitions_map: std.StringHashMap(f64),
        emission_probs: ?std.StringHashMap(f64),
        state_names: []const []const u8,
        track: *const Track,
    ) !void {
        std.debug.assert(name.len > 0);
        if (self.name.len > 0) allocator.free(self.name);
        self.name = try allocator.dupe(u8, name);

        if (germline_nuc) |gn| {
            if (self.germline_nuc.len > 0) allocator.free(self.germline_nuc);
            self.germline_nuc = try allocator.dupe(u8, gn);
        }
        if (ambiguous_emission_prob) |aep| {
            self.ambiguous_emission_logprob = log(aep);
        }
        if (ambiguous_char) |ac| {
            if (self.ambiguous_char.len > 0) allocator.free(self.ambiguous_char);
            self.ambiguous_char = try allocator.dupe(u8, ac);
        }

        // Parse transitions
        var total: f64 = 0.0;
        var it = transitions_map.iterator();
        while (it.next()) |entry| {
            const to_state = entry.key_ptr.*;
            const prob = entry.value_ptr.*;
            // Validate: must be "end" or a known state name
            if (!std.mem.eql(u8, to_state, "end")) {
                var found = false;
                for (state_names) |sn| {
                    if (std.mem.eql(u8, sn, to_state)) {
                        found = true;
                        break;
                    }
                }
                if (!found) return error.UnknownTransitionTarget;
            }
            total += prob;
            const trans = try allocator.create(Transition);
            trans.* = try Transition.init(allocator, to_state, prob);
            if (std.mem.eql(u8, to_state, "end")) {
                self.trans_to_end = trans;
            } else {
                try self.transitions.append(allocator, trans);
            }
        }
        if (@abs(total - 1.0) >= EPS) {
            return error.BadNormalization;
        }

        // Parse emissions (skip for "init" state)
        if (std.mem.eql(u8, name, "init")) return;

        const ep = emission_probs orelse return error.MissingEmissions;
        try self.emission.parse(allocator, ep, track);
    }

    /// Return the state's name.
    pub fn getName(self: *const State) []const u8 {
        return self.name;
    }

    /// First character of name (C++ `abbreviation()`).
    pub fn abbreviation(self: *const State) u8 {
        return if (self.name.len > 0) self.name[0] else '?';
    }

    /// Emission log-probability for a single digitized character.
    /// Handles ambiguous character specially.
    /// Corresponds to C++ `State::EmissionLogprob(uint8_t ch)`.
    pub fn emissionLogprob(self: *const State, ch: u8) f64 {
        if (self.ambiguous_char.len > 0 and ch == track_mod.AMBIGUOUS_INDEX) {
            return self.ambiguous_emission_logprob;
        }
        return self.emission.scoreIndex(ch);
    }

    /// Emission log-probability for a position in a Sequences object.
    /// Accumulates log-probs for each sequence at that position.
    /// Corresponds to C++ `State::EmissionLogprob(Sequences*, size_t pos)`.
    pub fn emissionLogprobSeqs(self: *const State, seqs: *const Sequences, pos: usize) f64 {
        var logprob: f64 = 0.0;
        for (0..seqs.nSeqs()) |iseq| {
            const seq = seqs.get(iseq);
            logprob = mathutils.addWithMinusInfinities(logprob, self.emissionLogprob(seq.seqq[pos]));
        }
        return logprob;
    }

    /// Transition log-probability to state at index `to_state` (post-reorder).
    /// Returns -inf if no transition exists to that state.
    pub inline fn transitionLogprob(self: *const State, to_state: usize) f64 {
        if (to_state >= self.transitions.items.len) return mathutils.NEG_INF;
        const maybe_t = self.transitions.items[to_state];
        return if (maybe_t) |t| t.log_prob else mathutils.NEG_INF;
    }

    /// Log-probability for the transition to "end" (-inf if no such transition).
    /// Corresponds to C++ `State::end_transition_logprob()`.
    pub fn endTransitionLogprob(self: *const State) f64 {
        if (self.trans_to_end) |t| return t.log_prob;
        return mathutils.NEG_INF;
    }

    /// Mark that this state transitions to `st`.
    /// Corresponds to C++ `State::AddToState(State*)`.
    pub fn addToState(self: *State, st: *const State) void {
        self.to_states.set(st.index);
    }

    /// Mark that `st` transitions to this state.
    /// Corresponds to C++ `State::AddFromState(State*)`.
    pub fn addFromState(self: *State, st: *const State) void {
        self.from_states.set(st.index);
    }

    /// Set this state's index.
    pub fn setIndex(self: *State, val: usize) void {
        self.index = val;
    }

    /// Reorder `transitions` to indexed order matching `state_indices` map.
    /// Resizes transitions to `n_states - 1` with null sentinels for absent transitions.
    /// Corresponds to C++ `State::ReorderTransitions(map<string,State*>&)`.
    pub fn reorderTransitions(
        self: *State,
        allocator: std.mem.Allocator,
        state_indices: *const std.StringHashMap(*State),
    ) !void {
        const n_states = state_indices.count();
        // Collect existing transitions before clearing
        var old_list = self.transitions;
        self.transitions = .{};
        defer old_list.deinit(allocator);

        // Build null-filled fixed array of size n_states-1
        const fixed_len = n_states - 1;
        var fixed = try allocator.alloc(?*Transition, fixed_len);
        defer allocator.free(fixed);
        for (fixed) |*p| p.* = null;

        for (old_list.items) |maybe_t| {
            const t = maybe_t orelse continue;
            const to_state = state_indices.get(t.to_state_name) orelse return error.UnknownTransitionTarget;
            fixed[to_state.index] = t;
        }

        // Populate new transitions list
        for (fixed) |maybe_t| {
            try self.transitions.append(allocator, maybe_t);
        }
    }

    /// Build `from_state_indices` from `from_states` bitset.
    /// Corresponds to C++ `State::SetFromStateIndices()`.
    pub fn setFromStateIndices(self: *State, allocator: std.mem.Allocator) !void {
        self.from_state_indices.clearRetainingCapacity();
        for (0..STATE_MAX) |i| {
            if (self.from_states.isSet(i)) {
                try self.from_state_indices.append(allocator, i);
            }
        }
    }

    /// Build `to_state_indices` from `to_states` bitset.
    pub fn setToStateIndices(self: *State, allocator: std.mem.Allocator) !void {
        self.to_state_indices.clearRetainingCapacity();
        for (0..STATE_MAX) |i| {
            if (self.to_states.isSet(i)) {
                try self.to_state_indices.append(allocator, i);
            }
        }
    }

    /// Rescale emission probabilities for overall mutation frequency.
    /// Corresponds to C++ `State::RescaleOverallMuteFreq(double factor)`.
    pub fn rescaleOverallMuteFreq(self: *State, allocator: std.mem.Allocator, factor: f64) !void {
        // Skip if germline is ambiguous or unset
        if (self.germline_nuc.len == 0 or
            (self.ambiguous_char.len > 0 and std.mem.eql(u8, self.germline_nuc, self.ambiguous_char)))
            return;

        const trk = self.emission.getTrack() orelse return;
        const alpha_size = trk.alphabetSize();
        var new_log_probs = try self.emission.logProbs(allocator);
        defer allocator.free(new_log_probs);
        std.debug.assert(new_log_probs.len == alpha_size);

        // Find old mutation frequency
        var old_mute_freq: f64 = -1.0;
        for (0..alpha_size) |ip| {
            const is_germ = std.mem.eql(u8, trk.symbol(ip), self.germline_nuc);
            if (is_germ) {
                old_mute_freq = 1.0 - exp(new_log_probs[ip]);
            }
        }
        std.debug.assert(old_mute_freq >= 0.0);
        var new_mute_freq = @min(0.95, factor * old_mute_freq);
        new_mute_freq = @max(EPS, new_mute_freq);
        if (new_mute_freq <= 0.0 or new_mute_freq >= 1.0) return error.BadMuteFreq;

        for (0..alpha_size) |ip| {
            const is_germ = std.mem.eql(u8, trk.symbol(ip), self.germline_nuc);
            if (is_germ) {
                new_log_probs[ip] = log(1.0 - new_mute_freq);
            } else {
                new_log_probs[ip] = log(exp(new_log_probs[ip]) * new_mute_freq / old_mute_freq);
            }
        }
        try self.emission.replaceLogProbs(allocator, new_log_probs);
    }

    /// Revert emission rescaling.
    /// Corresponds to C++ `State::UnRescaleOverallMuteFreq()`.
    pub fn unRescaleOverallMuteFreq(self: *State, allocator: std.mem.Allocator) void {
        if (self.germline_nuc.len == 0 or
            (self.ambiguous_char.len > 0 and std.mem.eql(u8, self.germline_nuc, self.ambiguous_char)))
            return;
        self.emission.unReplaceLogProbs(allocator);
    }

    /// Print this state to stderr.
    /// Corresponds to C++ `State::Print()`.
    pub fn print(self: *const State) void {
        std.debug.print("state: {s}", .{self.name});
        if (self.germline_nuc.len > 0) std.debug.print(" ({s})", .{self.germline_nuc});
        std.debug.print("\n", .{});
        std.debug.print("  transitions:\n", .{});
        for (self.transitions.items) |maybe_t| {
            if (maybe_t) |t| t.print();
        }
        if (self.trans_to_end) |t| t.print();
        if (std.mem.eql(u8, self.name, "init")) return;
        std.debug.print("  emissions:\n", .{});
        self.emission.print();
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "State: init state parse (no emissions)" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var state = State.init();
    defer state.deinit(allocator);

    var trans = std.StringHashMap(f64).init(allocator);
    defer trans.deinit();
    try trans.put("state1", 0.7);
    try trans.put("end", 0.3);

    const snames = [_][]const u8{"state1"};
    try state.parse(allocator, "init", null, null, null, trans, null, &snames, &track);

    try std.testing.expectEqualStrings("init", state.name);
    try std.testing.expectEqual(@as(usize, 1), state.transitions.items.len);
    try std.testing.expect(state.trans_to_end != null);
}

test "State: normal state parse with emissions" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var state = State.init();
    defer state.deinit(allocator);

    var trans = std.StringHashMap(f64).init(allocator);
    defer trans.deinit();
    try trans.put("end", 1.0);

    var em_probs = std.StringHashMap(f64).init(allocator);
    defer em_probs.deinit();
    try em_probs.put("A", 0.25);
    try em_probs.put("C", 0.25);
    try em_probs.put("G", 0.25);
    try em_probs.put("T", 0.25);

    const snames = [_][]const u8{};
    try state.parse(allocator, "state0", "A", null, "N", trans, em_probs, &snames, &track);

    try std.testing.expectEqualStrings("state0", state.name);
    try std.testing.expectEqualStrings("A", state.germline_nuc);
    try std.testing.expectApproxEqAbs(@log(0.25), state.emissionLogprob(0), 1e-9);
    // end transition stored in trans_to_end
    try std.testing.expectApproxEqAbs(@log(1.0), state.endTransitionLogprob(), 1e-9);
}

test "State: endTransitionLogprob absent returns -inf" {
    const allocator = std.testing.allocator;
    const sym_list = [_][]const u8{ "A", "C", "G", "T" };
    var track = try Track.initWithSymbols(allocator, "nucs", &sym_list, "N");
    defer track.deinit(allocator);

    var state = State.init();
    defer state.deinit(allocator);

    var trans = std.StringHashMap(f64).init(allocator);
    defer trans.deinit();
    try trans.put("state_x", 0.5);
    try trans.put("state_y", 0.5);

    var em_probs = std.StringHashMap(f64).init(allocator);
    defer em_probs.deinit();
    try em_probs.put("A", 0.25);
    try em_probs.put("C", 0.25);
    try em_probs.put("G", 0.25);
    try em_probs.put("T", 0.25);

    const snames = [_][]const u8{ "state_x", "state_y" };
    try state.parse(allocator, "mystate", null, null, null, trans, em_probs, &snames, &track);

    try std.testing.expect(std.math.isNegativeInf(state.endTransitionLogprob()));
}
