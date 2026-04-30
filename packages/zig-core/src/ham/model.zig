/// ham/model.zig — Zig port of ham/src/model.cc + ham/include/model.h
///
/// HMM model: loads from a YAML file, holds all states, finalizes connectivity.
///
/// C++ source: packages/ham/src/model.cc, packages/ham/include/model.h
/// C++ author: psathyrella/ham

const std = @import("std");
const Track = @import("track.zig").Track;
const track_mod = @import("track.zig");
const State = @import("state.zig").State;
const Transition = @import("transitions.zig").Transition;
const hmm_yaml = @import("hmm_yaml.zig");
const mathutils = @import("mathutils.zig");

/// EPS for normalization checks (-DEPS=1e-6).
const EPS: f64 = 1e-6;

/// HMM model: all states, a track, and finalized connectivity.
/// Corresponds to C++ `ham::Model`.
pub const Model = struct {
    name: []u8,
    /// Overall probability of this gene/HMM (from YAML extras.gene_prob).
    overall_prob: f64,
    /// Original mean mutation frequency (from YAML extras.overall_mute_freq).
    original_overall_mute_freq: f64,
    /// Ratio by which we've rescaled emissions (-inf if not rescaled).
    rescale_ratio: f64,
    /// Ambiguous character (e.g. "N"). Empty if not set.
    ambiguous_char: []u8,

    /// Track (alphabet). Owned by this Model.
    track: ?*Track,

    /// All non-init states in model order. Each pointer is heap-allocated and owned.
    states: std.ArrayListUnmanaged(*State),
    /// Map from state name to State pointer (pointers owned via `states` or `initial`).
    states_by_name: std.StringHashMap(*State),
    /// The "init" state. Owned.
    initial: ?*State,
    /// The synthetic "end" state (no emissions, only from-states). Owned.
    ending: *State,

    finalized: bool,

    /// Create an empty Model (not yet parsed).
    /// Corresponds to C++ `Model()`.
    pub fn init(allocator: std.mem.Allocator) !Model {
        const ending = try allocator.create(State);
        ending.* = State.init();
        return Model{
            .name = try allocator.dupe(u8, ""),
            .overall_prob = 0.0,
            .original_overall_mute_freq = 0.0,
            .rescale_ratio = mathutils.NEG_INF,
            .ambiguous_char = try allocator.dupe(u8, ""),
            .track = null,
            .states = .{},
            .states_by_name = std.StringHashMap(*State).init(allocator),
            .initial = null,
            .ending = ending,
            .finalized = false,
        };
    }

    pub fn deinit(self: *Model, allocator: std.mem.Allocator) void {
        allocator.free(self.name);
        allocator.free(self.ambiguous_char);
        // Free all states via states_by_name (owns all State pointers)
        var it = self.states_by_name.valueIterator();
        while (it.next()) |st_ptr| {
            st_ptr.*.deinit(allocator);
            allocator.destroy(st_ptr.*);
        }
        self.states_by_name.deinit();
        self.states.deinit(allocator);
        // ending is separate
        self.ending.deinit(allocator);
        allocator.destroy(self.ending);
        if (self.track) |t| {
            t.deinit(allocator);
            allocator.destroy(t);
        }
    }

    /// Parse a model from a YAML file.
    /// Corresponds to C++ `Model::Parse(string infname)`.
    pub fn parse(self: *Model, allocator: std.mem.Allocator, infname: []const u8) !void {
        var data = try hmm_yaml.parseFile(allocator, infname);
        defer data.deinit(allocator);

        // Model-wide info
        allocator.free(self.name);
        self.name = try allocator.dupe(u8, data.name);
        self.overall_prob = data.overall_prob;
        self.original_overall_mute_freq = data.original_overall_mute_freq;
        if (data.ambiguous_char.len > 0) {
            allocator.free(self.ambiguous_char);
            self.ambiguous_char = try allocator.dupe(u8, data.ambiguous_char);
        }

        // Build track
        if (data.tracks.items.len != 1) return error.ExpectedExactlyOneTrack;
        const td = &data.tracks.items[0];
        const trk = try allocator.create(Track);
        // NOTE: assign to self.track immediately so model.deinit handles cleanup on error.
        // Do NOT use an errdefer for trk — model.deinit owns it via self.track.
        trk.* = try Track.init(allocator, td.name);
        if (self.track) |old_t| {
            old_t.deinit(allocator);
            allocator.destroy(old_t);
        }
        self.track = trk;
        if (self.ambiguous_char.len > 0) {
            try trk.setAmbiguous(allocator, self.ambiguous_char);
        }
        for (td.symbols.items) |sym| try trk.addSymbol(allocator, sym);

        // Collect state names for transition validation
        var state_names: std.ArrayListUnmanaged([]const u8) = .{};
        defer state_names.deinit(allocator);
        for (data.states.items) |*sd| try state_names.append(allocator, sd.name);

        // Parse each state
        for (data.states.items) |*sd| {
            const st = try allocator.create(State);
            errdefer {
                st.deinit(allocator);
                allocator.destroy(st);
            }
            st.* = State.init();

            try st.parse(
                allocator,
                sd.name,
                sd.germline_nuc,
                sd.ambiguous_emission_prob,
                sd.ambiguous_char,
                sd.transitions,
                if (sd.emissions) |*em| em.probs else null,
                state_names.items,
                trk,
            );

            if (std.mem.eql(u8, sd.name, "init")) {
                self.initial = st;
            } else {
                if (self.states.items.len >= @import("state.zig").STATE_MAX) return error.TooManyStates;
                try self.states.append(allocator, st);
            }
            try self.states_by_name.put(st.name, st);
        }

        try self.finalize(allocator);
    }

    /// Post-process: set indices, compute to/from bitsets, reorder transitions.
    /// Corresponds to C++ `Model::Finalize()`.
    fn finalize(self: *Model, allocator: std.mem.Allocator) !void {
        std.debug.assert(!self.finalized);

        // Assign indices to all non-init states
        for (self.states.items, 0..) |st, i| st.setIndex(i);

        // Set to/from connectivity for each non-init state
        for (self.states.items) |st| try self.finalizeState(st);
        if (self.initial) |init_st| try self.finalizeState(init_st);

        // Reorder transitions in index order
        for (self.states.items) |st| try st.reorderTransitions(allocator, &self.states_by_name);
        if (self.initial) |init_st| try init_st.reorderTransitions(allocator, &self.states_by_name);

        // Check topology
        try self.checkTopology(allocator);

        // Build from_state_indices
        try self.addMaybeFasterFromStateStuff(allocator);

        // Materialize flat log_probs and free *Transition allocations. Must run
        // after every reader of to_state_name/to_state_ptr (finalizeState,
        // reorderTransitions, checkTopology, addMaybeFasterFromStateStuff).
        for (self.states.items) |st| try st.materializeTransitionLogProbs(allocator);
        if (self.initial) |init_st| try init_st.materializeTransitionLogProbs(allocator);

        self.finalized = true;
    }

    /// Set to/from bitsets for a state and update destination states' from-bitsets.
    /// Corresponds to C++ `Model::FinalizeState(State*)`.
    fn finalizeState(self: *Model, st: *State) !void {
        // Iterate over parse-order transitions (before reorder)
        for (st.transitions.items) |maybe_t| {
            const t = maybe_t orelse continue;
            const to_state = self.states_by_name.get(t.to_state_name) orelse return error.UnknownTransitionTarget;
            st.addToState(to_state);
            t.setToState(to_state);
            if (st != self.initial) to_state.addFromState(st);
        }
        if (st.trans_to_end) |_| {
            self.ending.addFromState(st);
        }
    }

    /// Build from_state_indices and to_state_indices for all states.
    /// Corresponds to C++ `Model::AddMaybeFasterFromStateStuff()`.
    fn addMaybeFasterFromStateStuff(self: *Model, allocator: std.mem.Allocator) !void {
        if (self.initial) |init_st| {
            try init_st.setFromStateIndices(allocator);
            try init_st.setToStateIndices(allocator);
        }
        for (self.states.items) |st| {
            try st.setFromStateIndices(allocator);
            try st.setToStateIndices(allocator);
        }
        try self.ending.setFromStateIndices(allocator);
        try self.ending.setToStateIndices(allocator);
    }

    /// Topology check: all states reachable from init, each has a non-self path.
    /// Corresponds to C++ `Model::CheckTopology()`.
    fn checkTopology(self: *Model, allocator: std.mem.Allocator) !void {
        const init_st = self.initial orelse return error.NoInitState;
        const n = self.states.items.len;

        var states_to_check: std.ArrayListUnmanaged(u16) = .{};
        defer states_to_check.deinit(allocator);
        try self.addToStateIndicesHelper(allocator, init_st, &states_to_check);

        var checked = try allocator.alloc(bool, n);
        defer allocator.free(checked);
        for (checked) |*c| c.* = false;

        while (states_to_check.items.len > 0) {
            const icheck = states_to_check.pop() orelse break;
            if (checked[icheck]) continue;

            var tmp: std.ArrayListUnmanaged(u16) = .{};
            defer tmp.deinit(allocator);
            try self.addToStateIndicesHelper(allocator, self.states.items[icheck], &tmp);
            const num_visited = tmp.items.len;

            if (num_visited == 0) {
                if (self.states.items[icheck].trans_to_end == null)
                    return error.StateHasNoTransitions;
            } else if (num_visited == 1 and tmp.items[0] == icheck) {
                if (self.states.items[icheck].trans_to_end == null)
                    return error.StateOnlySelfTransition;
            }
            checked[icheck] = true;
            for (tmp.items) |idx| {
                if (!checked[idx]) try states_to_check.append(allocator, idx);
            }
        }

        // At least one transition to end
        var found_end = false;
        for (self.states.items) |st| {
            if (st.trans_to_end != null) {
                found_end = true;
                break;
            }
        }
        if (!found_end) return error.NoTransitionToEnd;

        // All states reachable
        for (checked, 0..) |c, i| {
            if (!c) {
                std.debug.print("ERROR: state '{s}' not reachable from init\n", .{self.states.items[i].name});
                return error.StateNotReachable;
            }
        }
    }

    fn addToStateIndicesHelper(_: *const Model, allocator: std.mem.Allocator, st: *const State, visited: *std.ArrayListUnmanaged(u16)) !void {
        for (st.transitions.items) |maybe_t| {
            if (maybe_t) |t| {
                const to_st = t.to_state_ptr orelse continue;
                try visited.append(allocator, @intCast(to_st.index));
            }
        }
    }

    /// Rescale all state emissions by a factor derived from `overall_mute_freq`.
    /// Corresponds to C++ `Model::RescaleOverallMuteFreq(double)`.
    pub fn rescaleOverallMuteFreq(self: *Model, allocator: std.mem.Allocator, overall_mute_freq: f64) !void {
        std.debug.assert(!std.math.isNegativeInf(overall_mute_freq));
        if (self.original_overall_mute_freq == 0.0) return error.ZeroOriginalMuteFreq;
        const factor = @max(0.01, overall_mute_freq) / self.original_overall_mute_freq;
        for (self.states.items) |st| try st.rescaleOverallMuteFreq(allocator, factor);
    }

    /// Undo rescaling.
    /// Corresponds to C++ `Model::UnRescaleOverallMuteFreq()`.
    pub fn unRescaleOverallMuteFreq(self: *Model, allocator: std.mem.Allocator) void {
        for (self.states.items) |st| st.unRescaleOverallMuteFreq(allocator);
    }

    pub fn nStates(self: *const Model) usize {
        return self.states.items.len;
    }

    pub fn stateByName(self: *const Model, name: []const u8) ?*State {
        return self.states_by_name.get(name);
    }

    pub inline fn stateByIndex(self: *const Model, idx: usize) *State {
        return self.states.items[idx];
    }

    pub fn initState(self: *const Model) ?*State {
        return self.initial;
    }

    pub fn initialToStates(self: *const Model) ?*@import("state.zig").StateBitset {
        if (self.initial) |init_st| return &init_st.to_states;
        return null;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Model: parse simple 2-state HMM from string" {
    const allocator = std.testing.allocator;

    // Write a minimal YAML to a temp file
    const yaml = "!!python/object:python.hmmwriter.HMM\n" ++
        "extras:\n" ++
        "  gene_prob: 0.5\n" ++
        "  overall_mute_freq: 0.1\n" ++
        "name: test_gene\n" ++
        "tracks:\n" ++
        "  nukes:\n" ++
        "  - A\n" ++
        "  - C\n" ++
        "  - G\n" ++
        "  - T\n" ++
        "states:\n" ++
        "- !!python/object:python.hmmwriter.State\n" ++
        "  emissions: null\n" ++
        "  extras: {}\n" ++
        "  name: init\n" ++
        "  transitions:\n" ++
        "    state0: 1.0\n" ++
        "- !!python/object:python.hmmwriter.State\n" ++
        "  emissions:\n" ++
        "    probs:\n" ++
        "      A: 0.25\n" ++
        "      C: 0.25\n" ++
        "      G: 0.25\n" ++
        "      T: 0.25\n" ++
        "    track: nukes\n" ++
        "  extras:\n" ++
        "    ambiguous_char: N\n" ++
        "    ambiguous_emission_prob: 0.25\n" ++
        "    germline: A\n" ++
        "  name: state0\n" ++
        "  transitions:\n" ++
        "    end: 1.0\n";

    // Write to a temporary file
    const tmp_path = "/tmp/test_model.yaml";
    {
        const f = try std.fs.createFileAbsolute(tmp_path, .{});
        defer f.close();
        try f.writeAll(yaml);
    }
    defer std.fs.deleteFileAbsolute(tmp_path) catch {};

    var model = try Model.init(allocator);
    defer model.deinit(allocator);

    try model.parse(allocator, tmp_path);

    try std.testing.expectEqualStrings("test_gene", model.name);
    try std.testing.expectApproxEqAbs(0.5, model.overall_prob, 1e-9);
    try std.testing.expectEqual(@as(usize, 1), model.nStates());
    try std.testing.expect(model.initial != null);
    try std.testing.expect(model.stateByName("state0") != null);
    try std.testing.expect(model.stateByName("init") != null);
}
