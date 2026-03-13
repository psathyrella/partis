/// ham/hmm_yaml.zig — Minimal YAML reader for partis HMM parameter files.
///
/// Partis HMM YAML files use Python-tagged objects (`!!python/object:...`).
/// Standard YAML parsers reject these tags, so we implement a hand-rolled
/// line-oriented parser tailored exactly to the known HMM file structure.
///
/// Supported structure (from packages/ham/src/model.cc Parse()):
///   top-level scalar keys: name, extras.gene_prob, extras.overall_mute_freq,
///                          extras.ambiguous_char
///   tracks: map of track_name → list of symbols
///   states: sequence of state objects, each with:
///     name, extras (germline, ambiguous_char, ambiguous_emission_prob),
///     transitions (map of name → prob), emissions (probs map + track name)

const std = @import("std");

/// Parsed emission data for one state.
pub const EmissionData = struct {
    /// Map from symbol name to raw probability.
    /// Keys are owned allocator strings; must deinit keys separately.
    probs: std.StringHashMap(f64),
    /// Track name (e.g. "nukes"). Owned string.
    track_name: []u8,

    pub fn deinit(self: *EmissionData, allocator: std.mem.Allocator) void {
        var it = self.probs.keyIterator();
        while (it.next()) |k| allocator.free(k.*);
        self.probs.deinit();
        allocator.free(self.track_name);
    }
};

/// Parsed data for one HMM state.
pub const StateData = struct {
    name: []u8,
    germline_nuc: ?[]u8,
    ambiguous_char: ?[]u8,
    ambiguous_emission_prob: ?f64,
    /// Map from to_state_name to raw probability.
    /// Keys are owned allocator strings.
    transitions: std.StringHashMap(f64),
    /// null for "init" state or states with null emissions.
    emissions: ?EmissionData,

    pub fn deinit(self: *StateData, allocator: std.mem.Allocator) void {
        allocator.free(self.name);
        if (self.germline_nuc) |g| allocator.free(g);
        if (self.ambiguous_char) |a| allocator.free(a);
        var it = self.transitions.keyIterator();
        while (it.next()) |k| allocator.free(k.*);
        self.transitions.deinit();
        if (self.emissions) |*e| e.deinit(allocator);
    }
};

/// Parsed track data.
pub const TrackData = struct {
    name: []u8,
    symbols: std.ArrayListUnmanaged([]u8),

    pub fn deinit(self: *TrackData, allocator: std.mem.Allocator) void {
        allocator.free(self.name);
        for (self.symbols.items) |s| allocator.free(s);
        self.symbols.deinit(allocator);
    }
};

/// Full parsed HMM model data.
pub const HmmData = struct {
    name: []u8,
    overall_prob: f64,
    original_overall_mute_freq: f64,
    ambiguous_char: []u8,
    tracks: std.ArrayListUnmanaged(TrackData),
    states: std.ArrayListUnmanaged(StateData),

    pub fn deinit(self: *HmmData, allocator: std.mem.Allocator) void {
        allocator.free(self.name);
        allocator.free(self.ambiguous_char);
        for (self.tracks.items) |*t| t.deinit(allocator);
        self.tracks.deinit(allocator);
        for (self.states.items) |*s| s.deinit(allocator);
        self.states.deinit(allocator);
    }
};

/// Parse a partis HMM YAML file.
/// Returns an owned HmmData that the caller must deinit.
pub fn parseFile(allocator: std.mem.Allocator, path: []const u8) !HmmData {
    const content = try std.fs.cwd().readFileAlloc(allocator, path, 16 * 1024 * 1024);
    defer allocator.free(content);
    return parseString(allocator, content);
}

/// Context sections for the parser.
const Section = enum {
    top,
    top_extras,
    tracks,
    track_symbols,
    state_body,
    state_extras,
    state_transitions,
    state_emissions,
    state_emissions_probs,
};

fn countLeadingSpaces(line: []const u8) usize {
    var n: usize = 0;
    for (line) |c| {
        if (c == ' ') n += 1 else break;
    }
    return n;
}

fn trimLine(line: []const u8) []const u8 {
    return std.mem.trim(u8, line, " \t\r");
}

/// Flush current_state into result.states (if non-null) and reset.
fn flushState(allocator: std.mem.Allocator, result: *HmmData, current_state: *?StateData) !void {
    if (current_state.*) |cs| {
        try result.states.append(allocator, cs);
        current_state.* = null;
    }
}

/// Flush current_track into result.tracks (if non-null) and reset.
fn flushTrack(allocator: std.mem.Allocator, result: *HmmData, current_track: *?TrackData) !void {
    if (current_track.*) |ct| {
        try result.tracks.append(allocator, ct);
        current_track.* = null;
    }
}

/// Parse the HMM YAML string.
pub fn parseString(allocator: std.mem.Allocator, content: []const u8) !HmmData {
    var result = HmmData{
        .name = try allocator.dupe(u8, ""),
        .overall_prob = 0.0,
        .original_overall_mute_freq = 0.0,
        .ambiguous_char = try allocator.dupe(u8, ""),
        .tracks = .{},
        .states = .{},
    };
    errdefer result.deinit(allocator);

    var section: Section = .top;
    var current_track: ?TrackData = null;
    var current_state: ?StateData = null;

    var line_iter = std.mem.splitScalar(u8, content, '\n');
    while (line_iter.next()) |raw| {
        const line = if (raw.len > 0 and raw[raw.len - 1] == '\r') raw[0 .. raw.len - 1] else raw;

        // Detect start of a new state (list item with Python object tag)
        if (std.mem.startsWith(u8, line, "- !!python/object:")) {
            try flushState(allocator, &result, &current_state);
            current_state = StateData{
                .name = try allocator.dupe(u8, ""),
                .germline_nuc = null,
                .ambiguous_char = null,
                .ambiguous_emission_prob = null,
                .transitions = std.StringHashMap(f64).init(allocator),
                .emissions = null,
            };
            section = .state_body;
            continue;
        }

        // Skip top-level Python object tag and blank lines
        if (std.mem.startsWith(u8, line, "!!python/object:")) continue;
        if (line.len == 0) continue;

        const indent = countLeadingSpaces(line);
        const t = trimLine(line);

        // Universal top-level section transitions (can arrive from any section).
        // The `tracks:` key appears *after* all states in real partis YAML files.
        if (indent == 0 and std.mem.eql(u8, t, "tracks:")) {
            try flushState(allocator, &result, &current_state);
            section = .tracks;
            continue;
        }
        if (indent == 0 and std.mem.eql(u8, t, "states:")) {
            try flushTrack(allocator, &result, &current_track);
            section = .top; // states are started by "- !!python/object:" lines
            continue;
        }
        if (indent == 0 and std.mem.startsWith(u8, t, "name: ")) {
            allocator.free(result.name);
            result.name = try allocator.dupe(u8, t["name: ".len..]);
            continue;
        }

        switch (section) {
            .top => {
                if (std.mem.startsWith(u8, t, "name: ")) {
                    allocator.free(result.name);
                    result.name = try allocator.dupe(u8, t["name: ".len..]);
                } else if (std.mem.eql(u8, t, "extras:")) {
                    section = .top_extras;
                } else if (std.mem.eql(u8, t, "tracks:")) {
                    section = .tracks;
                }
                // "states:" is handled by the - !!python/object: detection
            },
            .top_extras => {
                if (indent == 0) {
                    section = .top;
                    // re-process
                    if (std.mem.startsWith(u8, t, "name: ")) {
                        allocator.free(result.name);
                        result.name = try allocator.dupe(u8, t["name: ".len..]);
                    } else if (std.mem.eql(u8, t, "tracks:")) {
                        section = .tracks;
                    }
                } else if (std.mem.startsWith(u8, t, "gene_prob: ")) {
                    result.overall_prob = try std.fmt.parseFloat(f64, t["gene_prob: ".len..]);
                } else if (std.mem.startsWith(u8, t, "overall_mute_freq: ")) {
                    result.original_overall_mute_freq = try std.fmt.parseFloat(f64, t["overall_mute_freq: ".len..]);
                } else if (std.mem.startsWith(u8, t, "ambiguous_char: ")) {
                    allocator.free(result.ambiguous_char);
                    result.ambiguous_char = try allocator.dupe(u8, t["ambiguous_char: ".len..]);
                }
            },
            .tracks => {
                if (indent == 0) {
                    try flushTrack(allocator, &result, &current_track);
                    section = .top;
                    if (std.mem.startsWith(u8, t, "name: ")) {
                        allocator.free(result.name);
                        result.name = try allocator.dupe(u8, t["name: ".len..]);
                    }
                } else if (indent == 2 and std.mem.endsWith(u8, t, ":") and !std.mem.startsWith(u8, t, "- ")) {
                    // New track: "  nukes:"
                    try flushTrack(allocator, &result, &current_track);
                    current_track = TrackData{
                        .name = try allocator.dupe(u8, t[0 .. t.len - 1]),
                        .symbols = .{},
                    };
                } else if (indent >= 2 and std.mem.startsWith(u8, t, "- ")) {
                    // Track symbol: "  - A"
                    if (current_track) |*ct| {
                        try ct.symbols.append(allocator, try allocator.dupe(u8, t[2..]));
                    }
                }
            },
            .track_symbols => unreachable,
            .state_body => {
                if (current_state == null) continue;
                if (std.mem.startsWith(u8, t, "name: ")) {
                    allocator.free(current_state.?.name);
                    current_state.?.name = try allocator.dupe(u8, t["name: ".len..]);
                } else if (std.mem.eql(u8, t, "extras:")) {
                    section = .state_extras;
                } else if (std.mem.eql(u8, t, "extras: {}")) {
                    // empty extras — do nothing
                } else if (std.mem.eql(u8, t, "transitions:")) {
                    section = .state_transitions;
                } else if (std.mem.eql(u8, t, "emissions: null")) {
                    // init state or null emissions — stay as null
                } else if (std.mem.eql(u8, t, "emissions:")) {
                    section = .state_emissions;
                    if (current_state.?.emissions == null) {
                        current_state.?.emissions = EmissionData{
                            .probs = std.StringHashMap(f64).init(allocator),
                            .track_name = try allocator.dupe(u8, ""),
                        };
                    }
                }
            },
            .state_extras => {
                if (current_state == null) continue;
                if (indent <= 2) {
                    section = .state_body;
                    // re-process at state_body level
                    if (std.mem.startsWith(u8, t, "name: ")) {
                        allocator.free(current_state.?.name);
                        current_state.?.name = try allocator.dupe(u8, t["name: ".len..]);
                    } else if (std.mem.eql(u8, t, "transitions:")) {
                        section = .state_transitions;
                    } else if (std.mem.eql(u8, t, "emissions: null")) {
                        // nothing
                    } else if (std.mem.eql(u8, t, "emissions:")) {
                        section = .state_emissions;
                        if (current_state.?.emissions == null) {
                            current_state.?.emissions = EmissionData{
                                .probs = std.StringHashMap(f64).init(allocator),
                                .track_name = try allocator.dupe(u8, ""),
                            };
                        }
                    }
                } else if (std.mem.startsWith(u8, t, "germline: ")) {
                    if (current_state.?.germline_nuc) |g| allocator.free(g);
                    current_state.?.germline_nuc = try allocator.dupe(u8, t["germline: ".len..]);
                } else if (std.mem.startsWith(u8, t, "ambiguous_char: ")) {
                    if (current_state.?.ambiguous_char) |a| allocator.free(a);
                    current_state.?.ambiguous_char = try allocator.dupe(u8, t["ambiguous_char: ".len..]);
                } else if (std.mem.startsWith(u8, t, "ambiguous_emission_prob: ")) {
                    current_state.?.ambiguous_emission_prob = try std.fmt.parseFloat(f64, t["ambiguous_emission_prob: ".len..]);
                }
            },
            .state_transitions => {
                if (current_state == null) continue;
                if (indent <= 2) {
                    section = .state_body;
                    if (std.mem.startsWith(u8, t, "name: ")) {
                        allocator.free(current_state.?.name);
                        current_state.?.name = try allocator.dupe(u8, t["name: ".len..]);
                    } else if (std.mem.eql(u8, t, "extras:")) {
                        section = .state_extras;
                    } else if (std.mem.eql(u8, t, "extras: {}")) {
                        // nothing
                    } else if (std.mem.eql(u8, t, "emissions: null")) {
                        // nothing
                    } else if (std.mem.eql(u8, t, "emissions:")) {
                        section = .state_emissions;
                        if (current_state.?.emissions == null) {
                            current_state.?.emissions = EmissionData{
                                .probs = std.StringHashMap(f64).init(allocator),
                                .track_name = try allocator.dupe(u8, ""),
                            };
                        }
                    }
                } else {
                    // "    state_name: 0.5"
                    if (std.mem.indexOf(u8, t, ": ")) |colon| {
                        const key = try allocator.dupe(u8, t[0..colon]);
                        const prob = try std.fmt.parseFloat(f64, t[colon + 2 ..]);
                        try current_state.?.transitions.put(key, prob);
                    }
                }
            },
            .state_emissions => {
                if (current_state == null) continue;
                if (indent <= 2) {
                    section = .state_body;
                    if (std.mem.startsWith(u8, t, "name: ")) {
                        allocator.free(current_state.?.name);
                        current_state.?.name = try allocator.dupe(u8, t["name: ".len..]);
                    } else if (std.mem.eql(u8, t, "extras:")) {
                        section = .state_extras;
                    } else if (std.mem.eql(u8, t, "extras: {}")) {
                        // nothing
                    } else if (std.mem.eql(u8, t, "transitions:")) {
                        section = .state_transitions;
                    }
                } else if (std.mem.startsWith(u8, t, "track: ")) {
                    if (current_state.?.emissions) |*em| {
                        allocator.free(em.track_name);
                        em.track_name = try allocator.dupe(u8, t["track: ".len..]);
                    }
                } else if (std.mem.eql(u8, t, "probs:")) {
                    section = .state_emissions_probs;
                }
            },
            .state_emissions_probs => {
                if (current_state == null) continue;
                // Probs lines have indent >= 6. Track line has indent == 4.
                if (indent <= 2) {
                    section = .state_body;
                    if (std.mem.startsWith(u8, t, "name: ")) {
                        allocator.free(current_state.?.name);
                        current_state.?.name = try allocator.dupe(u8, t["name: ".len..]);
                    } else if (std.mem.eql(u8, t, "extras:")) {
                        section = .state_extras;
                    } else if (std.mem.eql(u8, t, "extras: {}")) {
                        // nothing
                    } else if (std.mem.eql(u8, t, "transitions:")) {
                        section = .state_transitions;
                    }
                } else if (indent == 4) {
                    // Back to emissions level (e.g. "    track: nukes")
                    section = .state_emissions;
                    if (std.mem.startsWith(u8, t, "track: ")) {
                        if (current_state.?.emissions) |*em| {
                            allocator.free(em.track_name);
                            em.track_name = try allocator.dupe(u8, t["track: ".len..]);
                        }
                    }
                } else {
                    // "      A: 0.25"
                    if (std.mem.indexOf(u8, t, ": ")) |colon| {
                        const key = try allocator.dupe(u8, t[0..colon]);
                        const prob = try std.fmt.parseFloat(f64, t[colon + 2 ..]);
                        if (current_state.?.emissions) |*em| {
                            try em.probs.put(key, prob);
                        } else {
                            allocator.free(key);
                        }
                    }
                }
            },
        }
    }

    // Flush final state and track
    if (current_state) |cs| try result.states.append(allocator, cs);
    if (current_track) |ct| try result.tracks.append(allocator, ct);

    return result;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

test "hmm_yaml: parse simple HMM string" {
    const allocator = std.testing.allocator;

    // Use the exact same indentation as real partis YAML files.
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

    var data = try parseString(allocator, yaml);
    defer data.deinit(allocator);

    try std.testing.expectEqualStrings("test_gene", data.name);
    try std.testing.expectApproxEqAbs(0.5, data.overall_prob, 1e-9);
    try std.testing.expectApproxEqAbs(0.1, data.original_overall_mute_freq, 1e-9);
    try std.testing.expectEqual(@as(usize, 1), data.tracks.items.len);
    try std.testing.expectEqualStrings("nukes", data.tracks.items[0].name);
    try std.testing.expectEqual(@as(usize, 4), data.tracks.items[0].symbols.items.len);
    try std.testing.expectEqual(@as(usize, 2), data.states.items.len);

    // Check init state
    const init_state = &data.states.items[0];
    try std.testing.expectEqualStrings("init", init_state.name);
    try std.testing.expect(init_state.emissions == null);
    try std.testing.expectEqual(@as(usize, 1), init_state.transitions.count());

    // Check state0
    const s0 = &data.states.items[1];
    try std.testing.expectEqualStrings("state0", s0.name);
    try std.testing.expectEqualStrings("A", s0.germline_nuc.?);
    try std.testing.expect(s0.emissions != null);
    try std.testing.expectEqual(@as(usize, 4), s0.emissions.?.probs.count());
    try std.testing.expectApproxEqAbs(0.25, s0.emissions.?.probs.get("A").?, 1e-9);
}
