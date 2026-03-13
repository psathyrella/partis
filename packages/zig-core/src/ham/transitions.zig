/// ham/transitions.zig — Zig port of ham/src/transitions.cc + ham/include/transitions.h
///
/// Transition probability data structure used by State.
///
/// C++ source: packages/ham/src/transitions.cc, packages/ham/include/transitions.h
/// C++ author: psathyrella/ham

const std = @import("std");

/// Forward declaration: State is defined in state.zig but Transition needs a pointer to it.
/// We use `anyopaque` to break the circular dependency — callers cast as needed.
const StateOpaque = anyopaque;

/// Transition from one HMM state to another with an associated log-probability.
/// Corresponds to C++ `ham::Transition`.
pub const Transition = struct {
    /// Name of the destination state (owned string).
    to_state_name: []u8,
    /// Log transition probability (log of the raw probability passed to constructor).
    log_prob: f64,
    /// Pointer to the destination State (set by Model::FinalizeState; not owned).
    /// Cast to `*State` by callers that import state.zig.
    to_state_ptr: ?*StateOpaque,

    /// Create a Transition from a destination state name and a raw probability.
    /// Stores log(prob) to match C++ constructor behavior.
    /// Caller owns the returned Transition; call deinit when done.
    /// Corresponds to C++ `Transition(string to_state, double prob)`.
    pub fn init(allocator: std.mem.Allocator, to_state: []const u8, prob: f64) !Transition {
        return Transition{
            .to_state_name = try allocator.dupe(u8, to_state),
            .log_prob = @log(prob),
            .to_state_ptr = null,
        };
    }

    pub fn deinit(self: *Transition, allocator: std.mem.Allocator) void {
        allocator.free(self.to_state_name);
    }

    /// Set the destination state pointer (called by Model during finalization).
    /// Corresponds to C++ `void set_to_state(State*)`.
    pub fn setToState(self: *Transition, st: *StateOpaque) void {
        self.to_state_ptr = st;
    }

    /// Print the transition in C++ format: `%30s%14.3e\n`.
    /// Corresponds to C++ `Transition::Print()`.
    pub fn print(self: *const Transition) void {
        std.debug.print("{s:>30}{e:14.3}\n", .{ self.to_state_name, @exp(self.log_prob) });
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Transition init and log_prob" {
    const allocator = std.testing.allocator;
    var t = try Transition.init(allocator, "state_A", 0.5);
    defer t.deinit(allocator);
    try std.testing.expectEqualStrings("state_A", t.to_state_name);
    try std.testing.expectApproxEqAbs(@log(0.5), t.log_prob, 1e-10);
}

test "Transition: prob 1.0 gives log_prob 0.0" {
    const allocator = std.testing.allocator;
    var t = try Transition.init(allocator, "end", 1.0);
    defer t.deinit(allocator);
    try std.testing.expectApproxEqAbs(0.0, t.log_prob, 1e-10);
}
