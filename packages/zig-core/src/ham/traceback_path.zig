/// ham/traceback_path.zig — Zig port of ham/src/tracebackpath.cc + ham/include/tracebackpath.h
///
/// Viterbi traceback path: a sequence of state indices with an associated log-score.
///
/// C++ source: packages/ham/src/tracebackpath.cc, packages/ham/include/tracebackpath.h
/// C++ author: psathyrella/ham

const std = @import("std");
const Model = @import("model.zig").Model;

/// Viterbi traceback path: ordered sequence of state indices (in reverse order of traversal).
/// Corresponds to C++ `ham::TracebackPath`.
pub const TracebackPath = struct {
    /// Pointer to the HMM model (not owned; must outlive this path).
    hmm: ?*const Model,
    /// The path as a sequence of state indices.
    /// Stored in Viterbi fill order (index 0 = last state, end = first after init).
    /// Corresponds to C++ `vector<int> path_`.
    path: std.ArrayListUnmanaged(i32),
    /// Log-score of this path.
    score: f64,
    /// If true, use 1-char abbreviations for printing.
    abbreviate_flag: bool,

    /// Create a TracebackPath for a given model.
    /// Corresponds to C++ `TracebackPath(Model*)`.
    pub fn initWithModel(model: *const Model) TracebackPath {
        return TracebackPath{
            .hmm = model,
            .path = .{},
            .score = 0.0,
            .abbreviate_flag = false,
        };
    }

    /// Create an empty TracebackPath (no model).
    /// Corresponds to C++ `TracebackPath()`.
    pub fn init() TracebackPath {
        return TracebackPath{
            .hmm = null,
            .path = .{},
            .score = 0.0,
            .abbreviate_flag = false,
        };
    }

    pub fn deinit(self: *TracebackPath, allocator: std.mem.Allocator) void {
        self.path.deinit(allocator);
    }

    /// Append a state index to the path.
    /// Corresponds to C++ `void push_back(int state)`.
    pub fn pushBack(self: *TracebackPath, allocator: std.mem.Allocator, state: i32) !void {
        try self.path.append(allocator, state);
    }

    /// Clear all state indices.
    pub fn clear(self: *TracebackPath) void {
        self.path.clearRetainingCapacity();
    }

    pub fn size(self: *const TracebackPath) usize {
        return self.path.items.len;
    }

    pub fn setAbbreviate(self: *TracebackPath, abb: bool) void {
        self.abbreviate_flag = abb;
    }

    pub fn setModel(self: *TracebackPath, model: *const Model) void {
        self.hmm = model;
    }

    pub fn setScore(self: *TracebackPath, s: f64) void {
        self.score = s;
    }

    /// Access state index at position `val`.
    /// Corresponds to C++ `operator[](size_t val)`.
    pub fn at(self: *const TracebackPath, val: usize) i32 {
        return self.path.items[val];
    }

    /// Return the sequence of state names, in path order (reversed from storage order).
    /// Caller owns the returned ArrayList of owned strings.
    /// Corresponds to C++ `TracebackPath::name_vector()`.
    pub fn nameVector(self: *const TracebackPath, allocator: std.mem.Allocator) !std.ArrayListUnmanaged([]u8) {
        var result: std.ArrayListUnmanaged([]u8) = .{};
        const model = self.hmm orelse return result;
        const n = self.path.items.len;
        var k: usize = n;
        while (k > 0) {
            k -= 1;
            const idx: usize = @intCast(self.path.items[k]);
            const st = model.stateByIndex(idx);
            try result.append(allocator, try allocator.dupe(u8, st.name));
        }
        return result;
    }

    /// Print this path to stderr (in path order, reversed from storage order).
    /// Corresponds to C++ `operator<<(ostream&, const TracebackPath&)`.
    pub fn print(self: *const TracebackPath) void {
        const model = self.hmm orelse return;
        const n = self.path.items.len;
        var k: usize = n;
        while (k > 0) {
            k -= 1;
            const idx: usize = @intCast(self.path.items[k]);
            const st = model.stateByIndex(idx);
            if (self.abbreviate_flag) {
                std.debug.print("{c} ", .{st.abbreviation()});
            } else {
                std.debug.print("{s} ", .{st.name});
            }
        }
        std.debug.print("\n", .{});
    }

    /// Equality comparison (compares path arrays).
    pub fn eql(self: *const TracebackPath, other: *const TracebackPath) bool {
        if (self.path.items.len != other.path.items.len) return false;
        for (self.path.items, other.path.items) |a, b| {
            if (a != b) return false;
        }
        return true;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "TracebackPath: init, pushBack, at, size" {
    const allocator = std.testing.allocator;

    var path = TracebackPath.init();
    defer path.deinit(allocator);

    try path.pushBack(allocator, 5);
    try path.pushBack(allocator, 3);
    try path.pushBack(allocator, 7);

    try std.testing.expectEqual(@as(usize, 3), path.size());
    try std.testing.expectEqual(@as(i32, 5), path.at(0));
    try std.testing.expectEqual(@as(i32, 3), path.at(1));
    try std.testing.expectEqual(@as(i32, 7), path.at(2));
}

test "TracebackPath: score and clear" {
    const allocator = std.testing.allocator;

    var path = TracebackPath.init();
    defer path.deinit(allocator);

    path.setScore(-42.5);
    try std.testing.expectApproxEqAbs(-42.5, path.score, 1e-9);

    try path.pushBack(allocator, 1);
    try path.pushBack(allocator, 2);
    path.clear();
    try std.testing.expectEqual(@as(usize, 0), path.size());
}

test "TracebackPath: eql" {
    const allocator = std.testing.allocator;

    var p1 = TracebackPath.init();
    defer p1.deinit(allocator);
    var p2 = TracebackPath.init();
    defer p2.deinit(allocator);

    try p1.pushBack(allocator, 1);
    try p1.pushBack(allocator, 2);
    try p2.pushBack(allocator, 1);
    try p2.pushBack(allocator, 2);

    try std.testing.expect(p1.eql(&p2));

    try p2.pushBack(allocator, 3);
    try std.testing.expect(!p1.eql(&p2));
}
