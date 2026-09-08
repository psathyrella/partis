/// ham/cluster_path.zig — Zig port of ham/src/clusterpath.cc + ham/include/clusterpath.h
///
/// Sequence of gradually coalescing partitions with associated log-probabilities.
/// Used by the Glomerator for hierarchical clustering.
///
/// C++ source: packages/ham/src/clusterpath.cc, packages/ham/include/clusterpath.h
/// C++ author: psathyrella/ham

const std = @import("std");
const math = std.math;
const mathutils = @import("mathutils.zig");

/// A partition: a set of sequence IDs (strings) forming a cluster.
/// Corresponds to C++ `typedef set<string> Partition`.
/// We represent as a sorted ArrayList of owned strings for simplicity.
pub const Partition = std.ArrayListUnmanaged([]u8);

/// Sequence of partitions with associated log-probabilities, tracking the best.
/// Corresponds to C++ `ham::ClusterPath`.
pub const ClusterPath = struct {
    /// Sequence of partitions (front = oldest, back = most recent).
    partitions: std.ArrayListUnmanaged(Partition),
    /// Log-probabilities for each partition (parallel to partitions).
    logprobs: std.ArrayListUnmanaged(f64),
    /// Whether this path is finished (no more clustering steps).
    finished: bool,
    /// Index of the partition in the batch that gave rise to this path.
    initial_path_index: i32,
    /// Log-prob of the best partition seen so far.
    max_log_prob_of_partition: f64,
    /// Index of the best partition.
    i_best: usize,

    /// Create an empty ClusterPath with no partitions.
    pub fn init() ClusterPath {
        return ClusterPath{
            .partitions = .{},
            .logprobs = .{},
            .finished = false,
            .initial_path_index = 0,
            .max_log_prob_of_partition = mathutils.NEG_INF,
            .i_best = 0,
        };
    }

    /// Create a ClusterPath with an initial partition.
    /// Takes ownership of `initial_partition`.
    /// Corresponds to C++ `ClusterPath(Partition initial_partition, double initial_logprob)`.
    pub fn initWithPartition(
        allocator: std.mem.Allocator,
        initial_partition: Partition,
        initial_logprob: f64,
    ) !ClusterPath {
        var cp = ClusterPath.init();
        try cp.partitions.append(allocator, initial_partition);
        try cp.logprobs.append(allocator, initial_logprob);
        return cp;
    }

    pub fn deinit(self: *ClusterPath, allocator: std.mem.Allocator) void {
        for (self.partitions.items) |*p| {
            for (p.items) |s| allocator.free(s);
            p.deinit(allocator);
        }
        self.partitions.deinit(allocator);
        self.logprobs.deinit(allocator);
    }

    /// Update the log-probability for partition at index `il` and update best tracking.
    /// Corresponds to C++ `ClusterPath::set_logprob(size_t il, double logprob)`.
    pub fn setLogprob(self: *ClusterPath, il: usize, logprob: f64) void {
        self.logprobs.items[il] = logprob;
        if (self.max_log_prob_of_partition == mathutils.NEG_INF or logprob > self.max_log_prob_of_partition) {
            self.max_log_prob_of_partition = logprob;
            self.i_best = il;
        }
    }

    /// Append a new partition.
    /// Takes ownership of `partition`.
    /// Corresponds to C++ `ClusterPath::AddPartition(Partition, double, size_t n_max_partitions)`.
    pub fn addPartition(
        self: *ClusterPath,
        allocator: std.mem.Allocator,
        partition: Partition,
        logprob: f64,
        n_max_partitions: usize,
    ) !void {
        try self.partitions.append(allocator, partition);
        try self.logprobs.append(allocator, logprob);

        if (self.max_log_prob_of_partition == mathutils.NEG_INF or logprob > self.max_log_prob_of_partition) {
            self.max_log_prob_of_partition = logprob;
            self.i_best = self.logprobs.items.len - 1;
        }

        // Trim oldest entries if over the limit (n_max_partitions == 0 means unlimited)
        if (n_max_partitions > 0 and self.partitions.items.len > n_max_partitions) {
            var old = self.partitions.orderedRemove(0);
            for (old.items) |s| allocator.free(s);
            old.deinit(allocator);
            _ = self.logprobs.orderedRemove(0);
        }
    }

    /// Return the most recent (current) partition.
    /// Corresponds to C++ `ClusterPath::CurrentPartition()`.
    pub fn currentPartition(self: *const ClusterPath) *const Partition {
        return &self.partitions.items[self.partitions.items.len - 1];
    }

    /// Return the log-prob of the most recent partition.
    /// Corresponds to C++ `ClusterPath::CurrentLogProb()`.
    pub fn currentLogProb(self: *const ClusterPath) f64 {
        return self.logprobs.items[self.logprobs.items.len - 1];
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "ClusterPath: init empty" {
    const cp = ClusterPath.init();
    try std.testing.expectEqual(false, cp.finished);
    try std.testing.expectEqual(@as(usize, 0), cp.partitions.items.len);
}

test "ClusterPath: addPartition and best tracking" {
    const allocator = std.testing.allocator;
    var cp = ClusterPath.init();
    defer cp.deinit(allocator);

    // Add first partition with logprob -5.0
    var p1: Partition = .{};
    try p1.append(allocator, try allocator.dupe(u8, "seq1"));
    try cp.addPartition(allocator, p1, -5.0, 0);

    try std.testing.expectEqual(@as(f64, -5.0), cp.currentLogProb());
    try std.testing.expectEqual(@as(f64, -5.0), cp.max_log_prob_of_partition);
    try std.testing.expectEqual(@as(usize, 0), cp.i_best);

    // Add second partition with better logprob -2.0
    var p2: Partition = .{};
    try p2.append(allocator, try allocator.dupe(u8, "seq1"));
    try p2.append(allocator, try allocator.dupe(u8, "seq2"));
    try cp.addPartition(allocator, p2, -2.0, 0);

    try std.testing.expectEqual(@as(f64, -2.0), cp.currentLogProb());
    try std.testing.expectEqual(@as(f64, -2.0), cp.max_log_prob_of_partition);
    try std.testing.expectEqual(@as(usize, 1), cp.i_best);
}

test "ClusterPath: setLogprob updates best" {
    const allocator = std.testing.allocator;
    var cp = ClusterPath.init();
    defer cp.deinit(allocator);

    const p1: Partition = .{};
    try cp.addPartition(allocator, p1, -10.0, 0);

    cp.setLogprob(0, -3.0);
    try std.testing.expectEqual(@as(f64, -3.0), cp.max_log_prob_of_partition);
}
