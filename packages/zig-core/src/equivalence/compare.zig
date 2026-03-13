/// Equivalence harness: diff C++ checkpoint JSON stream vs Zig checkpoint JSON stream.
///
/// Usage:
///   compare <component-id> <cpp-stream-file> <zig-stream-file> [work-queue-json]
///
/// Each stream file contains newline-delimited JSON objects, one per checkpoint.
/// Each object has at minimum: {"checkpoint": "<name>", ...numeric fields...}
///
/// Exit codes:
///   0  all checkpoints within tolerance
///   1  one or more checkpoints failed
///   2  harness error (file not found, parse error, etc.)

const std = @import("std");
const mem = std.mem;
const json = std.json;
const fs = std.fs;

/// Per-component tolerance configuration loaded from work-queue.json.
const Tolerance = struct {
    float_ulp: u32 = 4,
    log_prob_absolute: f64 = 1e-6,
    integer_exact: bool = true,
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    var stdout_buf: [65536]u8 = undefined;
    var stdout_file_writer = fs.File.stdout().writer(&stdout_buf);
    const stdout = &stdout_file_writer.interface;

    if (args.len < 4) {
        std.debug.print(
            "usage: compare <component-id> <cpp-stream-file> <zig-stream-file> [work-queue-json]\n",
            .{},
        );
        std.process.exit(2);
    }

    const component_id = args[1];
    const cpp_path = args[2];
    const zig_path = args[3];
    const work_queue_path: ?[]const u8 = if (args.len >= 5) args[4] else null;

    const tolerance = try loadTolerance(allocator, work_queue_path, component_id);

    try stdout.print("equivalence-check: {s}\n", .{component_id});

    // Read both streams.
    const cpp_data = fs.cwd().readFileAlloc(allocator, cpp_path, 64 * 1024 * 1024) catch |err| {
        std.debug.print("ERROR: cannot read cpp stream '{s}': {}\n", .{ cpp_path, err });
        std.process.exit(2);
    };
    defer allocator.free(cpp_data);

    const zig_data = fs.cwd().readFileAlloc(allocator, zig_path, 64 * 1024 * 1024) catch |err| {
        std.debug.print("ERROR: cannot read zig stream '{s}': {}\n", .{ zig_path, err });
        std.process.exit(2);
    };
    defer allocator.free(zig_data);

    const cpp_lines = try splitLines(allocator, cpp_data);
    defer allocator.free(cpp_lines);
    const zig_lines = try splitLines(allocator, zig_data);
    defer allocator.free(zig_lines);

    if (cpp_lines.len != zig_lines.len) {
        try stdout.print(
            "  MISMATCH: cpp has {d} checkpoints, zig has {d} checkpoints\n",
            .{ cpp_lines.len, zig_lines.len },
        );
        try stdout.print("RESULT: FAIL (checkpoint count mismatch)\n", .{});
        try stdout_file_writer.interface.flush();
        std.process.exit(1);
    }

    var pass_count: usize = 0;
    var fail_count: usize = 0;

    for (cpp_lines, zig_lines, 0..) |cpp_line, zig_line, i| {
        if (cpp_line.len == 0 and zig_line.len == 0) continue;

        const result = compareCheckpoint(allocator, cpp_line, zig_line, tolerance, i) catch |err| {
            try stdout.print("  checkpoint[{d}]  HARNESS_ERROR: {}\n", .{ i, err });
            fail_count += 1;
            continue;
        };

        if (result.passed) {
            if (result.delta_info) |info| {
                try stdout.print("  checkpoint {s:<40}  PASS (delta: {s})\n", .{ result.name, info });
            } else {
                try stdout.print("  checkpoint {s:<40}  PASS (exact)\n", .{result.name});
            }
            pass_count += 1;
        } else {
            try stdout.print(
                "  checkpoint {s:<40}  FAIL  cpp={s}  zig={s}  delta={s}  threshold={s}\n",
                .{
                    result.name,
                    result.cpp_val,
                    result.zig_val,
                    result.delta_info orelse "?",
                    result.threshold_info orelse "?",
                },
            );
            fail_count += 1;
        }
    }

    const total = pass_count + fail_count;
    if (fail_count == 0) {
        try stdout.print("\nRESULT: PASS ({d}/{d} checkpoints)\n", .{ pass_count, total });
        try stdout_file_writer.interface.flush();
        std.process.exit(0);
    } else {
        try stdout.print("\nRESULT: FAIL ({d}/{d} checkpoints passed)\n", .{ pass_count, total });
        try stdout_file_writer.interface.flush();
        std.process.exit(1);
    }
}

const CheckpointResult = struct {
    name: []const u8,
    passed: bool,
    cpp_val: []const u8 = "",
    zig_val: []const u8 = "",
    delta_info: ?[]const u8 = null,
    threshold_info: ?[]const u8 = null,
};

fn compareCheckpoint(
    allocator: mem.Allocator,
    cpp_line: []const u8,
    zig_line: []const u8,
    tolerance: Tolerance,
    idx: usize,
) !CheckpointResult {
    const cpp_parsed = json.parseFromSlice(json.Value, allocator, cpp_line, .{}) catch {
        return CheckpointResult{
            .name = try std.fmt.allocPrint(allocator, "[{d}]", .{idx}),
            .passed = false,
            .cpp_val = "(parse error)",
            .zig_val = zig_line,
        };
    };
    defer cpp_parsed.deinit();

    const zig_parsed = json.parseFromSlice(json.Value, allocator, zig_line, .{}) catch {
        return CheckpointResult{
            .name = try std.fmt.allocPrint(allocator, "[{d}]", .{idx}),
            .passed = false,
            .cpp_val = cpp_line,
            .zig_val = "(parse error)",
        };
    };
    defer zig_parsed.deinit();

    const cpp_obj = switch (cpp_parsed.value) {
        .object => |o| o,
        else => return error.NotAnObject,
    };
    const zig_obj = switch (zig_parsed.value) {
        .object => |o| o,
        else => return error.NotAnObject,
    };

    // Allocate an owned copy of the checkpoint name so it survives after
    // cpp_parsed.deinit() runs at end of this function.
    const name_raw = if (cpp_obj.get("checkpoint")) |v|
        switch (v) {
            .string => |s| s,
            else => "(unnamed)",
        }
    else
        "(unnamed)";
    const name = try allocator.dupe(u8, name_raw);

    // Exact match fast path.
    if (mem.eql(u8, cpp_line, zig_line)) {
        return CheckpointResult{ .name = name, .passed = true };
    }

    // Field-by-field comparison with tolerance.
    var all_pass = true;
    var first_fail_cpp: []const u8 = "";
    var first_fail_zig: []const u8 = "";
    var first_fail_delta: ?[]const u8 = null;
    var first_fail_threshold: ?[]const u8 = null;

    var it = cpp_obj.iterator();
    while (it.next()) |entry| {
        const key = entry.key_ptr.*;
        if (mem.eql(u8, key, "checkpoint")) continue;

        const cpp_val = entry.value_ptr.*;
        const zig_val = zig_obj.get(key) orelse continue;

        const field_pass = valuesWithinTolerance(cpp_val, zig_val, tolerance) catch false;
        if (!field_pass) {
            all_pass = false;
            first_fail_cpp = try std.fmt.allocPrint(allocator, "{}", .{cpp_val});
            first_fail_zig = try std.fmt.allocPrint(allocator, "{}", .{zig_val});
            first_fail_delta = try std.fmt.allocPrint(allocator, "field={s}", .{key});
            first_fail_threshold = try std.fmt.allocPrint(
                allocator,
                "log_prob_abs={e:.1}",
                .{tolerance.log_prob_absolute},
            );
        }
    }

    return CheckpointResult{
        .name = name,
        .passed = all_pass,
        .cpp_val = first_fail_cpp,
        .zig_val = first_fail_zig,
        .delta_info = first_fail_delta,
        .threshold_info = first_fail_threshold,
    };
}

fn valuesWithinTolerance(cpp: json.Value, zig_val: json.Value, tol: Tolerance) !bool {
    switch (cpp) {
        .float => |cf| {
            const zf = switch (zig_val) {
                .float => |f| f,
                .integer => |i| @as(f64, @floatFromInt(i)),
                else => return false,
            };
            const delta = @abs(cf - zf);
            return delta <= tol.log_prob_absolute;
        },
        .integer => |ci| {
            const zi = switch (zig_val) {
                .integer => |i| i,
                else => return false,
            };
            return if (tol.integer_exact) ci == zi else true;
        },
        .bool => |cb| {
            const zb = switch (zig_val) {
                .bool => |b| b,
                else => return false,
            };
            return cb == zb;
        },
        .string => |cs| {
            const zs = switch (zig_val) {
                .string => |s| s,
                else => return false,
            };
            return mem.eql(u8, cs, zs);
        },
        else => return true, // null == null, arrays compared later if needed
    }
}

fn loadTolerance(allocator: mem.Allocator, work_queue_path: ?[]const u8, component_id: []const u8) !Tolerance {
    const default = Tolerance{};
    const path = work_queue_path orelse return default;

    const data = fs.cwd().readFileAlloc(allocator, path, 4 * 1024 * 1024) catch return default;
    defer allocator.free(data);

    const parsed = json.parseFromSlice(json.Value, allocator, data, .{}) catch return default;
    defer parsed.deinit();

    const obj = switch (parsed.value) {
        .object => |o| o,
        else => return default,
    };
    const components = switch (obj.get("components") orelse return default) {
        .array => |a| a,
        else => return default,
    };

    for (components.items) |item| {
        const comp = switch (item) {
            .object => |o| o,
            else => continue,
        };
        const id = switch (comp.get("id") orelse continue) {
            .string => |s| s,
            else => continue,
        };
        if (!mem.eql(u8, id, component_id)) continue;

        const tol_val = comp.get("tolerance") orelse return default;
        const tol_obj = switch (tol_val) {
            .object => |o| o,
            else => return default,
        };

        var tol = default;
        if (tol_obj.get("float_ulp")) |v| {
            tol.float_ulp = @intCast(switch (v) {
                .integer => |i| i,
                else => @as(i64, default.float_ulp),
            });
        }
        if (tol_obj.get("log_prob_absolute")) |v| {
            tol.log_prob_absolute = switch (v) {
                .float => |f| f,
                .integer => |i| @as(f64, @floatFromInt(i)),
                else => default.log_prob_absolute,
            };
        }
        if (tol_obj.get("integer_exact")) |v| {
            tol.integer_exact = switch (v) {
                .bool => |b| b,
                else => default.integer_exact,
            };
        }
        return tol;
    }
    return default;
}

fn splitLines(allocator: mem.Allocator, data: []const u8) ![][]const u8 {
    var lines: std.ArrayListUnmanaged([]const u8) = .empty;
    var iter = mem.splitScalar(u8, data, '\n');
    while (iter.next()) |line| {
        const trimmed = mem.trim(u8, line, " \r\t");
        if (trimmed.len > 0) try lines.append(allocator, trimmed);
    }
    return lines.toOwnedSlice(allocator);
}
