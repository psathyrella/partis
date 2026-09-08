/// igsw/main.zig — CLI entry point for partis-zig-igsw
///
/// Drop-in replacement for the C ig-sw binary.  Parses the same flags that
/// partis passes via get_ig_sw_cmd_str, constructs germline FASTA paths by
/// replicating GetFileName() from ig_align_main.cpp, then delegates to the
/// Zig ig_align wrapper.
///
/// C++ source: packages/ig-sw/src/ig_align/ig_align_main.cpp

const std = @import("std");
const ig_align = @import("ig_align.zig");

// ─── path construction (replicates C++ GetFileName) ──────────────────────────

const Locus = enum { IGH, IGK, IGL, TRA, TRB, TRG, TRD, DJ };

fn parseLocus(s: []const u8) ?Locus {
    if (std.mem.eql(u8, s, "IGH")) return .IGH;
    if (std.mem.eql(u8, s, "IGK")) return .IGK;
    if (std.mem.eql(u8, s, "IGL")) return .IGL;
    if (std.mem.eql(u8, s, "TRA")) return .TRA;
    if (std.mem.eql(u8, s, "TRB")) return .TRB;
    if (std.mem.eql(u8, s, "TRG")) return .TRG;
    if (std.mem.eql(u8, s, "TRD")) return .TRD;
    if (std.mem.eql(u8, s, "DJ"))  return .DJ;
    return null;
}

/// Build the list of germline FASTA paths from (vdj_dir, locus).
/// Returns paths as owned null-terminated strings; caller must free via gpa.
/// paths[0] = ref_path, paths[1..] = extra_ref_paths (D, J).
fn getFileNames(
    allocator: std.mem.Allocator,
    vdj_dir: []const u8,
    locus: Locus,
) !std.ArrayListUnmanaged([:0]u8) {
    var paths: std.ArrayListUnmanaged([:0]u8) = .{};
    errdefer {
        for (paths.items) |p| allocator.free(p);
        paths.deinit(allocator);
    }

    // gene letters and locus prefix follow C++ GetFileName() exactly
    const genes: []const u8 = switch (locus) {
        .IGH, .TRB, .TRD => "vdj",
        .IGK, .IGL, .TRA, .TRG => "vj",
        .DJ => "dj",
    };
    // DJ uses the "igh" prefix for the files
    const lower_locus: []const u8 = switch (locus) {
        .IGH, .DJ => "igh",
        .IGK => "igk",
        .IGL => "igl",
        .TRA => "tra",
        .TRB => "trb",
        .TRG => "trg",
        .TRD => "trd",
    };

    for (genes) |gene| {
        const path_str = try std.fmt.allocPrint(allocator, "{s}{s}{c}.fasta", .{ vdj_dir, lower_locus, gene });
        const path_z = try allocator.dupeZ(u8, path_str);
        allocator.free(path_str);
        try paths.append(allocator, path_z);
    }
    return paths;
}

// ─── CLI parsing ─────────────────────────────────────────────────────────────

const Args = struct {
    locus: []const u8 = "IGH",
    vdj_dir: []const u8 = "",
    match: i32 = 2,
    mismatch: i32 = 2,
    gap_o: i32 = 3,
    gap_e: i32 = 1,
    max_drop: u32 = 1000,
    min_score: i32 = 0,
    bandwidth: u32 = 150,
    n_threads: u8 = 1,
    query_fasta: []const u8 = "",
    output_sam: []const u8 = "",
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const argv = std.os.argv;
    var args: Args = .{};
    var positionals: u8 = 0;

    var i: usize = 1;
    while (i < argv.len) : (i += 1) {
        const flag = std.mem.span(argv[i]);
        // single-char flags that take a value
        if (flag.len == 2 and flag[0] == '-') {
            const ch = flag[1];
            switch (ch) {
                'l', 'd', 'm', 'u', 'o', 'p', 'e', 'b', 's', 'j' => {
                    if (i + 1 >= argv.len) {
                        std.debug.print("error: flag -{c} requires a value\n", .{ch});
                        std.process.exit(1);
                    }
                    i += 1;
                    const val = std.mem.span(argv[i]);
                    switch (ch) {
                        'l' => args.locus = val,
                        'd' => args.max_drop = try std.fmt.parseInt(u32, val, 10),
                        'm' => args.match = try std.fmt.parseInt(i32, val, 10),
                        'u' => args.mismatch = try std.fmt.parseInt(i32, val, 10),
                        'o' => args.gap_o = try std.fmt.parseInt(i32, val, 10),
                        'p' => args.vdj_dir = val,
                        'e' => args.gap_e = try std.fmt.parseInt(i32, val, 10),
                        'b' => args.bandwidth = try std.fmt.parseInt(u32, val, 10),
                        's' => args.min_score = try std.fmt.parseInt(i32, val, 10),
                        'j' => args.n_threads = try std.fmt.parseInt(u8, val, 10),
                        else => unreachable,
                    }
                },
                else => {
                    // unknown single-char flag: skip its value if present
                    if (i + 1 < argv.len) {
                        const next = std.mem.span(argv[i + 1]);
                        if (next.len == 0 or next[0] != '-') i += 1;
                    }
                },
            }
        } else if (flag.len > 2 and flag[0] == '-' and flag[1] == '-') {
            // long flag: skip it and its value
            if (i + 1 < argv.len) {
                const next = std.mem.span(argv[i + 1]);
                if (next.len == 0 or next[0] != '-') i += 1;
            }
        } else {
            // positional
            if (positionals == 0) {
                args.query_fasta = flag;
            } else if (positionals == 1) {
                args.output_sam = flag;
            }
            positionals += 1;
        }
    }

    if (positionals < 2 or args.vdj_dir.len == 0) {
        std.debug.print(
            \\Usage: partis-zig-igsw -l <LOCUS> -p <vdj_dir/> [-d <max_drop>] [-m <match>]
            \\                       [-u <mismatch>] [-o <gap_open>] [-e <gap_extend>]
            \\                       [-b <bandwidth>] [-s <min_score>] [-j <threads>]
            \\                       <query_fasta> <output_sam>
            \\
        , .{});
        std.process.exit(1);
    }

    const locus = parseLocus(args.locus) orelse {
        std.debug.print("error: unknown locus '{s}'\n", .{args.locus});
        std.process.exit(1);
    };

    // Build null-terminated path strings for the C API
    var paths = try getFileNames(allocator, args.vdj_dir, locus);
    defer {
        for (paths.items) |p| allocator.free(p);
        paths.deinit(allocator);
    }

    const ref_path: [*:0]const u8 = paths.items[0].ptr;

    // Build extra_ref_paths array (stack-allocated; max 2 entries)
    var extra_buf: [2][*:0]const u8 = undefined;
    const n_extra = paths.items.len - 1;
    for (paths.items[1..], 0..) |p, idx| extra_buf[idx] = p.ptr;
    const extra_ref_paths: []const [*:0]const u8 = extra_buf[0..n_extra];

    const qry_z = try allocator.dupeZ(u8, args.query_fasta);
    defer allocator.free(qry_z);
    const out_z = try allocator.dupeZ(u8, args.output_sam);
    defer allocator.free(out_z);

    ig_align.alignReads(
        ref_path,
        extra_ref_paths,
        qry_z.ptr,
        out_z.ptr,
        args.match,
        args.mismatch,
        args.gap_o,
        args.gap_e,
        args.max_drop,
        @intCast(args.min_score),
        args.bandwidth,
        args.n_threads,
        null,
        null,
    );
}
