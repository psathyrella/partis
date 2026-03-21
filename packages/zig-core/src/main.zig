/// main.zig — partis-zig-core executable entry point.
///
/// Delegates all HMM computation to bcrham.run(), which accepts the same
/// CLI arguments as the C++ bcrham binary and produces compatible CSV output.
///
/// Usage:
///   partis-zig-core --hmmdir <dir> --datadir <dir> --infile <file> \
///                   --outfile <file> --algorithm <viterbi|forward> \
///                   --locus <igh|igk|igl|tra|trb|trg|trd> [optional-args...]
///
/// Unknown arguments are silently ignored (bcrham compat).

const std = @import("std");
const bcrham = @import("ham/bcrham.zig");

pub fn main() !void {
    const allocator = if (@import("builtin").mode == .Debug) blk: {
        // In debug builds, use GPA for leak detection and use-after-free checks.
        break :blk gpa_instance.allocator();
    } else
        // In release builds, use the C allocator — GPA's per-allocation tracking
        // metadata adds massive overhead (observed 11 GB vs 278 MB on a 5k-sequence
        // partition workload).
        std.heap.c_allocator;
    defer if (@import("builtin").mode == .Debug) {
        _ = gpa_instance.deinit();
    };
    try bcrham.run(allocator, std.os.argv[1..]);
}

var gpa_instance = std.heap.GeneralPurposeAllocator(.{}){};
