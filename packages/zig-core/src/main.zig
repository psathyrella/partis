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
const builtin = @import("builtin");
const bcrham = @import("ham/bcrham.zig");

pub fn main() !void {
    const is_safe = builtin.mode == .Debug or builtin.mode == .ReleaseSafe;
    const allocator = if (is_safe) blk: {
        // In debug and release-safe builds, use GPA for leak detection and safety checks.
        break :blk gpa_instance.allocator();
    } else
        // In release-fast builds, use the C allocator — GPA's per-allocation tracking
        // metadata adds overhead on large workloads.
        std.heap.c_allocator;
    defer if (is_safe) {
        _ = gpa_instance.deinit();
    };
    try bcrham.run(allocator, std.os.argv[1..]);
}

var gpa_instance = std.heap.GeneralPurposeAllocator(.{}){};

// In test mode, force the analyzer to walk all decls reachable from the
// executable so that tests in transitively-imported files get discovered by
// `zig build test`. Without this, only tests in main.zig itself would run —
// and main.zig has no tests. (`refAllDeclsRecursive` is a no-op outside tests.)
//
// The eval-branch quota is bumped because `refAllDeclsRecursive` walks the
// full decl tree of `bcrham` (and everything it imports — args, model, state,
// dp_handler, glomerator, bcr_utils, etc.); the default 1000 is exceeded.
// 20000 was the smallest round number that compiled with current sources;
// raise it (don't lower) if a future module addition trips the quota again.
test {
    @setEvalBranchQuota(20000);
    std.testing.refAllDeclsRecursive(bcrham);
}
