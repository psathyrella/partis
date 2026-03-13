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
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    try bcrham.run(allocator, std.os.argv[1..]);
}
