/// igsw/ig_align.zig — Zig wrapper for ig_align C library
///
/// Wraps the ig_align C API via @cImport.  The C sources (ig_align.c, ksw.c,
/// kstring.c) are compiled as C objects by build.zig and linked into the
/// partis-zig-core library.  This Zig module re-exports the public entry
/// point and provides a Zig-friendly wrapper.
///
/// C source: packages/ig-sw/src/ig_align/ig_align.c  (+ ksw.c, kstring.c)
/// C author: psathyrella/ig-sw (based on lh3/ksw / Heng Li)

const std = @import("std");

// ksw.zig exports C-ABI ksw_align/ksw_global/ksw_qinit/ksw_extend, replacing
// the C ksw.c which uses SSE2 intrinsics. Importing it here compiles it into
// the same artifact so the linker resolves ig_align.c's references to ksw_*.
comptime {
    _ = @import("ksw.zig");
}

pub const c = @cImport({
    @cInclude("ig_align.h");
});

/// Wrapper around the C `ig_align_reads` function.
///
/// Aligns query FASTA sequences against one primary reference FASTA and zero
/// or more extra references (D, J genes) using Smith-Waterman, writing SAM
/// output to `output_path`.
///
/// Parameters mirror the C function exactly.  Default values used by partis:
///   match=2, mismatch=2, gap_o=3, gap_e=1, max_drop=1000, min_score=0,
///   bandwidth=1000, n_threads=1.
pub fn alignReads(
    ref_path: [*:0]const u8,
    extra_ref_paths: []const [*:0]const u8,
    qry_path: [*:0]const u8,
    output_path: [*:0]const u8,
    match: i32,
    mismatch: i32,
    gap_o: i32,
    gap_e: i32,
    max_drop: c_uint,
    min_score: c_int,
    bandwidth: c_uint,
    n_threads: u8,
    read_group: ?[*:0]const u8,
    read_group_id: ?[*:0]const u8,
) void {
    const n_extra: u8 = @intCast(extra_ref_paths.len);
    // The C function takes `const char **extra_ref_paths`.
    // We need a pointer to the first element (or null if empty).
    const extra_ptr: ?[*]const [*:0]const u8 = if (n_extra > 0) extra_ref_paths.ptr else null;
    c.ig_align_reads(
        ref_path,
        n_extra,
        if (extra_ptr) |p| @ptrCast(@constCast(p)) else null,
        qry_path,
        output_path,
        match,
        mismatch,
        gap_o,
        gap_e,
        max_drop,
        min_score,
        bandwidth,
        n_threads,
        read_group orelse null,
        read_group_id orelse null,
    );
}
