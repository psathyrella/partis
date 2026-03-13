/// C ABI surface for partis-zig-core.
///
/// This is the ONLY public C-facing interface. All internal Zig modules are
/// private. Python loads this shared library and calls these functions instead
/// of spawning a bcrham subprocess.
///
/// Stub bodies panic with a clear message until the corresponding component is
/// ported. The conversion-loop agent fills these in one component at a time.
///
/// Function signatures are derived from:
///   - python/partis.py  (bcrham subprocess invocation arguments)
///   - ham/src/main.cc   (bcrham entry-point argument parsing)
///
/// During first-run bootstrap (T040), the conversion-loop agent reads those
/// files and populates the full signature list here.

const std = @import("std");

// ─── internal ham modules (imported so their tests run with `zig build test`) ─
pub const ham_text = @import("ham/text.zig");
pub const ham_mathutils = @import("ham/mathutils.zig");
pub const ham_transitions = @import("ham/transitions.zig");
pub const ham_track = @import("ham/track.zig");
pub const ham_cluster_path = @import("ham/cluster_path.zig");
pub const ham_sequences = @import("ham/sequences.zig");
pub const ham_lexical_table = @import("ham/lexical_table.zig");
pub const ham_emission = @import("ham/emission.zig");
pub const ham_state = @import("ham/state.zig");
pub const ham_hmm_yaml = @import("ham/hmm_yaml.zig");
pub const ham_model = @import("ham/model.zig");
pub const ham_traceback_path = @import("ham/traceback_path.zig");
pub const ham_trellis = @import("ham/trellis.zig");
pub const ham_args = @import("ham/args.zig");
pub const ham_bcr_utils = @import("ham/bcr_utils/root.zig");
pub const ham_dp_handler = @import("ham/dp_handler.zig");
pub const ham_glomerator = @import("ham/glomerator/root.zig");
pub const ham_bcrham = @import("ham/bcrham.zig");
pub const igsw_ig_align = @import("igsw/ig_align.zig");

// ─── bcrham entry points ──────────────────────────────────────────────────────

/// Run the full bcrham pipeline for a set of sequences.
/// Corresponds to the bcrham subprocess invocation in python/partis.py.
/// Arguments mirror bcrham's --parameter-dir, --seqs-file, --outfile, etc.
export fn bcrham_run(
    parameter_dir: [*:0]const u8,
    seqs_file: [*:0]const u8,
    outfile: [*:0]const u8,
    algorithm: [*:0]const u8,
    n_max_per_batch: c_int,
) c_int {
    _ = parameter_dir;
    _ = seqs_file;
    _ = outfile;
    _ = algorithm;
    _ = n_max_per_batch;
    @panic("not yet implemented: bcrham_run");
}

/// Run the Forward algorithm only and return log-probability.
export fn bcrham_forward(
    parameter_dir: [*:0]const u8,
    seqs_file: [*:0]const u8,
    outfile: [*:0]const u8,
) c_int {
    _ = parameter_dir;
    _ = seqs_file;
    _ = outfile;
    @panic("not yet implemented: bcrham_forward");
}

/// Run the Viterbi algorithm and return the most probable path annotation.
export fn bcrham_viterbi(
    parameter_dir: [*:0]const u8,
    seqs_file: [*:0]const u8,
    outfile: [*:0]const u8,
) c_int {
    _ = parameter_dir;
    _ = seqs_file;
    _ = outfile;
    @panic("not yet implemented: bcrham_viterbi");
}

// ─── ig-sw entry points ───────────────────────────────────────────────────────

/// Run Smith-Waterman alignment for V/D/J gene segment matching.
export fn igsw_align(
    query_seq: [*:0]const u8,
    germline_dir: [*:0]const u8,
    outfile: [*:0]const u8,
) c_int {
    _ = query_seq;
    _ = germline_dir;
    _ = outfile;
    @panic("not yet implemented: igsw_align");
}

// ─── library version / health check ──────────────────────────────────────────

/// Return a null-terminated version string. Safe to call before any porting is done.
export fn partis_zig_core_version() [*:0]const u8 {
    return "0.0.1-stub";
}
