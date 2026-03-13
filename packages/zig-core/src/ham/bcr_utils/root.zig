/// ham/bcr_utils/root.zig — Zig port of ham::bcrutils (directory group)
///
/// Re-exports all sub-modules of the bcr_utils directory group.
/// C++ source: packages/ham/src/bcrutils.cc, packages/ham/include/bcrutils.h
/// C++ author: psathyrella/ham

pub const support_pair = @import("support_pair.zig");
pub const insertions = @import("insertions.zig");
pub const term_colors = @import("term_colors.zig");
pub const germ_lines = @import("germ_lines.zig");
pub const k_bounds = @import("k_bounds.zig");
pub const reco_event = @import("reco_event.zig");
pub const hmm_holder = @import("hmm_holder.zig");
pub const result = @import("result.zig");
pub const utils = @import("utils.zig");

// Convenience re-exports of the most-used types
pub const SupportPair = support_pair.SupportPair;
pub const Insertions = insertions.Insertions;
pub const TermColors = term_colors.TermColors;
pub const GermLines = germ_lines.GermLines;
pub const hasDGene = germ_lines.hasDGene;
pub const KSet = k_bounds.KSet;
pub const KBounds = k_bounds.KBounds;
pub const RecoEvent = reco_event.RecoEvent;
pub const HMMHolder = hmm_holder.HMMHolder;
pub const Result = result.Result;
