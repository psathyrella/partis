/// ham/glomerator/root.zig — Zig port of ham::Glomerator (directory group)
///
/// Re-exports all sub-modules of the glomerator directory group.
/// C++ source: packages/ham/src/glomerator.cc, packages/ham/include/glomerator.h
/// C++ author: psathyrella/ham

pub const query = @import("query.zig");
pub const glomerator = @import("glomerator.zig");

pub const Query = query.Query;
pub const Glomerator = glomerator.Glomerator;
