const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // ── partis-zig-core executable: bcrham-compatible CLI ────────────────
    // `zig build` produces this by default. No external library deps.
    const exe = b.addExecutable(.{
        .name = "partis-zig-core",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    b.installArtifact(exe);

    // ── compare: equivalence harness binary ──────────────────────────────
    const compare = b.addExecutable(.{
        .name = "compare",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/equivalence/compare.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    b.installArtifact(compare);

    // ── lib: C ABI shared library (requires ig-sw sources from partis) ───
    // Not built by default. Use:
    //   zig build lib -Digsw-src=<partis>/packages/ig-sw/src/ig_align
    const igsw_src = b.option([]const u8, "igsw-src",
        "Path to ig-sw C sources (partis/packages/ig-sw/src/ig_align)") orelse "";

    const cabi_mod = b.createModule(.{
        .root_source_file = b.path("src/cabi.zig"),
        .target = target,
        .optimize = optimize,
    });
    if (igsw_src.len > 0) {
        cabi_mod.addIncludePath(.{ .cwd_relative = igsw_src });
    }

    const lib = b.addLibrary(.{
        .name = "partis-zig-core",
        .root_module = cabi_mod,
        .linkage = .dynamic,
    });
    if (igsw_src.len > 0) {
        const flags = &.{ "-std=gnu99", "-O2" };
        lib.addCSourceFile(.{ .file = .{ .cwd_relative = b.pathJoin(&.{ igsw_src, "ig_align.c" }) }, .flags = flags });
        lib.addCSourceFile(.{ .file = .{ .cwd_relative = b.pathJoin(&.{ igsw_src, "ksw.c" }) }, .flags = flags });
        lib.addCSourceFile(.{ .file = .{ .cwd_relative = b.pathJoin(&.{ igsw_src, "kstring.c" }) }, .flags = flags });
        lib.addIncludePath(.{ .cwd_relative = igsw_src });
    }
    lib.linkSystemLibrary("z");
    lib.linkLibC();

    const lib_step = b.step("lib", "Build C ABI shared library (requires -Digsw-src=...)");
    lib_step.dependOn(&b.addInstallArtifact(lib, .{}).step);

    // ── partis-zig-igsw executable: ig-sw-compatible CLI (vendored C sources) ─
    // `zig build` produces this by default alongside partis-zig-core.
    const igsw_exe_mod = b.createModule(.{
        .root_source_file = b.path("src/igsw/main.zig"),
        .target = target,
        .optimize = optimize,
    });
    igsw_exe_mod.addIncludePath(b.path("src/igsw/c"));
    const igsw_exe = b.addExecutable(.{
        .name = "partis-zig-igsw",
        .root_module = igsw_exe_mod,
    });
    // -fno-sanitize=undefined: ig_align.c has a pre-existing shift-overflow UB
    // that is harmless in practice but trips Zig's clang UBSan in Debug builds.
    // ksw.c is NOT listed here — it is replaced by src/igsw/ksw.zig (cross-platform SIMD).
    const igsw_c_flags = &.{ "-std=gnu99", "-O2", "-fno-sanitize=undefined" };
    igsw_exe.addCSourceFile(.{ .file = b.path("src/igsw/c/ig_align.c"), .flags = igsw_c_flags });
    igsw_exe.addCSourceFile(.{ .file = b.path("src/igsw/c/kstring.c"), .flags = igsw_c_flags });
    igsw_exe.addIncludePath(b.path("src/igsw/c"));
    igsw_exe.linkSystemLibrary("z");
    igsw_exe.linkLibC();
    b.installArtifact(igsw_exe);

    // ── test: unit test runner (ham modules only, no ig-sw needed) ────────
    const ham_test_mod = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });
    const unit_tests = b.addTest(.{ .root_module = ham_test_mod });
    const run_unit_tests = b.addRunArtifact(unit_tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_unit_tests.step);
}
