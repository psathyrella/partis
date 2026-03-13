/// igsw/ksw.zig — Zig port of ksw.c (Farrar striped Smith-Waterman)
///
/// Exports C-ABI ksw_align, ksw_global, ksw_qinit, ksw_extend using
/// @Vector(16, u8) / @Vector(8, i16) instead of SSE2 intrinsics.
/// Compiles and runs correctly on any platform (ARM64, x86, RISC-V, …).
///
/// C source: packages/ig-sw/src/ig_align/ksw.c (Attractive Chaos, MIT)
/// Algorithm: Farrar (2007) striped Smith-Waterman

const std = @import("std");

// ─── vector types ─────────────────────────────────────────────────────────────

/// 16 unsigned bytes  ≡  __m128i in the u8 path
const V16u8 = @Vector(16, u8);
/// 8 signed 16-bit integers  ≡  __m128i in the i16 path
const V8i16 = @Vector(8, i16);
/// 8 unsigned 16-bit — for _mm_subs_epu16 emulation
const V8u16 = @Vector(8, u16);

// ─── helper: _mm_slli_si128 emulation ─────────────────────────────────────────
//
// _mm_slli_si128(v, 1 byte) shifts the 128-bit register left by 1 byte:
//   result[0] = 0,  result[i] = v[i-1]   (i >= 1)
// In @Vector terms: elements move toward higher indices; 0 fills index 0.
//
// For V16u8:  shift by 1 element = 1 byte
// For V8i16:  shift by 1 element = 2 bytes (≡ slli_si128 by 2)

inline fn slli1_u8(v: V16u8) V16u8 {
    const z: V16u8 = @splat(0);
    return @shuffle(u8, v, z, [16]i32{ -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 });
}

inline fn slli1_i16(v: V8i16) V8i16 {
    const z: V8i16 = @splat(0);
    return @shuffle(i16, v, z, [8]i32{ -1, 0, 1, 2, 3, 4, 5, 6 });
}

// ─── helper: _mm_subs_epu16 emulation ─────────────────────────────────────────
// Unsigned saturating subtract on 16-bit elements.
inline fn subsEpu16(a: V8i16, b: V8i16) V8i16 {
    const au: V8u16 = @bitCast(a);
    const bu: V8u16 = @bitCast(b);
    return @bitCast(au -| bu);
}

// ─── kswr_t  (must match C struct layout exactly) ────────────────────────────

pub const Kswr = extern struct {
    score: c_int = 0,
    te: c_int = -1,
    qe: c_int = -1,
    score2: c_int = -1,
    te2: c_int = -1,
    tb: c_int = -1,
    qb: c_int = -1,
};

const g_defr = Kswr{};

// KSW flags
const KSW_XBYTE: c_int = 0x10000;
const KSW_XSTOP: c_int = 0x20000;
const KSW_XSUBO: c_int = 0x40000;
const KSW_XSTART: c_int = 0x80000;

// ─── kswq_t internal struct ───────────────────────────────────────────────────
//
// Completely opaque to ig_align.c — it only ever holds a kswq_t* pointer.
// Single malloc block: [Kswq][padding][qp vectors][H0..Hmax vectors]
// All pointer fields point into the same block so free(q) releases everything.

const Kswq = struct {
    qlen: c_int,
    slen: c_int,
    shift: u8,
    mdiff: u8,
    max: u8,
    size: u8, // 1 = u8 path, 2 = i16 path
    // These point into the malloc block past this struct (16-byte aligned)
    qp: [*]align(16) [16]u8,
    H0: [*]align(16) [16]u8,
    H1: [*]align(16) [16]u8,
    E: [*]align(16) [16]u8,
    Hmax: [*]align(16) [16]u8,
};

// ─── ksw_qinit ────────────────────────────────────────────────────────────────

pub export fn ksw_qinit(
    size: c_int,
    qlen: c_int,
    query: [*]const u8,
    m: c_int,
    mat: [*]const i8,
) callconv(.c) ?*Kswq {
    const sz: usize = if (size > 1) 2 else 1;
    const p: usize = 8 * (3 - sz); // 16 for u8, 8 for i16
    const slen: usize = (@as(usize, @intCast(qlen)) + p - 1) / p;
    const mu: usize = @intCast(m);

    const total: usize = @sizeOf(Kswq) + 256 + 16 * slen * (mu + 4);
    const raw = std.c.malloc(total) orelse return null;

    // Zero the entire block so DP arrays start at zero
    @memset(@as([*]u8, @ptrCast(raw))[0..total], 0);

    const q: *Kswq = @ptrCast(@alignCast(raw));

    // Align vector area to 16 bytes
    const arr_start: usize = (@intFromPtr(raw) + @sizeOf(Kswq) + 15) & ~@as(usize, 15);
    const vec_base: [*]align(16) [16]u8 = @ptrFromInt(arr_start);

    q.qlen = qlen;
    q.slen = @intCast(slen);
    q.size = @intCast(sz);
    q.qp = vec_base;
    q.H0 = vec_base + slen * mu;
    q.H1 = q.H0 + slen;
    q.E = q.H1 + slen;
    q.Hmax = q.E + slen;

    // Compute shift (minimum score) and mdiff — identical for both paths
    const mat_len: usize = mu * mu;
    var min_s: i32 = 127;
    var max_s: i32 = 0;
    for (0..mat_len) |i| {
        const v: i32 = @intCast(mat[i]);
        if (v < min_s) min_s = v;
        if (v > max_s) max_s = v;
    }
    q.max = @intCast(max_s);
    // q->shift = 256 - min_score (as uint8, may wrap)
    q.shift = @truncate(@as(u32, @intCast(256 - min_s)));
    // q->mdiff = max_score + shift
    q.mdiff = @intCast(max_s + @as(i32, q.shift));

    const qlen_u: usize = @intCast(qlen);

    if (sz == 1) {
        // u8 path: 16 elements per segment; scores stored with +shift bias
        const t: [*]u8 = @ptrCast(q.qp);
        var ti: usize = 0;
        for (0..mu) |a| {
            const nlen = slen * p;
            const ma: [*]const i8 = mat + a * mu;
            for (0..slen) |i| {
                var k: usize = i;
                while (k < nlen) : (k += slen) {
                    t[ti] = if (k >= qlen_u) 0 else @intCast(@as(i32, @intCast(ma[query[k]])) + @as(i32, q.shift));
                    ti += 1;
                }
            }
        }
    } else {
        // i16 path: 8 elements per segment; raw scores (no bias)
        const t: [*]i16 = @ptrCast(@alignCast(q.qp));
        var ti: usize = 0;
        for (0..mu) |a| {
            const nlen = slen * p;
            const ma: [*]const i8 = mat + a * mu;
            for (0..slen) |i| {
                var k: usize = i;
                while (k < nlen) : (k += slen) {
                    t[ti] = if (k >= qlen_u) 0 else @intCast(ma[query[k]]);
                    ti += 1;
                }
            }
        }
    }
    return q;
}

// ─── ksw_u8 (internal) ────────────────────────────────────────────────────────

fn ksw_u8(q: *Kswq, tlen: c_int, target: [*]const u8, _gapo: c_int, _gape: c_int, xtra: c_int) Kswr {
    const slen: usize = @intCast(q.slen);
    if (slen == 0) return g_defr;
    var r = g_defr;

    const minsc: i32 = if ((xtra & KSW_XSUBO) != 0) xtra & 0xffff else 0x10000;
    const endsc: i32 = if ((xtra & KSW_XSTOP) != 0) xtra & 0xffff else 0x10000;

    const zero: V16u8 = @splat(0);
    const gapoe: V16u8 = @splat(@intCast(_gapo + _gape));
    const gape_v: V16u8 = @splat(@intCast(_gape));
    const shift_v: V16u8 = @splat(q.shift);

    var H0: [*]V16u8 = @ptrCast(q.H0);
    var H1: [*]V16u8 = @ptrCast(q.H1);
    const E: [*]V16u8 = @ptrCast(q.E);
    const Hmax: [*]V16u8 = @ptrCast(q.Hmax);

    // Zero DP arrays
    for (0..slen) |i| {
        H0[i] = zero;
        H1[i] = zero;
        E[i] = zero;
        Hmax[i] = zero;
    }

    var m_b: usize = 0;
    var n_b: usize = 0;
    var b: ?[*]u64 = null;
    defer if (b) |ptr| std.c.free(ptr);

    var gmax: i32 = 0;
    var te: i32 = -1;

    const tlen_u: usize = @intCast(tlen);
    const qp: [*]V16u8 = @ptrCast(q.qp);

    var i: usize = 0;
    while (i < tlen_u) : (i += 1) {
        // Load last H0 segment, shift left by 1 byte
        var h = slli1_u8(H0[slen - 1]);
        var f: V16u8 = zero;
        var max_v: V16u8 = zero;
        const S: [*]V16u8 = qp + target[i] * slen;

        var j: usize = 0;
        while (j < slen) : (j += 1) {
            // H' = max(H(i-1,j-1) + S(i,j) - shift, E(i,j), F(i,j))
            h = (h +| S[j]) -| shift_v;
            const e = E[j];
            h = @max(h, e);
            h = @max(h, f);
            max_v = @max(max_v, h);
            H1[j] = h;
            // E(i+1,j) = max(H(i,j)-gapo, E(i,j)) - gape → stored as primed
            h = h -| gapoe;
            E[j] = @max(h, e -| gape_v);
            // F(i,j+1)
            f = @max(h, f -| gape_v);
            // advance h
            h = H0[j];
        }

        // Lazy-F loop: propagate F that may have been updated
        var k: usize = 0;
        while (k < 16) : (k += 1) {
            f = slli1_u8(f);
            var done = true;
            var jj: usize = 0;
            while (jj < slen) : (jj += 1) {
                h = H1[jj];
                h = @max(h, f);
                H1[jj] = h;
                h = h -| gapoe;
                f = f -| gape_v; // lazy-F: only decrement, no max with h
                const diff = f -| h;
                if (@reduce(.Or, diff != zero)) done = false;
            }
            if (done) break;
        }

        // Horizontal max of max_v
        const imax: i32 = @reduce(.Max, max_v);

        if (imax >= minsc) {
            if (n_b == 0 or @as(i32, @truncate(@as(i64, @bitCast(b.?[n_b - 1])))) + 1 != @as(i32, @intCast(i))) {
                if (n_b == m_b) {
                    m_b = if (m_b != 0) m_b << 1 else 8;
                    b = @ptrCast(@alignCast(std.c.realloc(if (b) |p| p else null, 8 * m_b)));
                }
                b.?[n_b] = @as(u64, @intCast(imax)) << 32 | i;
                n_b += 1;
            } else if (@as(i32, @intCast(b.?[n_b - 1] >> 32)) < imax) {
                b.?[n_b - 1] = @as(u64, @intCast(imax)) << 32 | i;
            }
        }

        if (imax > gmax) {
            gmax = imax;
            te = @intCast(i);
            for (0..slen) |jj| Hmax[jj] = H1[jj];
            if (gmax + @as(i32, q.shift) >= 255 or gmax >= endsc) break;
        }

        // Swap H0 ↔ H1
        const tmp_H = H0;
        H0 = H1;
        H1 = tmp_H;
    }

    r.score = if (gmax + @as(i32, q.shift) < 255) gmax else 255;
    r.te = te;

    if (r.score != 255) {
        var best: i32 = -1;
        const qlen_loop: usize = slen * 16;
        const hmax_bytes: [*]u8 = @ptrCast(q.Hmax);
        var qi: usize = 0;
        while (qi < qlen_loop) : (qi += 1) {
            const v: i32 = @intCast(hmax_bytes[qi]);
            const pos: i32 = @intCast(qi / 16 + qi % 16 * slen);
            if (v > best) {
                best = v;
                r.qe = pos;
            } else if (v == best and pos < r.qe) {
                r.qe = pos;
            }
        }
        if (b) |bp| {
            const margin: i32 = @divTrunc(r.score + @as(i32, q.max) - 1, @as(i32, q.max));
            const low: i32 = te - margin;
            const high: i32 = te + margin;
            var bi: usize = 0;
            while (bi < n_b) : (bi += 1) {
                const e: i32 = @truncate(@as(i64, @bitCast(bp[bi])));
                if ((e < low or e > high) and @as(i32, @intCast(bp[bi] >> 32)) > r.score2) {
                    r.score2 = @intCast(bp[bi] >> 32);
                    r.te2 = e;
                }
            }
        }
    }
    return r;
}

// ─── ksw_i16 (internal) ───────────────────────────────────────────────────────

fn ksw_i16(q: *Kswq, tlen: c_int, target: [*]const u8, _gapo: c_int, _gape: c_int, xtra: c_int) Kswr {
    const slen: usize = @intCast(q.slen);
    if (slen == 0) return g_defr;
    var r = g_defr;

    const minsc: i32 = if ((xtra & KSW_XSUBO) != 0) xtra & 0xffff else 0x10000;
    const endsc: i32 = if ((xtra & KSW_XSTOP) != 0) xtra & 0xffff else 0x10000;

    const zero: V8i16 = @splat(0);
    const gapoe: V8i16 = @splat(@intCast(_gapo + _gape));
    const gape_v: V8i16 = @splat(@intCast(_gape));

    var H0: [*]V8i16 = @ptrCast(q.H0);
    var H1: [*]V8i16 = @ptrCast(q.H1);
    const E: [*]V8i16 = @ptrCast(q.E);
    const Hmax: [*]V8i16 = @ptrCast(q.Hmax);
    const qp: [*]V8i16 = @ptrCast(q.qp);

    for (0..slen) |i| {
        H0[i] = zero;
        H1[i] = zero;
        E[i] = zero;
        Hmax[i] = zero;
    }

    var m_b: usize = 0;
    var n_b: usize = 0;
    var b: ?[*]u64 = null;
    defer if (b) |ptr| std.c.free(ptr);

    var gmax: i32 = 0;
    var te: i32 = -1;
    const tlen_u: usize = @intCast(tlen);

    var i: usize = 0;
    while (i < tlen_u) : (i += 1) {
        var h = slli1_i16(H0[slen - 1]);
        var f: V8i16 = zero;
        var max_v: V8i16 = zero;
        const S: [*]V8i16 = qp + target[i] * slen;

        var j: usize = 0;
        while (j < slen) : (j += 1) {
            h = h +| S[j];
            const e = E[j];
            h = @max(h, e);
            h = @max(h, f);
            max_v = @max(max_v, h);
            H1[j] = h;
            h = subsEpu16(h, gapoe);
            E[j] = @max(h, subsEpu16(e, gape_v));
            f = @max(h, subsEpu16(f, gape_v));
            h = H0[j];
        }

        var k: usize = 0;
        while (k < 16) : (k += 1) {
            f = slli1_i16(f);
            var done = true;
            var jj: usize = 0;
            while (jj < slen) : (jj += 1) {
                h = H1[jj];
                h = @max(h, f);
                H1[jj] = h;
                h = subsEpu16(h, gapoe);
                f = subsEpu16(f, gape_v); // lazy-F: only decrement, no max with h
                if (@reduce(.Or, f > h)) done = false;
            }
            if (done) break;
        }

        const imax: i32 = @reduce(.Max, max_v);

        if (imax >= minsc) {
            if (n_b == 0 or @as(i32, @truncate(@as(i64, @bitCast(b.?[n_b - 1])))) + 1 != @as(i32, @intCast(i))) {
                if (n_b == m_b) {
                    m_b = if (m_b != 0) m_b << 1 else 8;
                    b = @ptrCast(@alignCast(std.c.realloc(if (b) |p| p else null, 8 * m_b)));
                }
                b.?[n_b] = @as(u64, @intCast(imax)) << 32 | i;
                n_b += 1;
            } else if (@as(i32, @intCast(b.?[n_b - 1] >> 32)) < imax) {
                b.?[n_b - 1] = @as(u64, @intCast(imax)) << 32 | i;
            }
        }

        if (imax > gmax) {
            gmax = imax;
            te = @intCast(i);
            for (0..slen) |jj| Hmax[jj] = H1[jj];
            if (gmax >= endsc) break;
        }

        // Swap H0 ↔ H1 (local vars, like ksw_u8)
        const tmp_H = H0;
        H0 = H1;
        H1 = tmp_H;
    }

    r.score = gmax;
    r.te = te;

    {
        var best: i32 = -1;
        const qlen_loop: usize = slen * 8;
        const hmax_shorts: [*]i16 = @ptrCast(@alignCast(q.Hmax));
        var qi: usize = 0;
        r.qe = -1;
        while (qi < qlen_loop) : (qi += 1) {
            const v: i32 = @intCast(hmax_shorts[qi]);
            const pos: i32 = @intCast(qi / 8 + qi % 8 * slen);
            if (v > best) {
                best = v;
                r.qe = pos;
            } else if (v == best and pos < r.qe) {
                r.qe = pos;
            }
        }
        if (b) |bp| {
            const margin: i32 = @divTrunc(r.score + @as(i32, q.max) - 1, @as(i32, q.max));
            const low: i32 = te - margin;
            const high: i32 = te + margin;
            var bi: usize = 0;
            while (bi < n_b) : (bi += 1) {
                const e: i32 = @truncate(@as(i64, @bitCast(bp[bi])));
                if ((e < low or e > high) and @as(i32, @intCast(bp[bi] >> 32)) > r.score2) {
                    r.score2 = @intCast(bp[bi] >> 32);
                    r.te2 = e;
                }
            }
        }
    }
    return r;
}

// ─── revseq ───────────────────────────────────────────────────────────────────

fn revseq(l: usize, s: [*]u8) void {
    var i: usize = 0;
    while (i < l >> 1) : (i += 1) {
        const tmp = s[i];
        s[i] = s[l - 1 - i];
        s[l - 1 - i] = tmp;
    }
}

// ─── ksw_align ────────────────────────────────────────────────────────────────

pub export fn ksw_align(
    qlen: c_int,
    query: [*]u8,
    tlen: c_int,
    target: [*]u8,
    m: c_int,
    mat: [*]const i8,
    gapo: c_int,
    gape: c_int,
    xtra: c_int,
    qry: ?*?*Kswq,
) callconv(.c) Kswr {
    // Determine size (1 = u8, 2 = i16)
    const use_byte = (xtra & KSW_XBYTE) != 0;

    // Get or allocate query profile
    const q: *Kswq = blk: {
        if (qry) |qp| {
            if (qp.*) |existing| {
                break :blk existing;
            } else {
                const nq = ksw_qinit(if (use_byte) 1 else 2, qlen, query, m, mat) orelse return g_defr;
                qp.* = nq;
                break :blk nq;
            }
        } else {
            break :blk ksw_qinit(if (use_byte) 1 else 2, qlen, query, m, mat) orelse return g_defr;
        }
    };
    const own_q = (qry == null); // we own q if caller passed null for qry

    const size = q.size;
    const func: *const fn (*Kswq, c_int, [*]const u8, c_int, c_int, c_int) Kswr =
        if (size == 2) &ksw_i16 else &ksw_u8;

    var r = func(q, tlen, target, gapo, gape, xtra);

    if (own_q) std.c.free(q);

    if ((xtra & KSW_XSTART) == 0 or ((xtra & KSW_XSUBO) != 0 and r.score < (xtra & 0xffff))) return r;

    // Find start positions by reversing and re-aligning
    const qe: usize = @intCast(r.qe + 1);
    const te_end: usize = @intCast(r.te + 1);
    revseq(qe, query);
    revseq(te_end, target);
    const q2 = ksw_qinit(@intCast(size), r.qe + 1, query, m, mat) orelse {
        revseq(qe, query);
        revseq(te_end, target);
        return r;
    };
    const rr = func(q2, tlen, target, gapo, gape, KSW_XSTOP | r.score);
    std.c.free(q2);
    revseq(qe, query);
    revseq(te_end, target);

    if (r.score == rr.score) {
        r.tb = r.te - rr.te;
        r.qb = r.qe - rr.qe;
    }
    return r;
}

// ─── ksw_extend (scalar, no SIMD) ────────────────────────────────────────────

pub export fn ksw_extend(
    qlen: c_int,
    query: [*]const u8,
    tlen: c_int,
    target: [*]const u8,
    m: c_int,
    mat: [*]const i8,
    gapo: c_int,
    gape: c_int,
    w: c_int,
    end_bonus: c_int,
    zdrop: c_int,
    h0: c_int,
    _qle: ?*c_int,
    _tle: ?*c_int,
    _gtle: ?*c_int,
    _gscore: ?*c_int,
    _max_off: ?*c_int,
) callconv(.c) c_int {
    const qlen_u: usize = @intCast(qlen);
    const tlen_u: usize = @intCast(tlen);
    const mu: usize = @intCast(m);
    const gapoe: i32 = gapo + gape;
    const h0v: i32 = if (h0 < 0) 0 else h0;

    // query profile
    const qp: [*]i8 = @ptrCast(std.c.malloc(qlen_u * mu) orelse return 0);
    defer std.c.free(qp);
    // eh array: {H, E} pairs
    const Eh = extern struct { h: i32, e: i32 };
    const eh: [*]Eh = @ptrCast(@alignCast(std.c.malloc((qlen_u + 1) * @sizeOf(Eh)) orelse return 0));
    defer std.c.free(eh);
    @memset(@as([*]u8, @ptrCast(eh))[0..(qlen_u + 1) * @sizeOf(Eh)], 0);

    // fill query profile
    var k: usize = 0;
    var idx: usize = 0;
    while (k < mu) : (k += 1) {
        const p: [*]const i8 = mat + k * mu;
        var j: usize = 0;
        while (j < qlen_u) : (j += 1) {
            qp[idx] = p[query[j]];
            idx += 1;
        }
    }

    // fill first row
    eh[0].h = h0v;
    eh[1].h = if (h0v > gapoe) h0v - gapoe else 0;
    {
        var j: usize = 2;
        while (j <= qlen_u and eh[j - 1].h > gape) : (j += 1) {
            eh[j].h = eh[j - 1].h - gape;
        }
    }

    // bandwidth
    var max_score: i32 = 0;
    {
        var ki: usize = 0;
        while (ki < mu * mu) : (ki += 1) {
            if (@as(i32, @intCast(mat[ki])) > max_score) max_score = @intCast(mat[ki]);
        }
    }
    var max_gap: i32 = @intFromFloat(@as(f64, @floatFromInt(qlen_u)) * @as(f64, @floatFromInt(max_score)) + @as(f64, @floatFromInt(end_bonus)) - @as(f64, @floatFromInt(gapo)));
    max_gap = @divFloor(max_gap, gape) + 1;
    if (max_gap < 1) max_gap = 1;
    const eff_w: i32 = if (w < max_gap) w else max_gap;

    var max: i32 = h0v;
    var max_i: i32 = -1;
    var max_j: i32 = -1;
    var max_ie: i32 = -1;
    var gscore: i32 = -1;
    var max_off: i32 = 0;
    var beg: i32 = 0;
    var end_: i32 = @intCast(qlen_u);

    var i: i32 = 0;
    while (i < @as(i32, @intCast(tlen_u))) : (i += 1) {
        const target_i: i8 = @intCast(target[@intCast(i)]);
        var f: i32 = 0;
        var h1: i32 = h0v - (gapo + gape * (i + 1));
        if (h1 < 0) h1 = 0;
        const q_row: [*]const i8 = qp + @as(usize, @intCast(target_i)) * qlen_u;
        var cur_m: i32 = 0;
        var mj: i32 = -1;
        if (beg < i - eff_w) beg = i - eff_w;
        if (end_ > i + eff_w + 1) end_ = i + eff_w + 1;
        if (end_ > @as(i32, @intCast(qlen_u))) end_ = @intCast(qlen_u);
        var j: i32 = beg;
        while (j < end_) : (j += 1) {
            const ju: usize = @intCast(j);
            var h = eh[ju].h;
            var e = eh[ju].e;
            eh[ju].h = h1;
            h += @intCast(q_row[ju]);
            h = @max(h, e);
            h = @max(h, f);
            h1 = h;
            if (h > cur_m) { cur_m = h; mj = j; }
            h -= gapoe;
            if (h < 0) h = 0;
            e -= gape;
            eh[ju].e = if (e > h) e else h;
            f -= gape;
            f = if (f > h) f else h;
        }
        eh[@intCast(end_)].h = h1;
        eh[@intCast(end_)].e = 0;
        if (j == @as(i32, @intCast(qlen_u))) {
            if (h1 > gscore) { gscore = h1; max_ie = i; }
        }
        if (cur_m == 0 or (zdrop > 0 and @as(i32, @intCast(@abs(i - max_i))) + @as(i32, @intCast(@abs(mj - max_j))) > @divFloor(max - cur_m, gape) + zdrop)) break;
        if (cur_m > max) {
            max = cur_m;
            max_i = i;
            max_j = mj;
            const off: i32 = @intCast(@abs(mj - i));
            if (off > max_off) max_off = off;
        }
        // update beg/end for next round
        j = mj;
        while (j >= beg and eh[@intCast(j)].h != 0) j -= 1;
        beg = j + 1;
        j = mj + 2;
        while (j <= end_ and eh[@intCast(j)].h != 0) j += 1;
        end_ = j;
    }

    if (_qle) |p| p.* = max_j + 1;
    if (_tle) |p| p.* = max_i + 1;
    if (_gtle) |p| p.* = max_ie + 1;
    if (_gscore) |p| p.* = gscore;
    if (_max_off) |p| p.* = max_off;
    return @intCast(max);
}

// ─── ksw_global (scalar, no SIMD) ────────────────────────────────────────────

pub export fn ksw_global(
    qlen: c_int,
    query: [*]const u8,
    tlen: c_int,
    target: [*]const u8,
    m: c_int,
    mat: [*]const i8,
    gapo: c_int,
    gape: c_int,
    w: c_int,
    n_cigar_: ?*c_int,
    cigar_: ?*?[*]u32,
) callconv(.c) c_int {
    const qlen_u: usize = @intCast(qlen);
    const tlen_u: usize = @intCast(tlen);
    const mu: usize = @intCast(m);
    const gapoe: i32 = gapo + gape;
    const MINUS_INF: i32 = -0x40000000;

    if (n_cigar_) |p| p.* = 0;
    if (qlen_u == 0 or tlen_u == 0) return 0;

    const n_col: usize = if (@as(usize, @intCast(qlen)) < @as(usize, @intCast(2 * w + 1)))
        qlen_u
    else
        @intCast(2 * w + 1);

    // z: backtrack matrix n_col * tlen
    const z: [*]u8 = @ptrCast(std.c.malloc(n_col * tlen_u) orelse return 0);
    defer std.c.free(z);
    const Eh = extern struct { h: i32, e: i32 };
    const qp: [*]i8 = @ptrCast(std.c.malloc(qlen_u * mu) orelse return 0);
    defer std.c.free(qp);
    const eh: [*]Eh = @ptrCast(@alignCast(std.c.calloc(qlen_u + 1, @sizeOf(Eh)) orelse return 0));
    defer std.c.free(eh);

    // query profile
    {
        var k: usize = 0;
        var ti: usize = 0;
        while (k < mu) : (k += 1) {
            const p: [*]const i8 = mat + k * mu;
            var j: usize = 0;
            while (j < qlen_u) : (j += 1) {
                qp[ti] = p[query[j]];
                ti += 1;
            }
        }
    }

    // fill first row
    eh[0].h = 0;
    eh[0].e = MINUS_INF;
    {
        var j: usize = 1;
        while (j <= qlen_u and @as(i32, @intCast(j)) <= w) : (j += 1) {
            eh[j].h = -(gapo + gape * @as(i32, @intCast(j)));
            eh[j].e = MINUS_INF;
        }
        while (j <= qlen_u) : (j += 1) {
            eh[j].h = MINUS_INF;
            eh[j].e = MINUS_INF;
        }
    }

    // DP
    var i: i32 = 0;
    while (i < @as(i32, @intCast(tlen_u))) : (i += 1) {
        const iu: usize = @intCast(i);
        var f: i32 = MINUS_INF;
        const beg: i32 = if (i > w) i - w else 0;
        const end_: i32 = if (i + w + 1 < @as(i32, @intCast(qlen_u))) i + w + 1 else @intCast(qlen_u);
        const q_row: [*]const i8 = qp + @as(usize, @intCast(target[iu])) * qlen_u;
        const zi: [*]u8 = z + iu * n_col;
        var h1: i32 = if (beg == 0) -(gapo + gape * (i + 1)) else MINUS_INF;
        var j: i32 = beg;
        while (j < end_) : (j += 1) {
            const ju: usize = @intCast(j);
            var h = eh[ju].h;
            var e = eh[ju].e;
            var d: u8 = undefined;
            eh[ju].h = h1;
            h += @intCast(q_row[ju]);
            d = if (h >= e) 0 else 1;
            h = if (h >= e) h else e;
            d = if (h >= f) d else 2;
            h = if (h >= f) h else f;
            h1 = h;
            h -= gapoe;
            e -= gape;
            d |= if (e > h) @as(u8, 1) << 2 else 0;
            e = if (e > h) e else h;
            eh[ju].e = e;
            f -= gape;
            d |= if (f > h) @as(u8, 2) << 4 else 0;
            f = if (f > h) f else h;
            zi[@intCast(j - beg)] = d;
        }
        eh[@intCast(end_)].h = h1;
        eh[@intCast(end_)].e = MINUS_INF;
    }

    const score = eh[qlen_u].h;

    // backtrack
    if (n_cigar_ != null and cigar_ != null) {
        var n_cigar: c_int = 0;
        var m_cigar: c_int = 0;
        var cigar: ?[*]u32 = null;
        var which: i32 = 0;
        var ti: i32 = @intCast(tlen_u - 1);
        var k: i32 = @intCast((if (ti + w + 1 < @as(i32, @intCast(qlen_u))) ti + w + 1 else @as(i32, @intCast(qlen_u))) - 1);
        while (ti >= 0 and k >= 0) {
            const beg_i: i32 = if (ti > w) ti - w else 0;
            const d = z[@as(usize, @intCast(ti)) * n_col + @as(usize, @intCast(k - beg_i))];
            which = (d >> @intCast(which * 2)) & 3;
            if (which == 0) {
                cigar = pushCigar(&n_cigar, &m_cigar, cigar, 0, 1);
                ti -= 1;
                k -= 1;
            } else if (which == 1) {
                cigar = pushCigar(&n_cigar, &m_cigar, cigar, 2, 1);
                ti -= 1;
            } else {
                cigar = pushCigar(&n_cigar, &m_cigar, cigar, 1, 1);
                k -= 1;
            }
        }
        if (ti >= 0) cigar = pushCigar(&n_cigar, &m_cigar, cigar, 2, @intCast(ti + 1));
        if (k >= 0) cigar = pushCigar(&n_cigar, &m_cigar, cigar, 1, @intCast(k + 1));
        // reverse
        var lo: usize = 0;
        var hi: usize = @intCast(n_cigar - 1);
        while (lo < hi) : ({ lo += 1; hi -= 1; }) {
            const tmp = cigar.?[lo];
            cigar.?[lo] = cigar.?[hi];
            cigar.?[hi] = tmp;
        }
        n_cigar_.?.* = n_cigar;
        cigar_.?.* = cigar;
    }

    return @intCast(score);
}

fn pushCigar(n_cigar: *c_int, m_cigar: *c_int, cigar: ?[*]u32, op: u32, len: u32) ?[*]u32 {
    if (n_cigar.* == 0 or op != (cigar.?[@intCast(n_cigar.* - 1)] & 0xf)) {
        var result = cigar;
        if (n_cigar.* == m_cigar.*) {
            m_cigar.* = if (m_cigar.* != 0) m_cigar.* << 1 else 4;
            result = @ptrCast(@alignCast(std.c.realloc(if (cigar) |p| @ptrCast(p) else null, @intCast(m_cigar.* * 4))));
        }
        result.?[@intCast(n_cigar.*)] = len << 4 | op;
        n_cigar.* += 1;
        return result;
    } else {
        cigar.?[@intCast(n_cigar.* - 1)] += len << 4;
        return cigar;
    }
}
