/* fast_math.c — fast log-sum-exp via tabulated softplus.
 *
 * Drop-in replacement for the `log(1 + exp(d))` work in AddInLogSpace
 * (mathutils.h:98). Replaces 1 log + 1 exp glibc call (~5–17 ns combined,
 * arch-dependent) with a single 65 K-entry linear-interpolation lookup
 * (~3.5–10 ns per call).
 *
 * Pre-flight measurements at issue #366 item 3.2:
 *   - Zen 5 (glibc 2.39): glibc 9.6 ns → LUT64K 7.3 ns per addInLogSpace
 *   - Broadwell:          glibc 16.7 ns → LUT64K 10.0 ns
 *   - LUT64K accuracy: max abs error 6.5e-9 over [-30, 0] (5 orders under
 *     the partis-test 1e-5 gate).
 *
 * Why a 65 K linear-interp table:
 *   The Cephes scalar polynomial alternative (~7-15 ns) loses to modern
 *   glibc's FMA-optimized scalar log/exp on every CPU we tested. LUT
 *   wins because it replaces ~10 ALU ops + branches with a single load
 *   + 1 multiply + 1 fma. Smaller LUTs (256, 4 K, 16 K) fail the worst-
 *   case error budget; LUT64K passes by 4×. Cubic interp at 1 K entries
 *   would also pass speed-wise but not accuracy-wise.
 *
 * The linear-interp arithmetic dominates over the table size in cycle
 * count (LUT256 ≈ LUT64K on every benchmarked CPU), so going larger is
 * effectively free as long as the hot region stays cache-resident — and
 * it does, because real DP traffic concentrates near d = 0.
 *
 * Built with -O3 -ffp-contract=off so that the linear-interp arithmetic
 * is bit-deterministic across g++/clang/zig from the same C source.
 * (Without that flag, cross-compiler drift would be 1–2 ULP — far under
 * the partis-test gate, but bit-equality between backends is cheap to
 * keep so we keep it.)
 *
 * Mirrored on the Zig side at packages/zig-core/src/ham/fast_math.c.
 * Both sides MUST stay in sync; the table is data so cross-backend
 * parity is trivial as long as the same source compiles to both.
 */
#include <math.h>

/* Range of d for which the table covers softplus(d) = log(1 + exp(d))
 * directly. Outside this range, we use the closed-form fall-throughs:
 *   d < LUT_LO:  exp(d) is well under double precision, softplus ≈ 0.
 *   d > 0:       use the symmetry softplus(d) = d + softplus(-d).
 *
 * Within addInLogSpace's standard form (see mathutils.h AddInLogSpace),
 * the argument is always (smaller - larger) ≤ 0. The d > 0 branch is
 * defensive — should not be exercised in practice, but keeps the
 * function safe to call independently.
 */
#define LUT_N  65536
#define LUT_LO -30.0
#define LUT_HI 0.0

static double softplus_lut[LUT_N + 1];

/* (LUT_N / |LUT_LO|) precomputed for the index calculation. Using a
 * macro keeps it a compile-time constant, which the compiler can fold
 * into the multiply at the call site. */
#define LUT_SCALE ((double) LUT_N / -(LUT_LO))

__attribute__((constructor))
static void init_softplus_lut(void) {
    /* log1p(exp(d)) is the numerically stable form of log(1 + exp(d))
     * for the table-build phase. Build is one-shot at process start. */
    for (int i = 0; i <= LUT_N; ++i) {
        double d = LUT_LO + (LUT_HI - LUT_LO) * (double) i / (double) LUT_N;
        softplus_lut[i] = log1p(exp(d));
    }
}

double fast_softplus(double d) {
    if (d < LUT_LO) return 0.0;
    if (d > 0.0)    return d + fast_softplus(-d);   /* symmetry */
    /* d ∈ [LUT_LO, 0] — index in [0, LUT_N] */
    double t = (d - LUT_LO) * LUT_SCALE;
    int i = (int) t;
    /* Clamp `i` so `i+1` stays in bounds. At d == 0 (or d so close to 0
     * that float rounding gives t == LUT_N), the unclamped i would equal
     * LUT_N and the `softplus_lut[i+1]` read would be one past the array.
     * The production result was correct only because `frac == 0` masked
     * the OOB read; Debug-mode bounds checks fired. With the clamp,
     * frac becomes 1.0 at the boundary and y0 + 1.0*(y1-y0) = y1 =
     * softplus_lut[LUT_N], identical to the unclamped value. */
    if (i >= LUT_N) i = LUT_N - 1;
    double frac = t - (double) i;
    double y0 = softplus_lut[i];
    double y1 = softplus_lut[i + 1];
    return y0 + frac * (y1 - y0);
}
