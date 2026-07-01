/* glibc_math.c — Force glibc log/exp, bypassing Zig's compiler-rt.
 *
 * Zig's compiler-rt bundles its own log/exp as local text symbols ('t')
 * that may differ from glibc by ~1 ULP on some CPUs.  C++ bcrham always
 * uses glibc, so we must match.
 *
 * We dlopen libm.so.6 and resolve log/exp from it directly, bypassing
 * any compiler-rt symbols in the binary.
 *
 * macOS has no libm.so.6 (no glibc), so there we just use the system libm's
 * log/exp.  Results won't be bit-identical to the Linux glibc reference, but
 * that's fine: the C++ backend (ig-sw/bcrham) can't be built on macOS anyway
 * (SSE2/emmintrin), so there is no Linux-parity requirement on that platform.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__

#include <math.h>

double glibc_log(double x) { return log(x); }
double glibc_exp(double x) { return exp(x); }

#else

#include <dlfcn.h>

typedef double (*math_fn)(double);

static math_fn real_log;
static math_fn real_exp;

__attribute__((constructor))
static void init_glibc_math(void) {
    void *libm = dlopen("libm.so.6", RTLD_LAZY);
    if (!libm) {
        fprintf(stderr, "glibc_math: dlopen libm.so.6 failed: %s\n", dlerror());
        abort();
    }
    real_log = (math_fn)dlsym(libm, "log");
    real_exp = (math_fn)dlsym(libm, "exp");
    if (!real_log || !real_exp) {
        fprintf(stderr, "glibc_math: dlsym failed: %s\n", dlerror());
        abort();
    }
    /* don't dlclose — keep libm loaded for the process lifetime */
}

double glibc_log(double x) { return real_log(x); }
double glibc_exp(double x) { return real_exp(x); }

#endif
