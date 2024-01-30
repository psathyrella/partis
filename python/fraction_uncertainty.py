from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
from .cache import cached_uncertainties
from scipy.stats import beta
from io import open

# ----------------------------------------------------------------------------------------
def err(obs, total, use_cache=True):
    """ Return uncertainty on the ratio obs / total """
    assert obs <= total
    assert obs >= 0
    assert total >= 0

    if total == 0.0:
        return (0.0, 0.0)

    key = str(obs) + '/' + str(total)
    if use_cache and key in cached_uncertainties.errs:
        return cached_uncertainties.errs[key]

    frac = float(obs) / total
    vol = 2./3  # +/- 1 sigma
    cpr = 1.  # constant prior
    def call_ppf(lower_tail_prob, eff_obs, total):
        return beta.ppf(lower_tail_prob, cpr + eff_obs, cpr + total - eff_obs)

    if obs == 0:
        lo = 0.
        eff_lo = call_ppf((1. - vol)/2, 1, total)  # take the width from obs of 1
        eff_hi = call_ppf((1. + vol)/2, 1, total)
        hi = eff_hi - eff_lo
    else:
        lo = call_ppf((1. - vol)/2, obs, total)
        if frac < lo:
            lo = 0.
            hi = call_ppf(vol, obs, total)  # one-sided c.i. with 2/3 (0.95) the mass (probably doesn't happen much any more, now that obs == 0 is handled separately
        else:
            hi = call_ppf((1. + vol)/2, obs, total)

    if frac > hi:  # same deal if obs/total very large (probably one)
        lo = call_ppf(1. - vol, obs, total)  # I would imagine that this, also, needs to be changed (same as 0) but I don't want to do it right now
        hi = 1.

    assert lo < frac or frac == 0.0
    assert frac < hi or frac == 1.0
    return (lo, hi)

# ----------------------------------------------------------------------------------------
def chk():  # check current version output against cached values
    iline = 0
    eps = 1e-6
    for key in cached_uncertainties.errs:
        n, d = [int(v) for v in key.split('/')]
        cached = cached_uncertainties.errs[key]
        new = err(n, d, use_cache=False)
        for iv in range(2):
            diff = new[iv] - cached[iv]
            if diff > eps:
                print('   %20s   %s: diff %.0e    (%f --> %f)' % (key, ['lo', 'hi'][iv], diff, cached[iv], new[iv]))
        # print ' %8.5f  %8.5f   %s' % (cached[0], new[0])
        iline += 1
        # if iline % 1000 == 0:
        #     print '  %d / %d  ok' % (iline, len(cached_uncertainties.errs))
        if iline > 200:
            break

# ----------------------------------------------------------------------------------------
def write():
    with open('tmpcache.py', 'w') as cachefile:
        cachefile.write('errs = {\n')
        # iline = 0
        for key in sorted(cached_uncertainties.errs):
            obs, total = [int(v) for v in key.split('/')]
            lo, hi = err(obs, total, use_cache=False)
            cachefile.write('\'%s\' : (%.12f, %.12f),\n' % (key, lo, hi))
            # iline += 1
            # if iline > 100:
            #     break
        cachefile.write('}\n')
