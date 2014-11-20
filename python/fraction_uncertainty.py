import cached_uncertainties
from scipy.stats import beta

# #!/usr/bin/env python

# import utils

# print '{'

# for i in range(100):
#     for j in range(1000):
#         if i < j:
#             lo, hi = fraction_uncertainty(i, j, use_beta=True, use_cache=False)
#             print '\'%d/%d\' : (%f, %f),' % (i, j, lo, hi)

# print '}'

# ----------------------------------------------------------------------------------------
def err(obs, total, use_beta=True, use_cache=True):
    """ Return [some hack-job approximation of] uncertainty on the ratio n / total """
    assert obs <= total
    if total == 0.0:
        return 0.0  # TODO wait, shouldn't this be a tuple?

    key = str(obs) + '/' + str(total)
    if use_cache and key in cached_uncertainties.errs:
        return cached_uncertainties.errs[key]

    frac = float(obs) / total
    if use_beta or frac == 0.0 or frac == 1.0:  # still need to use beta for 0 and 1 'cause the sqrt thing below gives garbage for those cases
        lo = beta.ppf(1./6, 1 + obs, 1 + total - obs)
        hi = beta.ppf(1. - 1./6, 1 + obs, 1 + total - obs)
        if frac < lo:  # if k/n very small (probably zero), take a one-sided c.i. with 2/3 the mass
            lo = 0.
            hi = beta.ppf(2./3, 1 + obs, 1 + total - obs)
        if frac > hi:  # same deal if k/n very large (probably one)
            lo = beta.ppf(1./3, 1 + obs, 1 + total - obs)
            hi = 1.
    else:  # square root shenaniganery
        err = (1./(total*total)) * (math.sqrt(obs)*total - obs*math.sqrt(total))
        lo = frac - err
        hi = frac + err

    assert lo < frac or frac == 0.0
    assert frac < hi or frac == 1.0
    return (lo,hi)
