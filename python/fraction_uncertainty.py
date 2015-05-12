import cached_uncertainties
from scipy.stats import beta

# ----------------------------------------------------------------------------------------
def err(obs, total, use_beta=True, use_cache=True, for_paper=False):
    """ Return uncertainty on the ratio n / total """
    assert obs <= total
    if total == 0.0:
        return (0.0, 0.0, False)

    key = str(obs) + '/' + str(total)
    if use_cache and key in cached_uncertainties.errs:
        return cached_uncertainties.errs[key] + (True, )

    frac = float(obs) / total
    if use_beta or frac == 0.0 or frac == 1.0:  # still need to use beta for 0 and 1 'cause the sqrt thing below gives garbage for those cases
        if for_paper:  # total volume of confidence interval
            vol = 0.95  # use 95% for paper
            cpr = 0.5  # constant from prior (jeffreys for paper)
        else:
            vol = 2./3  # otherwise +/- 1 sigma
            cpr = 1.  # constant prior
        lo = beta.ppf((1. - vol)/2, cpr + obs, cpr + total - obs)
        hi = beta.ppf((1. + vol)/2, cpr + obs, cpr + total - obs)
        if frac < lo:  # if k/n very small (probably zero), take a one-sided c.i. with 2/3 (0.95) the mass
            lo = 0.
            hi = beta.ppf(vol, cpr + obs, cpr + total - obs)
        if frac > hi:  # same deal if k/n very large (probably one)
            lo = beta.ppf(1. - vol, cpr + obs, cpr + total - obs)
            hi = 1.
    else:  # square root shenaniganery
        err = (1./(total*total)) * (math.sqrt(obs)*total - obs*math.sqrt(total))
        lo = frac - err
        hi = frac + err

    assert lo < frac or frac == 0.0
    assert frac < hi or frac == 1.0
    return (lo, hi) + (False, )
