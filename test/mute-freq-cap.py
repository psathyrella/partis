#!/usr/bin/env python
# Checks the cap that paramutils.read_mute_freqs_with_weights() puts on the substitution it makes for
# a position with zero observed mutations (see psathyrella/partis#411). Builds a parameter dir with
# chosen per-base counts rather than running partis, so it takes a second and doesn't need binaries.
#   ./test/mute-freq-cap.py
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import csv
import os
import shutil
import sys
import tempfile

sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)) + '/..')
from partis import paramutils
from partis import fraction_uncertainty
from partis import utils

GENE = 'IGHV3-23*01'
N_POSITIONS = 20

# ----------------------------------------------------------------------------------------
def write_parameter_dir(pdir, n_obs, sample_mute_freq, mute_freq=0.):
    # one gene, <N_POSITIONS> positions, each with <n_obs> observations and measured <mute_freq>
    utils.mkdir(pdir + '/mute-freqs')
    lo_err, hi_err = fraction_uncertainty.err(int(round(mute_freq * n_obs)), n_obs)
    germline_obs = n_obs - int(round(mute_freq * n_obs))
    with open('%s/mute-freqs/%s.csv' % (pdir, utils.sanitize_name(GENE)), utils.csv_wmode()) as cfile:
        cols = ['position', 'mute_freq', 'lo_err', 'hi_err'] + [c for n in utils.nukes for c in (n, n+'_obs', n+'_lo_err', n+'_hi_err')]
        writer = csv.DictWriter(cfile, cols)
        writer.writeheader()
        for pos in range(N_POSITIONS):
            row = {'position' : pos, 'mute_freq' : mute_freq, 'lo_err' : lo_err, 'hi_err' : hi_err}
            for inuke, nuke in enumerate(utils.nukes):
                obs = germline_obs if inuke == 0 else (n_obs - germline_obs if inuke == 1 else 0)
                row.update({nuke : obs / float(n_obs), nuke+'_obs' : obs, nuke+'_lo_err' : 0., nuke+'_hi_err' : 0.})
            writer.writerow(row)
    # the sample-wide hist read_mute_freqs_with_weights() takes its cap from: one bin holding
    # everything, centred on <sample_mute_freq>, so Hist.get_mean() returns that value
    half = 0.5 * sample_mute_freq
    with open('%s/all-mean-mute-freqs.csv' % pdir, utils.csv_wmode()) as hfile:
        writer = csv.DictWriter(hfile, ['bin_low_edge', 'contents', 'binlabel', 'error'])
        writer.writeheader()
        for low_edge, contents in [(-3 * half, 0.), (-half, 0.), (half, 1.), (3 * half, 0.)]:
            writer.writerow({'bin_low_edge' : low_edge, 'contents' : contents, 'binlabel' : '', 'error' : 0.})

# ----------------------------------------------------------------------------------------
def mute_freq_at(n_obs, sample_mute_freq, mute_freq=0.):
    workdir = tempfile.mkdtemp(prefix='mute-freq-cap-')
    try:
        write_parameter_dir(workdir, n_obs, sample_mute_freq, mute_freq=mute_freq)
        return paramutils.read_mute_freqs_with_weights(workdir, [GENE])[0]
    finally:
        shutil.rmtree(workdir)

# ----------------------------------------------------------------------------------------
def check(label, condition, detail):
    print('    %s  %-58s %s' % (utils.color('green', 'ok ') if condition else utils.color('red', 'FAIL'), label, detail))
    return condition

# ----------------------------------------------------------------------------------------
def run_tests():
    sample_mute_freq = 0.05
    passed = []

    # the substitution for zero observed mutations is half the width of the obs=1 interval, ~1.2/n,
    # so it grows as observations run out. it must never come back above the sample's mute freq.
    print('  zero observed mutations, sample mute freq %.3f:' % sample_mute_freq)
    print('    %-4s %12s %12s' % ('n', 'uncapped', 'returned'))
    for n_obs in [2, 5, 10, 20, 50, 120, 500, 2000]:
        lo_err, hi_err = fraction_uncertainty.err(0, n_obs)
        uncapped = 0.5 * (lo_err + hi_err)
        returned = mute_freq_at(n_obs, sample_mute_freq)
        print('    %-4d %12.5f %12.5f' % (n_obs, uncapped, returned))
        passed.append(check('n=%d: at or below the sample mute freq' % n_obs, returned <= sample_mute_freq + utils.eps,
                            '%.5f <= %.5f' % (returned, sample_mute_freq)))
        if uncapped > sample_mute_freq:  # cap binds
            passed.append(check('n=%d: cap binds, so the sample mute freq comes back' % n_obs, abs(returned - sample_mute_freq) < utils.eps,
                                '%.5f' % returned))
        else:  # cap doesn't bind, so the substitution is untouched
            passed.append(check('n=%d: cap does not bind, substitution untouched' % n_obs, abs(returned - uncapped) < utils.eps,
                                '%.5f' % returned))

    # a position with mutations measured is not a zero-count position, so the cap must not touch it,
    # even when the measured value is above the sample's mute freq
    print('  measured mutations, which the cap must leave alone:')
    for n_obs, mute_freq in [(100, 0.2), (100, 0.02), (20, 0.5)]:
        returned = mute_freq_at(n_obs, sample_mute_freq, mute_freq=mute_freq)
        passed.append(check('n=%d, measured %.2f: passes through' % (n_obs, mute_freq), abs(returned - mute_freq) < utils.eps,
                            '%.5f' % returned))

    # with the cap, a smaller sample can no longer push a zero-count position higher without limit
    print('  monotonicity in n:')
    returned = [mute_freq_at(n, sample_mute_freq) for n in [2, 10, 50, 500]]
    passed.append(check('never exceeds the sample mute freq at any n', max(returned) <= sample_mute_freq + utils.eps,
                        'max %.5f' % max(returned)))
    passed.append(check('non-increasing as n grows', all(a >= b - utils.eps for a, b in zip(returned, returned[1:])),
                        ' '.join('%.5f' % f for f in returned)))

    print('  %d/%d checks passed' % (sum(passed), len(passed)))
    return all(passed)

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    print('testing the zero-count mute freq cap (%s)' % GENE)
    sys.exit(0 if run_tests() else 1)
