#!/usr/bin/env python
# Checks the cap that paramutils.read_mute_freqs_with_weights() puts on the substitution it makes for
# a position with zero observed mutations (see psathyrella/partis#411). Builds a parameter dir with
# chosen per-base counts rather than running partis, so it takes a second and doesn't need binaries.
# Also checks the order in which utils.get_mean_mfreq() tries its sources.
#   ./test/mute-freq-cap.py
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import contextlib
import csv
import io
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
def write_parameter_dir(pdir, n_obs, sample_mute_freq, mute_freq=0., exact_mean=True):
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
    # exact mean <sample_mute_freq>, binned hist (one bin) at twice that, so the cap shows which was used
    # without <exact_mean>, only the binned hist, at <sample_mute_freq> (an older parameter dir)
    if exact_mean:
        with open('%s/%s' % (pdir, utils.mute_freq_sum_fname), utils.csv_wmode()) as sfile:
            writer = csv.DictWriter(sfile, ['region', 'sum', 'n_seqs'])
            writer.writeheader()
            writer.writerow({'region' : 'all', 'sum' : 1000 * sample_mute_freq, 'n_seqs' : 1000})
    half = sample_mute_freq if exact_mean else 0.5 * sample_mute_freq
    with open('%s/all-mean-mute-freqs.csv' % pdir, utils.csv_wmode()) as hfile:
        writer = csv.DictWriter(hfile, ['bin_low_edge', 'contents', 'binlabel', 'error'])
        writer.writeheader()
        for low_edge, contents in [(-3 * half, 0.), (-half, 0.), (half, 1.), (3 * half, 0.)]:
            writer.writerow({'bin_low_edge' : low_edge, 'contents' : contents, 'binlabel' : '', 'error' : 0.})

# ----------------------------------------------------------------------------------------
def mute_freq_at(n_obs, sample_mute_freq, mute_freq=0., exact_mean=True):
    workdir = tempfile.mkdtemp(prefix='mute-freq-cap-')
    try:
        write_parameter_dir(workdir, n_obs, sample_mute_freq, mute_freq=mute_freq, exact_mean=exact_mean)
        return paramutils.read_mute_freqs_with_weights(workdir, [GENE])[0]
    finally:
        shutil.rmtree(workdir)

# ----------------------------------------------------------------------------------------
def write_hist(fname, low_edges, contents):
    # <low_edges> and <contents> include the under- and overflow rows
    with open(fname, utils.csv_wmode()) as hfile:
        writer = csv.DictWriter(hfile, ['bin_low_edge', 'contents', 'binlabel', 'error'])
        writer.writeheader()
        for low_edge, content in zip(low_edges, contents):
            writer.writerow({'bin_low_edge' : low_edge, 'contents' : content, 'binlabel' : '', 'error' : 0.})

# ----------------------------------------------------------------------------------------
def mean_mfreq_from(exact=None, mfreq_contents=None, n_muted=None, n_calls=1, sums=None, **kwargs):
    # get_mean_mfreq() on a dir holding only the given sources, plus how many times it warned over <n_calls>
    # <sums>: raw contents of the sums file
    workdir = tempfile.mkdtemp(prefix='mute-freq-cap-')
    try:
        if exact is not None:
            with open('%s/%s' % (workdir, utils.mute_freq_sum_fname), utils.csv_wmode()) as sfile:
                writer = csv.DictWriter(sfile, ['region', 'sum', 'n_seqs'])
                writer.writeheader()
                writer.writerow({'region' : 'all', 'sum' : 1000 * exact, 'n_seqs' : 1000})
        if sums is not None:
            with open('%s/%s' % (workdir, utils.mute_freq_sum_fname), utils.csv_wmode()) as sfile:
                sfile.write(sums)
        if mfreq_contents is not None:  # three bins of width 0.04
            write_hist(workdir + '/all-mean-mute-freqs.csv', [-0.04, 0., 0.04, 0.08, 0.12], [0.] + mfreq_contents + [0.])
        if n_muted is not None:  # one bin for each of 0 to 10 mutations
            write_hist(workdir + '/all-mean-n-muted.csv', [n - 1.5 for n in range(13)], [0.] + [n_muted.get(n, 0.) for n in range(11)] + [0.])
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            returned = [utils.get_mean_mfreq(workdir, **kwargs) for _ in range(n_calls)]
        assert len(set(returned)) == 1
        return returned[0], out.getvalue().count('using binned mean')
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

    # a parameter dir without the exact mean's files falls back to the binned hist mean
    print('  older parameter dir, binned hist only:')
    returned = mute_freq_at(2, sample_mute_freq, exact_mean=False)  # n where the cap binds
    passed.append(check('cap falls back to the binned mean', abs(returned - sample_mute_freq) < utils.eps, '%.5f' % returned))

    # get_mean_mfreq() tries the exact sum, then the estimate from the two hists, then the binned mean
    # half the sequences have 1 of 100 positions mutated, half have 10: exact mean 0.055, binned mean 0.06
    print('  sources for the mean mute freq:')
    mfreq_contents, n_muted = [0.5, 0., 0.5], {1 : 0.5, 10 : 0.5}
    returned, n_warn = mean_mfreq_from(exact=0.03, mfreq_contents=mfreq_contents, n_muted=n_muted)
    passed.append(check('exact sum wins over the hists', abs(returned - 0.03) < utils.eps and n_warn == 0, '%.5f, %d warnings' % (returned, n_warn)))
    returned, n_warn = mean_mfreq_from(mfreq_contents=mfreq_contents, n_muted=n_muted)
    passed.append(check('no exact sum: estimate from the hists', abs(returned - 0.055) < 1e-5 and n_warn == 0, '%.5f, %d warnings' % (returned, n_warn)))
    returned, n_warn = mean_mfreq_from(mfreq_contents=mfreq_contents, n_calls=3)
    passed.append(check('no n-muted hist: binned mean, one warning', abs(returned - 0.06) < utils.eps and n_warn == 1, '%.5f, %d warnings in 3 calls' % (returned, n_warn)))
    returned, n_warn = mean_mfreq_from(mfreq_contents=[1., 0., 0.], n_muted={1 : 1.}, n_calls=3)  # nothing above the first bin
    passed.append(check('hists cannot tell: binned mean, one warning', abs(returned - 0.02) < utils.eps and n_warn == 1, '%.5f, %d warnings in 3 calls' % (returned, n_warn)))
    returned, n_warn = mean_mfreq_from(mfreq_contents=mfreq_contents, binned_fallback=False)
    passed.append(check('binned fallback off: None, no warning', returned is None and n_warn == 0, '%s, %d warnings' % (returned, n_warn)))

    # a sums file that cannot be read, or has no row for the region, raises
    print('  sums file that exists but cannot be used:')
    header = 'region,sum,n_seqs\n'
    for label, sums, region in [('empty', '', 'all'), ('header only', header, 'all'), ('no row for the region', header + 'v,30.0,1000\n', 'all'), ('row for all, v asked for', header + 'all,30.0,1000\n', 'v'),
                                ('no region column', 'sum,n_seqs\n30.0,1000\n', 'all'), ('no n_seqs column', 'region,sum\nall,30.0\n', 'all'), ('row cut off', header + 'all,30.0', 'all'),
                                ('sum not a number', header + 'all,abc,1000\n', 'all'), ('n_seqs not a number', header + 'all,30.0,xyz\n', 'all'), ('bad row for another region', header + 'all,30.0,1000\nv,abc,1000\n', 'all')]:
        try:
            mean_mfreq_from(sums=sums, mfreq_contents=mfreq_contents, n_muted=n_muted, region=region)
            errstr = 'no exception'
        except Exception as err:
            errstr = '%s: %s' % (type(err).__name__, err)
        passed.append(check('%s: raises, naming the file' % label, errstr.startswith('Exception: ') and utils.mute_freq_sum_fname in errstr, errstr.split('/')[-1][:60]))
    returned, n_warn = mean_mfreq_from(sums=header + 'all,0.0,0\n', mfreq_contents=mfreq_contents, n_muted=n_muted)
    passed.append(check('n_seqs 0: estimate from the hists', abs(returned - 0.055) < 1e-5 and n_warn == 0, '%.5f, %d warnings' % (returned, n_warn)))

    print('  %d/%d checks passed' % (sum(passed), len(passed)))
    return all(passed)

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    print('testing the zero-count mute freq cap (%s)' % GENE)
    sys.exit(0 if run_tests() else 1)
