#!/usr/bin/env python
# Checks the prior that hmmwriter.HmmWriter.process_mutation_info() uses for positions with too few
# observations and for positions missing from the counts file: it should be the sample's own mean
# mutation frequency, held within the per-position bounds, rather than a fixed 0.1. Copies the test
# parameter dir (only as a source of well-formed input files), replaces one gene's per-position counts
# and the sample-wide mean with chosen values, and builds that gene's HmmWriter, so it takes a few
# seconds and doesn't need binaries. Also checks that the light chain dummy d keeps its fixed 0.1.
#   ./test/mute-freq-prior.py
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import csv
import os
import shutil
import sys
import tempfile

partis_dir = os.path.dirname(os.path.abspath(__file__)).replace('/test', '')
sys.path.insert(1, partis_dir)
from partis import glutils
from partis import fraction_uncertainty
from partis import hmmwriter
from partis import utils

LOCUS = 'igh'
GENE = 'IGHV3-23*01'
SOURCE_PDIR = partis_dir + '/test/ref-results/test/parameters/data/hmm'
LIGHT_LOCUS = 'igk'
LIGHT_SOURCE_PDIR = partis_dir + '/test/paired/ref-results/test/parameters/data/igk/hmm'
DUMMY_D_MUTE_FREQ = 0.1  # what the dummy d has always had
MIN_OBS = 20  # default --min-observations-per-gene
THIN_POSITIONS = list(range(100, 120))  # one observation, no mutations
MISSING_POSITIONS = list(range(120, 140))  # no row in the counts file
N_WELL_COVERED, WELL_COVERED_MUTE_FREQ = 100, 0.03  # every other position

# ----------------------------------------------------------------------------------------
def write_counts(pdir, germline_seq, sample_mute_freq):
    with open('%s/mute-freqs/%s.csv' % (pdir, utils.sanitize_name(GENE)), utils.csv_wmode()) as cfile:
        cols = ['position', 'mute_freq', 'lo_err', 'hi_err'] + [c for n in utils.nukes for c in (n, n+'_obs', n+'_lo_err', n+'_hi_err')]
        writer = csv.DictWriter(cfile, cols)
        writer.writeheader()
        for pos, gl_nuke in enumerate(germline_seq):
            if pos in MISSING_POSITIONS:
                continue
            n_obs, mute_freq = (1, 0.) if pos in THIN_POSITIONS else (N_WELL_COVERED, WELL_COVERED_MUTE_FREQ)
            n_mutated = int(round(mute_freq * n_obs))
            lo_err, hi_err = fraction_uncertainty.err(n_mutated, n_obs)
            mut_nuke = [n for n in utils.nukes if n != gl_nuke][0]
            row = {'position' : pos, 'mute_freq' : mute_freq, 'lo_err' : lo_err, 'hi_err' : hi_err}
            for nuke in utils.nukes:
                obs = n_obs - n_mutated if nuke == gl_nuke else (n_mutated if nuke == mut_nuke else 0)
                row.update({nuke : obs / float(n_obs), nuke+'_obs' : obs, nuke+'_lo_err' : 0., nuke+'_hi_err' : 0.})
            writer.writerow(row)
    write_sample_mean(pdir, sample_mute_freq)

# ----------------------------------------------------------------------------------------
def write_sample_mean(pdir, sample_mute_freq):
    # sample-wide mean: one bin holding everything, centred on <sample_mute_freq>, so Hist.get_mean() returns that value
    half = 0.5 * sample_mute_freq
    with open('%s/all-mean-mute-freqs.csv' % pdir, utils.csv_wmode()) as hfile:
        writer = csv.DictWriter(hfile, ['bin_low_edge', 'contents', 'binlabel', 'error'])
        writer.writeheader()
        for low_edge, contents in [(-3 * half, 0.), (-half, 0.), (half, 1.), (3 * half, 0.)]:
            writer.writerow({'bin_low_edge' : low_edge, 'contents' : contents, 'binlabel' : '', 'error' : 0.})

# ----------------------------------------------------------------------------------------
def build_writer(sample_mute_freq):
    workdir = tempfile.mkdtemp(prefix='mute-freq-prior-')
    try:
        pdir = workdir + '/hmm'
        shutil.copytree(SOURCE_PDIR, pdir)
        glfo = glutils.read_glfo(pdir + '/germline-sets', LOCUS)
        write_counts(pdir, glfo['seqs']['v'][GENE], sample_mute_freq)
        args = argparse.Namespace(locus=LOCUS, min_observations_per_gene=MIN_OBS, allow_conserved_codon_deletion=False, no_per_base_mfreqs=False)
        return hmmwriter.HmmWriter(pdir, workdir + '/hmms', GENE, glfo, args)
    finally:
        shutil.rmtree(workdir)

# ----------------------------------------------------------------------------------------
def build_dummy_d_writer(sample_mute_freq):  # light chain dummy d, whose single base isn't a real germline position
    workdir = tempfile.mkdtemp(prefix='mute-freq-prior-')
    try:
        pdir = workdir + '/hmm'
        shutil.copytree(LIGHT_SOURCE_PDIR, pdir)
        write_sample_mean(pdir, sample_mute_freq)
        glfo = glutils.read_glfo(pdir + '/germline-sets', LIGHT_LOCUS)
        args = argparse.Namespace(locus=LIGHT_LOCUS, min_observations_per_gene=MIN_OBS, allow_conserved_codon_deletion=False, no_per_base_mfreqs=False)
        return hmmwriter.HmmWriter(pdir, workdir + '/hmms', glutils.dummy_d_genes[LIGHT_LOCUS], glfo, args)
    finally:
        shutil.rmtree(workdir)

# ----------------------------------------------------------------------------------------
def check(label, condition, detail):
    print('    %s  %-62s %s' % (utils.color('green', 'ok ') if condition else utils.color('red', 'FAIL'), label, detail))
    return condition

# ----------------------------------------------------------------------------------------
def run_tests():
    passed = []
    # naive (below the lower bound), typical, heavily mutated, and above the upper bound
    for sample_mute_freq in [0.001, 0.005, 0.05, 0.2, 0.5]:
        writer = build_writer(sample_mute_freq)
        bounds = writer.mute_freq_bounds
        prior = min(bounds['hi'], max(bounds['lo'], sample_mute_freq))
        # a thin position's own value is the #411 substitution for zero mutations in one observation, capped at the sample mean, then
        # held within the bounds; it gets weight 1 of <MIN_OBS> in the blend, and the prior gets the rest
        lo_err, hi_err = fraction_uncertainty.err(0, 1)
        own = min(bounds['hi'], max(bounds['lo'], min(0.5 * (lo_err + hi_err), sample_mute_freq)))
        expected = {'thin' : ((MIN_OBS - 1) * prior + own) / float(MIN_OBS), 'missing' : prior}
        print('  sample mean mute freq %.3f (expected prior %.3f):' % (sample_mute_freq, prior))
        for label, positions in [('thin', THIN_POSITIONS), ('missing', MISSING_POSITIONS)]:
            vals = [writer.mute_freqs[p] for p in positions]
            passed.append(check('%s positions: %.5f' % (label, expected[label]), all(abs(v - expected[label]) < utils.eps for v in vals),
                                'min %.5f max %.5f' % (min(vals), max(vals))))
        well_covered = [writer.mute_freqs[p] for p in range(60, 80)]
        passed.append(check('well-covered positions keep their measured value', all(abs(v - WELL_COVERED_MUTE_FREQ) < utils.eps for v in well_covered),
                            'min %.5f max %.5f' % (min(well_covered), max(well_covered))))
        passed.append(check('hmm still records the unclamped sample mean', abs(writer.hmm.extras['overall_mute_freq'] - sample_mute_freq) < utils.eps,
                            '%.5f' % writer.hmm.extras['overall_mute_freq']))

    # the light chain dummy d keeps its fixed value, since its single base is a placeholder, not something the sample's mean describes
    print('  light chain dummy d (%s):' % glutils.dummy_d_genes[LIGHT_LOCUS])
    for sample_mute_freq in [0.001, 0.005, 0.05, 0.2, 0.5]:
        writer = build_dummy_d_writer(sample_mute_freq)
        passed.append(check('sample mean %.3f: dummy d stays at %.2f' % (sample_mute_freq, DUMMY_D_MUTE_FREQ), abs(writer.mute_freqs[0] - DUMMY_D_MUTE_FREQ) < utils.eps,
                            '%.5f' % writer.mute_freqs[0]))

    print('  %d/%d checks passed' % (sum(passed), len(passed)))
    return all(passed)

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    print('testing the prior mute freq for thin and missing positions (%s)' % GENE)
    sys.exit(0 if run_tests() else 1)
