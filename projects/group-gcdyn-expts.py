#!/usr/bin/env python2
import argparse
import os
import sys
import colored_traceback.always
import numpy
import csv

sys.path.insert(1, './python')
import utils

# group deep learning predictions from --test-file (each of which correspond to the prediction on one tree) into "experiments" of
#  size --n-trees-per-expt, i.e. mimicking the case where we average over N trees in an experiment.
#  Note that this *ignores* values that are leftover after grouping into groups of --n-trees-per-expt

parser = argparse.ArgumentParser()
parser.add_argument('--test-file', required=True)
parser.add_argument('--outfile', required=True)
parser.add_argument('--n-trees-per-expt', required=True, type=int)
args = parser.parse_args()

clines = utils.csvlines(args.test_file)
true_vals = {}
for cln in clines:
    tval, pval = [float(cln[k]) for k in ['Truth', 'Predicted']]
    if tval not in true_vals:
        true_vals[tval] = []
    true_vals[tval].append(pval)
print '    read %d lines with %d different true values: %s' % (len(clines), len(true_vals), ' '.join('%.2f'%v for v in sorted(true_vals)))

final_vals = []
for tval, tlist in true_vals.items():
    for ival, istart in enumerate(range(0, len(tlist), args.n_trees_per_expt)):
        subvals = tlist[istart : istart + args.n_trees_per_expt]
        fline = {'' : ival, 'Truth' : tval, 'Predicted' : numpy.mean(subvals)}
        final_vals.append(fline)
print '      grouped into %d final lines with value counts: %s' % (len(final_vals), '  '.join('%.2f %d'%(v, len([l for l in final_vals if l['Truth']==v])) for v in true_vals))

print '    writing %d lines to %s' % (len(final_vals), args.outfile)
utils.mkdir(args.outfile, isfile=True)
with open(args.outfile, 'w') as ofile:
    writer = csv.DictWriter(ofile, final_vals[0].keys())
    writer.writeheader()
    for fline in final_vals:
        writer.writerow(fline)

