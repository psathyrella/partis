#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import sys
import csv
import argparse
import operator
import argparse
import yaml
import colored_traceback.always
from io import open

import partis.utils as utils

# ----------------------------------------------------------------------------------------
def is_acceptable(scol, acceptable_values, lval):
    if lval in acceptable_values:
        return True
    if args.any_allele and '_gene' in scol and any(utils.are_alleles(g, lval) for g in acceptable_values):
        return True
    return False

class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter)
parser.add_argument('--infname', default='test/ref-results/test/parameters/data/hmm/all-probs.csv', help='input all-probs.csv file from a previously-inferred partis parameter directory, for instance: test/reference-results/test/parameters/data/hmm/all-probs.csv')
parser.add_argument('--config-fname', help='yaml file with info on columns for which we want to specify particular values (and skip others). See default/example set below. To create a yaml config file to start from, uncomment the yaml.dump() line below and rerun with no arguments.')
parser.add_argument('--outfname')
parser.add_argument('--any-allele', action='store_true', help='if set, also include any other alleles of any of the genes specified in \'skip_column_vals\' (note: can also set it in the cfg file).')
parser.add_argument('--debug', action='store_true', default=True)  # it's kind of confusing without the debug printout
args = parser.parse_args()

non_summed_column = None
if args.config_fname is None:
    non_summed_column = 'v_gene'
    skip_column_vals = {  # to input your own dict on the command line, just convert with str() and quote it
        # 'cdr3_length' : ['33', '36', '39', '42', '45', '48'],  # <value> is list of acceptable values NOTE need to all be strings, otherwise you have to worry about converting the values in the csv file
        'v_gene' : ['IGHV1-2*02+G35A', 'IGHV1-2*02+T147C', 'IGHV1-2*02'],
        # 'd_gene' : ['IGHD3-22*01'],
        'j_gene' : ['IGHJ4*02'],
        'cdr3_length' : ['66',],
    }
    print('%s using default skip column/non-summed column values (which probably don\'t correspond to what you\'re actually interested in)' % utils.color('red', 'note'))
    # # uncomment to create a yaml file to start from:
    # with open('tmp.yaml', 'w') as tfile:
    #     yaml.dump({'non_summed_column' : non_summed_column, 'skip_column_vals' : skip_column_vals}, tfile)
else:
    with open(args.config_fname) as yamlfile:
        yamlfo = yaml.load(yamlfile, Loader=yaml.Loader)
    if 'non_summed_column' in yamlfo:
        non_summed_column = yamlfo['non_summed_column']
    skip_column_vals = yamlfo['skip_column_vals']
    for scol in skip_column_vals:
        skip_column_vals[scol] = [str(v) for v in skip_column_vals[scol]]  # yaml.load() converts to integers, which is usually nice, but here we don't want it to since we're not converting when reading all-probs.csv (I think there's options to yaml.load to change this, I just don't want to figure it out now)
    if 'any_allele' in yamlfo:
        if args.any_allele and not yamlfo['any_allele']:  # if it's set to true on the command line, but false in the file
            print(' %s overwriting --any-allele with value from cfg file %s' % (utils.color('red', 'warning'), args.config_fname))
        args.any_allele = yamlfo['any_allele']

info = {}
lines_skipped, lines_used = 0, 0
counts_skipped, counts_used = 0, 0
print('  reading probs from %s' % args.infname)
with open(args.infname) as csvfile:
    reader = csv.DictReader(csvfile)
    # if args.debug:
    #     print '  all columns in file: %s' % ' '.join(reader.fieldnames)
    if len(set(skip_column_vals) - set(reader.fieldnames)) > 0:
        raise Exception('keys in --skip-column-fname not in file: %s' % ' '.join(set(skip_column_vals) - set(reader.fieldnames)))
    for line in reader:
        skip_this_line = False
        for scol, acceptable_values in skip_column_vals.items():
            if not is_acceptable(scol, acceptable_values, line[scol]):
                skip_this_line = True
                lines_skipped += 1
                counts_skipped += int(line['count'])
                break
        if skip_this_line:
            continue

        if non_summed_column is not None:
            if line[non_summed_column] not in info:
                info[line[non_summed_column]] = 0
            info[line[non_summed_column]] += int(line['count'])

        lines_used += 1
        counts_used += int(line['count'])

# ----------------------------------------------------------------------------------------
import partis.fraction_uncertainty as fraction_uncertainty
def frac_err(obs, total):
    lo, hi = fraction_uncertainty.err(obs, total)
    return 0.5 * (hi - lo)
count_fraction = counts_used / float(counts_used + counts_skipped)

if args.debug:
    print('  applied restrictions:%s' % ('     (including all alleles of these genes)' if args.any_allele else ''))
    for scol, acceptable_values in skip_column_vals.items():
        print('      %15s in %s' % (scol, acceptable_values))
    print('   used:')
    print('     %6d / %-6d = %.3f  lines'  % (lines_used, lines_used + lines_skipped, lines_used / float(lines_used + lines_skipped)))
    print('     %6d / %-6d = %.3f +/- %.3f counts'  % (counts_used, counts_used + counts_skipped, count_fraction, frac_err(counts_used, counts_used + counts_skipped)))

    if non_summed_column is not None:
        print('    %18s      count      / %d = fraction' % (non_summed_column, counts_used))
        for val, count in sorted(list(info.items()), key=operator.itemgetter(1), reverse=True):  # sort by counts
        # for val, count in sorted(info.items()):  # sort by column value (e.g. cdr3 length)
            print('   %18s   %6d          %.3f +/- %.3f' % (val, count, count / float(counts_used), frac_err(count, counts_used)))

if args.outfname is not None:
    if args.debug:
        print('  writing total counts (plus %d info entries) to %s' % (len(info), args.outfname))
    with open(args.outfname, 'w') as outfile:
        yamlfo = {'counts' : counts_used,
                  'total' : counts_used + counts_skipped,
                  'fraction' : count_fraction,
                  'frac_err' : frac_err(counts_used, counts_used + counts_skipped),
                  'info' : info}
        yaml.dump(yamlfo, outfile, width=150)
