#!/usr/bin/env python
import csv
import argparse
import operator

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils

# infile = '/fh/fast/matsen_e/processed-data/partis/bf520-synth/v17/BF520-g-merged/hmm/all-probs.csv'
infile = '_output/test_reference-results_test_simu/multi-hmm/all-probs.csv'
outfile = 'processed-all-probs.csv'
non_summed_column = 'v_gene'
skip_columns = {  # columns for which we want to specify particular values (and skip others)
    # 'cdr3_length' : ['39',],  # <value> is list of acceptable values NOTE need to all be strings, otherwise you have to worry about converting the values in the csv file
    # 'v_gene' : ['IGHV1-2*02',],
    'v_gene' : ['IGHV1-2*02+G35A', 'IGHV1-2*02+T147C', 'IGHV1-2*02'],
    # IGHV1-2*02+G35A,2216
    # IGHV1-2*02+T147C,746
    # IGHV1-2*02,336
    # 'd_gene' : ['IGHD3-22*01'],
    'j_gene' : ['IGHJ4*02'],
    'cdr3_length' : ['66',],  #  TGTGCGAGAGGGCCATTCCCGAATTACTATGGTCCGGGGAGTTATTGGGGGGGTTTTGACCACTGG
}


info = {}
lines_skipped, lines_used = 0, 0
counts_skipped, counts_used = 0, 0
with open(infile) as csvfile:
    reader = csv.DictReader(csvfile)
    print '  all columns in file: %s' % ' '.join(reader.fieldnames)
    for line in reader:
        skip_this_line = False
        for scol, acceptable_values in skip_columns.items():
            if line[scol] not in acceptable_values:
                skip_this_line = True
                lines_skipped += 1
                counts_skipped += int(line['count'])
                break
        if skip_this_line:
            continue

        if line[non_summed_column] not in info:
            info[line[non_summed_column]] = 0
        info[line[non_summed_column]] += int(line['count'])
        lines_used += 1
        counts_used += int(line['count'])

print '  applied restrictions:'
for scol, acceptable_values in skip_columns.items():
    print '      %15s in %s' % (scol, acceptable_values)
print '   used:\n     %6d / %-6d  lines\n     %6d / %-6d  counts'  % (lines_used, lines_used + lines_skipped, counts_used, counts_used + counts_skipped)

print '  writing %d lines to %s' % (len(info), outfile)
with open(outfile, 'w') as csvfile:
    writer = csv.DictWriter(csvfile, (non_summed_column, 'count'))
    writer.writeheader()
    for val, count in sorted(info.items(), key=operator.itemgetter(1), reverse=True):
        writer.writerow({non_summed_column : val, 'count' : count})
