#!/usr/bin/env python
import argparse
import sys
import csv

sys.path.insert(1, './python')
import utils
from opener import opener

# ----------------------------------------------------------------------------------------
modulo = 100  # for the moment hard-code that we want every hundredth seq to go in each file
# ----------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument('--infname', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--start-indices', required=True)  # colon-separated list of start indices. E.g. with '0:1:2' we will write three output files. The first seq line in <infname> goes to 0, the next to 1, the third to 2, and then we skip 97 seqs, then yadda yadda
args = parser.parse_args()
args.start_indices = utils.get_arg_list(args.start_indices, intify=True)

infile = opener('r')(args.infname)
reader = csv.DictReader(infile, delimiter='\t')

utils.prep_dir(args.outdir, '*.tsv.bz2')
outfiles, writers = {}, {}
for iout in args.start_indices:
    outfname = args.outdir + ('/every-hundredth-subset-%d.tsv.bz2' % iout)
    outfiles[iout] = opener('w')(outfname)
    writers[iout] = csv.DictWriter(outfiles[iout], reader.fieldnames, delimiter='\t')
    writers[iout].writeheader()

iline = 0
for line in reader:
    if iline % 100000 == 0:
        print 'iline %d' % iline
    remainder = iline % modulo
    if remainder in args.start_indices:
        # print 'write %d to %d' % (iline, remainder)
        writers[remainder].writerow(line)
    # else:
        # print 'skip %d' % iline
    iline += 1
    # if iline > 100:
    #     sys.exit()

for iout in outfiles:
    outfiles[iout].close()
