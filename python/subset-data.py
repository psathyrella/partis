#!/usr/bin/env python
import argparse
import sys
import csv

sys.path.insert(1, './python')
import utils
from seqfileopener import get_seqfile_info
from opener import opener

parser = argparse.ArgumentParser()
parser.add_argument('--infname', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--start-indices', required=True)  # colon-separated list of start indices. E.g. with '0:1:2' we will write three output files. The first seq line in <infname> goes to 0, the next to 1, the third to 2, and then we skip 97 seqs, then yadda yadda
parser.add_argument('--modulo', type=int, default=100)
args = parser.parse_args()
args.start_indices = utils.get_arg_list(args.start_indices, intify=True)

print 'subsetting %s: every %d th sequence' % (args.infname, args.modulo)

infile = opener('r')(args.infname)
input_info, _ = get_seqfile_info(args.infname, is_data=True)  #, n_max_queries=1000)
for key, d in input_info.items():  # get field names (they should be the same for each row, this just grabs the first one)
    fieldnames = d.keys()
    break

utils.prep_dir(args.outdir)  #, '*.bz2')
outfiles, writers = {}, {}
for iout in args.start_indices:
    outfname = args.outdir + ('/every-' + str(args.modulo) + '-subset-%d.csv.bz2' % iout)
    outfiles[iout] = opener('w')(outfname)
    writers[iout] = csv.DictWriter(outfiles[iout], fieldnames, delimiter=',')
    writers[iout].writeheader()

iline = 0
n_written = 0
for line in input_info.values():
    if iline % 100000 == 0:
        print 'iline %d' % iline
    remainder = iline % args.modulo
    if remainder in args.start_indices:
        # print 'write %d to %d' % (iline, remainder)
        writers[remainder].writerow(line)
        n_written += 1
    # else:
        # print 'skip %d' % iline
    iline += 1
    # if iline > 100:
    #     sys.exit()

print 'wrote %d / %d lines' % (n_written, iline)
for iout in outfiles:
    outfiles[iout].close()
