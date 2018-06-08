#!/usr/bin/env python
import argparse
import sys
import os
import csv

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')
import utils
from hist import Hist
import plotting
import glutils

parser = argparse.ArgumentParser()
parser.add_argument('infile')
parser.add_argument('plotdir')
args = parser.parse_args()

def gk(uids):
    return ':'.join(uids)

glfo = glutils.read_glfo(args.infile.replace('.csv', '-glfo'), locus='igh')

annotations = {}
with open(args.infile) as csvfile:
    reader = csv.DictReader(csvfile)
    for line in reader:
        if line['v_gene'] == '':  # failed (i.e. couldn't find an annotation)
            continue
        utils.process_input_line(line)  # converts strings in the csv file to floats/ints/dicts/etc.
        utils.add_implicit_info(glfo, line)  # add stuff to <line> that's useful, isn't written to the csv since it's redundant
        annotations[gk(line['unique_ids'])] = line

chfo = {uid : utils.get_chimera_max_abs_diff(annotations[uid], iseq=0) for uid in annotations}
biggest_adiffs = sorted(chfo, key=lambda q: chfo[q][1], reverse=True)
for uid in biggest_adiffs[:10]:
    print chfo[uid]
    utils.print_reco_event(annotations[uid])

htmp = Hist(45, 0., 0.65)
for uid in annotations:
    htmp.fill(chfo[uid][1])
utils.prep_dir(args.plotdir, wildlings=['*.svg', '*.csv'])
plotname = 'mfreq-diff'
plotting.draw_no_root(htmp, plotdir=args.plotdir, plotname=plotname, shift_overflows=True, xtitle='abs mfreq diff', ytitle='seqs')
plotting.draw_no_root(htmp, plotdir=args.plotdir, plotname=plotname + '-log', shift_overflows=True, log='y', xtitle='abs mfreq diff', ytitle='seqs')
print 'writing to %s' % args.plotdir
htmp.write('%s/%s.csv' % (args.plotdir, plotname))
