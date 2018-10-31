#!/usr/bin/env python
import collections
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
parser.add_argument('--chunk-len', default=75, type=int)
parser.add_argument('--title')
args = parser.parse_args()
if args.title == 'good':
    args.title = 'none'
elif args.title == 'chimeras':
    args.title = 'all chimeras'

def gk(uids):
    return ':'.join(uids)


glfo, annotation_list, _ = utils.read_output(args.infile)
annotations = collections.OrderedDict((line['unique_ids'][0], line) for line in annotation_list)

chfo = {uid : {k : v for k, v in zip(('imax', 'max_abs_diff'), utils.get_chimera_max_abs_diff(annotations[uid], iseq=0, chunk_len=args.chunk_len))} for uid in annotations}
biggest_adiffs = sorted(chfo, key=lambda q: chfo[q]['max_abs_diff'], reverse=True)
for uid in biggest_adiffs[:5]:
    print '%-3d  %6.3f' % (chfo[uid]['imax'], chfo[uid]['max_abs_diff'])
    utils.print_reco_event(annotations[uid])

hmaxval = Hist(45, 0., 0.65)
for uid in annotations:
    hmaxval.fill(chfo[uid]['max_abs_diff'])
himax = Hist(75, 0., 400)
for uid in annotations:
    himax.fill(chfo[uid]['imax'])

utils.prep_dir(args.plotdir, wildlings=['*.svg', '*.csv'])

import matplotlib
from matplotlib import pyplot as plt
fig, ax = plotting.mpl_init()
xvals, yvals = zip(*[(v['imax'], v['max_abs_diff']) for v in chfo.values()])
plt.scatter(xvals, yvals, alpha=0.4)

print 'writing to %s' % args.plotdir
plotting.mpl_finish(ax, args.plotdir, 'hexbin', title=args.title, xlabel='break point', ylabel='abs mfreq diff')

plotting.draw_no_root(hmaxval, plotdir=args.plotdir, plotname='mfreq-diff', shift_overflows=True, xtitle='abs mfreq diff', ytitle='seqs')
hmaxval.write('%s/%s.csv' % (args.plotdir, 'mfreq-diff'))

plotting.draw_no_root(himax, plotdir=args.plotdir, plotname='imax', shift_overflows=True, xtitle='break point', ytitle='seqs')
himax.write('%s/%s.csv' % (args.plotdir, 'imax'))
