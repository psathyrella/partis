#!/usr/bin/env python
from subprocess import check_call, Popen
import os

samplesizes = (10000, 5000, 1000, 500, 100, 50)  #, 25)
# samplesizes = (50, 25)
baseplotdir = os.getenv('www') + '/partis'
for size in samplesizes:
    label = 'sample-size-' + str(size)
    cmd = './bin/run-driver.py --label ' + label + ' --plotdir ' + baseplotdir + '/' + label + ' --action plot-performance --simfname _output/test/simu.csv --n-queries ' + str(size)
    # Popen(cmd.split(' '))

plotdirs = [ baseplotdir + '/sample-size-' + str(s) + '/hmm/performance' for s in samplesizes ]
names = [ str(s) for s in samplesizes ]
colors = (595, 807, 834, 632, 596, 1)
linestyles = (1, 1, 1, 1, 1, 2)
plot_cmd = './bin/compare.py --plotdirs ' + ':'.join(plotdirs)
plot_cmd += ' --names ' + ':'.join(names) + ' --outdir ' + baseplotdir + '/check-sample-size'
# plot_cmd += ' --leaves-per-tree ' + ':'.join(['5' for _ in range(len(names))])
plot_cmd += ' --colors ' + ':'.join([str(c) for c in colors])
plot_cmd += ' --linestyles ' + ':'.join([str(l) for l in linestyles])
plot_cmd += ' --normalize --markersize 0 --no-errors --graphify --linewidths 6:5:4:3:2:2'
check_call(plot_cmd.split())
