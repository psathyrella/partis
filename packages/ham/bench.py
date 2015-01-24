#!/usr/bin/env python

from subprocess import check_output
import numpy
import math
from scipy import stats
nruns = 30
times, mems = [], []
for _ in range(nruns):
    outlist = check_output('memtime ./bench --hmmfname examples/casino.yaml --seqs test-10m.txt 2>&1', shell=True).split()
    time = float(outlist[outlist.index('sec') - 1])
    mem = float(outlist[outlist.index('kB') - 1])
    print '    %.2f sec %.2f MB' % (time, 0.001*mem)  # convert mem to MB
    times.append(time)
    mems.append(0.001*mem)

print '%.2f (+/- %.2f) sec %.2f (+/- %.2f) MB' % (numpy.mean(times), stats.sem(times), numpy.mean(mems), stats.sem(mems))
