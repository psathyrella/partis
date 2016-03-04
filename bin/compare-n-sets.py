#!/usr/bin/env python
#!/usr/bin/env python
import csv
import time
import sys
from subprocess import check_output, check_call, Popen
sys.path.insert(1, './python')
import random
import glob
import matplotlib.pyplot as plt
import scipy.special as sps
import numpy as np
import os

from hist import Hist

n_set_list = [1, 2, 5, 7, 10, 15, 20, 25]
baseplotdir = os.getenv('www') + '/partis/n-sets'

def divide_simulation():
    allinfo = {}
    headers = None
    with open('25-leaves.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for line in reader:
            if headers is None:
                headers = line.keys()
            if line['reco_id'] not in allinfo:
                allinfo[line['reco_id']] = []
            allinfo[line['reco_id']].append(line)
    
    for n_to_keep in n_set_list:
        outfile = open('%d-sub-leaves.csv' % n_to_keep, 'w')
        writer = csv.DictWriter(outfile, headers)
        writer.writeheader()
        iline, ievent, iseq = 0, 0, 0
        for recoid in allinfo:
            already_used = set()
            for _ in range(n_to_keep):
                iseq = None
                while iseq is None or iseq in already_used:
                    iseq = random.randint(0, len(allinfo[recoid]) - 1)
                already_used.add(iseq)
                writer.writerow(allinfo[recoid][iseq])

def run_inference():
    base_cmd = './bin/partis.py --action run-viterbi --parameter-dir /fh/fast/matsen_e/dralph/work/partis-dev/_output/021-018/simu-2.3-leaves-1.0-mutate-zipf/hmm --plot-performance --is-simu --slurm'
    for n_set in n_set_list:
        if n_set == 5:
            continue
        n_procs = max(5, int(20. * (n_set / 25.)))
        cmd = base_cmd + ' --n-sets ' + str(n_set) + ' --seqfile ' + str(n_set) + '-sub-leaves.csv --plotdir ' + baseplotdir + '/' + str(n_set) + ' --n-procs ' + str(n_procs) + ' --workdir /fh/fast/matsen_e/dralph/work/partis-dev/_tmp/' + str(random.randint(0, 99999))
        print cmd
        Popen(cmd.split())
        time.sleep(0.5)


# divide_simulation()
# run_inference()
# sys.exit()

# hall = Hist(n_set_list[-1], n_set_list[0] - 0.5, n_set_list[-1] + 0.5)
means = []
for n_set in n_set_list:
    plotdir = baseplotdir + '/' + str(n_set)
    hist = Hist(fname=plotdir + '/hmm/hamming_to_true_naive.csv')
    print '%2d   %.2f' % (n_set, hist.get_mean()),
    # hall.set_ibin(hall.find_bin(n_set), hist.get_mean())
    means.append(hist.get_mean())

import plotting
fig, ax = plotting.mpl_init()
# hall.mpl_plot(ax)
ax.plot(n_set_list, means)
plotting.mpl_finish(ax, baseplotdir, 'means', xlabel='N simultaneous seqs', ylabel='hamming to true naive')
