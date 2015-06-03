#!/usr/bin/env python
import sys
import time
import os
from subprocess import Popen, PIPE, check_call
import random

fs = '/fh/fast/matsen_e/dralph/work/partis-dev'
# basesimdir = fs + '/_output/moofit-branch-lengths-full-vdj'
basecmd = './bin/partis.py --action partition --rescale-emissions --initial-cachefname hmm_cached_info-full-vdj.csv'
basecmd += ' --seqfile ' + fs + '/_output/moofit-branch-lengths/simu-full-vdj.csv'  #basesimdir + '/simu.csv'
basecmd += ' --parameter-dir ' + fs + '/_output/moofit-branch-lengths/simu-full-vdj/hmm'  # + basesimdir + '/simu/hmm'

print '\nusing old imgt datadir\n'

procs, rand_ints = [], []
for n_queries in (100,):
# for n_queries in (10, 25, 100, 200, 300, 400, 500, 750, 1000, 2000):
    rand_int = random.randint(0, 99999)
    rand_ints.append(rand_int)
    n_procs = max(1, int(n_queries / 20.))
    cmd = basecmd + ' --outfname ' + fs + '/_output/sample-sizes-full-vdj/' + str(n_queries) + '-' + str(rand_int) + '-queries.csv'
    cmd += ' --datadir data/imgt-old-version'
    cmd += ' --n-max-queries ' + str(n_queries) + ' --n-procs ' + str(n_procs) + ':' + str(max(1, int(n_procs/10.)))
    if n_procs > 2:
        cmd += ' --slurm --workdir ' + fs + '/_tmp/' + str(rand_int)
    print n_queries, n_procs, str(max(1, int(n_procs/10.)))
    cmd = 'time ' + cmd
    print cmd
    sys.exit()
    proc = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
    procs.append(proc)
    time.sleep(0.1)

outdir = fs + '/_tmp/outstreams'
if not os.path.exists(outdir):
    os.makedirs(outdir)
for iproc in range(len(procs)):
    out, err = procs[iproc].communicate()
    with open(outdir + '/' + str(rand_ints[iproc]) + '.out', 'w') as outfile:
        outfile.write(out)
    with open(outdir + '/' + str(rand_ints[iproc]) + '.err', 'w') as errfile:
        errfile.write(err)
