#!/usr/bin/env python
# Checks that different --random-seed values give independent simulations, with one and with several
# processes, and that rerunning a seed reproduces it exactly (see psathyrella/partis#430). Runs
# partis simulate, so it needs the binaries; takes a minute or two.
#   ./test/simulate-seeds.py
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import shutil
import subprocess
import sys
import tempfile

partis_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, partis_dir)
from partis import utils

N_EVENTS = 20

# ----------------------------------------------------------------------------------------
def simulate(workdir, seed, n_procs, label):
    outfname = '%s/%s.yaml' % (workdir, label)
    cmd = '%s/bin/partis simulate --simulate-from-scratch --force-dont-generate-germline-set --n-sim-events %d --n-procs %d --random-seed %d --outfname %s --workdir %s/work-%s' % (partis_dir, N_EVENTS, n_procs, seed, outfname, workdir, label)
    with open('%s/%s.log' % (workdir, label), 'w') as logfile:
        subprocess.check_call(cmd.split(), stdout=logfile, stderr=subprocess.STDOUT)
    _, annotations, _ = utils.read_output(outfname)
    return annotations

# ----------------------------------------------------------------------------------------
def check(label, condition, detail):
    print('    %s  %-62s %s' % (utils.color('green', 'ok ') if condition else utils.color('red', 'FAIL'), label, detail))
    return condition

# ----------------------------------------------------------------------------------------
def run_tests():
    workdir = tempfile.mkdtemp(prefix='simulate-seeds-')
    passed = []
    try:
        for n_procs in [1, 4]:
            print('  --n-procs %d, %d events:' % (n_procs, N_EVENTS))
            naives = {s : [l['naive_seq'] for l in simulate(workdir, s, n_procs, 's%d-p%d' % (s, n_procs))] for s in [1, 2]}
            n_shared = len(set(naives[1]) & set(naives[2]))
            passed.append(check('seeds 1 and 2 share no naive sequences', n_shared == 0, '%d / %d shared' % (n_shared, N_EVENTS)))
        print('  rerun, --n-procs 4:')
        first, second = [simulate(workdir, 1, 4, 'rerun-%d' % i) for i in range(2)]
        same = [l['unique_ids'] for l in first] == [l['unique_ids'] for l in second] and [l['naive_seq'] for l in first] == [l['naive_seq'] for l in second]
        passed.append(check('same seed reproduces the same events and uids', same, ''))
    finally:
        shutil.rmtree(workdir)
    print('  %d/%d checks passed' % (sum(passed), len(passed)))
    return all(passed)

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    print('testing simulation random seeds')
    sys.exit(0 if run_tests() else 1)
