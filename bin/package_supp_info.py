#!/usr/bin/env python
import os
import glob
import time
from subprocess import check_call, Popen
import sys
sys.path.insert(1, './python')

from humans import humans

modulo = '10'

basedir = '_output'
tar_cmd = ['tar', 'cz', '-f']

for human in humans['both']:
    print human
    for subset in range(int(modulo)):
        dname = basedir + '/every-' + modulo + '-' + human + '-subset-' + str(subset) + '/data/hmm'
        if not os.path.exists(dname):
            raise Exception('ERROR missing subset ' + str(subset))
        # Popen(tar_cmd + [human + '-subset-' + str(subset) + '.tgz', dname])
        # check_call(['limit_procs', 'tar'])
        # time.sleep(0.2)
        # check_call(['du', '-hs', dname])
