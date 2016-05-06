#!/usr/bin/env python
import os
import glob
import time
from subprocess import check_call, Popen
import sys

current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
if not os.path.exists(current_script_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % current_script_dir
sys.path.insert(1, current_script_dir)

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
