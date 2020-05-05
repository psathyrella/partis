#!/usr/bin/env python
# I'm not sure if this gets used anywhere, but I'm too chicken to delete it since I just noticed it, but it's been ages since I worked on this stuff
import sys
import os
import argpase
import joblib

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')
import treeutils

parser = argparse.ArgumentParser()
parser.add_argument('fn')
parser.parse_args()
assert '.pickle' in args.fn

with open(args.fn) as picklefile:
    skmodels = joblib.load(picklefile)
write_pmml(args.fn.replace('.pickle', '.pmml'), skmodels, get_dtr_varnames(cg, dtr_cfgvals['vars']), tvar)
