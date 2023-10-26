#!/usr/bin/env python3
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always
import pandas as pd
import glob

# # if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/projects/spisak', '')
# sys.path.insert(1, partis_dir + '/python')
# import utils
sys.path.insert(1, partis_dir + '/packages/HILARy')
from utils import Compatible
from apriori import preprocess, Apriori
from inference import HILARy, CDR3Clustering

# ----------------------------------------------------------------------------------------
# ./projects/mobille-validation.py --version v1 # --actions simplot  --methods mobille:scoper:partis
# ./projects/mobille-validation.py --version trand --actions plot --n-random-seeds 10
parser = argparse.ArgumentParser()
parser.add_argument('--indir', default='/fh/fast/matsen_e/data/briney-great-rep/316188-test')
parser.add_argument('--outdir', required=True)
args = parser.parse_args()

compatible = Compatible()
incols = ['seq_id',
           'chain',
           'productive',
           'v_full',
           'j_full',
           'cdr3_nt',
           'v_start',
           'vdj_nt',
           'isotype']

dfs = []
globstr = '%s/consensus-cdr3nt-90_minimal/*.txt' % args.indir
for fn in glob.glob(globstr):
    print('  reading %s' % fn)
    df = pd.read_csv(fn, usecols=incols)
    dfs.append(compatible.df2airr(df))
if len(dfs) == 0:
    raise Exception('couldn\'t find any files with \'%s\'' % globstr)
df = pd.concat(dfs, ignore_index=True)
df['sequence_id'] = df.index

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
# TODO what's this file for?
idfn = '%s/ids.tsv.gz' % args.outdir
df[['seq_id', 'sequence_id']].to_csv(idfn, sep='\t', index=False)
df.drop('seq_id', axis=1, inplace=True)

tsvfn = '%s/processed-input.tsv.gz' % args.outdir
tsvcols = ['sequence_id',
           'v_call',
           'j_call',
           'junction',
           'v_sequence_alignment',
           'j_sequence_alignment',
           'v_germline_alignment',
           'j_germline_alignment']
df[tsvcols].to_csv(tsvfn, sep='\t', index=False)
df = pd.read_table(tsvfn, usecols=tsvcols)
df = preprocess(df)
ap = Apriori(df)
ap.get_histograms(df.loc[ap.productive])
ap.get_parameters()
ap.get_thresholds()

hilary = HILARy(ap)

prec = CDR3Clustering(ap.classes[hilary.group+['precise_threshold']])
sens = CDR3Clustering(ap.classes[hilary.group+['sensitive_threshold']])
df['precise_cluster'] = prec.infer(df.loc[ap.productive])
df['sensitive_cluster'] = sens.infer(df.loc[ap.productive])

hilary.to_do(df)
df['family'] = hilary.infer(df)


