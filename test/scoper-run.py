#!/usr/bin/env python
import numpy
import random
import glob
import time
import colored_traceback.always
import argparse
import subprocess
import sys
import os

sys.path.insert(1, './python')
import utils
import glutils
import paircluster

ig_or_tr = 'ig'

# ----------------------------------------------------------------------------------------
def getifn(args, locus):
    return '%s/%s/input.tsv' % (args.outdir, locus)

# ----------------------------------------------------------------------------------------
def getofn(args, locus):
    return '%s/%s/scoper-partitions.csv' % (args.outdir, locus)

# ----------------------------------------------------------------------------------------
def get_alignments(args, workdir):
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(utils.sub_loci(ig_or_tr))
    for locus in utils.sub_loci(ig_or_tr):
        swfn = '%s/%s/sw-cache.yaml' % (args.indir, locus)
        if not os.path.exists(swfn):
            n_total -= 1  # if no sw file, we probably only made one light locus (or i guess all input is missing, oh well)
            continue
        ofn = getifn(args, locus)
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue
        cmd = './bin/parse-output.py %s %s --airr-output' % (swfn, ofn)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : '%s/%s' % (args.outdir, locus),
            'workdir' : '%s/%s' % (args.outdir, locus),
        }]
    utils.run_scan_cmds(args, cmdfos, 'airr-convert-%s.log'%locus, n_total, n_already_there, ofn)

# ----------------------------------------------------------------------------------------
def run_scoper(args, workdir):
    ofn, cmdfos, n_already_there, n_total = None, [], 0, 0
    for locus in utils.sub_loci(ig_or_tr):
        if not os.path.exists(getifn(args, locus)):
            continue
        ofn = getofn(args, locus)
        if utils.output_exists(args, ofn):  # , offset=8):
            continue

        rcmds = ['library(scoper)']
        rcmds += ['db <- read.csv("%s", sep="\t")' % getifn(args, locus)]
        rcmds += ['db <- db[ , !(names(db) %in% c("clone_id"))]']  # it crashes if clone_id is already there
        rcmds += ['results <- spectralClones(db, method="novj", germline="germline_alignment_d_mask")']
        rcmds += ['df <- as.data.frame(results)']
        rcmds += ['write.csv(df[c("sequence_id", "clone_id")], "%s")' % ofn]

        lwkdir = '%s/%s' % (args.outdir, locus)
        utils.mkdir(lwkdir)
        utils.run_r(rcmds, lwkdir, logfname='%s/%s/scoper.log'%(args.outdir, locus), dryrun=args.dry)

    # ofn = paircluster.paired_fn(args.outdir, locus, actstr='partition')

# ----------------------------------------------------------------------------------------
def install():
    rcmds = ['install.packages(c("Biostrings", "GenomicAlignments", "IRanges", "alakazam", "shazam", "scoper"), repos="http://cran.rstudio.com/", dependencies=TRUE)']  # doesn't work on campanellabox (despite same R version as thneed -- maybe problem is ubuntu version?)
    workdir = '/tmp/%s/%d' % (os.getenv('USER'), random.randint(0,999999))
    os.makedirs(workdir)
    utils.run_r(rcmds, workdir)
    os.rmdir(workdir)

# install()
# sys.exit()

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--indir', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--n-max-procs', type=int)
args = parser.parse_args()

workdir = args.outdir + '/work'

outfname = get_alignments(args, workdir)
run_scoper(args, workdir)
