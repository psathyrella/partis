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
def wkdir(locus):
    return '%s/work/%s' % (args.outdir, locus)

# ----------------------------------------------------------------------------------------
def getifn(locus):
    return '%s/input.tsv' % wkdir(locus)

# ----------------------------------------------------------------------------------------
def scofn(locus, joint=False):
    return '%s/%spartition.tsv' % (wkdir(locus), 'joint-' if joint else '')

# ----------------------------------------------------------------------------------------
def simfn(locus):
    return paircluster.paired_fn(args.simdir, locus, suffix='.yaml')

# ----------------------------------------------------------------------------------------
def getofn(locus, joint=False):
    return paircluster.paired_fn(args.outdir, locus, single_chain=not joint, actstr='partition', suffix='.yaml')

# ----------------------------------------------------------------------------------------
def swfn(locus):
    return '%s/%s/sw-cache.yaml' % (args.indir, locus)

# ----------------------------------------------------------------------------------------
def gloci():  # if no sw file, we probably only made one light locus (or i guess all input is missing, oh well)
    return [l for l in utils.sub_loci(ig_or_tr) if os.path.exists(swfn(l))]

# ----------------------------------------------------------------------------------------
def get_alignments():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(gloci())
    for locus in gloci():
        ofn = getifn(locus)
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue
        utils.mkdir(wkdir(locus))
        cmd = './bin/parse-output.py %s %s --airr-output --skip-columns clone_id' % (swfn(locus), ofn)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : wkdir(locus),
            'workdir' : wkdir(locus),
        }]
    utils.run_scan_cmds(args, cmdfos, 'airr-convert.log', n_total, n_already_there, ofn)

# ----------------------------------------------------------------------------------------
def run_scoper():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, 0
    for locus in gloci():
        if not os.path.exists('%s/%s/sw-cache.yaml' % (args.indir, locus)):
            continue
        ofn = scofn(locus)
        if utils.output_exists(args, ofn): # and not args.dry:  # , offset=8):
            continue

        rcmds = ['library(scoper)']
        rcmds += ['db <- read.csv("%s", sep="\t")' % getifn(locus)]
        rcmds += ['db <- db[ , !(names(db) %in% c("clone_id"))]']  # it crashes if clone_id is already there UPDATE don't need this any more since i'm removing clone_id column above
        rcmds += ['results <- spectralClones(db)']  # , method="novj", germline="germline_alignment_d_mask"
        rcmds += ['df <- as.data.frame(results)']
        rcmds += ['write.table(df[c("sequence_id", "clone_id")], "%s", sep="\t", quote=FALSE, row.names=FALSE)' % ofn]

        utils.run_r(rcmds, wkdir(locus), logfname='%s/scoper.log'%wkdir(locus), print_time='%s scoper'%locus, dryrun=args.dry)

# ----------------------------------------------------------------------------------------
def run_joint_scoper():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, 0
    for locus in gloci():
        if not os.path.exists('%s/%s/sw-cache.yaml' % (args.indir, locus)):
            continue
        ofn = scofn(locus, joint=True)
        if utils.output_exists(args, ofn): # and not args.dry:  # , offset=8):
            continue

        cmd = 'Rscript packages/joint-scoper/scoperClones.R %s %s spectral HL vj 10' % (getifn(locus), ofn)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : wkdir(locus),
            'workdir' : wkdir(locus),
        }]
    utils.run_scan_cmds(args, cmdfos, 'joint-scoper.log', n_total, n_already_there, ofn)

# ----------------------------------------------------------------------------------------
def convert_output(joint=False):
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(gloci())
    for locus in gloci():
        pfn = getofn(locus, joint=joint)
        if utils.output_exists(args, pfn, debug=False):
            n_already_there += 1
            continue
        utils.mkdir(pfn, isfile=True)
        cmd = './bin/parse-output.py %s %s --airr-input --simfname %s' % (scofn(locus, joint=joint), pfn, simfn(locus))
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : pfn,
            'logdir' : wkdir(locus),
            'workdir' : wkdir(locus),
        }]
    utils.run_scan_cmds(args, cmdfos, 'partis-convert.log', n_total, n_already_there, pfn)  # not really much point in parallelizing this since the scoper step isn't parallelized but oh well

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
parser.add_argument('--simdir', required=True)
parser.add_argument('--joint', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--n-max-procs', type=int, help='NOT USED')
args = parser.parse_args()

get_alignments()
run_scoper()
convert_output()
run_joint_scoper()
convert_output(joint=True)
