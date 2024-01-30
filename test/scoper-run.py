#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import numpy
import random
import glob
import time
import colored_traceback.always
import argparse
import subprocess
import sys
import os
import glob

sys.path.insert(1, './python')
import utils
import glutils
import paircluster

ig_or_tr = 'ig'

# ----------------------------------------------------------------------------------------
def wkdir(locus):
    return '%s/work%s' % (args.outdir, '' if locus is None else '/'+locus)

# ----------------------------------------------------------------------------------------
def getifn(locus):
    return '%s/input.tsv' % wkdir(locus)

# ----------------------------------------------------------------------------------------
def scofn(locus):
    return '%s/partition.tsv' % wkdir(locus)

# ----------------------------------------------------------------------------------------
def glfd(locus):
    return '%s/germlines' % wkdir(locus)

# ----------------------------------------------------------------------------------------
def simfn(locus):
    if args.single_chain:
        return args.simdir
    return paircluster.paired_fn(args.simdir, locus, suffix='.yaml')

# ----------------------------------------------------------------------------------------
def getofn(locus, joint=False):
    if args.single_chain:
        return '%s/partition.yaml' % args.outdir
    return paircluster.paired_fn(args.outdir, locus, single_chain=not joint, actstr='partition', suffix='.yaml')

# ----------------------------------------------------------------------------------------
def swfn(locus):
    swd = '%s%s' % (args.indir, '' if locus is None else '/'+locus)
    if os.path.exists('%s/sw-cache.yaml'%swd):
        return '%s/sw-cache.yaml'%swd
    return utils.get_single_entry(glob.glob('%s/sw-*.yaml'%swd))

# ----------------------------------------------------------------------------------------
def gloci():  # if no sw file, we probably only made one light locus (or i guess all input is missing, oh well)
    if args.single_chain:
        return [None]
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

        glfo, _, _ = utils.read_output(swfn(locus), skip_annotations=True)
        glutils.write_glfo(glfd(locus), glfo, debug=True)

        cmd = './bin/parse-output.py %s %s --airr-output --skip-columns clone_id' % (swfn(locus), ofn)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : wkdir(locus),
            'workdir' : wkdir(locus),
        }]
    utils.run_scan_cmds(args, cmdfos, 'airr-convert.log', n_total, n_already_there, ofn)

# ----------------------------------------------------------------------------------------
def run_single_chain_scoper():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, 0
    for locus in gloci():
        ofn = scofn(locus)
        if utils.output_exists(args, ofn): # and not args.dry:  # , offset=8):
            continue

        rcmds = ['library(scoper)']
        rcmds += ['db <- read.csv("%s", sep="\t")' % getifn(locus)]
        rcmds += ['db <- db[ , !(names(db) %in% c("clone_id"))]']  # it crashes if clone_id is already there UPDATE don't need this any more since i'm removing clone_id column above
        rcmds += ['results <- spectralClones(db, method="%s")' % ('novj' if args.no_shared_mutation else 'vj')]  # , method="novj", germline="germline_alignment_d_mask"
        rcmds += ['df <- as.data.frame(results)']
        rcmds += ['write.table(df[c("sequence_id", "clone_id")], "%s", sep="\t", quote=FALSE, row.names=FALSE)' % ofn]

        utils.run_r(rcmds, wkdir(locus), logfname='%s/scoper.log'%wkdir(locus), print_time='%s scoper'%locus, dryrun=args.dry)

# ----------------------------------------------------------------------------------------
def run_joint_scoper():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, 0

    ofn = scofn('joint')
    if utils.output_exists(args, ofn): # and not args.dry:  # , offset=8):
        return
    utils.merge_csvs(getifn('joint'), [getifn(l) for l in gloci()])
    timecmd = 'time' #'/usr/bin/time -o %s/%s --append' % (wkdir('joint'), 'joint-scoper.log')  # this should append it to the same file, but i end up with two since it writes it before utils deals with stuff (NOTE gdit i got it to work in the mobille run script tho)
    cmd = '%s Rscript packages/joint-scoper/scoperClones.R %s %s spectral HL %s 10' % (timecmd, getifn('joint'), ofn, 'novj' if args.no_shared_mutation else 'vj')
    cmdfos += [{
        'cmd_str' : cmd,
        'outfname' : ofn,
        'logdir' : wkdir('joint'),
        'workdir' : wkdir('joint'),
    }]
    utils.run_scan_cmds(args, cmdfos, 'joint-scoper.log', n_total, n_already_there, ofn)

# ----------------------------------------------------------------------------------------
def convert_output(joint=False):
    # ----------------------------------------------------------------------------------------
    def get_joint_time():
        try:
            timeline, _ = utils.simplerun('grep maxresident %s/joint-scoper.log' % wkdir('joint'), return_out_err=True)
        except subprocess.CalledProcessError:
            print('    %s grep failed' % utils.color('red', 'error'))
            return -1
        user, system, elapsed, cpu, _, _ = timeline.strip().split()
        etstr = elapsed.replace('elapsed', '')
        assert etstr.count(':') == 1
        minutes, seconds = [float(vstr) for vstr in etstr.split(':')]
        etime = 60 * minutes + seconds
        return etime
    # ----------------------------------------------------------------------------------------
    if joint:
        print('    joint scoper time: %.1f' % get_joint_time())
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(gloci())
    for locus in gloci():
        pfn = getofn(locus, joint=joint)
        if utils.output_exists(args, pfn, debug=False):
            n_already_there += 1
            continue
        utils.mkdir(pfn, isfile=True)
        cmd = './bin/parse-output.py %s %s --airr-input' % (scofn('joint' if joint else locus), pfn)
        if args.simdir is not None:
            cmd += ' --simfname %s' % simfn(locus)
        if joint:
            cmd += ' --skip-other-locus --locus %s --glfo-dir %s' % (locus, glfd(locus))
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
parser.add_argument('--indir', required=True)  # should have partis sw file in it
parser.add_argument('--outdir', required=True)
parser.add_argument('--simdir')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--single-chain', action='store_true')
parser.add_argument('--no-shared-mutation', action='store_true')
# parser.add_argument('--infname')  # for use with --single-chain
parser.add_argument('--n-max-procs', type=int, help='NOT USED')
args = parser.parse_args()
if not args.single_chain:
    assert args.simdir is not None

get_alignments()
run_single_chain_scoper()  # note: this doesn't get used in any way for joint scoper, it's just in case we need the single chain results for something
convert_output()
if not args.single_chain:
    run_joint_scoper()
    convert_output(joint=True)
