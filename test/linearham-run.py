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

# sys.path.insert(1, './python')
# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/test', '')
sys.path.insert(1, partis_dir + '/python')
import utils
import glutils
import paircluster
from clusterpath import ClusterPath

docker_path = '/linearham/work'
ig_or_tr = 'ig'

# ----------------------------------------------------------------------------------------
def finalfn(locus):
    return paircluster.paired_fn(args.outdir, locus, single_chain=True, actstr='partition', suffix='.yaml')
# ----------------------------------------------------------------------------------------
def simfn(locus):
    return paircluster.paired_fn(args.simdir, locus, suffix='.yaml')
# ----------------------------------------------------------------------------------------
def basedir():
    assert args.partis_outdir.split('/')[-1] == 'partis'
    return '/'.join(args.partis_outdir.split('/')[:-1])  # i'm sure there's a better way to get the parent dir
# ----------------------------------------------------------------------------------------
def wkdir(locus, iclust=None):
    return '%s/work/%s%s' % (args.outdir, locus, '' if iclust is None else '/iclust-%d'%iclust)
# ----------------------------------------------------------------------------------------
def lhodir(locus, iclust=None):
    flist = glob.glob('%s/cluster-%d/mcmc*/burnin*' % (wkdir(locus, iclust=iclust), iclust))
    if len(flist) == 0:
        return '/x/y/z'
    return utils.get_single_entry(flist)
# ----------------------------------------------------------------------------------------
def lnhofn(locus, iclust=None):
    # cluster-0/mcmciter10000_mcmcthin10_tuneiter5000_tunethin100_numrates4_rngseed0/burninfrac0.1_subsampfrac0.05/aa_naive_seqs.dnamap
    return '%s/linearham_annotations_best.yaml' % lhodir(locus, iclust=iclust)
# ----------------------------------------------------------------------------------------
def ptnfn(locus):
    return paircluster.paired_fn(args.partis_outdir, locus, actstr='partition', suffix='.yaml')
# ----------------------------------------------------------------------------------------
def prmd(locus):
    return '%s/parameters/%s' % (args.partis_outdir, locus)
# ----------------------------------------------------------------------------------------
def dckr_trns(inpath):  # translate <inpath> to path within docker
    assert basedir() in inpath
    return inpath.replace(basedir(), docker_path)
# ----------------------------------------------------------------------------------------
def antn_plotdir(locus):
    return '%s/plots/%s/hmm' % (os.path.dirname(finalfn(locus)), locus)
# ----------------------------------------------------------------------------------------
def gloci():  # if no sw file, we probably only made one light locus (or i guess all input is missing, oh well)
    return [l for l in utils.sub_loci(ig_or_tr) if os.path.exists(simfn(l))]

# ----------------------------------------------------------------------------------------
def run_linearham():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, 0
    for locus in gloci():
        if not os.path.exists(simfn(locus)):
            continue
        if not os.path.exists(ptnfn(locus)):
            raise Exception('partition file doesn\'t exist (maybe forgot to run partis \'partition\' action first?) %s' % ptnfn(locus))
# ----------------------------------------------------------------------------------------
# TODO chown output from root
        for iclust in range(args.n_sim_events):
            n_total += 1
            ofn = lnhofn(locus, iclust=iclust)
            if utils.output_exists(args, ofn, debug=False): # and not args.dry:  # , offset=8):
                n_already_there += 1
                continue
            if iclust==0 and locus=='igh':
                print '    workdir: %s' % wkdir(locus, iclust=iclust)
            utils.mkdir(wkdir(locus, iclust=iclust))
            shlines = ['#!/bin/bash']
            # shlines += ['ls -ltrh %s %s' % (dckr_trns(ptnfn(locus)), dckr_trns(prmd(locus)))]
            shlines += ['grep get-linearham-info SConstruct']
            shlines += ['sed -i \'s/get-linearham-info/get-linearham-info --min-tree-metric-cluster-size 3/\' SConstruct']
            shlines += ['grep get-linearham-info SConstruct']
            shlines += ['scons --run-linearham --partis-yaml-file=%s --parameter-dir=%s --cluster-index=%d --outdir=%s' % (dckr_trns(ptnfn(locus)), dckr_trns(prmd(locus)), iclust, dckr_trns(wkdir(locus, iclust=iclust)))]
            bfn = '%s/run.sh' % wkdir(locus, iclust=iclust)  #  NOTE if i'd used utils.simplerun() i couldn've used its cmdfname arg
            with open(bfn, 'w') as bfile:
                for l in shlines:
                    bfile.write('%s\n'%l)
            utils.simplerun('chmod +x %s' % bfn, debug=False)
            cmd = 'sudo docker run -it --rm -v%s:%s quay.io/matsengrp/linearham %s' % (basedir(), docker_path, dckr_trns(bfn))
            cmdfos += [{
                'cmd_str' : cmd,
                'outfname' : ofn,
                'logdir' : wkdir(locus),
                'workdir' : wkdir(locus),
            }]
    utils.run_scan_cmds(args, cmdfos, 'linearham.log', n_total, n_already_there, ofn, dbstr='linearham run')

# ----------------------------------------------------------------------------------------
def processs_linearham_output():
    n_already_there, n_missing_iclusts, n_total_iclusts, n_total_out = 0, 0, 0, 0
    for locus in gloci():
        if not os.path.exists(simfn(locus)):
            continue
        ofn = finalfn(locus)
        n_total_out += 1
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue
# TODO don't overwrite by default
        antn_list = []
        for iclust in range(args.n_sim_events):
            n_total_iclusts += 1
            lhfn = lnhofn(locus, iclust=iclust)
            if not os.path.exists(lhfn):
                n_missing_iclusts += 1
                print '    missing %s' % lhfn
                continue
            glfo, iclust_antns, _ = utils.read_output(lhfn)
            antn_list += iclust_antns
        print '    writing %d clusters to %s' % (len(antn_list), ofn)
        utils.write_annotations(ofn, glfo, antn_list, utils.annotation_headers)

        cmd = './bin/parse-output.py %s %s/x.fa' % (ofn, wkdir(locus))
        cmd += ' --only-make-plots --simfname %s --plotdir %s --only-csv-plots --only-plot-performance' % (simfn(locus), antn_plotdir(locus))
        utils.simplerun(cmd, logfname='%s/plot-performance.log'%wkdir(locus)) #, dryrun=args.dry)

    if n_missing_iclusts > 0:
        print '  missing %d / %d' % (n_missing_iclusts, n_total_iclusts)
    if n_already_there > 0:
        print '  %d / %d final output files already there (e.g. %s' % (n_already_there, n_total_out, ofn)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--simdir', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--partis-outdir', required=True)
parser.add_argument('--n-sim-events', type=int, required=True)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--n-max-procs', type=int, help='NOT USED')
parser.add_argument('--n-procs', default=1, type=int)
args = parser.parse_args()

run_linearham()
processs_linearham_output()
