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
from clusterpath import ClusterPath

ig_or_tr = 'ig'

# ----------------------------------------------------------------------------------------
def wkdir(locus):
    return '%s/work/%s' % (args.outdir, locus)

# ----------------------------------------------------------------------------------------
def imgt_indir(locus):
    return '%s/imgt-input' % wkdir(locus)

# ----------------------------------------------------------------------------------------
def idlstr(locus, imgt=False):
    istr = '%s-%s' % (args.id_str, locus)
    if imgt:
        istr = istr.replace('-', '_')  # friggin imgt changes - to _
    return istr

# ----------------------------------------------------------------------------------------
def imgt_infname(locus):
    return '%s/%s.fa' % (imgt_indir(locus), idlstr(locus))

# ----------------------------------------------------------------------------------------
def imgt_outdir(locus):
    return '%s/%s' % (args.base_imgt_outdir, idlstr(locus, imgt=True)) #utils.getprefix(os.path.basename(imgt_infname(locus))).replace('-', '_'))

# ----------------------------------------------------------------------------------------
def imtxzfn(locus):
    return '%s.txz' % imgt_outdir(locus)

# ----------------------------------------------------------------------------------------
def mbofn(locus):
    return '%s/%s/%s_final_clusters_Fo.txt' % (wkdir(locus), idlstr(locus, imgt=True), idlstr(locus, imgt=True))

# ----------------------------------------------------------------------------------------
def simfn(locus):
    return paircluster.paired_fn(args.simdir, locus, suffix='.yaml')

# ----------------------------------------------------------------------------------------
def getofn(locus):
    return paircluster.paired_fn(args.outdir, locus, single_chain=True, actstr='partition', suffix='.yaml')

# ----------------------------------------------------------------------------------------
def gloci():  # if no sw file, we probably only made one light locus (or i guess all input is missing, oh well)
    return [l for l in utils.sub_loci(ig_or_tr) if os.path.exists(simfn(l))]

# ----------------------------------------------------------------------------------------
def prep_for_imgt():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(gloci())
    for locus in gloci():
        ofn = imgt_infname(locus)
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue
        utils.mkdir(wkdir(locus))
        cmd = './bin/parse-output.py %s %s' % (simfn(locus), ofn)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : wkdir(locus),
            'workdir' : wkdir(locus),
        }]
    utils.run_scan_cmds(args, cmdfos, 'imgt-input.log', n_total, n_already_there, ofn, dbstr='imgt conversion')

# ----------------------------------------------------------------------------------------
def run_mobille():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(gloci())
    for locus in gloci():
        if not os.path.exists(simfn(locus)):
            continue
        ofn = mbofn(locus)
        if utils.output_exists(args, ofn, debug=False): # and not args.dry:  # , offset=8):
            n_already_there += 1
            continue

        if not os.path.isdir(imgt_outdir(locus)):
            if os.path.exists(imtxzfn(locus)):
                utils.simplerun('mkdir -p %s && cd %s && tar xf %s' % (imgt_outdir(locus), imgt_outdir(locus), imtxzfn(locus)), shell=True, dryrun=args.dry)
            else:
                print '    missing imgt output txz %s' % imtxzfn(locus)

        # cmd = 'time bash run_MobiLLe.sh ../Data/Real_datasets/IMGT_highvquest_output/toy_dataset _output/tmp'
        shlines = [
            '#!/bin/bash',
            'cd packages/MobiLLe/src',
            '/usr/bin/time -f"mobille time: %%e" bash run_MobiLLe.sh %s %s 2>&1' % (imgt_outdir(locus), wkdir(locus)),  # it makes its own sub outdir
        ]
        bfn = '%s/run.sh' % wkdir(locus)
        with open(bfn, 'w') as bfile:
            for l in shlines:
                bfile.write('%s\n'%l)
        utils.simplerun('chmod +x %s' % bfn, debug=False)
        cmd = bfn
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : wkdir(locus),
            'workdir' : wkdir(locus),
        }]
    utils.run_scan_cmds(args, cmdfos, 'mobille.log', n_total, n_already_there, ofn, dbstr='mobille run')

# ----------------------------------------------------------------------------------------
def convert_output():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(gloci())
    for locus in gloci():
        pfn = getofn(locus)
        if utils.output_exists(args, pfn, debug=False):
            n_already_there += 1
            continue
        utils.mkdir(pfn, isfile=True)
        partition = []
        with open(mbofn(locus)) as mfile:
            for line in mfile:
                iclust, clstr = line.split('\t')
                iclust = int(iclust)
                partition.append([s.strip() for s in clstr.split()])
        _, _, true_cpath = utils.read_output(simfn(locus), skip_annotations=True)
        true_partition = true_cpath.best()
        plines = ClusterPath(partition=partition).get_partition_lines(true_partition=true_partition, calc_missing_values='best')
        print '    writing partition to %s' % pfn 
        utils.write_annotations(pfn, {}, [], utils.annotation_headers, partition_lines=plines)

# # ----------------------------------------------------------------------------------------
# def install():
# pip3 install --user python-Levenshtein scikit-bio tqdm palettable

# install()
# sys.exit()

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--simdir', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--id-str', required=True)
parser.add_argument('--base-imgt-outdir', required=True)
parser.add_argument('--prep', action='store_true', help='only prep for imgt')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--n-max-procs', type=int, help='NOT USED')
args = parser.parse_args()

prep_for_imgt()
if not args.prep:
    run_mobille()
    convert_output()
