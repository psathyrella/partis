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
from io import open

import utils
import glutils
import paircluster
from clusterpath import ClusterPath

ig_or_tr = 'ig'

imgt_by_hand_str = """
you need to get these with/from the imgt web site by hand, after first running just the first convert_file() call below
then running e.g. tar czf /fh/fast/matsen_e/dralph/partis/paired-loci/vs-shm/v3/imgt-input.tgz /fh/fast/matsen_e/dralph/partis/paired-loci/vs-shm/v3/seed-*/scratch-mute-freq-*/mobille/work/*/imgt-input/*.fa)
(or maybe just seeds 0 + 1)
then extract files on machine with browser, then mv to single dir: find . -name '*.fa' -exec mv {} . \;
then add them all individually to imgt/high v-quest search page

"""

# ----------------------------------------------------------------------------------------
def wkdir(locus):
    return '%s/work%s' % (args.outdir, '' if locus is None else '/'+locus)

# ----------------------------------------------------------------------------------------
def imgt_indir(locus):
    return '%s/imgt-input' % wkdir(locus)

# ----------------------------------------------------------------------------------------
def idlstr(locus, imgt=False):
    if args.single_chain:
        return args.id_str
    istr = '%s-%s' % (args.id_str, locus)
    if imgt:
        istr = istr.replace('-', '_').replace('.', '_') #+'_fa'  # friggin imgt changes - and . to _ UPDATE arg the new version doesn't have '_fa'
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
def infn(locus):
    if args.simdir is None:
        return args.inpath
    else:
        return simfn(locus)

# ----------------------------------------------------------------------------------------
def simfn(locus):
    return paircluster.paired_fn(args.simdir, locus, suffix='.yaml')

# ----------------------------------------------------------------------------------------
def getofn(locus):
    if args.single_chain:
        return '%s/partition.yaml' % args.outdir
    return paircluster.paired_fn(args.outdir, locus, single_chain=True, actstr='partition', suffix='.yaml')

# ----------------------------------------------------------------------------------------
def gloci():  # if no sw file, we probably only made one light locus (or i guess all input is missing, oh well)
    if args.single_chain:
        return [None]
    return [l for l in utils.sub_loci(ig_or_tr) if os.path.exists(simfn(l))]

# ----------------------------------------------------------------------------------------
def igbofn(locus):
    return '%s/igblast-out.tsv' % wkdir(locus)

# ----------------------------------------------------------------------------------------
def antn_plotdir(locus):
    return '%s/plots/%s/hmm' % (os.path.dirname(getofn(locus)), locus)

# ----------------------------------------------------------------------------------------
def convert_file(ifnfcn, ofnfcn, airr_input=False, igbl_gldir=False, plot_performance=False):
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(gloci())
    for locus in gloci():
        ofn = ofnfcn(locus)
        if utils.output_exists(args, imgt_outdir(locus), debug=True, outlabel='imgt output'):
            continue
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue
        utils.mkdir(wkdir(locus))
        cmd = './bin/parse-output.py %s %s' % (ifnfcn(locus), ofn)
        if airr_input:
            cmd += ' --airr-input'
        if igbl_gldir:
            cmd += ' --glfo-dir %s/igbl-db/fasta/partis-format-parsed --locus %s' % (args.igbdir, locus)  #  --template-glfo-dir %s/data/germlines/%s, utils.get_partis_dir(), args.species)  # just ran with the template set once so i could write the parsed glfo dirs
        if plot_performance:
            cmd += ' --simfname %s --plotdir %s --only-csv-plots --only-plot-performance' % (simfn(locus), antn_plotdir(locus))
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : wkdir(locus),
            'workdir' : wkdir(locus),
        }]
    utils.run_scan_cmds(args, cmdfos, 'imgt-input.log', n_total, n_already_there, ofn, dbstr='imgt conversion')

# # ----------------------------------------------------------------------------------------
# def initialize_germline_info(outdir):
#     if len(glob.glob('%s/%s/%s*-unaligned.fasta.nsq' % (args.igbdir, args.locus, args.locus))) == len(utils.regions):
#         return
# TODO use newer code from igblast 1.17.1 immcantation dir

#     cmds = ['#!/bin/bash']
#     cmds += ['export PATH=%s:$PATH' % args.condapath]
#     # cmds += ['igblastn -help']
#     cmds += ['cd %s/%s' % (args.igbdir, args.locus)]
#     for tmpreg in utils.regions:  # need a dummy d for igk/igl
#         cmds += ['perl ../edit_imgt_file.pl %s%s.fasta > %s%s-unaligned.fasta' % (args.locus, tmpreg, args.locus, tmpreg)]
#         cmds += ['makeblastdb -parse_seqids -dbtype nucl -in %s%s-unaligned.fasta' % (args.locus, tmpreg)]
#     utils.simplerun('\n'.join(cmds) + '\n', cmdfname=outdir + '/run.sh')

# ----------------------------------------------------------------------------------------
# # airr: --outfmt 19
# export IGDATA=~/share/igblast
# igblastn \
#     -germline_db_V ~/share/igblast/database/imgt_human_ig_v\
#     -germline_db_D ~/share/igblast/database/imgt_human_ig_d \
#     -germline_db_J ~/share/igblast/database/imgt_human_ig_v \
#     -auxiliary_data ~/share/igblast/optional_file/human_gl.aux \
#     -domain_system imgt -ig_seqtype Ig -organism human \
#     -outfmt '7 std qseq sseq btop' \
#     -query HD13M.fasta \
#     -out HD13M.fmt7
def run_igblast():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(gloci())
    for locus in gloci():
        if not os.path.exists(infn(locus)):
            continue
        ofn = igbofn(locus)
        if utils.output_exists(args, ofn, debug=False): # and not args.dry:  # , offset=8):
            n_already_there += 1
            continue
        shlines = ['#!/bin/bash']
        # shlines += ['cd %s' % args.igbdir]
        shlines += ['export IGDATA=%s' % args.igbdir]
        # shlines += ['export PATH=%s:$PATH' % args.condapath]
        shlines += ['%s/bin/igblastn' % args.igbdir]
        for tmpreg in utils.regions:
            shlines[-1] += ' -germline_db_%s %s/igbl-db/database/imgt_%s_ig_%s' % (tmpreg.upper(), args.igbdir, args.species, tmpreg)
        shlines[-1] += ' -min_D_match 5'  # default results in tons of seqs not having a D call, but it crashes if you set it less than 5 (help:    Required minimal consecutive nucleotide base matches for D genes)
        # maybe:
        # -extend_align5end
        #   Extend V gene alignment at 5' end
        # -extend_align3end
        #   Extend J gene alignment at 3' end
        shlines[-1] += ' -auxiliary_data %s/optional_file/%s_gl.aux' % (args.igbdir, args.species)
        shlines[-1] += ' -domain_system imgt -ig_seqtype Ig -organism %s -outfmt 19' % args.species #\'7 std qseq sseq btop\'' % args.species
        shlines[-1] += ' -num_threads %d' % args.n_procs
        shlines[-1] += ' -query ' + imgt_infname(locus) + ' -out ' + ofn
        bfn = '%s/run.sh' % wkdir(locus)  #  NOTE if i'd used utils.simplerun() i couldn've used its cmdfname arg
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
    utils.run_scan_cmds(args, cmdfos, 'igblast.log', n_total, n_already_there, ofn, dbstr='igblast run')

# ----------------------------------------------------------------------------------------
def run_mobille():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, len(gloci())
    for locus in gloci():
        if not os.path.exists(infn(locus)):
            continue
        ofn = mbofn(locus)
        if utils.output_exists(args, ofn, debug=False): # and not args.dry:  # , offset=8):
            n_already_there += 1
            continue

        if not os.path.isdir(imgt_outdir(locus)):
            if os.path.exists(imtxzfn(locus)):
                utils.simplerun('mkdir -p %s && cd %s && tar xf %s' % (imgt_outdir(locus), imgt_outdir(locus), imtxzfn(locus)), shell=True, dryrun=args.dry)
            else:
                raise Exception('    missing imgt output txz %s\n%s' % (imtxzfn(locus), imgt_by_hand_str))

        # cmd = 'time bash run_MobiLLe.sh ../Data/Real_datasets/IMGT_highvquest_output/toy_dataset _output/tmp'
        shlines = [
            '#!/bin/bash',
            'cd packages/MobiLLe/src',
            '/usr/bin/time -f"mobille time: %%e" bash run_MobiLLe.sh %s %s 2>&1' % (imgt_outdir(locus), wkdir(locus)),  # it makes its own sub outdir
        ]
        bfn = '%s/run.sh' % wkdir(locus)  #  NOTE if i'd used utils.simplerun() i couldn't used its cmdfname arg
        utils.mkdir(bfn, isfile=True)
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
def convert_mobille_output():
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
        true_partition = None
        if args.simdir is not None:
            _, _, true_cpath = utils.read_output(simfn(locus), skip_annotations=True)
            true_partition = true_cpath.best()
        plines = ClusterPath(partition=partition).get_partition_lines(true_partition=true_partition, calc_missing_values='best')
        print('    writing partition to %s' % pfn) 
        utils.write_annotations(pfn, {}, [], utils.annotation_headers, partition_lines=plines)

# # ----------------------------------------------------------------------------------------
# def install():
# pip3 install --user python-Levenshtein scikit-bio tqdm palettable

# install()
# sys.exit()

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('action')
parser.add_argument('--inpath')  # set this if 1) running on data or 2) running on single chain simulation (ok latter i think won't quite work yet, you'll need to fix)
parser.add_argument('--simdir')  # set this if you're running on paired simulation (as from cf-paired-loci.py)
parser.add_argument('--outdir', required=True)
parser.add_argument('--id-str', default='')
parser.add_argument('--base-imgt-outdir', required=True)
parser.add_argument('--prep', action='store_true', help='only prep for imgt')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--n-max-procs', type=int, help='NOT USED')
parser.add_argument('--igbdir', default='%s/packages/ncbi-igblast-1.17.1'%os.getcwd())
parser.add_argument('--species', default='human')
parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--single-chain', action='store_true')
args = parser.parse_args()
args.base_imgt_outdir = os.path.realpath(args.base_imgt_outdir)
args.outdir = os.path.realpath(args.outdir)
assert [args.inpath, args.simdir].count(None) == 1  # set exactly one of these
if args.inpath is not None:
    assert args.single_chain  # would just need to be implemented

convert_file(infn, imgt_infname)  # convert to fasta
if args.prep:
    sys.exit(0)
if args.action == 'mobille':
    run_mobille()
    convert_mobille_output()
elif args.action == 'igblast':
    run_igblast()
    convert_file(igbofn, getofn, airr_input=True, igbl_gldir=True, plot_performance=True)
else:
    raise Exception('unsuppored action %s' % args.action)
