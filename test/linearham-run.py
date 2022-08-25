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
import pwd

# sys.path.insert(1, './python')
# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/test', '')
sys.path.insert(1, partis_dir + '/python')
import utils
import glutils
import paircluster
from clusterpath import ClusterPath
import treeutils

docker_path = '/linearham/work'
ig_or_tr = 'ig'

# ----------------------------------------------------------------------------------------
def finalfn(locus, inferred_ancestors=False, itree=None):
    odir = args.outdir
    if inferred_ancestors:
        odir += '/with-inferred-ancestors'
    if itree is not None:
        odir += '/itree-%d' % itree
    return paircluster.paired_fn(odir, locus, single_chain=not inferred_ancestors, actstr='partition', suffix='.yaml')
# ----------------------------------------------------------------------------------------
def simfn(locus):
    return paircluster.paired_fn(args.simdir, locus, suffix='.yaml')
# ----------------------------------------------------------------------------------------
def basedir():
    if args.partis_outdir.split('/')[-1] == 'partis':  # paired clustering validation (i.e. from cf-paired-loci.py)
        return '/'.join(args.partis_outdir.split('/')[:-1])  # i'm sure there's a better way to get the parent dir
    else:
        return args.partis_outdir  # hmm, maybe?
# ----------------------------------------------------------------------------------------
def wkdir(locus, iclust=None):
    return '%s/work/%s%s' % (args.outdir, locus, '' if iclust is None else '/iclust-%d'%iclust)
# ----------------------------------------------------------------------------------------
def lhodir(locus, iclust=None):
    glstr = '%s/cluster-%d/mcmc*/burnin*' % (wkdir(locus, iclust=iclust), iclust)
    flist = glob.glob(glstr)
    if len(flist) == 0:
        return glstr
    return utils.get_single_entry(flist)
# ----------------------------------------------------------------------------------------
def lnhofn(locus, iclust=None, trees=False):
    # cluster-0/mcmciter10000_mcmcthin10_tuneiter5000_tunethin100_numrates4_rngseed0/burninfrac0.1_subsampfrac0.05/aa_naive_seqs.dnamap
    return '%s/linearham_%s' % (lhodir(locus, iclust=iclust), 'run.trees' if trees else 'annotations_best.yaml' )
# ----------------------------------------------------------------------------------------
def ptnfn(locus, for_work=False):
    pdir = wkdir(locus) if args.ignore_unmutated_seqs and for_work else args.partis_outdir
    return paircluster.paired_fn(pdir, locus, actstr='partition', suffix='.yaml')
# ----------------------------------------------------------------------------------------
def prmd(locus):
    return '%s/parameters/%s' % (args.partis_outdir, locus)
# ----------------------------------------------------------------------------------------
def dckr_trns(inpath):  # translate <inpath> to path within docker
    if not args.docker:
        return inpath
    assert basedir() in inpath
    return inpath.replace(basedir(), docker_path)
# ----------------------------------------------------------------------------------------
def antn_plotdir(locus):
    return '%s/plots/%s/hmm' % (os.path.dirname(finalfn(locus)), locus)
# ----------------------------------------------------------------------------------------
def gloci():  # if no sw file, we probably only made one light locus (or i guess all input is missing, oh well)
    return [l for l in utils.sub_loci(ig_or_tr) if os.path.exists(simfn(l) if args.simdir is not None else ptnfn(l))]

# ----------------------------------------------------------------------------------------
def get_clusters(locus):
    if args.n_sim_events is None:
        _, _, cpath = utils.read_output(ptnfn(locus), skip_annotations=True)
        return cpath.best()
    else:
        return [None for _ in range(args.n_sim_events)]

# ----------------------------------------------------------------------------------------
def run_linearham():
    existing_ofns, cmdfos, n_already_there, n_too_small, n_total = [], [], 0, 0, 0
    for locus in gloci():
        if args.simdir is not None and not os.path.exists(simfn(locus)):
            continue
        if not os.path.exists(ptnfn(locus)):
            raise Exception('partition file doesn\'t exist (maybe forgot to run partis \'partition\' action first?) %s' % ptnfn(locus))
        if args.ignore_unmutated_seqs:
            glfo, antn_list, cpath = utils.read_output(ptnfn(locus))
            new_atns = []
            for atn in antn_list:
               iseqs_to_keep = [i for i, n in enumerate(atn['n_mutations']) if n>0]
               new_atn = utils.get_non_implicit_copy(atn)
               if len(iseqs_to_keep) == 0:
                   if len(atn['unique_ids']) >= args.min_cluster_size:
                       print '  %s no mutated seqs in %s cluster with size %d, so keeping all seqs' % (utils.wrnstr(), locus, len(atn['unique_ids']))  #  (yeah yeah this is dumb but whatever, i just don\'t want the indices to get screwed up)
                   utils.add_implicit_info(glfo, new_atn)  # could probably just keep the original one, but maybe i'll want to modify them at some point
               else:
                   utils.restrict_to_iseqs(new_atn, iseqs_to_keep, glfo)
               new_atns.append(new_atn)
            utils.write_annotations(ptnfn(locus, for_work=True), glfo, new_atns, utils.annotation_headers)
        for iclust, tclust in enumerate(get_clusters(locus)):
            if tclust is not None and len(tclust) < args.min_cluster_size:
                n_too_small += 1
                continue
            n_total += 1
            ofn = lnhofn(locus, iclust=iclust)
            if utils.output_exists(args, ofn, debug=False): # and not args.dry:  # , offset=8):
                n_already_there += 1
                existing_ofns.append(ofn)
                continue
            if iclust==0 and locus=='igh':
                print '    workdir: %s' % wkdir(locus, iclust=iclust)
            utils.mkdir(wkdir(locus, iclust=iclust))
            shlines = ['#!/bin/bash']
            # shlines += ['ls -ltrh %s %s' % (dckr_trns(ptnfn(locus, for_work=True)), dckr_trns(prmd(locus)))]
            # it would be nice to also have an option to run on my local linearham install (rather than docker), which works fine, except that it requires the output dir to be a subdir of the linearham code dir, which would require some linking shenanigans or something so not doing it atm
            shlines += ['scons --run-linearham --partis-yaml-file=%s --parameter-dir=%s --cluster-index=%d --outdir=%s' % (dckr_trns(ptnfn(locus, for_work=True)), dckr_trns(prmd(locus)), iclust, dckr_trns(wkdir(locus, iclust=iclust)))]
            if args.fast:
                shlines[-1] += ' --mcmc-iter=1000 --tune-iter=500'
            bfn = '%s/run.sh' % wkdir(locus, iclust=iclust)  #  NOTE if i'd used utils.simplerun() i couldn've used its cmdfname arg
            with open(bfn, 'w') as bfile:
                for l in shlines:
                    bfile.write('%s\n'%l)
            utils.simplerun('chmod +x %s' % bfn, debug=False)
            if args.docker:
                cmd = 'sudo docker run -it --rm -v%s:%s %s %s' % (basedir(), docker_path, 'linearham-local' if args.local_docker_image else 'quay.io/matsengrp/linearham', dckr_trns(bfn))
            else:
                cmd = bfn
            cmdfos += [{
                'cmd_str' : cmd,
                'outfname' : ofn,
                'logdir' : wkdir(locus, iclust=iclust),
                'workdir' : wkdir(locus, iclust=iclust),
            }]
    if n_too_small > 0:
        print '    skipped %d clusters smaller than %d (leaving %d)' % (n_too_small, args.min_cluster_size, n_total)
    utils.run_scan_cmds(args, cmdfos, 'linearham.log', n_total, n_already_there, None, existing_ofns=existing_ofns, dbstr='linearham run')

# ----------------------------------------------------------------------------------------
def read_lh_trees(treefname, glfo, antn_list):
    # ----------------------------------------------------------------------------------------
    def fix_ambig_regions(input_atn, new_seqfos, itree=None):
        if len(set(len(s) for s in input_atn['seqs'])) > 1:
            print '  %s seqs not all the same length, so giving up on fixing ambiguous regions' % utils.wrnstr()
            return
        ambig_positions = []
        for ichar in range(len(input_atn['seqs'][0])):
            chars = [s[ichar] for s in input_atn['seqs']]
            if set(chars) == set([utils.ambig_base]):
                ambig_positions.append(ichar)
        if len(ambig_positions) > 0:
            if itree == 0:
                print '  %d entirely ambiguous positions (in %d seqs), so setting those positions to %s in inferred ancestors' % (len(ambig_positions), len(input_atn['seqs']), utils.ambig_base)
            for sfo in new_seqfos:
                assert len(sfo['seq']) == len(input_atn['seqs'][0])
                # utils.color_mutants(sfo['seq'], ''.join([utils.ambig_base if i in ambig_positions else c for i, c in enumerate(sfo['seq'])]), print_result=True)
                sfo['seq'] = ''.join([utils.ambig_base if i in ambig_positions else c for i, c in enumerate(sfo['seq'])])
    # ----------------------------------------------------------------------------------------
    input_antn = utils.get_single_entry(antn_list)  # should just be length 1 i think
    treestrs = treeutils.get_treestrs_from_file(treefname)
    new_antns = []
    for itree, treestr in enumerate(treestrs):
        dtree = treeutils.get_dendro_tree(treestr=treestr, debug=False)  # this is super slow because it's got to read all the sequences (although, really, why tf does that have to be slow)
        # print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=dtree))
        inferred_nodes = [n for n in dtree.preorder_node_iter() if n.taxon.label not in input_antn['unique_ids']]
        new_seqfos = [{'name' : n.taxon.label, 'seq' : n.annotations['ancestral'].value} for n in inferred_nodes]
        newatn = utils.get_non_implicit_copy(input_antn)
        fix_ambig_regions(input_antn, new_seqfos, itree=itree)
        utils.add_implicit_info(glfo, newatn)
        utils.add_seqs_to_line(newatn, new_seqfos, glfo, debug=False)
        for sfo in new_seqfos:
            iseq = newatn['unique_ids'].index(sfo['name'])
            if 'multiplicities' in newatn:
                assert newatn['multiplicities'][iseq] is None
                newatn['multiplicities'][iseq] = 1
        newatn['tree-info'] = {'lb' : {'tree' : dtree.as_string(schema='newick')}}
        new_antns.append(newatn)
    return new_antns

# ----------------------------------------------------------------------------------------
def processs_linearham_output():
    n_already_there, n_too_small, missing_iclusts, n_total_iclusts, n_total_out = 0, 0, [], 0, 0
    missing_icpaths = []
    for locus in gloci():
        if args.simdir is not None and  not os.path.exists(simfn(locus)):
            continue
        clusters = get_clusters(locus)

        if '/fh/fast/' not in wkdir(locus):
            pwstruct = pwd.getpwuid(os.getuid())
            utils.simplerun('sudo chown -R %s:%d %s' % (pwstruct.pw_name, pwstruct.pw_gid, wkdir(locus)), dryrun=args.dry)  # NOTE not really the right group

        ofn = finalfn(locus)
        n_total_out += 1
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            n_total_iclusts += utils.non_none([args.n_sim_events, len(clusters) - clusters.count(None)])
            continue

        glfo = None
        antn_list = []  # one linearham (best) annotation for each cluster
        anc_antns = []  # many annotations for each cluster (one for each sampled/inferred linearham tree), each with inferred intermediates added to the annotation
        for iclust, tclust in enumerate(clusters):
            if tclust is not None and len(tclust) < args.min_cluster_size:
                n_too_small += 1
                continue
            n_total_iclusts += 1
            lhfn = lnhofn(locus, iclust=iclust)
            if not os.path.exists(lhfn):
                missing_iclusts.append(iclust)
                missing_icpaths.append(lhfn)
                continue
            if args.dry:
                antn_list.append(None)  # just to print the right length
                continue
            glfo, iclust_antns, _ = utils.read_output(lhfn)
            antn_list.append(utils.get_single_entry(iclust_antns))
            anc_antns.append(read_lh_trees(lnhofn(locus, iclust=iclust, trees=True), glfo, iclust_antns))

        print '    %s %d clusters to %s' % ('would write' if args.dry else 'writing', len(antn_list), ofn)
        if not args.dry:
            utils.write_annotations(ofn, glfo, antn_list, utils.annotation_headers)
            n_sampled_trees = set([len(l) for l in anc_antns])
            if len(n_sampled_trees) > 1:
                print '  %s different number of sampled trees for different clusters (using smallest, i.e. discarding some): %s' % (utils.wrnstr(), ' '.join(str(n) for n in sorted(n_sampled_trees)))
            for itree in range(min(n_sampled_trees) if len(n_sampled_trees)>0 else 0):
                utils.write_annotations(finalfn(locus, inferred_ancestors=True, itree=itree), glfo, [alist[itree] for alist in anc_antns], utils.annotation_headers)

        if args.simdir is not None:
            cmd = './bin/parse-output.py %s %s/x.fa' % (ofn, wkdir(locus))
            cmd += ' --only-make-plots --simfname %s --plotdir %s --only-csv-plots --only-plot-performance' % (simfn(locus), antn_plotdir(locus))
            utils.simplerun(cmd, logfname='%s/plot-performance.log'%wkdir(locus), dryrun=args.dry)

    if n_too_small > 0:
        print '    skipped %d clusters smaller than %d (leaving %d)' % (n_too_small, args.min_cluster_size, n_total_iclusts)
    if len(missing_iclusts) > 0:
        print '  missing %d / %d: iclusts %s (e.g. %s)' % (len(missing_iclusts), n_total_iclusts, ' '.join(str(i) for i in missing_iclusts), missing_icpaths[0])
    if n_already_there > 0:
        print '  %d / %d final output files already there (e.g. %s' % (n_already_there, n_total_out, ofn)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--simdir')
parser.add_argument('--outdir', required=True)
parser.add_argument('--partis-outdir', required=True)
parser.add_argument('--n-sim-events', type=int)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--docker', action='store_true')
parser.add_argument('--local-docker-image', action='store_true')
parser.add_argument('--fast', action='store_true')
# parser.add_argument('--remove-duplicate-seqs', action='store_true')
parser.add_argument('--ignore-unmutated-seqs', action='store_true')
parser.add_argument('--n-max-procs', type=int, help='NOT USED')
parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--min-cluster-size', default=5, type=int)
args = parser.parse_args()
if not args.docker:
    raise Exception('will crash cause linearham needs output dir as subdir of code dir')

run_linearham()
processs_linearham_output()
