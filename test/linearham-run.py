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
import pwd
import collections
from io import open

import partis.utils as utils
import partis.glutils as glutils
import partis.paircluster as paircluster
from partis.clusterpath import ClusterPath
import partis.treeutils as treeutils

docker_path = '/linearham/work'
ig_or_tr = 'ig'

# ----------------------------------------------------------------------------------------
def finalfn(locus, nwk=False, sampled_trees=False, alt_antns=False, iclust=None, lineage_plot=False):
    odir = args.outdir
    if sampled_trees:
        assert iclust is not None
        return '%s/sampled-trees/%s/iclust-%d.nwk' % (odir, locus, iclust)  # this is kind of complicated, but the linearham output dirs are a total mess (tons of duplicated information and intermediate files, plus poor naming schemes), so it's nice to at least link to one thing we're pretty likely to need
    elif alt_antns:
        assert iclust is not None
        return '%s/alternative-annotations/iclust-%d/partition-%s.yaml' % (odir, iclust, locus)  # same comment as for sampled trees above
    elif nwk:
        return '%s/trees-%s.nwk' % (odir, locus)  # one tree for each cluster
    elif lineage_plot:
        return '%s/lineage-plots/iclust-%d-%s.png' % (odir, iclust, locus)
    else:
        return paircluster.paired_fn(odir, locus, single_chain=True, actstr='partition', suffix='.yaml')  # atm only writing single chain (relic of paired clustering paper, where we wanted linearahm to show up as a "single chain" method so it go put in the right plot), although eventually may want to also write joint
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
def wkdir(locus, iclust=None, for_lh_cmd=False):
    odir = args.outdir
    if not args.docker and for_lh_cmd:
        if args.linearham_dir+'/' not in odir:
            raise Exception()
        odir = odir.replace(args.linearham_dir+'/', '')
    return '%s/work/%s%s%s' % (odir, locus, '' if iclust is None else '/iclust-%d'%iclust, '' if args.seed_unique_id is None else '/lineage-%s'%args.seed_unique_id)
# ----------------------------------------------------------------------------------------
def lhodir(locus, iclust=None):  # NOTE i think the only reason this fcn exists (i.e. we don't just use wkdir() is that we used to need the glob stuff), i.e. eventually can probably remove this
    ldir = wkdir(locus, iclust=iclust)
    glstr = '%s/cluster-%d/mcmc*/burnin*' % (wkdir(locus, iclust=iclust), iclust)  # cruft for backwards compatibility
    flist = glob.glob(glstr)
    if len(flist) > 0:
        raise Exception('old-style mcmc*/burnin*/ subdirs exist, you may want to move them to where we currently expect output\n  %s\n  %s' % (' '.join(flist), ldir))
    return ldir
# ----------------------------------------------------------------------------------------
def lnhofn(locus, iclust=None, trees=False, logprobs=False, best=False):
    return '%s/linearham_%s' % (lhodir(locus, iclust=iclust), 'run.trees' if trees else ('run.log' if logprobs else ('annotations_%s.yaml'%('best' if best else 'all'))))
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
    def tfn(l): return simfn(l) if args.simdir is not None else ptnfn(l)
    tloci = [l for l in utils.sub_loci(ig_or_tr) if os.path.exists(tfn(l))]
    if len(tloci) == 0:
        print('  %s returning zero loci (files probably don\'t exist: %s)' % (utils.wrnstr(), ' '.join(tfn(l) for l in utils.sub_loci(ig_or_tr))))
    return tloci

# ----------------------------------------------------------------------------------------
def get_clusters(locus):
    if args.n_sim_events is None:
        _, _, cpath = utils.read_output(ptnfn(locus), skip_annotations=True)
        return cpath.best()  # NOTE do *not* sort here, since the indices get passed on cmd line to linearham so we can't change order
    else:
        return [None for _ in range(args.n_sim_events)]

# ----------------------------------------------------------------------------------------
def run_linearham():
    # ----------------------------------------------------------------------------------------
    def prep_cmd(cmdfos, locus, iclust, ofn):
        workd = wkdir(locus, iclust=iclust)
        if iclust==0 and locus=='igh':
            print('    workdir: %s' % workd)
        utils.mkdir(workd)
        if not os.path.exists(prmd(locus)):
            print('  %s parameter dir doesn\'t exist, maybe need to ln it in from somewhere else? %s' % (utils.wrnstr, prmd(locus)))
        shlines = ['#!/bin/bash']
        # shlines += ['ls -ltrh %s %s' % (dckr_trns(ptnfn(locus, for_work=True)), dckr_trns(prmd(locus)))]
        if not args.docker:
            shlines += ['cd %s' % args.linearham_dir]
        shlines += ['scons --run-linearham --partis-yaml-file=%s --parameter-dir=%s --no-nestly-subdirs --cluster-index=%d --outdir=%s' % (dckr_trns(ptnfn(locus, for_work=True)), dckr_trns(prmd(locus)), iclust, dckr_trns(wkdir(locus, iclust=iclust, for_lh_cmd=True)))]
        if args.fast:
            shlines[-1] += ' --mcmc-iter=%s --tune-iter=%s' % (linearham_defaults[0]['mcmciter'], linearham_defaults[0]['tuneiter'])
        if args.seed_unique_id is not None:
            shlines[-1] += ' --lineage-unique-ids=%s' % args.seed_unique_id  # NOTE linearham breaks if you actually pass it more than one at once
        if args.asr_pfilters is not None:
            shlines[-1] += ' --asr-pfilters=%s' % args.asr_pfilters
        bfn = '%s/run.sh' % workd  #  NOTE if i'd used utils.simplerun() i couldn've used its cmdfname arg
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
            'logdir' : workd,
            'workdir' : workd,
        }]
    # ----------------------------------------------------------------------------------------
    example_existing_ofn, cmdfos, n_already_there, n_too_small, n_non_seed, n_total = None, [], 0, 0, 0, 0
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
                       print('  %s no mutated seqs in %s cluster with size %d, so keeping all seqs' % (utils.wrnstr(), locus, len(atn['unique_ids'])))  #  (yeah yeah this is dumb but whatever, i just don\'t want the indices to get screwed up)
                   utils.add_implicit_info(glfo, new_atn)  # could probably just keep the original one, but maybe i'll want to modify them at some point
               else:
                   utils.restrict_to_iseqs(new_atn, iseqs_to_keep, glfo)
               new_atns.append(new_atn)
            utils.write_annotations(ptnfn(locus, for_work=True), glfo, new_atns, utils.annotation_headers)
        for iclust, tclust in enumerate(get_clusters(locus)):
            if tclust is not None and len(tclust) < args.min_cluster_size:
                n_too_small += 1
                continue
            if args.seed_unique_id is not None and args.seed_unique_id not in tclust:
                n_non_seed += 1
                continue
            n_total += 1
            ofn = lnhofn(locus, iclust=iclust)
            if utils.output_exists(args, ofn, outlabel='linearham'):
                n_already_there += 1
                if example_existing_ofn is None:
                    example_existing_ofn = ofn
                continue
            prep_cmd(cmdfos, locus, iclust, ofn)

    if n_too_small > 0:
        print('    skipped %d clusters smaller than %d (leaving %d)' % (n_too_small, args.min_cluster_size, n_total))
    if n_non_seed > 0:
        print('    skipped %d clusters that didn\'t contain the seed id (leaving %d)' % (n_non_seed, n_total))
    utils.run_scan_cmds(args, cmdfos, 'linearham.log', n_total, n_already_there, None, example_existing_ofn=example_existing_ofn, dbstr='linearham run')

# ----------------------------------------------------------------------------------------
def read_lh_trees(locus, iclust, best=False):  # read trees and inferred ancestral seqs (which are in the nwk strings)
    # ----------------------------------------------------------------------------------------
    def fix_ambig_regions(new_seqfos, tatn, dbg=False):
        if len(set(len(s) for s in tatn['seqs'])) > 1:
            print('  %s seqs not all the same length, so giving up on fixing ambiguous regions' % utils.wrnstr())
            return
        ambig_positions = []
        for ichar in range(len(tatn['seqs'][0])):
            chars = [s[ichar] for s in tatn['seqs']]
            if set(chars) == set([utils.ambig_base]):
                ambig_positions.append(ichar)
        if len(ambig_positions) > 0:
            if dbg:
                print('  %d entirely ambiguous positions (in %d seqs), so setting those positions to %s in inferred ancestors' % (len(ambig_positions), len(tatn['seqs']), utils.ambig_base))
            for sfo in new_seqfos:
                assert len(sfo['seq']) == len(tatn['seqs'][0])
                # utils.color_mutants(sfo['seq'], ''.join([utils.ambig_base if i in ambig_positions else c for i, c in enumerate(sfo['seq'])]), print_result=True)
                sfo['seq'] = ''.join([utils.ambig_base if i in ambig_positions else c for i, c in enumerate(sfo['seq'])])
    # ----------------------------------------------------------------------------------------
    def get_full_atn(lh_atn, itree):
        if 'tree-info' not in lh_atn or 'linearham' not in lh_atn['tree-info']:
            raise Exception('tree info not in linearham annotation, probably need to run updated version of linearham/scripts/write_lh_annotations.py on %s' % lnhofn(locus, iclust=iclust))
        fatn = utils.get_full_copy(lh_atn, glfo)
        treestr = lh_atn['tree-info']['linearham']['trees'][itree]  # they should all be ~equally likely (linearham samples N (e.g. 45) trees (+annotation for each), then collapses duplicate annotations, and appends all trees contributing to each unique annotation into this key)
        dtree = treeutils.get_dendro_tree(treestr=treestr, debug=False)  # this is super slow because it's got to read all the sequences (although, really, why tf does that have to be slow)
        inferred_nodes = [n for n in dtree.preorder_node_iter() if n.taxon.label not in lh_atn['unique_ids']]  # note that these inferred ancestors correspond to the (arbitrary) chosen tree, but not necessarily to the other trees
        new_seqfos = [{'name' : n.taxon.label, 'seq' : n.annotations['ancestral'].value} for n in inferred_nodes]
        fix_ambig_regions(new_seqfos, fatn, dbg=True)
        utils.add_seqs_to_line(fatn, new_seqfos, glfo, debug=False)  # don't think there's any reason not to modify the annotation we read from the file
        fatn['tree-info']['lb'] = {'tree' : dtree.as_string(schema='newick')}  # don't use the original tree str, since it has ancestral sequence annotations which are super slow to read (now that I'm adding the linearham trees to ['tree-info']['linearham']['trees'] maybe it'd be better not to also have it here, but atm everything assumes trees are in ['lb'])
        return fatn
    # ----------------------------------------------------------------------------------------
    glfo, lhalist, _ = utils.read_output(lnhofn(locus, iclust=iclust, best=best))  # NOTE that linearham initially makes N annotations for each of N sampled trees, but then collapses duplicate annotations, and sets a 'logprob' key in the annotation based just on the number of times each annotation was seen
    # treestrs = treeutils.get_treestrs_from_file(lnhofn(locus, iclust=iclust, trees=True))  # file with just one tree per line (don't need these any more now i'm adding the trees to their annotations in linearham)
    # logprobs = XXX lnhofn(locus, iclust=iclust, logprobs=True) 'linearham_run.log'  # tsv with log probs + stuff (key is LHLogLikelihood)
    final_atns = []
    for lh_atn in lhalist:
        for itree in range(len(lh_atn['tree-info']['linearham']['trees'])):  # have to make a new annotation for each tree in each annotation's list of alternative trees, since each tree may have different inferred ancestors
            final_atns.append(get_full_atn(lh_atn, itree))
    if best:
        return final_atns[0], glfo  # will be a list of length the number of trees contributing to the most likely annotation
    else:
        return final_atns, glfo

# ----------------------------------------------------------------------------------------
def processs_linearham_output():
    if args.original_outdir is not None and os.path.exists(args.outdir):  # rm the stupid link bullshit
        try:
            os.remove(args.outdir)
        except IsADirectoryError:
            print('  %s couldn\'t remove %s (it\'s a dir, not a link)' % (utils.wrnstr(), args.outdir))
        # os.rmdir(os.path.dirname(args.outdir))  # UPDATE nah can't do this, other processes might be sharing this parent dir OLD: maybe? it would be nice to keep going upwards and delete everything we made, but I'm not sure how to stop in exactly the right place
        args.outdir = args.original_outdir

    n_already_there, n_too_small, n_non_seed, missing_iclusts, n_total_iclusts, n_total_out = 0, 0, 0, [], 0, 0
    missing_icpaths = []
    for locus in gloci():
        if args.simdir is not None and not os.path.exists(simfn(locus)):
            continue
        clusters = get_clusters(locus)

        if args.docker and '/fh/fast/' not in wkdir(locus):  # NOTE not super sure this is right? adding first clause without checking anyway
            pwstruct = pwd.getpwuid(os.getuid())
            utils.simplerun('sudo chown -R %s:%d %s' % (pwstruct.pw_name, pwstruct.pw_gid, wkdir(locus)), dryrun=args.dry)  # NOTE not really the right group

        fofn = finalfn(locus)
        n_total_out += 1
        if utils.output_exists(args, fofn, outlabel='final annotation', debug=True):  # utils.all_outputs_exist(args, [fofn, finalfn(locus, inferred_ancestors=True, itree=0)], debug=False):  # would be nice to also include all the tree sampled subdirs? but maybe not worth it
            n_already_there += 1
            n_total_iclusts += utils.non_none([args.n_sim_events, len(clusters) - clusters.count(None)])
            continue

        # collect best linearham annotation for each cluster
        glfo = None
        best_antns, all_antns = [], []  # many annotations for each cluster (one for each sampled/inferred linearham tree), each with inferred intermediates added to the annotation
        for iclust, tclust in enumerate(clusters):
            if tclust is not None and len(tclust) < args.min_cluster_size:
                n_too_small += 1
                continue
            if args.seed_unique_id is not None and args.seed_unique_id not in tclust:
                n_non_seed += 1
                continue
            n_total_iclusts += 1
            lhfn = lnhofn(locus, iclust=iclust)
            if not os.path.exists(lhfn):
                missing_iclusts.append(iclust)
                missing_icpaths.append(lhfn)
                continue
            if args.dry:
                best_antns.append(None)  # just to print the right length
                continue
            tatn, glfo = read_lh_trees(locus, iclust, best=True)  # ancestors get added here (glfos should all be the same)
            best_antns.append(tatn)
            a_atns, _ = read_lh_trees(locus, iclust)
            all_antns.append(a_atns)

        if len(best_antns) == 0:
            print('  %s zero annotations to write, so exiting' % utils.wrnstr())
            sys.exit(0)
        print('    %s %d cluster%s to %s output file: %s' % ('would write' if args.dry else 'writing', len(best_antns), utils.plural(len(best_antns)), 'single chain' if os.path.basename(os.path.dirname(fofn))=='single-chain' else '???', fofn))  # ??? need updating if i decide to write the joint ones
        if not args.dry:
            utils.write_annotations(fofn, glfo, best_antns, utils.annotation_headers + ['logprob'])  # best annotation for each cluster, with inferred ancestral sequences added from one of the sampled trees for that annotation
            with open(finalfn(locus, nwk=True), 'w') as tfile:  # for convenience, write to fasta the tree from the best annotation from which we took inferred ancestors
                tfile.write('\n'.join(l['tree-info']['lb']['tree'] for l in best_antns))
            for iclust in range(len(clusters)):
                ftfn = finalfn(locus, sampled_trees=True, iclust=iclust)  # sampled tree nwk file (including ancestral sequences)
                utils.makelink(os.path.dirname(ftfn), os.path.abspath(lnhofn(locus, iclust=iclust, trees=True)), ftfn)
                alt_atn_fn = finalfn(locus, iclust=iclust, alt_antns=True)  # alternative (all) annotations, each with a list of all the trees that contributed to it
                utils.write_annotations(alt_atn_fn, glfo, all_antns[iclust], utils.annotation_headers + ['logprob'])
                for pngfname in glob.glob('%s/aa_lineage_seqs.pfilter*.png'%lhodir(locus, iclust=iclust)):  # will only be there if --seed-unique-id is set
                    lpfn = finalfn(locus, lineage_plot=True, iclust=iclust)
                    utils.makelink(os.path.dirname(lpfn), pngfname, lpfn)

        if args.simdir is not None:
            cmd = './bin/parse-output.py %s %s/x.fa' % (fofn, wkdir(locus))
            cmd += ' --only-make-plots --simfname %s --plotdir %s --only-csv-plots --only-plot-performance' % (simfn(locus), antn_plotdir(locus))
            utils.simplerun(cmd, logfname='%s/plot-performance.log'%wkdir(locus), dryrun=args.dry)

    if n_too_small > 0:
        print('    skipped %d clusters smaller than %d (leaving %d)' % (n_too_small, args.min_cluster_size, n_total_iclusts))
    if n_non_seed > 0:
        print('    skipped %d clusters that didn\'t contain the seed id (leaving %d)' % (n_non_seed, n_total))
    if len(missing_iclusts) > 0:
        print('  missing %d / %d linearham output files (e.g. %s)' % (len(missing_iclusts), n_total_iclusts, missing_icpaths[0])) # , ' '.join(str(i) for i in missing_iclusts)
    if n_already_there > 0:
        print('      %d / %d final output files already there (e.g. %s' % (n_already_there, n_total_out, fofn))

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--simdir')
parser.add_argument('--outdir', required=True, help='note that if not running with docker, this must be a subdir of --linearham-dir')
parser.add_argument('--partis-outdir', required=True)
parser.add_argument('--linearham-dir', help='if not running with docker, you must set this to the full path of the linearham code, so that it can be removed from the outdir path')
parser.add_argument('--n-sim-events', type=int)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--docker', action='store_true')
parser.add_argument('--local-docker-image', action='store_true')
parser.add_argument('--fast', action='store_true')
parser.add_argument('--asr-pfilters')  # make it a string to match linearham
# parser.add_argument('--remove-duplicate-seqs', action='store_true')
parser.add_argument('--ignore-unmutated-seqs', action='store_true')
parser.add_argument('--n-max-procs', type=int, help='NOT USED')
parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--min-cluster-size', default=5, type=int)
parser.add_argument('--seed-unique-id', help='colon-separated list of uids for which do make detailed linearham analyses (see linearham help).')
args = parser.parse_args()
if args.fast:
    linearham_defaults[0]['mcmciter'] = '1000'
    linearham_defaults[0]['tuneiter'] = '500'

# linearham (well, i think scons) gets confused if they don't have the full path
args.simdir = utils.fpath(args.simdir)
args.partis_outdir = utils.fpath(args.partis_outdir)
args.outdir = utils.fpath(args.outdir)

args.original_outdir = None
if not args.docker and (args.linearham_dir is None or args.linearham_dir not in args.outdir):
    args.original_outdir = args.outdir
    utils.mkdir(args.original_outdir)
    args.outdir = '%s/work/%s' % (args.linearham_dir, args.original_outdir.lstrip('/'))  # ok it's a little over verbose to use the full original path, but whatever
    if not os.path.exists(args.outdir):
        print('     --outdir is not a subdir of --linearham-dir, so %s it so it looks like one' % ('would link' if args.dry else 'linking'))
        if not args.dry:
            utils.mkdir(os.path.dirname(args.outdir))
            utils.makelink(os.path.dirname(args.outdir), args.original_outdir, args.outdir)

run_linearham()
processs_linearham_output()
