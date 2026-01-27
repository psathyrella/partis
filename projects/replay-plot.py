#!/usr/bin/env python3
from __future__ import absolute_import
from __future__ import print_function
import sys
import csv
from io import open
from six.moves import zip
from six.moves import range
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import glob
import argparse
import colored_traceback.always
import json
import itertools
import numpy
import collections
import matplotlib as mpl
mpl.rcParams["mathtext.default"] = 'regular'

import partis.utils as utils
import partis.plotting as plotting
import partis.hutils as hutils
from partis.hist import Hist
import partis.treeutils as treeutils
import partis.datautils as datautils

colors = {
    'bst-data' : '#cc0000', #'black',
    'gct-data' : '#cc0000',  # light red
    'gct-data-d15' : '#ea7979',  # darker red
    'gct-data-d20' : '#cc0000',
    'gct-data-w10' : '#2b65ec',  # light blue
    'bst-data-d20' : '#006600',
    'iqt-data' : '#a821c7',
    'iqt-data-d20' : '#a821c7',
    'iqt-data-d15' : '#a821c7',
    'iqt-data-w10' : '#a821c7',
    'simu' : '#006600', ##808080',
    'simu-iqtree' : '#006600',
}
linestyles = {
    'simu' : 'dashed',
    # 'data-d20' : 'dashed',
}
linewidths = {
    'simu' : 4,
    'simu-iqtree' : 2,
    # 'data-d20' : 4,
}
pltlabels = {
    'abundances' : 'abundances',
    'max-abundances' : 'max abundance in each GC',
    'csizes' : 'N leaves per tree',
    'leaf-muts-nuc' : 'N muts (leaf nodes)',
    'n_stops' : 'N stop codons',
}
tpstrs = ['affinity', 'muts-nuc', 'muts-aa']
labelstrs = ['affinity ($ -\Delta log_{10} K_D $)', 'N nuc muts', 'N AA muts']
all_timepoints = ['d20', 'w10', 'd15']
for tstr, lbstr in zip(tpstrs, labelstrs):
    for nstr in ['leaf', 'internal']:
        pltlabels['%s-%s'%(nstr, tstr)] = lbstr

# ----------------------------------------------------------------------------------------
def abfn(tlab, abtype='abundances'):
    return '%s/%s/%s.csv' % (args.outdir, tlab, abtype)

# ----------------------------------------------------------------------------------------
def abdn_hargs(hlist):
    xmax = max(h.xmax for h in hlist)
    xbounds = [0.5, xmax+0.5]
    xticks = list(range(1, int(xmax)+1, 2 if xmax<10 else 5))
    ymin = min(h.get_minimum(exclude_empty=True) for h in hlist)
    ymax = max(h.get_maximum() for h in hlist)
    ybounds = [0.9 * ymin, 1.1 * ymax]
    yticks, yticklabels = plotting.get_auto_y_ticks(ybounds[0], ybounds[1], log='y')
    ybounds = [yticks[0], yticks[-1]]
    return xbounds, ybounds, xticks, yticks, yticklabels

# ----------------------------------------------------------------------------------------
def read_gctree_csv(all_seqfos, label):
    gc_counts = {tk : set() for tk in ['all', 'skipped']}
    timestr = None
    if label.count('-') == 2:  # timepoint-spilt data e.g. gct-data-d20 (as opposed to e.g. gct-data)
        gstr, dstr, timestr = label.split('-')
        assert gstr == 'gct'
        assert dstr == 'data'
        assert timestr in all_timepoints
    print('      reading gctree node data from %s' % '%s/nextflow/results/merged-results/gctree-node-data.csv'%args.gcreplay_dir)
    with open('%s/nextflow/results/merged-results/gctree-node-data.csv'%args.gcreplay_dir) as cfile:
        reader = csv.DictReader(cfile)
        for line in reader:
            gcid = datautils.get_gcid(line['PR'].replace('PR', ''), line['HK_key_mouse'], line['HK_key_node'], line['HK_key_gc'])
            gc_counts['all'].add(gcid)
            if timestr is not None and rpmeta[gcid]['time'] != timestr:
                gc_counts['skipped'].add(gcid)
                continue
            if gcid not in all_seqfos:
                all_seqfos[gcid] = []
            affinity = line['delta_bind']  # _CGG_FVS_additive
            hdist = utils.hamming_distance(args.naive_seq, line['IgH_nt_sequence']+line['IgK_nt_sequence'])
            n_aa_muts = int(line['n_mutations_HC']) + int(line['n_mutations_LC'])
            # print('      %2d %s' % (n_aa_muts, utils.color_mutants(utils.ltranslate(args.naive_seq), utils.ltranslate(line['IgH_nt_sequence']+line['IgK_nt_sequence']), amino_acid=True)))
            abdn = int(line['abundance'])
            for iseq in range(max(1, abdn)):  # internal nodes have abundance 0, but want to add them once
                all_seqfos[gcid].append({'name' : '%s-%s-%d' % (gcid, line['name'], iseq),
                                         'base-name' : line['name'],
                                         'seq' : line['IgH_nt_sequence']+line['IgK_nt_sequence'],
                                         'n_muts' : hdist,
                                         'n_muts_aa' : n_aa_muts,
                                         'affinity' : None if affinity == '' else float(affinity),
                                         'abundance' : abdn,
                                         'sampled' : abdn > 0,
                                         })
    print('    kept %d / %d GCs%s (skipped %d)' % (len(all_seqfos), len(gc_counts['all']), '' if timestr is None else ' at time \'%s\''%timestr, len(gc_counts['skipped'])))

# ----------------------------------------------------------------------------------------
def nstr(node):
    return 'leaf' if node.is_leaf() else 'internal'

# ----------------------------------------------------------------------------------------
def isdata(label):
    return 'bst-data' in label or 'iqt-data' in label

# ----------------------------------------------------------------------------------------
def read_input_files(label):
    # ----------------------------------------------------------------------------------------
    def read_gcd_meta(sldir):
        # ----------------------------------------------------------------------------------------
        def mfname():  # have to look for old meta file name for backwards compatibility (copied form partis/bin/gcdyn-simu-run.py
            mfn = '%s/meta.csv' % sldir
            if os.path.exists(mfn):  # if new name exists, return that
                return mfn
            return '%s/leaf-meta.csv' % sldir  # otherwise return old name
        # ----------------------------------------------------------------------------------------
        mfos = {}
        with open(mfname()) as mfile:
            reader = csv.DictReader(mfile)
            for line in reader:
                itree = int(line['name'].split('-')[0])
                if not isdata(label) and args.n_max_simu_trees is not None and itree > args.n_max_simu_trees - 1:
                    print('    --n-max-simu-trees: breaking after reading leaf meta for %d trees' % itree)
                    break
                mfos[line['name']] = line
        return mfos
    # ----------------------------------------------------------------------------------------
    def check_seqfos_nodes():  # compare uids in seqfos vs those in dtree nodes
        # ----------------------------------------------------------------------------------------
        def check_gcn(gcn):
            sfo_ids = collections.Counter(s.get('base-name', s['name']) for s in all_seqfos[gcn])
            dtr_ids = collections.Counter(n.taxon.label for n in all_dtrees[gcn].preorder_node_iter())
            only_sfo = set(sfo_ids) - set(dtr_ids)
            if len(only_sfo) > 0:
                tcts['only-seqfo'].append(len(only_sfo))
            only_tree = set(dtr_ids) - set(sfo_ids)
            if len(only_tree) > 0:
                tcts['only-tree'].append(len(only_tree))
                internal_ids = [n.taxon.label for n in all_dtrees[gcn].preorder_internal_node_iter()]
                tcts['internal'].append(len(set(internal_ids) & only_tree))
            for cid in set(sfo_ids) & set(dtr_ids):
                if sfo_ids[cid] != dtr_ids[cid]:
                    tcts['diff-count'].append('%d %d' % (sfo_ids[cid], dtr_ids[cid]))
        # ----------------------------------------------------------------------------------------
        tcts = {'diff-count' : [], 'only-seqfo' : [], 'only-tree' : [], 'internal' : []}
        for gcn in all_seqfos:
            check_gcn(gcn)
        print('     checking tree vs seqfo ids for %d trees:' % len(all_seqfos))
        if len(tcts['only-seqfo']) > 0:
            print('        %d seqs (from %d trees) only in seqfos: %s' % (sum(tcts['only-seqfo']), len(tcts['only-seqfo']), ' '.join(str(c) for c in tcts['only-seqfo'])))
        if len(tcts['only-tree']) > 0:
            print('        %d seqs (from %d trees) only in tree (%d of these are internal nodes): %s' % (sum(tcts['only-tree']), len(tcts['only-tree']), sum(tcts['internal']), ' '.join(str(c) for c in tcts['only-tree'])))
        if len(tcts['diff-count']) > 0:
            def kfn(s): return int(s.split()[0])
            print('        %d seqs had different counts in sfo vs tree: %s' % (len(tcts['diff-count']), '   '.join(sorted(tcts['diff-count'], reverse=True, key=kfn))))
    # ----------------------------------------------------------------------------------------
    def read_gctree_data(all_seqfos, plotvals, n_missing, n_tot):
        print('    reading gctree data from %s' % args.gcreplay_dir)
        read_gctree_csv(all_seqfos, label)  # read seqs plus affinity and mutation info from csv file (still have to read trees below to get leaf/internal info)
        too_small = []
        for gcn in all_seqfos:
            if args.min_seqs_per_gc is not None and len(all_seqfos[gcn]) < args.min_seqs_per_gc:  # NOTE not subsampling with args.max_seqs_per_gc (can't be bothered, it seems more important for abundance stuff anyway)
                too_small.append(gcn)
                continue
            gpm = rpmeta[gcn]
            gcidstr = datautils.get_gcid(gpm['PR'], gpm['mouse'], gpm['node'], gpm['gc'])
            gctree_dir = '%s/nextflow/results/gctrees/%s' % (args.gcreplay_dir, gcidstr)
            dtree = treeutils.get_dendro_tree(treefname='%s/gctree.inference.1.nk'%gctree_dir)
            all_dtrees[gcn] = dtree
            nodefo = {n.taxon.label : n for n in dtree.preorder_node_iter()}
            for sfo in all_seqfos[gcn]:
                snode = nodefo[sfo['base-name']]
                n_tot[nstr(snode)].append(sfo['base-name'])
                if sfo['affinity'] is None:
                    n_missing[nstr(snode)].append(sfo['base-name'])
                    continue
                for tk in plotvals:  # atm, fills affinity, n_muts, n_muts_aa
                    plotvals[tk][nstr(snode)].append(sfo[tk])
                if snode is dtree.seed_node:  # check --naive-seq (should really just not have it as an arg)
                    if sfo['seq'] != args.naive_seq:
                        raise Exception('--naive seq doesn\'t match seq for root node from data tree')
        if sum(len(l) for l in n_missing.values()) > 0:
            print('      %s missing/none affinity values for: %d / %d leaves, %d / %d internal' % (utils.wrnstr(), len(n_missing['leaf']), len(n_tot['leaf']), len(n_missing['internal']), len(n_tot['internal'])))
        if len(too_small) > 0:
            for gcn in too_small:
                del all_seqfos[gcn]
            print('    --min-seqs-per-gc: skipped %d / %d gcs with fewer than %d seqs: %s' % (len(too_small), len(all_seqfos), args.min_seqs_per_gc, ' '.join(sorted(too_small))))
        n_trees = len(all_seqfos) - len(too_small)
        return n_trees
    # ----------------------------------------------------------------------------------------
    def read_simu_like_files(all_seqfos, plotvals, n_missing, n_tot):
        # ----------------------------------------------------------------------------------------
        def get_affy(plotvals, label, dendro_trees, mfos):
            n_trees = 0
            for itree, dtree in enumerate(dendro_trees):
                if 'data' not in label and args.n_max_simu_trees is not None and itree > args.n_max_simu_trees - 1:
                    print('    --n-max-simu-trees: breaking after reading affinity for %d trees' % itree)
                    break
                n_trees += 1
                for node in dtree.preorder_node_iter():
                    name = node.taxon.label
                    n_tot[nstr(node)].append(name)
                    if name not in mfos or mfos[name]['affinity'] is None:
                        n_missing[nstr(node)].append(name)
                        continue
                    for tk in plotvals:
                        cfn = float if tk=='affinity' else int
                        if tk not in mfos[name]:
                            continue
                        plotvals[tk][nstr(node)].append(cfn(mfos[name][tk]))
            if sum(len(l) for l in n_missing.values()) > 0:
                print('      %s missing/none affinity/n_muts values for: %d / %d leaves, %d / %d internal' % (utils.wrnstr(), len(n_missing['leaf']), len(n_tot['leaf']), len(n_missing['internal']), len(n_tot['internal'])))
            return n_trees
        # ----------------------------------------------------------------------------------------
        def scale_affinities(antn_list, mfos):
            naive_id_affys = [(u, a) for l in antn_list for u, a, n in zip(l['unique_ids'], l['affinities'], l['n_mutations']) if n==0]
            if len(naive_id_affys) == 0:
                print('  %s no unmutated seqs in annotation, using default naive affinity %f' % (utils.wrnstr(), args.default_naive_affinity))  # the naive seq seems to always be in there, i guess cause i'm sampling intermediate common ancestors
                naive_affy = args.default_naive_affinity
            else:
                _, naive_affy = naive_id_affys[0]
            affy_std = numpy.std([m['affinity'] for m in mfos.values()], ddof=1)
            print('    rescaling to naive mean %.3f and std 1 (dividing by current std %.4f)' % (naive_affy, affy_std))
            for mfo in mfos.values():
                mfo['affinity'] = (mfo['affinity'] - naive_affy) / affy_std
        # ----------------------------------------------------------------------------------------
        if isdata(label):
            bdirs = {'bst' : args.beast_dir, 'iqt' : args.iqtree_data_dir}
            tstr = 'combo' if label.count('-')==1 else label.split('-')[-1]
            assert tstr in ['combo'] + all_timepoints
            sldir = '%s/%s-trees' % (bdirs[label.split('-')[0]], tstr)
        elif 'simu' in label:
            sldir = args.simu_like_dir
            if '-iqtree' in label:
                sldir += '/iqtree'
        else:
            assert False
        print('  reading %s from %s' % (('%s data'%label.split('-')[0]) if isdata(label) else 'simulation', sldir))
        if args.bcr_phylo:
            glfo, antn_list, _ = utils.read_output('%s/fake-paired-annotations.yaml' % sldir, dont_add_implicit_info=True)
            mfos, tmp_seqfos = {}, []
            for itn, atn in enumerate(antn_list):
                if 'data' not in label and args.n_max_simu_trees is not None and itn > args.n_max_simu_trees - 1:
                        print('    --n-max-simu-trees: breaking after reading annotations for %d trees' % itn)
                        break
                if atn.get('is_fake_paired') == True:
                    atn['seqs'] = atn['input_seqs']  # UGH
                for iseq, (uid, n_muts, affy) in enumerate(zip(atn['unique_ids'], atn['n_mutations'], atn['affinities'])):
                    assert uid not in mfos  # jeez i really hope there aren't repeated uids in different trees
                    mfos[uid] = {'n_muts' : n_muts, 'n_muts_aa' : utils.shm_aa(atn, iseq=iseq), 'affinity' : affy}
                tnsfos = utils.seqfos_from_line(atn, use_input_seqs=True) #, extra_keys=['n_mutations'])
                for sfo in tnsfos:
                    sfo['gcn'] = itn
                tmp_seqfos += tnsfos
            dendro_trees = [treeutils.get_dendro_tree(treestr=l['tree']) for l in antn_list]
            scale_affinities(antn_list, mfos)
        else:
            mfos = read_gcd_meta(sldir)  # this applies args.n_max_simu_trees
            tmp_seqfos = utils.read_fastx('%s/seqs.fasta'%sldir, queries=None if isdata(label) or args.n_max_simu_trees is None else mfos.keys())
            dendro_trees = [treeutils.get_dendro_tree(treestr=s) for s in treeutils.get_treestrs_from_file('%s/trees.nwk'%sldir, n_max_trees=None if isdata(label) else args.n_max_simu_trees)]
            if isdata(label):
                gcids = [l['gcid'] for l in utils.csvlines('%s/gcids.csv'%sldir)]
        # loop through all seqfos, setting n_muts from meta info and adding to correct gcn in all_seqfos
        for sfo in tmp_seqfos:
            if 'naive' in sfo['name']:
                continue
            sfo['n_muts'] = int(mfos[sfo['name']]['n_muts'])
            if 'n_muts_aa' not in mfos[sfo['name']]:  # gcdyn simu we need to calculate it
                mfos[sfo['name']]['n_muts_aa'] = utils.hamming_distance(utils.ltranslate(sfo['seq']), args.naive_seq_aa, amino_acid=True)
            sfo['n_muts_aa'] = int(mfos[sfo['name']]['n_muts_aa'])
            sfo['n_stops'] = utils.is_there_a_stop_codon(sfo['seq'], '', '', 0, return_n_stops=True)
            if isdata(label):  # for beast + iqtree data, we write the gcid to the leaf meta file (it's got dashes, so too hard to add it to sequence id)
                gcn = mfos[sfo['name']]['gc']
            else:
                gcn = sfo['gcn'] if args.bcr_phylo else sfo['name'].split('-')[0]
            if gcn not in all_seqfos:
                all_seqfos[gcn] = []
            all_seqfos[gcn].append(sfo)
        # put trees into all_dtrees and fill 'sampled'
        for itr, dtree in enumerate(dendro_trees):
            tnode = list(dtree.leaf_node_iter())[0]
            gcn = gcids[itr] if isdata(label) else (itr if args.bcr_phylo else tnode.taxon.label.split('-')[0])
            if gcn not in all_seqfos:
                print('        gcn \'%s\' not in all_seqfos (probably skipped with --n-max-simu-trees), so not processing its tree' % gcn)
                continue
            all_dtrees[gcn] = dtree
            ndict = {n.taxon.label : n for n in dtree.preorder_node_iter()}
            for sfo in all_seqfos[gcn]:
                sfo['sampled'] = ndict[sfo['name']].is_leaf()  # in simulation we assume leaves are sampled, whereas in data (above) we use abundance (gctree output) (beast data, which we handle here, seems to always put sampled seqs as leaves [unlike gctree])
        n_trees = get_affy(plotvals, label, dendro_trees, mfos)
        return n_trees
    # ----------------------------------------------------------------------------------------
    def get_abundance_info(all_seqfos):
        abdnvals, max_abvals, n_leaf_fos = [], [], 0
        unique_partition = []
        for gcn in all_seqfos:
            gcabvals = []
            leaf_fos = [s for s in all_seqfos[gcn] if s['sampled']]
            def kfn(s): return s['seq']
            tclust = []
            for tseq, sgroup in itertools.groupby(sorted(leaf_fos, key=kfn), key=kfn):
                abundance = len(list(sgroup))
                gcabvals.append(abundance)
                tclust.append(tseq)
            unique_partition.append(tclust)
            abdnvals += gcabvals
            max_abvals.append(max(gcabvals))
            n_leaf_fos += len(leaf_fos)
        return abdnvals, max_abvals, n_leaf_fos, unique_partition
    # ----------------------------------------------------------------------------------------
    all_seqfos, all_dtrees = collections.OrderedDict(), collections.OrderedDict()
    plotvals = {t : {k : [] for k in ['leaf', 'internal']} for t in ['affinity', 'n_muts', 'n_muts_aa']}
    n_missing, n_tot = {'internal' : [], 'leaf' : []}, {'internal' : [], 'leaf' : []}  # per-seq (not per-gc) counts
    if 'data' in label and 'bst' not in label and 'iqt' not in label:
        n_trees = read_gctree_data(all_seqfos, plotvals, n_missing, n_tot)
    elif 'simu' in label or isdata(label):  # beast + iqtree data are formatted like simulation
        n_trees = read_simu_like_files(all_seqfos, plotvals, n_missing, n_tot)
    else:
        assert False

    # NOTE that trees from gctree are genotype-collapsed, whereas those from simulation + beast have different nodes for duplicate seqs
    #   so *always* use seqfos for anything to do with abundance/N muts, etc, only use tree nodes for checking leaf/internal, topology, etc

    partition = [[s['name'] for s in sfos if s['sampled']] for sfos in all_seqfos.values()]

    check_seqfos_nodes()

    hists = {}
    lblstr = plotting.legends.get(label, label)

    abdnvals, max_abvals, n_leaf_fos, unique_partition = get_abundance_info(all_seqfos)
    htmp = hutils.make_hist_from_list_of_values(abdnvals, 'int', 'abundances')
    htmp.title = lblstr if args.short_legends else '%s (%d nodes in %d trees)' % (lblstr, n_leaf_fos, n_trees)
    htmp.xtitle = pltlabels['abundances']
    hists['abundances'] = {'distr' : htmp}

    htmp = hutils.make_hist_from_list_of_values(max_abvals, 'int', 'max-new-abundances')
    htmp.title = '%s (%d trees)' % (lblstr, n_trees)
    htmp.xtitle = pltlabels['max-abundances']
    hists['max-abundances'] = {'distr' : htmp}

    for pkey, pvals in plotvals['affinity'].items():
        htmp = Hist(xmin=-15, xmax=10, n_bins=100, value_list=pvals, title=lblstr, xtitle=pltlabels['%s-affinity'%pkey])  # NOTE don't try to get clever about xmin/xmax since they need to be the same for all different hists
        if not args.short_legends:
            htmp.title += ' (%d nodes in %d trees)' % (len(pvals), n_trees)
        hists['%s-affinity'%pkey] = {'distr' : htmp}

    xmin, xmax = 1, 150  # have to be the same for all labels, so easier to just set hard coded ones
    n_bins = 16
    dx = int((xmax - xmin) / n_bins)
    xbins = [l-0.5 for l in range(xmin, xmax + dx, dx)]
    hists['csizes'] = {'distr' : plotting.make_csize_hist(partition, n_bins=len(xbins), xbins=xbins, xtitle='N leaves')}
    hists['csizes']['distr'].title = lblstr

    hists['unique_csizes'] = {'distr' : plotting.make_csize_hist(unique_partition, n_bins=len(xbins), xbins=xbins, xtitle='N unique leaves')}
    hists['unique_csizes']['distr'].title = lblstr

    for mut_type in ['n_muts', 'n_muts_aa']:
        for tstr in ['leaf', 'internal']:
            mstr = '%s-muts-%s' % (tstr, 'aa' if '_aa' in mut_type else 'nuc')
            htmp = hutils.make_hist_from_list_of_values(plotvals[mut_type][tstr], 'int', mstr)
            htmp.xtitle = pltlabels[mstr]
            htmp.title = lblstr if args.short_legends else '%s (%d nodes in %d trees)' % (lblstr, len(plotvals[mut_type][tstr]), n_trees)
            hists[mstr] = {'distr' : htmp}

    n_stop_list = [s['n_stops'] for sfos in all_seqfos.values() for s in sfos]
    hists['n_stops'] = {'distr' : hutils.make_hist_from_list_of_values(n_stop_list, 'int', 'n_stops')}
    hists['n_stops']['distr'].xtitle = pltlabels['n_stops']

    return hists

# ----------------------------------------------------------------------------------------
# NOTE may be better to eventually switch to optimal transport rather than this weighted average/center of mass approach
# NOTE may also/instead want to use log of y difference (The max abundance seems ok, but the high tail of the abundance distr is getting totally washed out/overwhelmed by the low end)
def hist_distance(h1, h2, dbgstr='hist', weighted=False, debug=False):
    if debug:
        print('    %s distance%s:' % (dbgstr, ' (weighted)' if weighted else ''))
        print('      xval     v1      v2    abs diff')
    xvals = sorted(set(x for h in [h1, h2] for x in h.get_bin_centers()))
    dvals = []
    for xval in xvals:
        ib1, ib2 = [h.find_bin(xval) for h in [h1, h2]]
        if not h1.is_overflow(ib1) and not h2.is_overflow(ib2):
            if h1.low_edges[ib1] != h2.low_edges[ib2]:
                raise Exception('h1 low edge %.3f (ibin %d) doesn\'t equal h2 low edge %.3f (ibin %d)' % (h1.low_edges[ib1], ib1, h2.low_edges[ib2], ib2))
        v1, v2 = [0. if h.is_overflow(i) else h.bin_contents[i] for h, i in zip([h1, h2], [ib1, ib2])]  # NOTE ignores contents of over and underflow bins, which isn't great, but there shouldn't be any?
        dvals.append((xval if weighted else 1) * abs(v1 - v2))
        if debug:
            def fstr(v): return utils.color('blue', '-', width=6) if v==0 else '%6.2f'%v
            print('      %3.0f  %s  %s  %s' % (xval, fstr(v1), fstr(v2), fstr(dvals[-1])))
    if debug:
        print('    %.1f' % sum(dvals))
    return sum(dvals)

# ----------------------------------------------------------------------------------------
def compare_plots(htype, plotdir, hists, labels, hname, diff_vals, log='', irow=-1):
    if hname == 'n_stops':
        log = 'y'
    ytitle, plottitle = hists[0].ytitle, ''
    if args.normalize:
        # print('  %s I\'m not really sure it makes sense to normalize the mean hists (maybe could just skip them)' % utils.wrnstr())
        for htmp in hists:
            htmp.normalize()
        if 'fraction of' in hists[0].ytitle:  # NOTE normalize() call sets the ytitle
            ytitle = 'fraction of total'
    xbounds, ybounds, xticks, yticks, yticklabels = abdn_hargs(hists) if hname=='abundances' and log=='y' else (None, None, None, None, None)  # seems not to need this? hutils.multi_hist_filled_bin_xbounds(hists)
    text_dict = None if 'affinity' not in hists[0].xtitle or not args.bcr_phylo else {'x' : 0.2, 'y' : 0.6, 'text' : 'simu scaled to\nnaive mean, std 1'}
    adjust = None
    if '\n' in ytitle:
        adjust = {'left' : 0.23 if hname=='abundances' else 0.18}
    if any(isdata(l) for l in labels) and 'muts' in hname:
        xbounds = [-0.5, 20.5]
    if 'affinity' in hname: # and log == '':
        xbounds = [-9 if log=='y' else -3, 3]
    if 'internal' in hname or 'leaf' in hname:
        plottitle = '%s nodes' % ('internal' if 'internal' in hname else 'leaf')
    shift_overflows = True  # 'muts' in hname or 'affinity' in hname
    make_legend_only_plot, no_legend = args.write_legend_only_plots, args.write_legend_only_plots
    if args.legend_plots is not None:
        if hname in args.legend_plots:
            make_legend_only_plot, no_legend = args.write_legend_only_plots, args.write_legend_only_plots
        else:
            make_legend_only_plot, no_legend = False, True
    fn = plotting.draw_no_root(None, plotdir=plotdir, plotname='%s-%s%s'%(htype, hname, '' if log=='' else '-log'), more_hists=hists, log=log, xtitle=hists[0].xtitle, ytitle=ytitle, plottitle=plottitle,
                               bounds=xbounds, ybounds=ybounds, xticks=xticks, yticks=yticks, yticklabels=yticklabels, errors=htype!='max', square_bins=htype=='max', linewidths=[linewidths.get(l, 3) for l in labels],
                               alphas=[0.6 for _ in hists], colors=[colors.get(l) for l in labels], linestyles=[linestyles.get(l, '-') for l in labels], translegend=[-0.65, 0.1] if affy_like(hname) and 'abundance' not in hname else [-0.3, 0.1], write_csv=True,
                               hfile_labels=labels, text_dict=text_dict, adjust=adjust, remove_empty_bins='csize' in hname, make_legend_only_plot=make_legend_only_plot, no_legend=no_legend, shift_overflows=shift_overflows)
    fnames[irow].append(fn)

    # hdict = {l : h for l, h in zip(labels, hists)}
    # if hname != 'max-abdn-shm':  # min/max aren't the same for all hists, so doesn't work yet
    #     dname = '%s-%s'%(htype, hname)
    #     diff_vals[dname] = hist_distance(hdict['simu'], hdict['data'], weighted='abundances' in hname, dbgstr=dname)

# ----------------------------------------------------------------------------------------
ustr = """
Compare various distributions (especially abundances) in simu and gcreplay data.
NOTE that there's other scripts that process gcreplay results for partis input here: partis/datascripts/meta/taraki-gctree-2021-20
  replay-plot.py --simu-like-dir <dir> --outdir <dir>
"""
parser = argparse.ArgumentParser(usage=ustr)
parser.add_argument('--gcreplay-dir', default='/fh/fast/matsen_e/shared/replay/gcreplay', help='dir with gctree results on gcreplay data from which we read seqs, affinity, mutation info, and trees)')  # old location: /fh/fast/matsen_e/data/taraki-gctree-2021-10/gcreplay
parser.add_argument('--beast-dir', default='/fh/fast/matsen_e/data/taraki-gctree-2021-10/beast-processed-data/v6', help='dir with beast results on gcreplay data (same format as simulation)')
parser.add_argument('--iqtree-version')
parser.add_argument('--iqtree-data-dir', default='/fh/fast/matsen_e/data/taraki-gctree-2021-10/iqtree-processed-data', help='dir with iqtree results on gcreplay data (from datascripts/taraki-gctree-2021-10/iqtree-run.py then projects/gcdyn/scripts/data-parse.py')
parser.add_argument('--simu-like-dir', help='Dir from which to read simulation results, either from gcdyn or bcr-phylo (if the latter, set --bcr-phylo)')
parser.add_argument('--outdir')
parser.add_argument('--min-seqs-per-gc', type=int, help='if set, skip families/gcs with fewer than this many seqs NOTE doesn\'t [yet] apply to affinity plots')
# parser.add_argument('--max-seqs-per-gc', type=int, help='if set, downsample any families/gcs with more than this many seqs NOTE doesn\'t [yet] apply to affinity plots')
parser.add_argument('--plot-labels', default='gct-data-d15:gct-data-d20:gct-data-w10:simu', help='which/both of data/simu to plot')
parser.add_argument('--max-gc-plots', type=int, default=0, help='only plot individual (per-GC) plots for this  many GCs')
parser.add_argument('--normalize', action='store_true')
parser.add_argument('--short-legends', action='store_true', help='if set, don\'t add N trees/N nodes stuff to legends')
parser.add_argument('--bcr-phylo', action='store_true', help='set this if you\'re using bcr-phylo (rather than gcdyn) simulation')
parser.add_argument('--write-legend-only-plots', action='store_true', help='if set, instead of putting legends on each plot, write a separate legend-only svg for each plot')
parser.add_argument('--legend-plots', help='colon-separated list of plots on which to include legends (if set to \'all\', add legends to all plots)')
parser.add_argument('--naive-seq', default="GAGGTGCAGCTTCAGGAGTCAGGACCTAGCCTCGTGAAACCTTCTCAGACTCTGTCCCTCACCTGTTCTGTCACTGGCGACTCCATCACCAGTGGTTACTGGAACTGGATCCGGAAATTCCCAGGGAATAAACTTGAGTACATGGGGTACATAAGCTACAGTGGTAGCACTTACTACAATCCATCTCTCAAAAGTCGAATCTCCATCACTCGAGACACATCCAAGAACCAGTACTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGGGACTTCGATGTCTGGGGCGCAGGGACCACGGTCACCGTCTCCTCAGACATTGTGATGACTCAGTCTCAAAAATTCATGTCCACATCAGTAGGAGACAGGGTCAGCGTCACCTGCAAGGCCAGTCAGAATGTGGGTACTAATGTAGCCTGGTATCAACAGAAACCAGGGCAATCTCCTAAAGCACTGATTTACTCGGCATCCTACAGGTACAGTGGAGTCCCTGATCGCTTCACAGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAATGTGCAGTCTGAAGACTTGGCAGAGTATTTCTGTCAGCAATATAACAGCTATCCTCTCACGTTCGGCTCGGGGACTAAGCTAGAAATAAAA")
parser.add_argument('--n-max-simu-trees', type=int, help='stop after reading this many trees from simulation')
parser.add_argument("--random-seed", type=int, default=1)
parser.add_argument("--default-naive-affinity", type=float, default=1./100, help="this is the default for bcr-phylo, so maybe be correct if we don\'t have an unmutated sequence")
args = parser.parse_args()
all_choices = ['gct-data', 'gct-data-d15', 'gct-data-d20', 'gct-data-w10', 'bst-data', 'bst-data-d15', 'bst-data-d20', 'iqt-data', 'iqt-data-d20', 'simu', 'simu-iqtree']
args.plot_labels = utils.get_arg_list(args.plot_labels, choices=all_choices)
args.legend_plots = utils.get_arg_list(args.legend_plots)
if len(args.plot_labels) > 4 and not args.write_legend_only_plots:
    print('  note: setting --write-legend-only-plots since --plot-labels is longer than 4')
    args.write_legend_only_plots = True
args.naive_seq_aa = utils.ltranslate(args.naive_seq)
if args.iqtree_version is not None:
    args.iqtree_data_dir += '/' + args.iqtree_version

if not os.path.exists(args.iqtree_data_dir):
    raise Exception('<--iqtree-data-dir>/<--iqtree-version> %s doesn\'t exist (probably need to set --iqtree-version)' % args.iqtree_data_dir)

numpy.random.seed(args.random_seed)

rpmeta = datautils.read_gcreplay_metadata(args.gcreplay_dir)

def affy_like(tp):  # ick ( plots that get filled similarly to how affinity plots get filled, i.e. not how abundance-like stuff gets filled)
    return 'affinity' in tp or tp in ['csizes', 'unique_csizes', 'leaf-muts', 'internal-muts', 'abundances', 'max-abundances', 'n_stops']
abrows = [
    ['abundances', 'max-abundances', 'leaf-affinity', 'internal-affinity'],
    [],
    ['leaf-muts-nuc', 'internal-muts-nuc', 'leaf-muts-aa', 'internal-muts-aa'],
    [],
    ['csizes', 'unique_csizes', 'n_stops'],
    [],
]
abtypes = [a for arow in abrows if len(arow)>0 for a in arow]
hclists = {t : {'distr' : [], 'max' : []} for t in abtypes}
for tlab in args.plot_labels:
    print('  %s' % utils.color('blue', tlab))
    lhists = read_input_files(tlab)  # puts affinity hists in lhists
    for abtype in abtypes:
        for hn in lhists[abtype]:
            hclists[abtype][hn].append(lhists[abtype][hn])

cfpdir = '%s/plots/comparisons' % args.outdir
utils.prep_dir(cfpdir, wildlings=['*.csv', '*.svg'])
print('  plotting comparisons to %s' % cfpdir)
diff_vals = {}
fnames = [[] for _ in abrows]
for irow, arow in enumerate(abrows):
    if len(arow) == 0:  # log rows
        continue
    for hname in arow:
        for htype, hlist in hclists[hname].items():
            if hlist.count(None) == len(hlist):
                continue
            compare_plots(htype, cfpdir, hlist, args.plot_labels, hname, diff_vals, irow=irow)
            if 'affinity' in hname or 'abundance' in hname:
                compare_plots(htype, cfpdir, hlist, args.plot_labels, hname, diff_vals, log='y', irow=irow+1)

for ifl in range(len(fnames)):
    fnames.append([])
    for fnm in fnames[ifl]:
        lfn = fnm.replace('.svg', '-legend.svg')
        if os.path.exists(lfn):
            fnames[-1].append(lfn)

plotting.make_html(args.outdir+'/plots', fnames=fnames)
dfn = '%s/diff-vals.yaml' % args.outdir
with open(dfn, 'w') as dfile:
    json.dump(diff_vals, dfile)
