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
from Bio.SeqUtils.CheckSum import seguid
import itertools
import pandas as pd
import numpy
import collections
import ast

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/projects', '')
sys.path.insert(1, partis_dir)

import python.utils as utils
import python.glutils as glutils
import python.plotting as plotting
import python.lbplotting as lbplotting
import python.hutils as hutils
from python.hist import Hist
import python.treeutils as treeutils
import python.paircluster as paircluster

colors = {
        'data' :  plotting.default_colors[1],
        'data-d15' : '#006600',
        'data-d20' : plotting.default_colors[1],
        'data-w10' : '#2b65ec',
        'simu' : '#808080',
    }
pltlabels = {'hdists' : 'root-tip dist', 'max-abdn-shm' : 'SHM in most\nabundant seq'}  # 'median SHM of seqs\nw/max abundance'}

# ----------------------------------------------------------------------------------------
def abfn(tlab, abtype='abundances'):
    return '%s/%s/%s.csv' % (args.outdir, tlab, abtype)

# ----------------------------------------------------------------------------------------
def get_gcid(prn, mouse, gcn):
    return 'pr-%s-m-%s-gc-%s' % (prn, mouse, gcn)

# ----------------------------------------------------------------------------------------
def read_gcreplay_metadata():
    # ----------------------------------------------------------------------------------------
    def readcfn(prn):
        with open('%s/metadata.PR%s.csv' % (args.gcreplay_dir, prn)) as cfile:
            reader = csv.DictReader(cfile)
            for line in reader:
                line['pr'] = prn
                if prn == '1':
                    line['time'] = line['imm_duration']
                    del line['imm_duration']
                elif prn == '2':
                    line['time'] = 'd15'
                    line['strain'] = 'wt'
                else:
                    assert False
                gcid = get_gcid(line['pr'], line['mouse'], line['gc'])
                assert gcid not in rpmeta
                rpmeta[gcid] = line
    # ----------------------------------------------------------------------------------------
    print('  reading gcreplay meta info from %s' % args.gcreplay_dir)
    rpmeta = {}
    for prn in ['1', '2']:
        readcfn(prn)
    return rpmeta

# ----------------------------------------------------------------------------------------
def write_abdn_csv(label, all_seqfos):  # summarize abundance (and other) info into csv files, for later reading (kind of weird to separate it like this, it's because this used to be in another script)
    # ----------------------------------------------------------------------------------------
    def process_family(fam_fos):
        n_seqs = len(fam_fos)
        counters['total'] += 1
        init_sizes.append(n_seqs)
        if args.min_seqs_per_gc is not None and n_seqs < args.min_seqs_per_gc:
            counters['too-small'] += 1
            return
        if args.max_seqs_per_gc is not None and n_seqs > args.max_seqs_per_gc:
            i_to_keep = numpy.random.choice(list(range(n_seqs)), size=args.max_seqs_per_gc)
            fam_fos = [fam_fos[i] for i in i_to_keep]
            final_sizes.append(len(fam_fos))  # don't actually use this number any more, but we use the len of final_sizes

        # This dictionary will map sequence checksum to the list of squence ids that have that
        # sequence checksum.
        ids_by_checksum = collections.defaultdict(list)
        hdvals, hd_dict = (
            [],
            {},
        )  # list of hamming distance to naive for each sequence (hd_dict is just for max_abdn below)
        sequence_count = 0
        for sfo in fam_fos:
            sequence_count = sequence_count + 1
            ids_by_checksum[seguid(sfo['seq'])].append(sfo['name'])
            hdvals.append(sfo['n_muts'])
            hd_dict[sfo['name']] = hdvals[-1]

        abundance_distribution = collections.defaultdict(int)
        for id_list in ids_by_checksum.values():
            id_count = len(id_list)
            abundance_distribution[id_count] = abundance_distribution[id_count] + 1
        assert sequence_count == sum(k * v for k, v in abundance_distribution.items())
        for amax, max_abdn_idlists in itertools.groupby(
            sorted(list(ids_by_checksum.values()), key=len, reverse=True), key=len
        ):
            # print(amax, [len(l) for l in max_abdn_idlists])
            break  # just want the first one (max abundance)

        base = hash(''.join(s['seq'] for s in fam_fos))
        abundances[base] = pd.Series(
            list(abundance_distribution.values()), index=list(abundance_distribution.keys())
        )

        fdicts["hdists"][base] = hdvals
        fdicts["max-abdn-shm"][base] = [
            int(numpy.median([hd_dict[u] for x in max_abdn_idlists for u in x]))
        ]

    # ----------------------------------------------------------------------------------------
    counters = {'too-small' : 0, 'total' : 0}
    init_sizes, final_sizes = [], []
    abundances = {}
    fdicts = {"hdists": {}, "max-abdn-shm": {}}
    for gcn, fam_fos in all_seqfos.items():
        process_family(fam_fos)

    print('    %d initial families with sizes: %s' % (len(init_sizes), " ".join(str(s) for s in sorted(init_sizes, reverse=True))))
    if counters['too-small'] > 0:
        print("      skipped %d / %d files with fewer than %d seqs" % (counters['too-small'], counters['total'], args.min_seqs_per_gc))
    if counters['too-small'] == counters['total']:
        raise Exception('skipped all families (see previous line)')
    if len(final_sizes) > 0:
        print("      downsampled %d samples to %d" % (len(final_sizes), args.max_seqs_per_gc))

    print("    writing %s to %s" % (label, os.path.dirname(abfn(label))))
    if not os.path.exists(os.path.dirname(abfn(label))):
        os.makedirs(os.path.dirname(abfn(label)))

    to_write = pd.DataFrame(abundances).fillna(0).astype(int)
    to_write.to_csv(abfn(label, 'abundances'))

    for fstr, fdct in fdicts.items():
        with open(abfn(label, abtype=fstr), 'w') as cfile:
            writer = csv.DictWriter(cfile, ["fbase", "vlist"])
            writer.writeheader()
            for fbase, vlist in fdicts[fstr].items():
                writer.writerow({"fbase": fbase, "vlist": ":".join(str(v) for v in vlist)})

# ----------------------------------------------------------------------------------------
def abdn_hargs(hlist):
    xmax = max(h.xmax for h in hlist)
    xbounds = [0.5, xmax+0.5]
    xticks = list(range(1, int(xmax)+1, 2 if xmax<10 else 5))
    ymin = min(h.get_minimum(exclude_empty=True) for h in hlist)
    ymax = max(h.get_maximum() for h in hlist)
    ybounds = [0.9 * ymin, 1.1 * ymax]
    yticks, yticklabels = plotting.get_auto_y_ticks(ybounds[0], ybounds[1])
    ybounds = [yticks[0], yticks[-1]]
    return xbounds, ybounds, xticks, yticks, yticklabels

# ----------------------------------------------------------------------------------------
def plot_abdn_stuff(lhists, plotdir, label, abtype):
    # read values from csv file
    with open(abfn(label, abtype=abtype)) as afile:
        reader = csv.DictReader(afile)
        if abtype == 'abundances':
            plotvals = {k : {} for k in reader.fieldnames if k!=''}
            for line in reader:
                abn = int(line[''])  # i can't figure out how to set this column label in the other script
                for bn in plotvals:  # <bn> is GC name
                    plotvals[bn][abn] = int(line[bn])
        else:
            plotvals = {}
            for line in reader:
                plotvals[line['fbase']] = [int(s) for s in line['vlist'].split(':')]

    # collect values from each GC
    max_vals = []
    distr_hists = []
    for bn in plotvals:  # bn is the base name of the fasta file (i.e. gc name)
        if abtype == 'abundances':
            max_vals.append(max([a for a, n in plotvals[bn].items() if n>0]))
            htmp = hutils.make_hist_from_dict_of_counts(plotvals[bn], 'int', bn)
        else:
            htmp = hutils.make_hist_from_list_of_values(plotvals[bn], 'int', bn) #, xmin_force=-0.5 if abtype=='max-abdn-shm' else None)
        if len(distr_hists) < args.max_gc_plots:
            nstxt = '%d seqs'%htmp.integral(True, multiply_by_bin_center=abtype=='abundances')
            mvtxt = 'mean %.1f' % htmp.get_mean()
            fn = htmp.fullplot(plotdir, '%s-distr-gc-%s'%(abtype, bn), pargs={'square_bins' : True, 'errors' : False, 'color' : colors[label]},
                               fargs={'title' : bn, 'xlabel' : pltlabels.get(abtype, abtype), 'ylabel' : 'counts', 'log' : 'y', 'title' : '%s: %s'%(label, bn)}, texts=[[0.6, 0.8, nstxt], [0.6, 0.75, mvtxt]])
            lbplotting.add_fn(fnames, fn, new_row=len(distr_hists)==0)
        distr_hists.append(htmp)

    hmax = None
    if abtype == 'abundances':
        hmax = hutils.make_hist_from_list_of_values(max_vals, 'int', 'max-abdn')
        hmax.title = '%s (%d trees)' % (label, len(distr_hists))
        hmax.xtitle = 'max abundance in GC'

    # plot mean distribution over GCs
    mean_hdistr = plotting.make_mean_hist(distr_hists)
    mean_hdistr.title = '%s (%d trees)' % (label, len(distr_hists))
    mean_hdistr.xtitle = pltlabels.get(abtype, abtype)
    mean_hdistr.ytitle = 'N seqs in bin\nmean+/-std (over GCs)'

    lhists[abtype] = {'distr' : mean_hdistr, 'max' : hmax}

# ----------------------------------------------------------------------------------------
def read_input_files(label):
    # ----------------------------------------------------------------------------------------
    def nstr(node):
        return 'leaf' if node.is_leaf() else 'internal'
    # ----------------------------------------------------------------------------------------
    def get_simu_affy(label, dendro_trees, affy_vals):
        for dtree in dendro_trees:
            for node in dtree.preorder_node_iter():
                name = node.taxon.label
                n_tot[nstr(node)].append(name)
                if affy_vals.get(name) is None:
                    n_missing[nstr(node)].append(name)
                    continue
                plotvals[nstr(node)].append(affy_vals[name])
        if sum(len(l) for l in n_missing.values()) > 0:
            print('      %s missing/none affinity values for: %d / %d leaves, %d / %d internal' % (utils.wrnstr(), len(n_missing['leaf']), len(n_tot['leaf']), len(n_missing['internal']), len(n_tot['internal'])))
        return plotvals
    # ----------------------------------------------------------------------------------------
    def read_data_csv(all_seqfos, label):
        gc_counts = {tk : set() for tk in ['all', 'skipped']}
        timestr = None
        if '-' in label:
            dstr, timestr = label.split('-')
            assert dstr == 'data'
            assert timestr in ['d20', 'w10', 'd15']
        print('      reading node data from %s' % '%s/nextflow/results/merged-results/gctree-node-data.csv'%args.gcreplay_dir)
        with open('%s/nextflow/results/merged-results/gctree-node-data.csv'%args.gcreplay_dir) as cfile:
            reader = csv.DictReader(cfile)
            for line in reader:
                prstr, prsuffix = line['PR'].split('.')
                assert prstr[:2] == 'PR'
                prn = prstr[2:]
                gcid = get_gcid(prn, line['HK_key_mouse'], line['HK_key_gc'])
                gc_counts['all'].add(gcid)
                if timestr is not None and rpmeta[gcid]['time'] != timestr:
                    gc_counts['skipped'].add(gcid)
                    continue
                if gcid not in all_seqfos:
                    all_seqfos[gcid] = []
                affinity = line['delta_bind']  # _CGG_FVS_additive
                hdist = utils.hamming_distance(args.naive_seq, line['IgH_nt_sequence']+line['IgK_nt_sequence'])  # this is different to nmuts in the next line, not yet sure why (UPDATE: one is nuc/other is AA, docs/column names maybe will be changed in gcreplay)
                # nmuts = int(line['n_mutations_HC']) + int(line['n_mutations_LC'])
                # print '  %2d  %2d  %s' % (hdist, nmuts, utils.color('red', '<--') if hdist!=nmuts else '')
                # print '      %s' % utils.color_mutants(NAIVE_SEQUENCE, line['IgH_nt_sequence']+line['IgK_nt_sequence'])
                for iseq in range(int(line['abundance'])):
                    all_seqfos[gcid].append({'name' : '%s-%s-%d' % (gcid, line['name'], iseq),
                                            'seq' : line['IgH_nt_sequence']+line['IgK_nt_sequence'],
                                            # 'n_muts' : nmuts,
                                            'n_muts' : hdist,
                                            'affinity' : None if affinity == '' else float(affinity),
                                            })
        print('    kept %d / %d GCs%s (skipped %d)' % (len(all_seqfos), len(gc_counts['all']), '' if timestr is None else ' at time \'%s\''%timestr, len(gc_counts['skipped'])))
    # ----------------------------------------------------------------------------------------
    def read_gcd_meta():
        mfos = {}
        with open('%s/leaf-meta.csv'%args.simu_dir) as mfile:
            reader = csv.DictReader(mfile)
            for line in reader:
                itree = int(line['name'].split('-')[0])
                if args.n_max_simu_trees is not None and itree > args.n_max_simu_trees - 1:
                        print('    --n-max-simu-trees: breaking after reading leaf meta for %d trees' % itree)
                        break
                mfos[line['name']] = line
        return mfos
    # ----------------------------------------------------------------------------------------
    def scale_affinities(antn_list, mfos):
        naive_id_affys = [(u, a) for l in antn_list for u, a, n in zip(l['unique_ids'], l['affinities'], l['n_mutations']) if n==0]
        if len(naive_id_affys) == 0:
            print('  %s no unmutated seqs in annotation, using default naive affinity %f' % (utils.wrnstr(), args.default_naive_affinity))  # the naive seq seems to always be in there, i guess cause i'm sampling intermediate common ancestors
            naive_affy = args.default_naive_affinity
        else:
            _, naive_affy = naive_id_affys[0]
        affy_std = numpy.std([m['affinity'] for m in mfos.values()], ddof=1)
        print('    rescaling to (naive) mean %.3f std %.4f' % (naive_affy, affy_std))
        for mfo in mfos.values():
            mfo['affinity'] = (mfo['affinity'] - naive_affy) / affy_std
    # ----------------------------------------------------------------------------------------
    all_seqfos = collections.OrderedDict()
    plotvals = {k : [] for k in ['leaf', 'internal']}
    n_missing, n_tot = {'internal' : [], 'leaf' : []}, {'internal' : [], 'leaf' : []}  # per-seq (not per-gc) counts
    if 'data' in label:
        print('    reading replay data from %s' % args.gcreplay_dir)
        read_data_csv(all_seqfos, label)  # read seqs plus affinity and mutation info from csv file (still have to read trees below to get leaf/internal info)
        n_too_small = 0
        for gcn in all_seqfos:
            if args.min_seqs_per_gc is not None and len(all_seqfos[gcn]) < args.min_seqs_per_gc:  # NOTE not subsampling with args.max_seqs_per_gc as in write_abdn_csv() (can't be bothered, it seems more important for abundance stuff anyway)
                n_too_small += 1
                continue
            gctree_dir = utils.get_single_entry(glob.glob('%s/nextflow/results/gctrees/PR%s*-%s-%s-%s-GC'%(args.gcreplay_dir, rpmeta[gcn]['pr'], rpmeta[gcn]['mouse'], rpmeta[gcn]['node'], rpmeta[gcn]['gc'])))
            dtree = treeutils.get_dendro_tree(treefname='%s/gctree.inference.1.nk'%gctree_dir)
            nodefo = {n.taxon.label : n for n in dtree.preorder_node_iter()}
            for sfo in all_seqfos[gcn]:
                n_gcid_dashes = get_gcid('a', 'b', 'c').count('-')
                # gcstr = sfo['name'].split('-')[ : n_gcid_dashes + 1]
                nname, iseq = sfo['name'].split('-')[n_gcid_dashes + 1 : ]
                n_tot[nstr(nodefo[nname])].append(nname)
                if sfo['affinity'] is None:
                    n_missing[nstr(nodefo[nname])].append(nname)
                    continue
                plotvals[nstr(nodefo[nname])].append(sfo['affinity'])
        if sum(len(l) for l in n_missing.values()) > 0:
            print('      %s missing/none affinity values for: %d / %d leaves, %d / %d internal' % (utils.wrnstr(), len(n_missing['leaf']), len(n_tot['leaf']), len(n_missing['internal']), len(n_tot['internal'])))
        if n_too_small > 0:
            print('    skipped %d / %d gcs with fewer than %d seqs' % (n_too_small, len(all_seqfos), args.min_seqs_per_gc))
        n_trees = len(all_seqfos) - n_too_small
    elif label == 'simu':
        print('  reading simulation from %s' % args.simu_dir)
        if args.bcr_phylo:
            glfo, antn_list, _ = utils.read_output('%s/fake-paired-annotations.yaml' % args.simu_dir, dont_add_implicit_info=True)
            mfos, tmp_seqfos = {}, []
            for itn, atn in enumerate(antn_list):
                if args.n_max_simu_trees is not None and itn > args.n_max_simu_trees - 1:
                        print('    --n-max-simu-trees: breaking after reading annotations for %d trees' % itn)
                        break
                for iseq, (uid, n_muts, affy) in enumerate(zip(atn['unique_ids'], atn['n_mutations'], atn['affinities'])):
                    assert uid not in mfos  # jeez i really hope there aren't repeated uids in different trees
                    mfos[uid] = {'n_muts' : n_muts, 'affinity' : affy}
                tnsfos = utils.seqfos_from_line(atn, use_input_seqs=True) #, extra_keys=['n_mutations'])
                for sfo in tnsfos:
                    sfo['gcn'] = itn
                tmp_seqfos += tnsfos
            dendro_trees = [treeutils.get_dendro_tree(treestr=l['tree']) for l in antn_list]
            scale_affinities(antn_list, mfos)
        else:
            mfos = read_gcd_meta()
            tmp_seqfos = utils.read_fastx('%s/seqs.fasta'%args.simu_dir, queries=None if args.n_max_simu_trees is None else mfos.keys())
            dendro_trees = [treeutils.get_dendro_tree(treestr=s) for s in treeutils.get_treestrs_from_file('%s/trees.nwk'%args.simu_dir, n_max_trees=args.n_max_simu_trees)]
        for sfo in tmp_seqfos:
            if 'naive' in sfo['name']:
                continue
            sfo['n_muts'] = int(mfos[sfo['name']]['n_muts'])
            gcn = sfo['gcn'] if args.bcr_phylo else sfo['name'].split('-')[0]
            if gcn not in all_seqfos:
                all_seqfos[gcn] = []
            all_seqfos[gcn].append(sfo)
        plotvals = get_simu_affy(label, dendro_trees, {u : float(mfos[u]['affinity']) for u in mfos})
        n_trees = len(dendro_trees)
    else:
        assert False

    write_abdn_csv(label, all_seqfos)

    hists = {}
    for pkey, pvals in plotvals.items():
        htmp = Hist(xmin=-15, xmax=10, n_bins=30, value_list=pvals, title=label, xtitle='%s affinity'%pkey)  # NOTE don't try to get clever about xmin/xmax since they need to be the same for all different hists
        htmp.title += ' (%d nodes in %d trees)' % (len(pvals), n_trees)
        hists['%s-affinity'%pkey] = {'distr' : htmp}
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
def compare_plots(hname, plotdir, hists, labels, abtype, diff_vals):
    ytitle = hists[0].ytitle
    if not args.dont_normalize:
        for htmp in hists:
            htmp.normalize()
        if 'fraction of' in hists[0].ytitle:
            ytitle = 'fraction of total'
    xbounds, ybounds, xticks, yticks, yticklabels = abdn_hargs(hists) if abtype=='abundances' else (None, None, None, None, None)  # seems not to need this? hutils.multi_hist_filled_bin_xbounds(hists)
    text_dict = None if 'affinity' not in hists[0].xtitle or not args.bcr_phylo else {'x' : 0.2, 'y' : 0.6, 'text' : 'simu scaled to\nnaive mean, std 1'}
    fn = plotting.draw_no_root(None, plotdir=plotdir, plotname='%s-%s'%(hname, abtype), more_hists=hists, log='y' if abtype=='abundances' else '', xtitle=hists[0].xtitle, ytitle=ytitle,
                               bounds=xbounds, ybounds=ybounds, xticks=xticks, yticks=yticks, yticklabels=yticklabels, errors=hname!='max', square_bins=hname=='max', linewidths=[4, 3],
                               plottitle='mean distr. over GCs' if 'N seqs in bin' in ytitle else '',  # this is a shitty way to identify the mean_hdistr hists, but best i can come up with atm
                               alphas=[0.6 for _ in hists], colors=[colors[l] for l in labels], linestyles=['--' if l=='simu' else '-' for l in labels], translegend=[-0.65, 0] if 'affinity' in abtype else [-0.2, 0], write_csv=True,
                               hfile_labels=labels, text_dict=text_dict)
    fnames[0].append(fn)

    # hdict = {l : h for l, h in zip(labels, hists)}
    # if abtype != 'max-abdn-shm':  # min/max aren't the same for all hists, so doesn't work yet
    #     dname = '%s-%s'%(hname, abtype)
    #     diff_vals[dname] = hist_distance(hdict['simu'], hdict['data'], weighted='abundances' in abtype, dbgstr=dname)

# ----------------------------------------------------------------------------------------
ustr = """
Compare various distributions (especially abundances) in simu and gcreplay data.
NOTE that there's other scripts that process gcreplay results for partis input here: partis/datascripts/meta/taraki-gctree-2021-20
  replay-plot.py --simu-dir <dir> --outdir <dir>
"""
parser = argparse.ArgumentParser(usage=ustr)
parser.add_argument('--gcreplay-dir', default='/fh/fast/matsen_e/data/taraki-gctree-2021-10/gcreplay', help='dir with gctree results on gcreplay data from which we read seqs, affinity, mutation info, and trees)')
parser.add_argument('--simu-dir', help='Dir from which to read simulation results, either from gcdyn or bcr-phylo (if the latter, set --bcr-phylo)')
parser.add_argument('--outdir')
parser.add_argument('--min-seqs-per-gc', type=int, help='if set, skip families/gcs with fewer than this many seqs NOTE doesn\'t [yet] apply to affinity plots')
parser.add_argument('--max-seqs-per-gc', type=int, help='if set, downsample any families/gcs with more than this many seqs NOTE doesn\'t [yet] apply to affinity plots')
parser.add_argument('--plot-labels', default='data-d15:data-d20:data-w10:simu', help='which/both of data/simu to plot')
parser.add_argument('--max-gc-plots', type=int, default=0, help='only plot individual (per-GC) plots for this  many GCs')
parser.add_argument('--dont-normalize', action='store_true')
parser.add_argument('--bcr-phylo', action='store_true', help='set this if you\'re using bcr-phylo (rather than gcdyn) simulation')
parser.add_argument('--naive-seq', default="GAGGTGCAGCTTCAGGAGTCAGGACCTAGCCTCGTGAAACCTTCTCAGACTCTGTCCCTCACCTGTTCTGTCACTGGCGACTCCATCACCAGTGGTTACTGGAACTGGATCCGGAAATTCCCAGGGAATAAACTTGAGTACATGGGGTACATAAGCTACAGTGGTAGCACTTACTACAATCCATCTCTCAAAAGTCGAATCTCCATCACTCGAGACACATCCAAGAACCAGTACTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGGGACTTCGATGTCTGGGGCGCAGGGACCACGGTCACCGTCTCCTCAGACATTGTGATGACTCAGTCTCAAAAATTCATGTCCACATCAGTAGGAGACAGGGTCAGCGTCACCTGCAAGGCCAGTCAGAATGTGGGTACTAATGTAGCCTGGTATCAACAGAAACCAGGGCAATCTCCTAAAGCACTGATTTACTCGGCATCCTACAGGTACAGTGGAGTCCCTGATCGCTTCACAGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAATGTGCAGTCTGAAGACTTGGCAGAGTATTTCTGTCAGCAATATAACAGCTATCCTCTCACGTTCGGCTCGGGGACTAAGCTAGAAATAAAA")
parser.add_argument('--n-max-simu-trees', type=int, help='stop after reading this many trees from simulation')
parser.add_argument("--random-seed", type=int, default=1, help="random seed for subsampling")
parser.add_argument("--default-naive-affinity", type=float, default=1./100, help="this is the default for bcr-phylo, so maybe be correct if we don\'t have an unmutated sequence")
args = parser.parse_args()
args.plot_labels = utils.get_arg_list(args.plot_labels, choices=['data', 'data-d15', 'data-d20', 'data-w10', 'simu'])

numpy.random.seed(args.random_seed)

rpmeta = read_gcreplay_metadata()

abtypes = ['leaf-affinity', 'internal-affinity', 'abundances', 'hdists', 'max-abdn-shm']
hclists = {t : {'distr' : [], 'max' : []} for t in abtypes}
fnames = [[], [], []]
for tlab in args.plot_labels:
    print('  %s' % utils.color('blue', tlab))
    utils.prep_dir('%s/plots/%s'%(args.outdir, tlab), wildlings=['*.csv', '*.svg'])
    lhists = read_input_files(tlab)  # puts affinity hists in lhists
    for abtype in abtypes:
        if 'affinity' not in abtype:
            plot_abdn_stuff(lhists, '%s/plots/%s'%(args.outdir, tlab), tlab, abtype)  # adds abdn hists to lhists
        for hn in lhists[abtype]:
            hclists[abtype][hn].append(lhists[abtype][hn])

cfpdir = '%s/plots/comparisons' % args.outdir
utils.prep_dir(cfpdir, wildlings=['*.csv', '*.svg'])
print('  plotting comparisons to %s' % cfpdir)
diff_vals = {}
for abtype in abtypes:
    for hname, hlist in hclists[abtype].items():
        if hlist.count(None) == len(hlist):
            continue
        compare_plots(hname, cfpdir, hlist, args.plot_labels, abtype, diff_vals)

if len(fnames[0]) > 4:
    fnames = [fnames[0][:4], fnames[0][4:]] + fnames[1:]
plotting.make_html(args.outdir+'/plots', fnames=fnames)
dfn = '%s/diff-vals.yaml' % args.outdir
with open(dfn, 'w') as dfile:
    json.dump(diff_vals, dfile)
