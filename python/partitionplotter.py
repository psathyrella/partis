import os
import math
import numpy
import itertools
import sys
import time
import collections

from hist import Hist
import hutils
import utils
from clusterpath import ClusterPath
import mds
import treeutils

# ----------------------------------------------------------------------------------------
class PartitionPlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, args):
        self.args = args
        import plotting
        self.plotting = sys.modules['plotting']

        self.n_clusters_per_joy_plot = 50 if self.args.meta_info_key_to_color is None else 30
        self.n_max_joy_plots = 12
        self.n_max_mutations = 65
        self.n_joyplots_in_html = {'shm-vs-size' : self.n_max_joy_plots, 'overview' : 2}  # the rest of them are still in the dir, they're just not displayed in the html (note this is just
        self.min_high_mutation_cluster_size = 1
        self.n_biggest_to_plot = 24  # this functions as the number of plots for mds, sfs, and laplacian spectra NOTE this is overridden by self.args.queries_to_include (i.e. those queries always get plotted), but *not* by self.args.meta_info_to_emphasize
        self.n_plots_per_row = 4

        self.size_vs_shm_min_cluster_size = 3  # don't plot singletons and pairs for really big repertoires
        self.mds_max_cluster_size = 50000  # it's way tf too slow NOTE also max_cluster_size in make_mds_plots() (should combine them or at least put them in the same place)
        self.laplacian_spectra_min_clusters_size = 4
        self.max_clusters_to_apply_size_vs_shm_min_cluster_size = 500  # don't apply the previous thing unless the repertoire's actually pretty large

        self.n_mds_components = 2

    # ----------------------------------------------------------------------------------------
    def init_subd(self, subd, base_plotdir):
        plotdir = base_plotdir + '/' + subd
        utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])
        return subd, plotdir

    # ----------------------------------------------------------------------------------------
    def get_cdr3_title(self, annotation):
        naive_cdr3_seq, _ = utils.subset_iseq(annotation, 0, restrict_to_region='cdr3')
        title = ''
        if len(naive_cdr3_seq) % 3 != 0:
            # print '  out of frame: adding %s' % ((3 - len(naive_cdr3_seq) % 3) * 'N')
            naive_cdr3_seq += (3 - len(naive_cdr3_seq) % 3) * 'N'
            title += ' (out of frame)'
        title = utils.ltranslate(naive_cdr3_seq) + title
        return title

    # ----------------------------------------------------------------------------------------
    def make_single_hexbin_size_vs_shm_plot(self, sorted_clusters, annotations, repertoire_size, plotdir, plotname, log_cluster_size=False, repfrac_ylabel=True, debug=False):  # NOTE not using <repertoire_size> any more, but don't remember if there was a reason I should leave it
        # ----------------------------------------------------------------------------------------
        def add_emph_clusters(sorted_clusters):
            for cluster in sorted_clusters:
                if len(cluster) < self.size_vs_shm_min_cluster_size:
                    continue
                tqtis = set()  # queries to emphasize in this cluster
                if self.args.queries_to_include is not None:
                    tqtis |= set(cluster) & set(self.args.queries_to_include)
                if self.args.meta_info_to_emphasize is not None and self.any_meta_emph(annotations, cluster):
                    key, val = self.args.meta_info_to_emphasize.items()[0]
                    tqtis.add(utils.meta_emph_str(key, val, formats=self.args.meta_emph_formats))
                if len(tqtis) == 0:
                    continue
                xval, yval = numpy.mean(getnmutelist(cluster)), len(cluster)  # it kind of sucks to recalulate the x and y vals here, but by default (i.e. no emphasis) we don't keep track of which cluster goes with which x val, and it's nice to keep that simple
                if log_cluster_size:
                    yval = math.log(yval)
                ax.plot([xval], [yval], color='red', marker='.', markersize=5)
                ax.text(xval + 0.1, yval + (0.01 if log_cluster_size else 0.1), ' '.join(tqtis), color='red', fontsize=8)

        # ----------------------------------------------------------------------------------------
        import matplotlib.pyplot as plt
        def getnmutelist(cluster):
            return annotations[':'.join(cluster)]['n_mutations']

        fig, ax = self.plotting.mpl_init()

        clusters_to_use = [cluster for cluster in sorted_clusters if numpy.mean(getnmutelist(cluster)) < self.n_max_mutations]  # have to do it as a separate line so the zip/* don't crash if no clusters pass the criterion
        skipped_small_clusters = False
        if len(clusters_to_use) > self.max_clusters_to_apply_size_vs_shm_min_cluster_size:  # if repertoire is really big, ignore smaller clusters to keep the plots from being huge
            clusters_to_use = [cluster for cluster in clusters_to_use if len(cluster) >= self.size_vs_shm_min_cluster_size]
            skipped_small_clusters = True
        if len(clusters_to_use) == 0:
            print '  %s no clusters to plot for hexbin size vs shm' % utils.color('yellow', 'warning')
            return
        xvals, yvals = zip(*[[numpy.mean(getnmutelist(cluster)), len(cluster)] for cluster in clusters_to_use])
        if log_cluster_size:
            yvals = [math.log(yv) for yv in yvals]
        # hb = ax.hexbin(xvals, yvals, gridsize=self.n_max_mutations, cmap=plt.cm.Blues, bins='log')
        hb = ax.scatter(xvals, yvals, alpha=0.4)

        nticks = 5
        ymin, ymax = yvals[-1], yvals[0]
        ymin = 1  # make it 1, even if we aren't plotting small clusters, to make it more obvious that we skipped them
        yticks = [ymax + itick * (ymin - ymax) / float(nticks - 1) for itick in range(nticks)]
        if repfrac_ylabel:
            ytlfcn = lambda yt: utils.get_repfracstr(yt, repertoire_size)
        else:
            ytlfcn = lambda yt: ('%.0f' % yt) if yt > 5 else ('%.1f' % yt)
        yticklabels = [ytlfcn(math.exp(yt) if log_cluster_size else yt) for yt in yticks]

        # if necessary, add red labels to clusters
        if self.args.queries_to_include is not None or self.args.meta_info_to_emphasize is not None:
            add_emph_clusters(sorted_clusters)

        ylabel = ('family size\n(frac. of %d)' % repertoire_size) if repfrac_ylabel else 'clonal family size'
        if log_cluster_size:
            ylabel = '(log) ' + ylabel
            plotname += '-log'
        if skipped_small_clusters:
            fig.text(0.8, 0.25, 'skipping clusters\nsmaller than %d' % self.size_vs_shm_min_cluster_size, color='green', fontsize=8)
        self.plotting.mpl_finish(ax, plotdir, plotname, xlabel='mean N mutations', ylabel=ylabel, xbounds=(0, self.n_max_mutations), ybounds=(ymin, 1.05 * ymax), yticks=yticks, yticklabels=yticklabels)

    # ----------------------------------------------------------------------------------------
    def addfname(self, fnames, fname, force_new_row=False):
        fname += '.svg'
        if force_new_row or len(fnames[-1]) >= self.n_plots_per_row:
            fnames.append([fname])
        else:
            fnames[-1].append(fname)

    # ----------------------------------------------------------------------------------------
    def meta_emph(self, annotations, cluster, uid):  # return True if <uid> from <cluster> satisfies criteria in self.args.meta_info_to_emphasize
        antn = annotations[':'.join(cluster)]
        key, val = self.args.meta_info_to_emphasize.items()[0]
        return key in antn and utils.meta_info_equal(key, val, utils.per_seq_val(antn, key, uid), formats=self.args.meta_emph_formats)

    # ----------------------------------------------------------------------------------------
    def any_meta_emph(self, annotations, cluster):
        return any(self.meta_emph(annotations, cluster, u) for u in cluster)

    # ----------------------------------------------------------------------------------------
    def plot_this_cluster(self, sorted_clusters, iclust, annotations, plottype=None):
        if len(sorted_clusters[iclust]) == 1:
            return False
        if plottype == 'mds' and len(sorted_clusters[iclust]) > self.mds_max_cluster_size:
            print '     skipping mds plots for cluster with size %d > %d' % (len(sorted_clusters[iclust]), self.mds_max_cluster_size)
            return False
        if self.args.cluster_indices is not None and iclust not in self.args.cluster_indices:
            return False
        if iclust < self.n_biggest_to_plot:
            return True
        if self.args.queries_to_include is not None and len(set(self.args.queries_to_include) & set(sorted_clusters[iclust])) > 0:  # seed is added to <args.queries_to_include> in bin/partis
            return True
        return False  # falls through if <iclust> is too big, or if there's no --queries-to-include (which includes the seed)

    # ----------------------------------------------------------------------------------------
    def make_shm_vs_cluster_size_plots(self, sorted_clusters, annotations, base_plotdir, debug=False):
        def get_fname(iclustergroup=None, high_mutation=False, hexbin=False):
            if iclustergroup is not None:  # index of this group of clusters
                return 'size-vs-shm-%d' % iclustergroup
            elif high_mutation:
                return 'size-vs-shm-high-mutation'
            elif hexbin:
                return 'size-vs-shm-hexbin'
            else:
                assert False
        subd, plotdir = self.init_subd('shm-vs-size', base_plotdir)

        repertoire_size = sum([len(c) for c in sorted_clusters])
        cluster_indices = {':'.join(sorted_clusters[i]) : i for i in range(len(sorted_clusters))}  # index over all clusters, in the order that the mds plots will appear (compare to the two other indices I need within plotting.make_single_joyplot())

        # size vs shm joy plots
        iclustergroup = 0
        fnd = {'joy' : [], 'hex' : []}
        high_mutation_clusters = []
        sorted_cluster_groups = [sorted_clusters[i : i + self.n_clusters_per_joy_plot] for i in range(0, len(sorted_clusters), self.n_clusters_per_joy_plot)]
        if debug:
            print 'divided repertoire of size %d with %d clusters into %d cluster groups' % (repertoire_size, len(sorted_clusters), len(sorted_cluster_groups))
        for subclusters in sorted_cluster_groups:
            if iclustergroup > self.n_max_joy_plots:  # note that when this is activated, the high mutation plot is no longer guaranteed to have every high mutation cluster (but it should have every high mutation cluster that was bigger than the cluster size when we started skipping here)
                continue
            title = 'per-family SHM (%d / %d)' % (iclustergroup + 1, len(sorted_cluster_groups))  # NOTE it's important that this denominator is still right even when we don't make plots for all the clusters (which it is, now)
            high_mutation_clusters += self.plotting.make_single_joyplot(subclusters, annotations, repertoire_size, plotdir, get_fname(iclustergroup=iclustergroup), cluster_indices=cluster_indices, title=title, high_x_val=self.n_max_mutations,
                                                                        queries_to_include=self.args.queries_to_include, meta_info_to_emphasize=self.args.meta_info_to_emphasize, meta_info_key_to_color=self.args.meta_info_key_to_color, meta_emph_formats=self.args.meta_emph_formats,
                                                                        make_legend=self.args.meta_info_key_to_color is not None, debug=debug)  # have to make legend for every plot
            if len(fnd['joy']) < self.n_joyplots_in_html['shm-vs-size']:
                fnd['joy'].append(get_fname(iclustergroup=iclustergroup))
            iclustergroup += 1
        if self.args.meta_info_key_to_color is not None:
            fnd['leg'] = [get_fname(iclustergroup=0)+'-legend']
        if len(high_mutation_clusters) > 0 and len(high_mutation_clusters[0]) > self.min_high_mutation_cluster_size:
            high_mutation_clusters = [cluster for cluster in high_mutation_clusters if len(cluster) > self.min_high_mutation_cluster_size]
            self.plotting.make_single_joyplot(high_mutation_clusters, annotations, repertoire_size, plotdir, get_fname(high_mutation=True), plot_high_x=True, cluster_indices=cluster_indices, title='families with mean > %d mutations' % self.n_max_mutations,
                                              high_x_val=self.n_max_mutations, queries_to_include=self.args.queries_to_include, meta_info_to_emphasize=self.args.meta_info_to_emphasize, meta_info_key_to_color=self.args.meta_info_key_to_color, meta_emph_formats=self.args.meta_emph_formats,
                                              make_legend=self.args.meta_info_key_to_color is not None, debug=debug)
            fnd['high'] = [get_fname(high_mutation=True)]

        # size vs shm hexbin plots
        self.make_single_hexbin_size_vs_shm_plot(sorted_clusters, annotations, repertoire_size, plotdir, get_fname(hexbin=True))
        fnd['hex'].append(get_fname(hexbin=True))
        self.make_single_hexbin_size_vs_shm_plot(sorted_clusters, annotations, repertoire_size, plotdir, get_fname(hexbin=True), log_cluster_size=True)
        fnd['hex'].append(get_fname(hexbin=True) + '-log')

        fnames, rfnames = [[]], [[]]
        if 'leg' in fnd:
            self.addfname(fnames, fnd['leg'][0])
            self.addfname(rfnames, subd + '/' + fnd['leg'][0])
        for ifn, fn in enumerate(fnd['joy']):
            self.addfname(fnames, fn)
            if ifn < self.n_joyplots_in_html['overview']:
                self.addfname(rfnames, subd + '/' + fn)
        if 'high' in fnd:
            self.addfname(fnames, fnd['high'][0])
            self.addfname(rfnames, subd + '/' + fnd['high'][0])
        for ifn, fn in enumerate(fnd['hex']):
            self.addfname(fnames, fn, force_new_row=ifn==0)
            self.addfname(rfnames, subd + '/' + fn, force_new_row=ifn==0)

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames, new_table_each_row=True, extra_links=[('all shm-vs-size plots', subd)])

        return rfnames

    # ----------------------------------------------------------------------------------------
    def make_mds_plots(self, sorted_clusters, annotations, base_plotdir, max_cluster_size=10000, reco_info=None, color_rule=None, run_in_parallel=False, debug=False):
        debug = True
        # ----------------------------------------------------------------------------------------
        def get_fname(ic):
            return 'icluster-%d' % ic
        # ----------------------------------------------------------------------------------------
        def get_cluster_info(full_cluster, iclust):
            # ----------------------------------------------------------------------------------------
            def addseq(name, seq):
                found = False
                for sfo in seqfos:  # mds barfs if we have duplicate sequences, so if the sequence is already in there with a different name we just rename it (ick)
                    if sfo['seq'] == seq:
                        found = True
                        if sfo['name'] in tqtis:
                            name += '@(%s)' % tqtis[sfo['name']].lstrip('_')
                        sfo['name'] = name
                        break
                if not found:
                    seqfos.append({'name' : name, 'seq' : seq})
                color_scale_vals[name] = 0
                tqtis[name] = name  # i think if we're calling this fcn we always want to label with its actual name
            # ----------------------------------------------------------------------------------------
            def addqtis(tqtis, new_queries):  # NOTE that this can't do anything to prevent two points being on top of each other
                for uid, val in new_queries.items():
                    if uid in tqtis:
                        tqtis[uid] += '@%s' % val  # doesn't nicely condense duplicates like the joy plots do, but they should be much much less likely here (also note it can't use ',' since it gets passed on the command line)
                    else:
                        tqtis[uid] = '_%s' % val
            # ----------------------------------------------------------------------------------------
            full_info = annotations[':'.join(full_cluster)]
            # title = '%s   (size: %d)' % (self.get_cdr3_title(full_info), len(full_cluster))
            title = 'MDS comp. (index %d size %d)' % (iclust, len(full_cluster))

            kept_indices = list(range(len(full_cluster)))
            if len(kept_indices) > max_cluster_size:
                uids_to_choose_from = set([full_cluster[i] for i in kept_indices])  # note similarity to code in seqfileopener.post_process()
                if self.args.queries_to_include is not None:
                    uids_to_choose_from -= set(self.args.queries_to_include)
                if self.args.meta_info_to_emphasize is not None:
                    uids_to_choose_from -= set([u for u in uids_to_choose_from if self.meta_emph(annotations, full_cluster, u)])
                n_to_remove = len(kept_indices) - max_cluster_size
                if n_to_remove >= len(uids_to_choose_from):  # i.e. if we'd have to start removing queries that are in <self.args.queries_to_include>
                    removed_uids = uids_to_choose_from
                else:
                    removed_uids = numpy.random.choice(list(uids_to_choose_from), n_to_remove, replace=False)  # i think this'll still crash if len(uids_to_choose_from) is zero, but, meh
                kept_indices = sorted(set(kept_indices) - set([full_cluster.index(uid) for uid in removed_uids]))
                title += ' (subset: %d / %d)' % (len(kept_indices), len(full_cluster))

            seqfos = [{'name' : full_info['unique_ids'][iseq], 'seq' : full_info['seqs'][iseq]} for iseq in kept_indices]
            color_scale_vals = {full_cluster[iseq] : full_info['n_mutations'][iseq] for iseq in kept_indices}

            tqtis = {}
            if self.args.queries_to_include is not None:
                addqtis(tqtis, {u : u for u in self.args.queries_to_include})
            if self.args.meta_info_to_emphasize is not None:
                key, val = self.args.meta_info_to_emphasize.items()[0]
                addqtis(tqtis, {full_cluster[i] : utils.meta_emph_str(key, val, formats=self.args.meta_emph_formats) for i in kept_indices if self.meta_emph(annotations, full_cluster, full_cluster[i])})  # leading '_' is so dot doesn't cover up label
            addseq('_naive', full_info['naive_seq'])  # note that if any naive sequences that were removed above are in self.args.queries_to_include, they won't be labeled in the plot (but, screw it, who's going to ask to specifically label a sequence that's already specifically labeled?)
            addseq('_consensus', utils.cons_seq_of_line(full_info))  # leading underscore is 'cause the mds will crash if there's another sequence with the same name, and e.g. christian's simulation spits out the naive sequence with name 'naive'. No, this is not a good long term fix

            # remove duplicates, since they crash mds
            new_seqfos, all_seqs = [], {}
            for sfo in seqfos:
                # ----------------------------------------------------------------------------------------
                def add_to_label(nid, oid):  # add label for new id <nid> to any label for old id <oid>
                    if oid in tqtis:
                        tqtis[oid] += '@'
                        tqtis[nid] = tqtis[nid].lstrip('_')
                    else:
                        tqtis[oid] = ''
                    tqtis[oid] += tqtis[nid]
                    del tqtis[nid]
                # ----------------------------------------------------------------------------------------
                if sfo['seq'] in all_seqs:  # if we already added this seq with a different uid, we skip it
                    if sfo['name'] in tqtis:  # *and* if we need this seq to be labeled, then we have to add its label under the previous/other uid
                        add_to_label(sfo['name'], all_seqs[sfo['seq']])
                    continue
                new_seqfos.append(sfo)
                all_seqs[sfo['seq']] = sfo['name']
            seqfos = new_seqfos

            return seqfos, color_scale_vals, tqtis, title

        # ----------------------------------------------------------------------------------------
        def get_labels_for_coloring(full_cluster, color_rule):
            full_info = annotations[':'.join(full_cluster)]
            if color_rule == 'nearest-target':  # color by the index of the nearest cluster index (bcr-phylo simulation only)
                if 'target_seqs' not in reco_info[full_cluster[0]]:
                    return
                labels = {uid : str(reco_info[uid]['nearest_target_indices'][0]) for uid in full_cluster}
                labels['_naive'] = 'foop'
            elif color_rule == 'wtf':
                labels = {uid : uid.split('@')[1] for uid in full_cluster}
                labels['_naive'] = 'foop'
            else:
                assert False

            return labels

        # ----------------------------------------------------------------------------------------
        def prep_cmdfo(iclust, seqfos, tqtis, color_scale_vals, title):
            subworkdir = '%s/mds-%d' % (self.args.workdir, iclust)
            utils.prep_dir(subworkdir)
            tmpfname = '%s/seqs.fa' % subworkdir
            with open(tmpfname, 'w') as tmpfile:
                for sfo in seqfos:
                    csval = None
                    if sfo['name'] in color_scale_vals:
                        csval = color_scale_vals[sfo['name']]
                    tmpfile.write('>%s%s\n%s\n' % (sfo['name'], (' %d' % csval) if csval is not None else '' , sfo['seq']))
            cmdstr = '%s/bin/mds-run.py %s --aligned --plotdir %s --plotname %s --workdir %s --seed %d' % (utils.get_partis_dir(), tmpfname, plotdir, get_fname(iclust), subworkdir, self.args.random_seed)
            if tqtis is not None:
                cmdstr += ' --queries-to-include %s' % ':'.join(','.join([u, l]) for u, l in tqtis.items())
            if title is not None:
                cmdstr += ' --title=%s' % title.replace(' ', '@')
            return {'cmd_str' : cmdstr, 'workdir' : subworkdir, 'outfname' : '%s/%s.svg' % (plotdir, get_fname(iclust)), 'workfnames' : [tmpfname]}

        # ----------------------------------------------------------------------------------------
        subd, plotdir = self.init_subd('mds', base_plotdir)

        start = time.time()
        if debug:
            if not run_in_parallel:
                print '    making mds plots starting with %d clusters' % len(sorted_clusters)
                print '       size (+naive)   mds    plot   total'
        plotted_cluster_lengths, skipped_cluster_lengths = [], []
        fnames = [[]]
        cmdfos = []
        for iclust in range(len(sorted_clusters)):
            if not self.plot_this_cluster(sorted_clusters, iclust, annotations, plottype='mds'):
                skipped_cluster_lengths.append(len(sorted_clusters[iclust]))
                continue
            plotted_cluster_lengths.append(len(sorted_clusters[iclust]))

            seqfos, color_scale_vals, tqtis, title = get_cluster_info(sorted_clusters[iclust], iclust)

            labels = None
            if color_rule is not None:
                labels = get_labels_for_coloring(sorted_clusters[iclust], color_rule)
                # print '   %s setting color_scale_vals to None so we can use colors for nearest target seq index' % utils.color('red', 'note')
                color_scale_vals = None  # not sure this is really the best way to do this

            if debug and not run_in_parallel:
                substart = time.time()
                subset_str = '' if len(sorted_clusters[iclust]) <= max_cluster_size else utils.color('red', '/%d' % len(sorted_clusters[iclust]), width=6, padside='right')  # -1 is for the added naive seq
                tmpfo = annotations[':'.join(sorted_clusters[iclust])]
                # n_naive_in_cluster = len([iseq for iseq in range(len(sorted_clusters[iclust])) if tmpfo['n_mutations'][iseq] == 0])  # work out if there was a sequence already in the cluster that was the same as the naive sequence
                # print '      %4d%6s' % (len(seqfos) - 1 + n_naive_in_cluster, subset_str),
                print '      %4d%6s' % (len(seqfos), subset_str),

            if run_in_parallel:
                assert labels is None  # would need to implement this (or just switch to non-parallel version if you need to run with labels set)
                cmdfos.append(prep_cmdfo(iclust, seqfos, tqtis, color_scale_vals, title))
            else:
                mds.run_bios2mds(self.n_mds_components, None, seqfos, self.args.workdir, self.args.random_seed,
                                 aligned=True, plotdir=plotdir, plotname=get_fname(iclust),
                                 queries_to_include=tqtis, color_scale_vals=color_scale_vals, labels=labels, title=title)
                if debug:
                    print '  %5.1f' % (time.time() - substart)
            self.addfname(fnames, '%s' % get_fname(iclust))

        if run_in_parallel and len(cmdfos) > 0:
            utils.run_cmds(cmdfos, clean_on_success=True)  #, debug='print')

        if len(skipped_cluster_lengths) > 0:
            print '    mds: skipped %d clusters with lengths: %s' % (len(skipped_cluster_lengths), utils.cluster_size_str(skipped_cluster_lengths, only_passing_lengths=True))

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        print '      made %d mds plots (%.1fs) with sizes %s' % (sum(len(x) for x in fnames), time.time() - start, utils.cluster_size_str(plotted_cluster_lengths, only_passing_lengths=True))

        return [[subd + '/' + fn for fn in fnames[i]] for i in range(min(2, len(fnames)))]

    # ----------------------------------------------------------------------------------------
    def make_laplacian_spectra_plots(self, sorted_clusters, annotations, plotdir, cpath=None, debug=False):  # NOTE it's kind of weird to have this here, but all the other tree-dependent plotting in treeutils, but it's because this is for comparing clusters, whereas the stuff in treeutils is all about lb values, which are mostly useful within clusters
        subd, plotdir = self.init_subd('laplacian-spectra', plotdir)

        fnames = [[]]
        for iclust in range(len(sorted_clusters)):
            if not self.plot_this_cluster(sorted_clusters, iclust, annotations):
                continue
            annotation = annotations[':'.join(sorted_clusters[iclust])]
            if len(annotation['unique_ids']) < self.laplacian_spectra_min_clusters_size:
                continue
            if len(set(sorted_clusters[iclust])) < len(sorted_clusters[iclust]):
                repeated_uids = [u for u, count in collections.Counter(sorted_clusters[iclust]).items() if count > 1]
                print '  skipping laplacian spectra plotting for cluster with %d duplicate uids (%s)' % (len(repeated_uids), ' '.join(repeated_uids))
                continue
            if 'tree-info' in annotation and 'lb' in annotation['tree-info']:
                treestr = annotation['tree-info']['lb']['tree']
            else:  # if this is simulation, and calculate_tree_metrics() was called with use_true_clusters=True, then we probably have to get our own trees here for the actual clusters in the best partition
                treefo = treeutils.get_tree_for_inf_line(annotation, cpath=cpath, annotations=annotations, debug=debug)
                print '  %s no tree in annotation, so getting new tree from/with \'%s\'' % (utils.color('yellow', 'warning'), treefo['origin'])
                treestr = treefo['tree'].as_string(schema='newick').strip()
            treeutils.run_laplacian_spectra(treestr, plotdir=plotdir, plotname='icluster-%d' % iclust, title='size %d' % len(annotation['unique_ids']))
            if len(fnames[-1]) < self.n_plots_per_row:
                self.addfname(fnames, 'icluster-%d' % iclust)

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_sfs_plots(self, sorted_clusters, annotations, base_plotdir, restrict_to_region=None, debug=False):
        def addplot(oindexlist, ofracslist, n_seqs, fname, cdr3titlestr, red_text=None):
            hist = Hist(30, 0., 1.)
            for ofracs in ofracslist:
                hist.fill(ofracs)
            fig, ax = self.plotting.mpl_init()
            hist.mpl_plot(ax, remove_empty_bins=True)
            ax.text(0.65, 0.8 * ax.get_ylim()[1], 'h: %+.1f' % utils.fay_wu_h(line=None, restrict_to_region=restrict_to_region, occurence_indices=oindexlist, n_seqs=n_seqs), fontsize=17, fontweight='bold')
            if red_text is not None:
                ax.text(0.65, 0.7 * ax.get_ylim()[1], red_text, fontsize=17, color='red', fontweight='bold')
            titlestr = '%s (size: %d)' % (cdr3titlestr, n_seqs)

            regionstr = restrict_to_region + ' ' if restrict_to_region is not None else ''
            self.plotting.mpl_finish(ax, plotdir, fname, title=titlestr, xlabel=regionstr + 'mutation frequency', ylabel=regionstr + 'density of mutations', xticks=[0, 1], log='')  # xticks=[min(occurence_fractions), max(occurence_fractions)], 
            self.addfname(fnames, fname)

        subd, plotdir = self.init_subd('sfs', base_plotdir)

        fnames = [[]]
        for iclust in range(len(sorted_clusters)):
            if not self.plot_this_cluster(sorted_clusters, iclust, annotations):
                continue
            annotation = annotations[':'.join(sorted_clusters[iclust])]
            occurence_indices, occurence_fractions = utils.get_sfs_occurence_info(annotation, restrict_to_region=restrict_to_region)
            red_text = None
            assert args.meta_info_to_emphasize is None  # would need to be implemented
            if self.args.queries_to_include is not None and len(set(self.args.queries_to_include) & set(sorted_clusters[iclust])) > 0:
                red_text = '%s' % ' '.join(set(self.args.queries_to_include) & set(sorted_clusters[iclust]))
            addplot(occurence_indices, occurence_fractions, len(sorted_clusters[iclust]), 'icluster-%d' % iclust, self.get_cdr3_title(annotation), red_text=red_text)

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_cluster_size_distribution(self, base_plotdir, sorted_clusters, annotations):
        subd, plotdir = self.init_subd('sizes', base_plotdir)
        fname = 'cluster-sizes'
        hcolors = None
        cslist = [len(c) for c in sorted_clusters]
        # is_log_x = len(cslist) > 100 # and len([s for s in cslist if s>50]) > 30
        csize_hists = {'best' : hutils.make_hist_from_list_of_values(cslist, 'int', fname, is_log_x=True)}  # seems kind of wasteful to make a bin for every integer (as here), but it's not going to be *that* many, and we want to be able to sample from them, and it's always a hassle getting the bins you want UPDATE ok now sometimes using log bins, for aesthetic plotting reasons, but maybe also ok for sampling?
        self.plotting.plot_cluster_size_hists(plotdir, fname, csize_hists)
        fnlist = [subd + '/' + fname + '.svg']
        ytitle = None
        if self.args.meta_info_key_to_color is not None:  # plot mean fraction of cluster that's X for each cluster size
            fname = 'cluster-size-fractions'
            mekey = self.args.meta_info_key_to_color
            all_emph_vals, emph_colors = self.plotting.meta_emph_init(mekey, sorted_clusters, annotations, formats=self.args.meta_emph_formats)
            hcolors = {v : c for v, c in emph_colors}
            plotvals = {v : [] for v in all_emph_vals}  # for each possible value, a list of (cluster size, fraction of seqs in cluster with that val) for clusters that contain seqs with that value
            for csize, cluster_group in itertools.groupby(sorted_clusters, key=len):
                for tclust in cluster_group:
                    antn = annotations.get(':'.join(tclust))
                    if antn is None:
                        continue
                    def gfcn(x): return utils.meta_emph_str(mekey, utils.per_seq_val(antn, mekey, x, use_default=True), formats=self.args.meta_emph_formats)
                    vgroups = utils.group_seqs_by_value(tclust, gfcn, return_values=True)
                    emph_fracs = {v : len(grp) / float(csize) for v, grp in vgroups}
                    emph_fracs.update({v : 0. for v in all_emph_vals - set(emph_fracs)})  # need to include families with no seqs with a given value (otherwise all the lines go to 1 as cluster size goes to 1)
                    for v, frac in emph_fracs.items():
                        plotvals[v].append((csize, frac))
            bhist = csize_hists['best']
            csize_hists.update({v : Hist(template_hist=bhist) for v in all_emph_vals})  # for each possible value, a list of (cluster size, fraction of seqs in cluster with that val) for clusters that contain seqs with that value
            del csize_hists['best']
            for e_val, cvals in plotvals.items():
                ehist = csize_hists[e_val] #utils.meta_emph_str(mekey, e_val, formats=self.args.meta_emph_formats)]
                for ibin in ehist.ibiniter(include_overflows=True):
                    ib_vals = [f for s, f in cvals if ehist.find_bin(s)==ibin]  # fracs whose cluster sizes fall in this bin (should all be quite similar in size if our bins are sensible, so shouldn't need to do an average weighted for cluster size)
                    if len(ib_vals) == 0:
                        continue
                    mval = numpy.mean(ib_vals)
                    err = mval / math.sqrt(2) if len(ib_vals) == 1 else numpy.std(ib_vals, ddof=1) / math.sqrt(len(ib_vals))  # that isn't really right for len 1, but whatever
                    ehist.set_ibin(ibin, mval, err)
            ytitle = 'mean fraction of each cluster'

        self.plotting.plot_cluster_size_hists(plotdir, fname, csize_hists, hcolors=hcolors, ytitle=ytitle)
        for hname, thist in csize_hists.items():
            thist.write('%s/%s%s.csv' % (plotdir, fname, '' if hname=='best' else '-'+hname))
        fnlist.append(subd + '/' + fname + '.svg')
        return [fnlist]

    # ----------------------------------------------------------------------------------------
    def remove_failed_clusters(self, partition, annotations):
        # remove clusters with failed annotations
        failed_clusters = []
        for cluster in partition:
            if ':'.join(cluster) not in annotations:
                print '    %s cluster %s not in annotations' % (utils.color('red', 'warning'), ':'.join(cluster))
                failed_clusters.append(cluster)
        for fclust in failed_clusters:
            partition.remove(fclust)

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, partition, annotations, reco_info=None, cpath=None, no_mds_plots=False):
        if self.args.only_csv_plots:
            print '  --only-csv-plots not implemented for partition plots, so returning without plotting'
            return
        if not utils.check_cmd('R', options=['--slave', '--version'], return_bool=True):
            no_mds_plots = True
            print '  note: R does not seem to be installed, so skipping mds partition plots'
        print '  plotting partitions'
        sys.stdout.flush()
        start = time.time()

        fnames = []
        self.remove_failed_clusters(partition, annotations)
        sorted_clusters = sorted(partition, key=lambda c: len(c), reverse=True)
        fnames += self.make_shm_vs_cluster_size_plots(sorted_clusters, annotations, plotdir)
        if not no_mds_plots:
            fnames += self.make_mds_plots(sorted_clusters, annotations, plotdir, reco_info=reco_info, run_in_parallel=True) #, color_rule='wtf')
        # fnames += self.make_laplacian_spectra_plots(sorted_clusters, annotations, plotdir, cpath=cpath)
        # fnames += self.make_sfs_plots(sorted_clusters, annotations, plotdir)
        csfns = self.make_cluster_size_distribution(plotdir, sorted_clusters, annotations)
        fnames[0] += csfns[0]

        subdirs = ['shm-vs-size', 'mds'] #, 'laplacian-spectra']  # , 'sfs
        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames, new_table_each_row=True, htmlfname=plotdir + '/overview.html', extra_links=[(subd, '%s.html'%subd) for subd in subdirs])

        print '    partition plotting time: %.1f sec' % (time.time()-start)
