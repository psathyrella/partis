import sys
import time
import os
import string
import scipy
import numpy
import itertools
import subprocess
import warnings

import utils

# ----------------------------------------------------------------------------------------
def read_kmeans_clusterfile(clusterfname, seqfos, debug=False):

    # holy crap the need for this function [and its consequent form] make me angry

    all_uids = set([sfo['name'] for sfo in seqfos])
    partition = []
    with open(clusterfname) as clusterfile:
        lines = [l.strip() for l in clusterfile.readlines()]
        iline = -1
        while iline < len(lines) - 1:
            iline += 1

            clidline = lines[iline]
            if debug:
                print 'clid  ', clidline
            if clidline[0] != '$' or int(clidline.lstrip('$').strip('`')) != len(partition) + 1:
                raise Exception('couldn\'t convert %s to the expected cluster id %d' % (clidline, len(partition) + 1))
            partition.append([])

            while True:
                if iline + 2 >= len(lines):
                    break
                iline += 1
                uidline = lines[iline]
                if debug:
                    print 'uid   ', uidline
                iline += 1
                floatline = lines[iline]  # some info about the kmean cluster quality i think? don't care a.t.m.
                if debug:
                    print 'float ', floatline

                uids = set([u for u in uidline.split()])
                if len(uids - all_uids) > 0:
                    raise Exception('read unexpected uid[s] \'%s\' from %s' % (' '.join(uids - all_uids), clusterfname))
                all_uids -= uids
                partition[-1] += list(uids)

                floats = [float(istr) for istr in floatline.split()]
                if len(floats) != len(uids):
                    raise Exception('uid line %d and floats line %d have different lengths:\n  %s\n  %s' % (len(uids), len(floats), uidline, floatline))

                if lines[iline + 1] == '':
                    iline += 1
                    break

            # if emptyline != '':
            #     raise Exception('expected empty line but got \'%s\'' % emptyline)

    if len(all_uids) > 0:
        raise Exception('didn\'t read %d expected queries from %s (%s)' % (len(all_uids), clusterfname, ' '.join(all_uids)))

    os.remove(clusterfname)
    return partition

# ----------------------------------------------------------------------------------------
def read_component_file(mdsfname, n_components, seqfos):
    if not os.path.exists(mdsfname):  # r failed (probably with complex eigenvalues)
        os.remove(os.path.dirname(mdsfname) + '/run.r')  # this should really be integrated with utils.run_r()
        return {sfo['name'] : [0. for _ in range(n_components)] for sfo in seqfos}

    pc_names = None
    pcvals = {}
    with open(mdsfname) as mdsfile:
        for line in mdsfile:
            if pc_names is None:  # should be the first line
                pc_names = [pcn.strip() for pcn in line.split()]
                expected_names = ['PC%d' % i for i in range(1, n_components + 1)]
                if pc_names != expected_names:
                    raise Exception('expected components (%s) don\'t match those read from %s (%s)' % (' '.join(expected_names), mdsfname, ' '.join(pc_names)))
                continue
            valstrs = [vs.strip() for vs in line.split()]
            if len(valstrs) != n_components + 1:
                raise Exception('wrong number of columns (expected one uid and %d components) in %s: %s' % (n_components, mdsfname, line))
            pcvals[valstrs[0]] = [float(v) for v in valstrs[1:]]

    expected_uids = set([sfo['name'] for sfo in seqfos])
    found_uids = set(pcvals)
    if found_uids != expected_uids:
        if len(found_uids - expected_uids) > 0:
            raise Exception('  extra queries read from component file %s:\n %s' % (mdsfname, ' '.join(found_uids - expected_uids)))
        if len(expected_uids - found_uids) > 0:
            raise Exception('  queries missing from component file %s:\n %s' % (mdsfname, ' '.join(expected_uids - found_uids)))
        assert False  # no, it's not possible to get here. why do you ask?

    os.remove(mdsfname)
    return pcvals

# ----------------------------------------------------------------------------------------
def plot_mds(n_components, pcvals, plotdir, plotname, labels=None, partition=None, queries_to_include=None, gridsize=65, color_scale_vals=None, hexbin=False, title=None):
    import plotting
    # TODO switch to mpl_init/mpl_finalize
    import matplotlib
    from matplotlib import pyplot as plt
    # plt.rcParams['axes.facecolor'] = '#f1efef'
    colors = ['blue', 'forestgreen', 'red', 'grey', 'orange', 'green', 'skyblue', 'maroon', 'salmon', 'chocolate', 'magenta']
    single_color = '#4b92e7'
    def plot_component_pair(ipair, svgfname, color_map):
        fig = plt.figure(1)
        ax = plt.axes([0., 0., 1., 1.])
        # ax.set_facecolor('black')
        if hexbin:
            xvals, yvals = zip(*[(v[ipair], v[ipair + 1]) for v in pcvals.values()])
            hb = ax.hexbin(xvals, yvals, gridsize=gridsize, cmap=plt.cm.Blues, bins='log')
        else:
            for uid, vals in pcvals.items():
                plt.scatter(vals[ipair], vals[ipair + 1], color=color_map.get(uid, single_color), alpha=0.4 if labels is not None else 1.)

        if queries_to_include is not None:
            queries_to_include_in_this_cluster = set(pcvals) & set(queries_to_include)
            for uid in queries_to_include_in_this_cluster:
                xval, yval = pcvals[uid]
                ax.plot([xval], [yval], color='red', marker='.', markersize=10)
                ax.text(xval, yval, uid, color='red', fontsize=8)

        # smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        # smap.set_array([])
        # fig.colorbar(smap)
        if title is not None:
            # plt.title(title, fontweight='bold')  # wtf, doesn't work
            ax.text(ax.get_xlim()[0] + 0.5 * (ax.get_xlim()[1] - ax.get_xlim()[0]), 0.9 * ax.get_ylim()[1], title, color='red', fontsize=12)

        plt.savefig(svgfname)
        plt.close()

    if n_components % 2 != 0:
        print '%s odd number of components' % utils.color('red', 'warning')

    # set values in <color_map>
    color_map = {}
    if labels is not None or partition is not None:
        assert labels is None or partition is None  # should only specify one of them
        assert color_scale_vals is None
        if partition is None:
            def keyfunc(q):  # should really integrate this with utils.collapse_naive_seqs()/utils.split_partition_with_criterion()
                return labels[q]
            partition = [list(group) for _, group in itertools.groupby(sorted(pcvals, key=keyfunc), key=keyfunc)]
        if len(partition) > len(colors):
            raise Exception('more clusters/labels %d than colors %d' % (len(partition), len(colors)))
        color_map = {uid : colors[iclust] for iclust in range(len(partition)) for uid in partition[iclust]}  # just for coloring the plot
    elif color_scale_vals is not None:  # map with a number for each sequence (e.g. number of mutations) that we use to make a new color scale
        cmap, norm = plotting.get_normalized_cmap_and_norm([v for k, v in color_scale_vals.items() if k != '_naive'])
        color_map = {uid : cmap(norm(color_scale_vals[uid])) for uid in pcvals}

    for ipair in range(0, n_components - 1, 2):
        pcstr = '' if n_components == 2 else ('-pc-%d-vs-%d' % (ipair, ipair + 1))
        plot_component_pair(ipair, '%s/%s%s.svg' % (plotdir, plotname, pcstr), color_map)

# ----------------------------------------------------------------------------------------
def run_bios2mds(n_components, n_clusters, seqfos, base_workdir, seed, aligned=False, reco_info=None, region=None,
                 max_runs=100, max_iterations=1000, method='euclidean',
                 plotdir=None, plotname='mds', queries_to_include=None, color_scale_vals=None, labels=None, title=None, debug=False):
    workdir = base_workdir + '/mds'
    msafname = workdir + '/msa.fa'
    mdsfname = workdir + '/components.txt'
    clusterfname = workdir + '/clusters.txt'
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    if len(set([sfo['seq'] for sfo in seqfos])) < len(seqfos):  # it'll just crash when it's running mds later, but this is faster
        raise Exception('duplicate sequences in seqfos')

    if aligned:  # NOTE unlike the sklearn version below, this doesn't modify <seqfos>
        with open(msafname, 'w') as fastafile:
            for sfo in seqfos:
                fastafile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))
    else:
        utils.align_many_seqs(seqfos, outfname=msafname)

    # build the R cmd file
    cmdlines = [
        'options(rgl.useNULL=TRUE)',
        'require(bios2mds, quietly=TRUE)',
        'set.seed(%d)' % seed,
        'human <- import.fasta("%s")' % msafname,
        'active <- mat.dif(human, human)',  # mat.dif or mat.dis?
    ]

    if n_components is not None:
        cmdlines += ['mmds_active <- mmds(active, pc=%d)' % n_components]
        cmdlines += ['capture.output(mmds_active$coord, file="%s")' % mdsfname]
    else:
        raise Exception('need to implement')

    if n_clusters is not None:
        cmdlines += [
            'kmeans.run1 <- kmeans.run(mmds_active$coord, nb.clus=%d, nb.run=%d, iter.max=%d, method="%s")' % (n_clusters, max_runs, max_iterations, method),
            # 'kmeans.run1$clusters',
            # 'kmeans.run1$elements',
            'options(width=10000)',
            'capture.output(kmeans.run1$clusters, file="%s")' % clusterfname,
            # sil.score(mat, nb.clus = c(2:13), nb.run = 100, iter.max = 1000,  # run for every possible number of clusters (?)
            #               method = "euclidean")
            # random.msa  # builds a random [...]
        ]

    rstart = time.time()
    try:
        utils.run_r(cmdlines, workdir)  #, print_time='kmeans')
    except subprocess.CalledProcessError as e:  # typically happens because of complex eigenvalues
        print e
        print '   mds failed on cluster'  # NOTE will still crash in read_kmeans_clusterfile(), but I'm not using that a.t.m.
        title = (title if title is not None else '') + ' mds failed'
    pcvals = read_component_file(mdsfname, n_components, seqfos)
    partition = read_kmeans_clusterfile(clusterfname, seqfos) if n_clusters is not None else None
    rstop = time.time()

    os.remove(msafname)
    os.rmdir(workdir)

    plotstart = time.time()
    if plotdir is not None:
        # utils.prep_dir(plotdir, wildlings=['*.svg'])
        plot_mds(n_components, pcvals, plotdir, plotname, partition=partition if n_clusters is not None else None, queries_to_include=queries_to_include, color_scale_vals=color_scale_vals, labels=labels, title=title)
        if reco_info is not None:
            labels = {uid : reco_info[uid][region + '_gene'] for uid in pcvals}
            plot_mds(n_components, pcvals, plotdir, 'true-genes', labels=labels, queries_to_include=queries_to_include, color_scale_vals=color_scale_vals, title=title)
    print '    %5.1f  %5.1f' % (rstop - rstart, time.time() - plotstart),

    return partition

# ----------------------------------------------------------------------------------------
def run_sklearn_mds(n_components, n_clusters, seqfos, seed, reco_info=None, region=None, aligned=False, n_init=4, max_iter=300, eps=1e-3, n_jobs=-1, plotdir=None):
    assert n_clusters is not None
    if 'sklearn' not in sys.modules:
        from sklearn import manifold  # these are both slow af to import, even on local ssd
        from sklearn.cluster import KMeans

    if len(set(sfo['name'] for sfo in seqfos)) != len(seqfos):
        raise Exception('duplicate sequence ids in <seqfos>')

    print 'align'
    if not aligned:  # NOTE unlike the bios2mds version above, this modifies <seqfos>
        seqfos = utils.align_many_seqs(seqfos)

    print '  distances'
    # translations = string.maketrans('ACGT-', '01234')
    # def convert(seq):
    #     return [int(c) for c in seq.translate(translations)]
    # converted_seqs = [convert(x['seq']) for x in seqfos]
    # similarities = scipy.spatial.distance.pdist(converted_seqs, 'hamming')
    # similarities = scipy.spatial.distance.squareform(similarities)
    similarities = scipy.spatial.distance.squareform([utils.hamming_fraction(seqfos[i]['seq'], seqfos[j]['seq']) for i in range(len(seqfos)) for j in range(i + 1, len(seqfos))])

    print '  mds'
    random_state = numpy.random.RandomState(seed=seed)
    mds = sys.modules['sklearn'].manifold.MDS(n_components=n_components, n_init=n_init, max_iter=max_iter, eps=eps, random_state=random_state, dissimilarity="precomputed", n_jobs=n_jobs)
    pos = mds.fit_transform(similarities)
    # pos = mds.fit(similarities).embedding_

    print '  kmeans'
    kmeans = sys.modules['sklearn'].cluster.KMeans(n_clusters=n_clusters, random_state=random_state).fit(pos)
    pcvals = {seqfos[iseq]['name'] : pos[iseq] for iseq in range(len(seqfos))}
    labels = {seqfos[iseq]['name'] : kmeans.labels_[iseq] for iseq in range(len(seqfos))}
    def keyfunc(q):  # should really integrate this with utils.collapse_naive_seqs()/utils.split_partition_with_criterion()
        return labels[q]
    partition = [list(group) for _, group in itertools.groupby(sorted(pcvals, key=keyfunc), key=keyfunc)]

    if plotdir is not None:
        utils.prep_dir(plotdir, wildlings=['*.svg'])
        print '  plot'
        plot_mds(n_components, pcvals, plotdir, 'mds', partition=partition)

        if reco_info is not None:
            labels = {uid : reco_info[uid][region + '_gene'] for uid in pcvals}
            plot_mds(n_components, pcvals, plotdir, 'true-genes', labels=labels)

    return partition
