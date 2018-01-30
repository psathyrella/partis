import sys
import os
import string
from matplotlib import pyplot as plt
import scipy
import numpy
from sklearn import manifold
from sklearn.metrics import euclidean_distances
# from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

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
colors = ['red', 'blue', 'forestgreen', 'grey', 'orange', 'green', 'skyblue4', 'maroon', 'salmon', 'chocolate4', 'magenta']

# ----------------------------------------------------------------------------------------
def kmeans_cluster(n_clusters, seqfos, all_qr_seqs, base_workdir, seed, reco_info=None, region=None, n_components=None, max_iterations=1000, max_runs=10, debug=False):
    # NOTE duplication in plotting fcn
    workdir = base_workdir + '/mds'
    msafname = workdir + '/msa.fa'
    mdsfname = workdir + '/components.txt'
    clusterfname = workdir + '/clusters.txt'
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    utils.align_many_seqs(seqfos, outfname=msafname)

    # build the R cmd file
    cmdlines = [
        'require(bios2mds, quietly=TRUE)',
        'set.seed(%d)' % seed,
        'human <- import.fasta("%s")' % msafname,
        'active <- mat.dif(human, human)',  # mat.dif or mat.dis?
    ]
    if n_components is not None:  # TODO now that you moved mmds to here, you need to move the plotting stuff into here (or, really, get it to output the PCA components and stop doing the plotting in R)
        cmdlines += ['mmds_active <- mmds(active, pc=%d)' % n_components]
        cmdlines += ['capture.output(mmds_active$coord, file="%s")' % mdsfname]
    else:
        raise Exception('need to implement')
    cmdlines += [  # TODO iter.max (max iteratoins and nb.run (max runs), also maybe add method (default "euclidean")
        'kmeans.run1 <- kmeans.run(mmds_active$coord, nb.clus=%d, iter.max=%d, nb.run=%d)' % (n_clusters, max_iterations, max_runs),
        # 'kmeans.run1$clusters',
        # 'kmeans.run1$elements',
        'options(width=10000)',
        'capture.output(kmeans.run1$clusters, file="%s")' % clusterfname,

        # sil.score(mat, nb.clus = c(2:13), nb.run = 100, iter.max = 1000,  # run for every possible number of clusters (?)
        #               method = "euclidean")
        # random.msa  # builds a random [...]
    ]

    utils.run_r(cmdlines, workdir, print_time='kmeans')
    pcvals = read_component_file(mdsfname, n_components, seqfos)
    partition = read_kmeans_clusterfile(clusterfname, seqfos)

    clusterfos = []
    for cluster in partition:
        cfo = {
            'seqfos' : [{'name' : uid, 'seq' : all_qr_seqs[uid]} for uid in cluster],
            # 'centroid'  # placeholder to remind you that vsearch clustering adds this, but I think it isn't subsequently used
        }
        cfo['cons_seq'] = utils.cons_seq(0.1, unaligned_seqfos=cfo['seqfos'])
        clusterfos.append(cfo)

        # debug = True
        # if debug and reco_info is not None:
        #     print len(cluster)
        #     for uid in cluster:
        #         print '    %s' % utils.color_gene(reco_info[uid][XXX region + '_gene'])  # need to pass in <region> if I want to uncomment this

    os.remove(msafname)
    os.rmdir(workdir)


# ----------------------------------------------------------------------------------------
    print '  plot'
    if n_components % 2 != 0:
        print '%s odd number of components' % utils.color('red', 'warning')

    if reco_info is not None:
        all_genes = list(set([reco_info[seqfo['name']][region + '_gene'] for seqfo in seqfos]))
        if len(all_genes) > len(colors):
            raise Exception('more genes %d than colors %d' % (len(all_genes), len(colors)))
        gene_colors = {all_genes[ig] : colors[ig] for ig in range(len(all_genes))}

    def plot_component_pair(ipair, plotname):
        fig = plt.figure(1)
        ax = plt.axes([0., 0., 1., 1.])
        for uid, vals in pcvals.items():
            plt.scatter(vals[ipair], vals[ipair + 1], color=colors[cluster_indices[uid]])
        # plt.scatter(pos[:, 0], pos[:, 1], color='forestgreen', lw=0, label='MDS')
        # plt.legend(scatterpoints=1, loc='best', shadow=False)
        plt.savefig(plotname)

    def plot_simu_component_pair(ipair, plotname):
        fig = plt.figure(1)
        ax = plt.axes([0., 0., 1., 1.])
        for seqfo in seqfos:
            vals = pcvals[seqfo['name']]
            gene = reco_info[seqfo['name']][region + '_gene']
            plt.scatter(vals[ipair], vals[ipair + 1], color=gene_colors[gene])
        plt.savefig(plotname)

    cluster_indices = {uid : iclust for iclust in range(len(partition)) for uid in partition[iclust]}  # just for coloring the plot
    for ipair in range(0, n_components - 1, 2):
        print '  %d' % ipair
        plot_component_pair(ipair, 'tmp-%d.svg' % ipair)
        plot_simu_component_pair(ipair, 'tmp-simu-%d.svg' % ipair)
# ----------------------------------------------------------------------------------------

    return clusterfos

# ----------------------------------------------------------------------------------------
def run_sklearn_mds(seqfos, n_components, n_clusters, seed, aligned=False, max_iter=1000, eps=1e-9, n_jobs=-1):
    print 'align'
    if not aligned:
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
    mds = manifold.MDS(n_components=n_components, max_iter=max_iter, eps=eps, random_state=random_state, dissimilarity="precomputed", n_jobs=n_jobs)
    pos = mds.fit_transform(similarities)  # hm, should this be mds.fit(similarities).embedding_?

    print '  kmeans'
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state).fit(pos)

    print '  plot'
    if n_components % 2 != 0:
        print '%s odd number of components' % utils.color('red', 'warning')

    def plot_component_pair(ipair, plotname):
        colors = ['forestgreen', 'red', 'blue', 'black', 'yellow', 'purple', 'orange']
        fig = plt.figure(1)
        ax = plt.axes([0., 0., 1., 1.])
        for iseq in range(len(seqfos)):
            plt.scatter(pos[iseq][ipair], pos[iseq][ipair + 1], color=colors[kmeans.labels_[iseq]])
        # plt.scatter(pos[:, 0], pos[:, 1], color='forestgreen', lw=0, label='MDS')
        # plt.legend(scatterpoints=1, loc='best', shadow=False)
        plt.savefig(plotname)

    for ipair in range(0, n_components - 1, 2):
        print '  %d' % ipair
        plot_component_pair(ipair, 'tmp-%d.svg' % ipair)

