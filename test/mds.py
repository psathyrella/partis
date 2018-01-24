#!/usr/bin/env python
import string
import sys
import scipy
import numpy as np
import argparse

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

sys.path.insert(0, 'python')
import utils

# ----------------------------------------------------------------------------------------
def seqs(n_components, n_clusters, seed):
    # seqfos = [
    #     {'name' : 'a', 'seq' : 'CCGGTTGAACGGTA', 'color' : 'forestgreen'},
    #     {'name' : 'b', 'seq' : 'ACGGTTGAACGGTA', 'color' : 'red'},
    #     {'name' : 'c', 'seq' : 'ACTGTTGTACGGTA', 'color' : 'red'},
    #     # {'name' : 'd', 'seq' : 'ATATT-GTACGGTA', 'color' : 'black'},
    #     # # {'name' : 'e', 'seq' : 'TCGGTTGAACGGTA', 'color' : 'black'},
    # ]

    msafname = 'msa.fa'
    seqfos = utils.read_fastx('data/germlines/human/igh/ighv.fasta')
    utils.align_many_seqs(seqfos, outfname=msafname)
    seqfos = utils.read_fastx(msafname)

    print '  distances'
    # translations = string.maketrans('ACGT-', '01234')
    # def convert(seq):
    #     return [int(c) for c in seq.translate(translations)]
    # converted_seqs = [convert(x['seq']) for x in seqfos]
    # similarities = scipy.spatial.distance.pdist(converted_seqs, 'hamming')
    # similarities = scipy.spatial.distance.squareform(similarities)
    similarities = scipy.spatial.distance.squareform([utils.hamming_fraction(seqfos[i]['seq'], seqfos[j]['seq']) for i in range(len(seqfos)) for j in range(i + 1, len(seqfos))])

    print '  mds'
    random_state = np.random.RandomState(seed=3)
    mds = manifold.MDS(n_components=n_components, max_iter=1000, eps=1e-9, random_state=random_state, dissimilarity="precomputed", n_jobs=-1)
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

# ----------------------------------------------------------------------------------------
def example():
    n_samples = 20
    seed = np.random.RandomState(seed=3)
    X_true = seed.randint(0, 20, 2 * n_samples).astype(np.float)
    X_true = X_true.reshape((n_samples, 2))
    X_true -= X_true.mean()
    
    similarities = euclidean_distances(X_true)

    # Add noise to the similarities
    noise = np.random.rand(n_samples, n_samples)
    noise = noise + noise.T
    noise[np.arange(noise.shape[0]), np.arange(noise.shape[0])] = 0
    similarities += noise

    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed, dissimilarity="precomputed", n_jobs=-1)
    pos = mds.fit_transform(similarities)

    # Rescale the data
    pos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((pos ** 2).sum())

    # Rotate the data
    clf = PCA(n_components=2)
    X_true = clf.fit_transform(X_true)
    pos = clf.fit_transform(pos)

    fig = plt.figure(1)
    ax = plt.axes([0., 0., 1., 1.])

    plt.scatter(X_true[:, 0], X_true[:, 1], color='navy', lw=0, label='True Position')
    plt.scatter(pos[:, 0], pos[:, 1], color='turquoise', lw=0, label='MDS')
    plt.legend(scatterpoints=1, loc='best', shadow=False)

    # # Plot the edges
    # start_idx, end_idx = np.where(pos)
    # # a sequence of (*line0*, *line1*, *line2*), where::
    # #            linen = (x0, y0), (x1, y1), ... (xm, ym)
    # segments = [[X_true[i, :], X_true[j, :]]
    #             for i in range(len(pos)) for j in range(len(pos))]
    # values = np.abs(similarities)
    # lc = LineCollection(segments,
    #                     zorder=0, cmap=plt.cm.Blues,
    #                     norm=plt.Normalize(0, values.max()))
    # lc.set_array(similarities.flatten())
    # lc.set_linewidths(0.5 * np.ones(len(segments)))
    # ax.add_collection(lc)

    plt.savefig('tmp.svg')

parser = argparse.ArgumentParser()
parser.add_argument('--n-clusters', type=int, required=True)
parser.add_argument('--n-components', type=int, default=2)
parser.add_argument('--seed', type=int, default=1)
args = parser.parse_args()

# example()
seqs(n_components=args.n_components, n_clusters=args.n_clusters, seed=args.seed)
