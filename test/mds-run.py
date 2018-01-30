#!/usr/bin/env python
import string
import sys
import scipy
import numpy as np
import argparse

from matplotlib import pyplot as plt

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

sys.path.insert(0, 'python')
import utils
import mds

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
    random_state = np.random.RandomState(seed=seed)
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

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--n-clusters', type=int, required=True)
parser.add_argument('--n-components', type=int, default=2)
parser.add_argument('--seed', type=int, default=1)
args = parser.parse_args()

seqfos = utils.read_fastx('v-qr.fa', n_max_queries=50)
run_sklearn_mds(seqfos, args.n_components, args.n_clusters, args.seed)
