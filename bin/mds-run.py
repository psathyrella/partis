#!/usr/bin/env python
import string
import sys
import scipy
import numpy as np
import argparse
import random

from matplotlib import pyplot as plt

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

sys.path.insert(0, 'python')
import utils
import mds

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--n-clusters', type=int, required=True)
parser.add_argument('--n-components', type=int, default=2)
parser.add_argument('--plotdir')
parser.add_argument('--workdir', default='/tmp/dralph/mds/' + str(random.randint(0, 999999)))
parser.add_argument('--seed', type=int, default=1)
args = parser.parse_args()

seqfos = utils.read_fastx('v-qr.fa', n_max_queries=500)
for iseq in range(len(seqfos)):
    seqfos[iseq]['name'] = str(iseq)

# mds.run_sklearn_mds(args.n_components, args.n_clusters, seqfos, args.seed, plotdir=args.plotdir)
mds.run_bios2mds(args.n_components, args.n_clusters, seqfos, args.workdir, args.seed, plotdir=args.plotdir)
