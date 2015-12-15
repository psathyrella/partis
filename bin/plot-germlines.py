#!/usr/bin/env python
import sys
sys.path.insert(1, './python')
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# import seaborn as sns
# sns.set_style('ticks')
import numpy
from subprocess import check_call
from collections import OrderedDict

import utils
import plotting

# ----------------------------------------------------------------------------------------
plotdir = os.getenv('www') + '/tmp'
datadir = 'data/imgt'
glfo = {}
glfo['seqs'] = utils.read_germlines(datadir)
glfo['aligned-v-genes'] = utils.read_germlines(datadir, only_region='v', aligned=True)
vgenes = glfo['aligned-v-genes']['v'].keys()
pversions = OrderedDict()
for vg in vgenes:
    pv = utils.primary_version(vg)
    if pv not in pversions:
        pversions[pv] = []
    pversions[pv].append(vg)

# ----------------------------------------------------------------------------------------
def getmatrix(genelist):
    smatrix = [[] for _ in range(len(genelist))]
    for iv in range(len(genelist)):
        for jv in range(len(genelist)):
            if jv < iv + 1:
                smatrix[iv].append(0.)
                continue
            s1, s2 = [glfo['aligned-v-genes']['v'][genelist[index]] for index in [iv, jv]]
            # utils.color_mutants(s1, s2, print_result=True)
            fraction, length = utils.hamming_fraction(s1, s2, return_len_excluding_ambig=True, extra_bases='.')
            distance = int(fraction * length)
            # print '%3.0f / %3d = %5.3f   %s   %s' % (fraction * length, length, fraction, utils.color_gene(vgenes[iv]), utils.color_gene(vgenes[jv]))
            smatrix[iv].append(fraction)
    return smatrix

# ----------------------------------------------------------------------------------------
def plotheatmap(genelist, plotname, title=''):
    smatrix = getmatrix(genelist)
    assert len(smatrix) == len(smatrix[0])  # uh, I think I need this to be true
    fig, ax = plotting.mpl_init(fontsize=7)  #plt.subplots()
    plt.gcf().subplots_adjust(bottom=0.14, left=0.18, right=0.95, top=0.92)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    data = numpy.array(smatrix)
    cmap = plt.cm.Blues  #cm.get_cmap('jet')
    cmap.set_under('w')
    heatmap = ax.pcolor(data, cmap=cmap, vmin=0., vmax=numpy.amax(smatrix))
    cbar = plt.colorbar(heatmap)

    xticklabels = [utils.summarize_gene_name(g) for g in genelist]
    ticks = [n - 0.5 for n in range(1, len(xticklabels) + 1, 1)]
    # xticklabels = [str(int(n + 0.5)) for n in ticks]
    yticklabels = xticklabels
    # if n_biggest_clusters > 20:
    #     modulo = 3
    #     ticks = [ticks[it] for it in range(0, len(ticks), modulo)]
    #     xticklabels = [b_cluster_lengths[it] for it in range(0, len(b_cluster_lengths), modulo)]
    #     yticklabels = [a_cluster_lengths[it] for it in range(0, len(a_cluster_lengths), modulo)]
    plt.xticks(ticks, xticklabels, rotation=90)
    plt.yticks(ticks, yticklabels)
    # plt.xlabel(legends.get(meth2, meth2) + ' cluster size')  # I don't know why it's reversed, it just is
    # plt.ylabel(legends.get(meth1, meth1) + ' cluster size')
    # ax.set_xlim(0, n_biggest_clusters)
    # ax.set_ylim(0, n_biggest_clusters)

    plt.title(title)

    if not os.path.exists(plotdir + '/plots'):
        os.makedirs(plotdir + '/plots')
    plt.savefig(plotdir + '/plots/' + plotname + '.svg')
    plt.close()

# ----------------------------------------------------------------------------------------
vgenes = vgenes[ : 5]
for pv in pversions:
    print pv
    plotheatmap(pversions[pv], utils.sanitize_name(pv), pv)
    # sys.exit()

check_call(['./bin/makeHtml', plotdir, '2', 'foop', 'svg'])
check_call(['./bin/permissify-www', plotdir])
