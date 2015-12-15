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

# remove primary versions that only have one gene
for pv in pversions:
    if len(pversions[pv]) == 1:
        print 'removing single-gene pv %s' % pv
        del pversions[pv]

# ----------------------------------------------------------------------------------------
def indel_difference_fraction(seq1, seq2):
    """ fraction of positions in the aligned sequences <seq1> <seq2> which are dots in exactly one of them """
    if len(seq1) != len(seq2):
        raise Exception('sequences different length:\n  %s\n  %s' % (seq1, seq2))

    n_differences = 0
    for ich in range(len(seq1)):
        ch1 = seq1[ich]
        ch2 = seq2[ich]
        if ch1 == '.' and ch2 != '.' or  ch2 == '.' and ch1 != '.':
            n_differences += 1

    return float(n_differences) / len(seq1)

# ----------------------------------------------------------------------------------------
def substitution_difference_fraction(seq1, seq2):
    """ fraction of positions in the aligned sequences <seq1> <seq2> (excluding positions which are dots in either or both sequences) which are different bases """
    if len(seq1) != len(seq2):
        raise Exception('sequences different length:\n  %s\n  %s' % (seq1, seq2))

    n_substitutions, n_total = 0, 0
    for ich in range(len(seq1)):
        ch1 = seq1[ich]
        ch2 = seq2[ich]
        if ch1 == '.' or ch2 == '.':
            continue

        n_total += 1
        if ch1 != ch2:
            n_substitutions += 1

    return float(n_substitutions) / n_total

# s1 = 'a.cdx.'
# s2 = 'a.cd..'
# utils.color_mutants(s1, s2, print_result=True)
# print substitution_difference_fraction(s1, s2)
# sys.exit()
# ----------------------------------------------------------------------------------------
def getmatrix(genelist, difftype):
    smatrix = [[] for _ in range(len(genelist))]
    for iv in range(len(genelist)):
        for jv in range(len(genelist)):
            if jv < iv + 1:
                smatrix[iv].append(0.)
                continue
            s1, s2 = [glfo['aligned-v-genes']['v'][genelist[index]] for index in [iv, jv]]
            # utils.color_mutants(s1, s2, print_result=True)
            if difftype == 'hamming':
                fraction, length = utils.hamming_fraction(s1, s2, return_len_excluding_ambig=True, extra_bases='.')
                # print '%3.0f / %3d = %5.3f   %s   %s' % (fraction * length, length, fraction, utils.color_gene(vgenes[iv]), utils.color_gene(vgenes[jv]))
            elif difftype == 'indels':
                fraction = indel_difference_fraction(s1, s2)
            elif difftype == 'subs':
                fraction = substitution_difference_fraction(s1, s2)
            else:
                raise Exception('unexpected difftype %s' % difftype)
            smatrix[iv].append(fraction)
    return smatrix

# ----------------------------------------------------------------------------------------
def plotheatmap(genelist, plotdir, plotname, difftype, title=''):
    smatrix = getmatrix(genelist, difftype)
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
baseplotdir = os.getenv('www') + '/tmp'
vgenes = vgenes[ : 5]
for difftype in ['indels', 'subs']:
    # # individual primary version plots
    # for pv in pversions:
    #     print pv
    #     plotheatmap(pversions[pv], baseplotdir + '/' + difftype, utils.sanitize_name(pv), difftype, title='primary version \"' + pv + '\"')
    # check_call(['./bin/makeHtml', baseplotdir + '/' + difftype, '2', 'foop', 'svg'])

    # plots comparing two different primary versions
    pvlist = pversions.keys()
    for ipv in range(len(pvlist)):
        for jpv in range(ipv + 1, len(pvlist)):
            print pvlist[ipv], pvlist[jpv]
        # plotheatmap(pversions[pv], baseplotdir + '/' + difftype, utils.sanitize_name(pv), difftype, title='primary version \"' + pv + '\"')


check_call(['./bin/permissify-www', baseplotdir])
