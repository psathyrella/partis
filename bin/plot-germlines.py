#!/usr/bin/env python
import sys
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# import seaborn as sns
# sns.set_style('ticks')
import numpy
from subprocess import check_call
import itertools
from collections import OrderedDict
current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
if not os.path.exists(current_script_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % current_script_dir
sys.path.insert(1, current_script_dir)

import utils
import plotting

# ----------------------------------------------------------------------------------------
datadir = 'data/imgt'
xtitles = {
    'indels' : 'fraction of positions indel\'d',
    'subs' : 'substitution fraction'
}
glfo = utils.read_glfo(datadir)
vgenes = glfo['aligned-genes']['v'].keys()
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
def get_gene_pair_matrix(genelist, difftype):
    """ return matrix comparing all pairs of genes in <genelist> """
    smatrix = [[] for _ in range(len(genelist))]
    for iv in range(len(genelist)):
        for jv in range(len(genelist)):
            if jv < iv + 1:
                smatrix[iv].append(0.)
                continue
            s1, s2 = [glfo['aligned-genes']['v'][genelist[index]] for index in [iv, jv]]
            # utils.color_mutants(s1, s2, print_result=True)
            if difftype == 'hamming':
                fraction, length = utils.hamming_fraction(s1, s2, return_len_excluding_ambig=True, extra_bases='.')
            elif difftype == 'indels':
                fraction = indel_difference_fraction(s1, s2)
            elif difftype == 'subs':
                fraction = substitution_difference_fraction(s1, s2)
            else:
                raise Exception('unexpected difftype %s' % difftype)
            smatrix[iv].append(fraction)
    return smatrix

# ----------------------------------------------------------------------------------------
def get_gene_set_mean_matrix(genesets, difftype):
    """ return matrix comparing the sets of genes in <genenames>, i.e. each entry is the average over all pairs of sequences in set 1 and set 2. """
    setnames, genenames = genesets.keys(), genesets.values()
    n_sets = len(genenames)
    smatrix = [[] for _ in range(n_sets)]
    for iv in range(n_sets):
        for jv in range(n_sets):
            # if setnames[iv] != '3/OR15' or setnames[jv] != '4/OR15':
            #     smatrix[iv].append(0.)
            #     continue
            if jv < iv + 1:
                smatrix[iv].append(0.)
                continue
            # print '  %s %s' % (setnames[iv], setnames[jv])
            seqs1 = [glfo['aligned-genes']['v'][g] for g in genenames[iv]]
            seqs2 = [glfo['aligned-genes']['v'][g] for g in genenames[jv]]

            total, nfractions = 0., 0
            for is1 in range(len(seqs1)):
                # print '   ', utils.color_gene(genenames[iv][is1])
                for is2 in range(len(seqs2)):
                    # print '     ', utils.color_gene(genenames[jv][is2]),
                    s1 = seqs1[is1]
                    s2 = seqs2[is2]
                    # utils.color_mutants(s1, s2, print_result=True, extra_str='    ')
                    if difftype == 'hamming':
                        fraction, length = utils.hamming_fraction(s1, s2, return_len_excluding_ambig=True, extra_bases='.')
                    elif difftype == 'indels':
                        fraction = indel_difference_fraction(s1, s2)
                    elif difftype == 'subs':
                        fraction = substitution_difference_fraction(s1, s2)
                    else:
                        raise Exception('unexpected difftype %s' % difftype)
                    # print '      %.3f' % fraction
                    total += fraction
                    nfractions += 1
            meanfraction = 0. if nfractions == 0 else float(total) / nfractions
            # print '   mean %.3f' % meanfraction
            smatrix[iv].append(meanfraction)

    return smatrix

# ----------------------------------------------------------------------------------------
def plotheatmap(plotdir, plotname, difftype, genelist=None, genesets=None, title='', xtitle=''):
    assert genelist is None or genesets is None
    if genelist is not None:
        smatrix = get_gene_pair_matrix(genelist, difftype)
        xticklabels = [utils.summarize_gene_name(g) for g in genelist]
    elif genesets is not None:
        smatrix = get_gene_set_mean_matrix(genesets, difftype)
        xticklabels = [sn for sn in genesets.keys()]
    else:
        raise Exception('no gene list specified')
    assert len(smatrix) == len(smatrix[0])  # uh, I think I need this to be true
    fig, ax = plotting.mpl_init()
    plt.tick_params(axis='both', which='major', labelsize=7)
    plt.gcf().subplots_adjust(bottom=0.14, left=0.18, right=0.95, top=0.92)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    data = numpy.array(smatrix)
    cmap = plt.cm.Blues  #cm.get_cmap('jet')
    cmap.set_under('w')
    heatmap = ax.pcolor(data, cmap=cmap, vmin=0., vmax=0.5)  #vmax=numpy.amax(smatrix))
    cbar = plt.colorbar(heatmap, shrink=0.9, pad=0.01)
    cbar.set_label(xtitle, rotation=270, labelpad=30)
    cbar.ax.tick_params(labelsize=10) 

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
for difftype in ['indels', 'subs']:
    print difftype
    # individual primary version plots
    for pv in pversions:
        print '   ', pv
        plotheatmap(baseplotdir + '/' + difftype, utils.sanitize_name(pv), difftype, genelist=pversions[pv], title='primary version \"' + pv + '\"', xtitle=xtitles[difftype])

    # plots comparing two different primary versions
    plotheatmap(baseplotdir + '/' + difftype,
                'compare-pvs',
                difftype,
                genesets=pversions, title='compare means over pairs of primary versions', xtitle=xtitles[difftype])

    check_call(['./bin/makeHtml', baseplotdir + '/' + difftype, '2', 'foop', 'svg'])

check_call(['./bin/permissify-www', baseplotdir])
