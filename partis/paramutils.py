from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import csv
import operator
import sys

from . import glutils
from . import utils
from io import open

# ----------------------------------------------------------------------------------------
def simplify_state_name(state_name):
    if state_name.find('IG') == 0 or state_name.find('TR') == 0:
        return state_name[state_name.rfind('_') + 1 : ]
    elif state_name == 'insert_left':
        return 'i_l'
    elif state_name == 'insert_right':
        return 'i_r'
    else:
        return state_name

# ----------------------------------------------------------------------------------------
def read_mute_counts(indir, gene, locus, extra_genes=None, debug=False):  # NOTE I'm adding the <extra_genes> arg in a hackish way because i need this to not crash in one specific instance (running bin/test-germline-inference.py) where the file for <gene> doesn't exist, but I don't remember/understand how this fcn and the following function work well enough to do this more sensibly
    # NOTE also that this new hack that allows a different gene's counts to be used might break something later on if the genes have different lengths? I have no idea
    # ----------------------------------------------------------------------------------------
    def read_single_file(gtmp):
        mfname = indir + '/mute-freqs/' + utils.sanitize_name(gtmp) + '.csv'
        if not os.path.exists(mfname):
            return None
        observed_counts = {}
        with open(mfname, 'r') as mutefile:
            reader = csv.DictReader(mutefile)
            for line in reader:
                pos = int(line['position'])
                assert pos not in observed_counts
                observed_counts[pos] = {n : int(line[n + '_obs']) for n in utils.nukes}
        if debug:
            print('    read %d per-base mute counts from %s' % (len(observed_counts), mfname))
        return observed_counts

    # ----------------------------------------------------------------------------------------
    if extra_genes is not None:  # I don't want to fix it cause it'd be kinda hard, and also I don't think it ever happens under normal circumstances -- it's only called with this arg from simulation, in which case you should always have parameters for the gene you're asking for
        print('%s Reading per-base mutation counts for genes (%s) in addition to the desired one (%s), which doesn\'t really make sense, since the counts will be wrong at the positions at which the genes differ.' % (utils.color('red', 'warning'), utils.color_genes(extra_genes), utils.color_gene(gene)))
        print('   This should only happen if you\'re doing something weird, probably running simulation asking for genes for which you don\'t have parameters.')
        print('   If this is the case and you only care that it doesn\'t crash, and not that the mutation model is particularly accurate, this is fine.')
    if gene == glutils.dummy_d_genes[locus]:
        return {}

    if extra_genes is None:
        approved_genes = [gene]
    else:
        assert gene not in extra_genes
        approved_genes = [gene] + extra_genes

    for gtmp in approved_genes:
        observed_counts = read_single_file(gtmp)
        if observed_counts is not None:  # HACK this just uses the first one that's there (in the vast majority of cases it'll just be <gene> -- i think the only way it can be missing is if you hard code a specific gene (e.g. in bin/test-germline-inference.py) and it isn't in the parameter directory you passed
            break
    return observed_counts  # raw per-{ACGT} counts for each position, summed over genes ("raw" as in not a weighted average over a bunch of genes as in read_mute_freqs_with_weights())

# ----------------------------------------------------------------------------------------
def read_mute_freqs_with_weights(indir, approved_genes, debug=False):  # it would be nice to eventually align the genes before combining
    # returns:
    #  - mute_freqs: inverse error-weighted average mute freq over all genes for each position
    #     - also includes weighted and unweigthed means over positions

    if len(approved_genes) == 0:
        raise Exception('no approved genes')

    if approved_genes[0] == glutils.dummy_d_genes[utils.get_locus(approved_genes[0])]:
        return {'overall_mean' : 0.5, 'unweighted_overall_mean' : 0.5}

    if debug:
        print('    reading mute freqs from %s for %d gene%s: %s' % (indir, len(approved_genes), utils.plural(len(approved_genes)), utils.color_genes(approved_genes)))

    # add an observation for each position, for each gene where we observed that position NOTE this would be more sensible if they were aligned first
    observed_freqs = {}
    for gene in approved_genes:
        mutefname = indir + '/mute-freqs/' + utils.sanitize_name(gene) + '.csv'
        if not os.path.exists(mutefname):
            continue
        with open(mutefname, 'r') as mutefile:
            reader = csv.DictReader(mutefile)
            for line in reader:
                pos = int(line['position'])
                freq = float(line['mute_freq'])
                lo_err = float(line['lo_err'])  # NOTE lo_err in the file is really the lower *bound*
                hi_err = float(line['hi_err'])  #   same deal
                assert freq >= 0.0 and lo_err >= 0.0 and hi_err >= 0.0  # you just can't be too careful

                if freq < utils.eps or abs(1.0 - freq) < utils.eps:  # if <freq> too close to 0 or 1, replace it with the midpoint of its uncertainty band
                    freq = 0.5 * (lo_err + hi_err)

                if pos not in observed_freqs:
                    observed_freqs[pos] = []

                observed_freqs[pos].append({'freq' : freq, 'err' : max(abs(freq-lo_err), abs(freq-hi_err))})  # append one for each gene

    # set final mute_freqs[pos] to the (inverse error-weighted) average over all the observations [i.e. genes] for each position
    mute_freqs = {}
    for pos in observed_freqs:
        total, sum_of_weights = 0.0, 0.0
        for obs in observed_freqs[pos]:  # loop over genes
            assert obs['err'] > 0.0
            weight = 1.0 / obs['err']
            total += weight * obs['freq']
            sum_of_weights += weight
        assert sum_of_weights > 0.0
        mean_freq = total / sum_of_weights
        mute_freqs[pos] = mean_freq

    # NOTE I'm sure that this weighting scheme makes sense for comparing differeing genes at the same position, but I'm less sure it makes sense for the overall mean. But, I don't want to track down all the places that changing it might affect right now
    mute_freqs['overall_mean'] = 0.
    weighted_denom = sum([1. / obs['err'] for pos in observed_freqs for obs in observed_freqs[pos]])
    if weighted_denom > 0.:
        mute_freqs['overall_mean'] = sum([obs['freq'] / obs['err'] for pos in observed_freqs for obs in observed_freqs[pos]]) / weighted_denom

    # I need the inverse-error-weighted numbers to sensibly combine genes, but then I also need unweigthed values that I can easily write to the yaml files for other people to use
    mute_freqs['unweighted_overall_mean'] = 0.
    unweighted_denom = sum([len(observed_freqs[pos]) for pos in observed_freqs])
    if unweighted_denom > 0.:
        mute_freqs['unweighted_overall_mean'] = sum([obs['freq'] for pos in observed_freqs for obs in observed_freqs[pos]]) / unweighted_denom

    if debug:
        iskipstart = 35  # i.e. for v genes skip the middle positions
        positions = sorted(observed_freqs)
        if len(positions) > 2 * iskipstart:
            print('      %s%s%s' % (' '.join([('%4d' % p) for p in positions[:iskipstart]]), utils.color('blue', ' [...] '), ' '.join([('%4d' % p) for p in positions[len(positions) - iskipstart :]])))
            print('      %s%s%s' % (' '.join([('%4.2f' % mute_freqs[p]) for p in positions[:iskipstart]]), utils.color('blue', ' [...] '), ' '.join([('%4.2f' % mute_freqs[p]) for p in positions[len(positions) - iskipstart :]])))
        else:
            print('      %s' % ' '.join([('%4d' % p) for p in positions]))
            print('      %s' % ' '.join([('%4.2f' % mute_freqs[p]) for p in positions]))
        print('        overall mean: %5.3f (unweighted %5.3f)' % (mute_freqs['overall_mean'], mute_freqs['unweighted_overall_mean']))

    return mute_freqs

# ----------------------------------------------------------------------------------------
def make_mutefreq_plot(plotdir, gene_name, positions, debug=False):
    from . import plotting
    """ NOTE shares a lot with make_transition_plot() in bin/plot-hmms.py. """
    nuke_colors = {'A' : 'red', 'C' : 'blue', 'G' : 'orange', 'T' : 'green'}
    fig, ax = plotting.mpl_init()
    fig.set_size_inches(plotting.plot_ratios[utils.get_region(gene_name)])

    ibin = 0
    if debug:
        print('  %s' % utils.color_gene(utils.unsanitize_name(gene_name)))
    legend_colors = set()
    for info in positions:
        posname = info['name']

        # make label below bin for position and germline nuke
        ax.text(-0.5 + ibin, -0.075, simplify_state_name(posname), rotation='vertical', size=8)
        ax.text(-0.5 + ibin, -0.15, info.get('gl_nuke', '?'), fontsize=10, fontweight='bold')
        sorted_nukes, _ = zip(*sorted(list(info['nuke_freqs'].items()), key=operator.itemgetter(1), reverse=True))
        if 'gl_nuke' in info and info['gl_nuke'] in info['nuke_freqs']:  # put the germline nuke first if we have it (second clause is for states with germline N))
            sorted_nukes = [info['gl_nuke']] + [n for n in sorted_nukes if n != info['gl_nuke']]

        total = 0.0
        alpha = 0.6
        for nuke in sorted_nukes:
            prob = info['nuke_freqs'][nuke]
            color = nuke_colors[nuke]

            label_to_use = None
            if color not in legend_colors:
                label_to_use = nuke
                legend_colors.add(color)

            # horizontal line at height total+prob
            ax.plot([-0.5 + ibin, 0.5 + ibin], [total + prob, total + prob], color=color, alpha=alpha, linewidth=3, label=label_to_use)

            # vertical line from total to total + prob
            ax.plot([ibin, ibin], [total + 0.01, total + prob], color=color, alpha=alpha, linewidth=3)

            # # write [ACGT] at midpoint between total and total+prob
            # midpoint = 0.5*(prob + 2*total)
            # ... *redacted*

            total += prob

        ibin += 1

    ax.get_xaxis().set_visible(False)
    plotting.mpl_finish(ax, plotdir, gene_name, ybounds=(-0.01, 1.01), xbounds=(-3, len(positions) + 3), leg_loc=(0.95, 0.1), adjust={'left' : 0.1, 'right' : 0.8}, leg_prop={'size' : 8})
