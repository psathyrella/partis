import os
import csv
import operator

import glutils
import utils

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
def read_mute_info(indir, this_gene, locus, approved_genes=None):  # NOTE this would probably be more accurate if we made some effort to align the genes before combining all the approved ones
    if this_gene == glutils.dummy_d_genes[locus]:
        return {'overall_mean' : 0.5}, {}

    if approved_genes is None:
        approved_genes = [this_gene, ]

    # add an observation for each position, for each gene where we observed that position NOTE this would be more sensible if they were aligned first
    observed_freqs, observed_counts = {}, {}
    total_counts = 0
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
                    observed_counts[pos] = {n : 0 for n in utils.nukes}

                observed_freqs[pos].append({'freq' : freq, 'err' : max(abs(freq-lo_err), abs(freq-hi_err))})  # append one for each gene
                for nuke in utils.nukes:
                    observed_counts[pos][nuke] += int(line[nuke + '_obs'])
                total_counts += int(line[nuke + '_obs'])

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

    observed_counts['total_counts'] = total_counts
    return mute_freqs, observed_counts

# ----------------------------------------------------------------------------------------
def make_mutefreq_plot(plotdir, gene_name, positions):
    import plotting
    """ NOTE shares a lot with make_transition_plot() in bin/plot-hmms.py. """
    nuke_colors = {'A' : 'red', 'C' : 'blue', 'G' : 'orange', 'T' : 'green'}
    fig, ax = plotting.mpl_init()
    fig.set_size_inches(plotting.plot_ratios[utils.get_region(gene_name)])

    ibin = 0
    print utils.color_gene(utils.unsanitize_name(gene_name))
    legend_colors = set()
    for info in positions:
        posname = info['name']

        # make label below bin
        ax.text(-0.5 + ibin, -0.075, simplify_state_name(posname), rotation='vertical', size=8)

        total = 0.0
        alpha = 0.6
        for nuke, prob in sorted(info['nuke_freqs'].items(), key=operator.itemgetter(1), reverse=True):
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
