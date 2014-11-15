import os
import csv

import utils
from opener import opener

# ----------------------------------------------------------------------------------------
def read_mute_info(indir, this_gene, approved_genes=None):
    if approved_genes == None:
        approved_genes = [this_gene,]
    observed_freqs = {}
    # add an observation for each position, for each gene where we observed that position
    for gene in approved_genes:
        mutefname = indir + '/mute-freqs/' + utils.sanitize_name(gene) + '.csv'
        if not os.path.exists(mutefname):
            continue
        with opener('r')(mutefname) as mutefile:
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
                observed_freqs[pos].append({'freq':freq, 'err':max(abs(freq-lo_err), abs(freq-hi_err))})

    # set final mute_freqs[pos] to the (inverse error-weighted) average over all the observations for each position
    mute_freqs = {}
    overall_total, overall_sum_of_weights = 0.0, 0.0  # also calculate the mean over all positions
    for pos in observed_freqs:
        total, sum_of_weights = 0.0, 0.0
        for obs in observed_freqs[pos]:
            assert obs['err'] > 0.0
            weight = 1.0 / obs['err']
            total += weight * obs['freq']
            sum_of_weights += weight
        assert sum_of_weights > 0.0
        mean_freq = total / sum_of_weights
        mute_freqs[pos] = mean_freq
        overall_total += total
        overall_sum_of_weights += sum_of_weights

    mute_freqs['overall_mean'] = overall_total / overall_sum_of_weights
    return mute_freqs
