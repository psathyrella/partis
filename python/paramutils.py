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

# ----------------------------------------------------------------------------------------
def make_mutefreq_plot(self, gene_name, states):
    nuke_colors = {'A':kRed+1, 'C':kBlue-7, 'G':kOrange-3, 'T':kGreen+2}

    ibin = 0
    drawn_name_texts, lines, vlines, texts = {}, {}, {}, {}
    for state in states:
        if len(state.emissions) == 0:
            assert state.name == 'init'
            continue

        # make label below bin
        drawn_name_texts[state.name] = TPaveText(-0.5 + ibin, -0.1, 0.5 + ibin, -0.05)
        drawn_name_texts[state.name].SetBorderSize(0)
        drawn_name_texts[state.name].SetFillColor(0)
        drawn_name_texts[state.name].SetFillStyle(0)
        drawn_name_texts[state.name].AddText(-0.5 + ibin, -0.075, simplify_state_name(state.name))

        total = 0.0
        lines[state.name], vlines[state.name], texts[state.name] = [], [], []
        for nuke, prob in sorted(state.emissions['probs'].items(), key=operator.itemgetter(1), reverse=True):
            # horizontal line at height total+prob
            lines[state.name].append(TLine(-0.5 + ibin, total + prob, 0.5 + ibin, total + prob))
            lines[state.name][-1].SetLineWidth(6)

            # vertical line from total to total+prob
            vlines[state.name].append(TLine(ibin, total, ibin, total + prob))
            vlines[state.name][-1].SetLineWidth(6)
            vlines[state.name][-1].SetLineColor(nuke_colors[nuke])

            # write [ACGT] at midpoint between total and total+prob
            midpoint = 0.5*(prob + 2*total)
            texts[state.name].append(TPaveText(-0.5 + ibin, midpoint-0.04, 0.5 + ibin, midpoint + 0.01))
            texts[state.name][-1].AddText(-0.5 + ibin, midpoint, nuke)
            texts[state.name][-1].SetBorderSize(0)
            texts[state.name][-1].SetFillColor(0)
            texts[state.name][-1].SetFillStyle(0)

            total += prob

        ibin += 1

    cvn = TCanvas('cvn', '', 1000, 300)
    n_bins = ibin
    hframe = TH1F(model.name + '-emission-frame', utils.unsanitize_name(model.name), n_bins, -0.5, n_bins - 0.5)
    hframe.SetNdivisions(202, 'y')
    hframe.SetNdivisions(0, 'x')
    hframe.Draw()

    for state_name in lines.keys():
        drawn_name_texts[state_name].Draw()
        for itrans in range(len(lines[state_name])):
            lines[state_name][itrans].Draw()
            vlines[state_name][itrans].Draw()
            # texts[state_name][itrans].Draw()

    cvn.SaveAs(self.base_plotdir + '/emissions/plots/' + gene_name + '.png')
    

