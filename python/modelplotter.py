#!/usr/bin/env python

import os
import argparse
from collections import OrderedDict
import glob
import operator
import yaml
import sys
from subprocess import check_call

import plotting
import paramutils
import utils
import modelplotter

# ----------------------------------------------------------------------------------------
def find_state_number(name):
    assert name.find('IGH') == 0
    state_number = int(name[ name.rfind('_') + 1 : ])
    assert state_number >= 0 and state_number <= 400
    return state_number

# ----------------------------------------------------------------------------------------
class ModelPlotter(object):
    def __init__(self, args, base_plotdir, skip_boring_states=''):
        raise Exception('needs to be converted off root')
        self.base_plotdir = base_plotdir
        self.skip_boring_states = skip_boring_states
        plot_types = ('transitions', 'emissions')
        for ptype in plot_types:
            plotdir = self.base_plotdir + '/' + ptype + '/plots'
            utils.prep_dir(plotdir, '*.png')

        if args.hmmdir != None:
            filelist = glob.glob(args.hmmdir + '/*.yaml')
        else:
            filelist = utils.get_arg_list(args.infiles)
        if len(filelist) == 0:
            print 'ERROR zero files passed to modelplotter'
            sys.exit()
        for infname in filelist:
            gene_name = os.path.basename(infname).replace('.yaml', '')  # the sanitized name, actually
# # ----------------------------------------------------------------------------------------
#             if utils.get_region(gene_name) == 'v' and 'IGHV4-39_star_' not in gene_name:
#                 continue
# # ----------------------------------------------------------------------------------------
            with open(infname) as infile:
                model = yaml.load(infile)
                self.make_transition_plot(gene_name, model)
                self.make_emission_plot(gene_name, model)

        for ptype in plot_types:
            check_call(['./bin/makeHtml', self.base_plotdir + '/' + ptype, '1', 'null', 'png'])
        check_call(['./bin/permissify-www', self.base_plotdir])

            # break

    # ----------------------------------------------------------------------------------------
    def make_transition_plot(self, gene_name, model):
        ibin = 0
        drawn_name_texts, lines, texts = {}, {}, {}
        for state in model.states:
            if utils.get_region(gene_name) in self.skip_boring_states:
                if state.name != 'init' and len(state.transitions) == 1:  # skip uninteresting states
                    to_state = state.transitions.keys()[0]  # skip states with only transitions to end
                    if to_state == 'end':
                        continue
                    if find_state_number(state.name) + 1 == find_state_number(to_state):  # skip states with only transitions to next state
                        continue

            drawn_name_texts[state.name] = TPaveText(-0.5 + ibin, -0.1, 0.5 + ibin, -0.05)
            drawn_name_texts[state.name].SetBorderSize(0)
            drawn_name_texts[state.name].SetFillColor(0)
            drawn_name_texts[state.name].SetFillStyle(0)
            drawn_name_texts[state.name].AddText(-0.5 + ibin, -0.075, paramutils.simplify_state_name(state.name))

            sorted_to_states = {}
            for name in state.transitions.keys():
                if name.find('IGH') == 0:
                    sorted_to_states[name] = int(paramutils.simplify_state_name(name))
                else:
                    sorted_to_states[name] = name
            sorted_to_states = sorted(sorted_to_states.items(), key=operator.itemgetter(1))

            total = 0.0
            lines[state.name], texts[state.name] = [], []
            for to_state, simple_to_state in sorted_to_states:

                prob = state.transitions[to_state]
                lines[state.name].append(TLine(-0.5 + ibin, total + prob, 0.5 + ibin, total + prob))
                lines[state.name][-1].SetLineColor(kGreen+2)
                lines[state.name][-1].SetLineWidth(6)

                midpoint = 0.5*(prob + 2*total)
                texts[state.name].append(TPaveText(-0.5 + ibin, midpoint-0.04, 0.5 + ibin, midpoint + 0.01))
                texts[state.name][-1].AddText(-0.5 + ibin, midpoint, paramutils.simplify_state_name(to_state))
                texts[state.name][-1].SetBorderSize(0)
                texts[state.name][-1].SetFillColor(0)
                texts[state.name][-1].SetFillStyle(0)

                total += prob
    
            ibin += 1

        cvn = TCanvas('mod-cvn', '', 1000, 400)
        n_bins = ibin
        hframe = TH1D(model.name + '-transition-frame', utils.unsanitize_name(model.name), n_bins, -0.5, n_bins - 0.5)
        if utils.get_region(gene_name) in self.skip_boring_states:
            hframe.SetTitle(hframe.GetTitle() + ' (skipped boring states)')
        hframe.SetNdivisions(202, 'y')
        hframe.SetNdivisions(0, 'x')
        hframe.Draw()

        for state_name in lines.keys():
            drawn_name_texts[state_name].Draw()
            for itrans in range(len(lines[state_name])):
                lines[state_name][itrans].Draw()
                texts[state_name][itrans].Draw()

        cvn.SaveAs(self.base_plotdir + '/transitions/plots/' + gene_name + '.png')

    # ----------------------------------------------------------------------------------------
    def make_emission_plot(self, gene_name, model):
        plotting_info = []
        for state in model.states:
            if state.emissions is None:
                assert state.name == 'init'
                continue

            plotting_info.append({})
            plotting_info[-1]['name'] = state.name
            plotting_info[-1]['nuke_freqs'] = state.emissions['probs']

        paramutils.make_mutefreq_plot(self.base_plotdir + '/emissions', gene_name, plotting_info)

    # # ----------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument('--hmmdir')
parser.add_argument('--infiles')
args = parser.parse_args()
assert args.hmmdir is not None or args.infiles is not None

if __name__ == '__main__':
    # hmmdir = os.getenv('HOME') + '/work/partis/caches/' + args.label + '/' + args.flavor + '_parameters/hmms'
    assert os.path.exists(os.getenv('www'))
    mplot = ModelPlotter(args, os.getenv('www') + '/modelplots/', skip_boring_states='')  #'v')
