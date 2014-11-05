#!/usr/bin/env python

import os
import glob
import operator
import yaml
import sys
from subprocess import check_call
sys.argv.append('-b')  # root just loves its stupid little splashes
from ROOT import TH1F, TCanvas, kRed, gROOT, TLine, TLegend, kBlue, kGreen, TPaveText, TStyle

import plotting
import utils

# ----------------------------------------------------------------------------------------
def simplify_state_name(state_name):
    if state_name.find('IGH') == 0:
        return state_name[state_name.rfind('_') + 1 : ]
    elif state_name == 'insert_left':
        return 'i_l'
    elif state_name == 'insert_right':
        return 'i_r'
    else:
        return state_name

# ----------------------------------------------------------------------------------------
def find_state_number(name):
    assert name.find('IGH') == 0
    state_number = int(name[ name.rfind('_') + 1 : ])
    assert state_number >= 0 and state_number <= 400
    return state_number

# ----------------------------------------------------------------------------------------
class ModelPlotter(object):
    def __init__(self, modeldir, base_plotdir):
        self.base_plotdir = base_plotdir
        self.cvn = TCanvas('cvn', '', 4000, 1000)
        plot_types = ('transitions', 'emissions', 'pair-emissions')
        for ptype in plot_types:
            plotdir = self.base_plotdir + '/' + ptype + '/plots'
            utils.prep_dir(plotdir, '*.png')

        filelist = glob.glob(modeldir + '/IGHJ*.yaml')
        for infname in filelist:
            gene_name = os.path.basename(infname).replace('.yaml', '')  # the sanitized name, actually
            with open(infname) as infile:
                model = yaml.load(infile)
                self.make_transition_plot(gene_name, model)
                # for state in model.states:
                #     self.plot_transitions(gene_name, state)
                #     # self.plot_emissions(gene_name, state)

            for ptype in plot_types:
                check_call(['makeHtml', self.base_plotdir + '/' + ptype, '1', 'null', 'png'])
            check_call(['./permissify-www', self.base_plotdir])

            # break

    # ----------------------------------------------------------------------------------------
    def make_transition_plot(self, gene_name, model):
        ibin = 0
        drawn_name_texts, lines, texts = {}, {}, {}
        for state in model.states:
            if len(state.transitions) == 1:
                to_state = state.transitions.keys()[0]
                if to_state == 'end':
                    continue
                if find_state_number(state.name) + 1 == find_state_number(to_state):
                    continue

            drawn_name_texts[state.name] = TPaveText(-0.5 + ibin, -0.1, 0.5 + ibin, -0.05)
            drawn_name_texts[state.name].SetBorderSize(0)
            drawn_name_texts[state.name].SetFillColor(0)
            drawn_name_texts[state.name].SetFillStyle(0)
            drawn_name_texts[state.name].AddText(-0.5 + ibin, -0.075, simplify_state_name(state.name))

            sorted_to_states = {}
            for name in state.transitions.keys():
                if name.find('IGH') == 0:
                    sorted_to_states[name] = int(simplify_state_name(name))
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
                texts[state.name].append(TPaveText(-0.5 + ibin, midpoint, 0.5 + ibin, midpoint + 0.05))
                texts[state.name][-1].AddText(-0.5 + ibin, midpoint, simplify_state_name(to_state))
                texts[state.name][-1].SetBorderSize(0)
                texts[state.name][-1].SetFillColor(0)
                texts[state.name][-1].SetFillStyle(0)

                total += prob
    
            ibin += 1

        n_bins = ibin
        hframe = TH1F(model.name + '-frame', '', n_bins, -0.5, n_bins - 0.5)
        hframe.SetNdivisions(202, 'y')
        hframe.SetNdivisions(0, 'x')
        hframe.Draw()

        for state_name in lines.keys():
            drawn_name_texts[state_name].Draw()
            for itrans in range(len(lines[state_name])):
                lines[state_name][itrans].Draw()
                texts[state_name][itrans].Draw()

        self.cvn.SaveAs(self.base_plotdir + '/transitions/plots/' + gene_name + '.png')

    # # ----------------------------------------------------------------------------------------
    # def plot_emissions(self, gene_name, state):
        

# plot transitions

# plot emissions

# plot pair emissions

if __name__ == '__main__':
    mplot = ModelPlotter(os.getenv('HOME') + '/work/partis/caches/fixed/hmm_parameters/hmms', os.getenv('www') + '/modelplots')
    # /IGHD1-14_star_01.yaml
