#!/usr/bin/env python

import os
import glob
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
        self.cvn = TCanvas('cvn', '', 10000, 2000)
        plot_types = ('transitions', 'emissions', 'pair-emissions')

        filelist = glob.glob(modeldir + '/IGHJ*.yaml')
        for infname in filelist:
            gene_name = os.path.basename(infname).replace('.yaml', '')  # the sanitized name, actually
            for ptype in plot_types:
                plotdir = self.base_plotdir + '/' + gene_name + '/' + ptype + '/plots'
                utils.prep_dir(plotdir, '*.svg')
            with open(infname) as infile:
                model = yaml.load(infile)
                self.make_transition_hist(model)
                sys.exit()
                for state in model.states:
                    self.plot_transitions(gene_name, state)
                    # self.plot_emissions(gene_name, state)

            for ptype in plot_types:
                check_call(['makeHtml', self.base_plotdir + '/' + gene_name + '/' + ptype, '3', 'null', 'svg'])
            check_call(['./permissify-www', self.base_plotdir])

            break

    # ----------------------------------------------------------------------------------------
    def plot_transitions(self, gene_name, state):
        hist = plotting.make_hist(state.transitions, 'string', state.name + '-transitions')
        if hist.GetNbinsX() == 1:
            if hist.GetXaxis().GetBinLabel(1) == 'end':
                # print '%s only has transitions to end' % state.name
                return
            if state.name.find('IGH') != 0 or hist.GetXaxis().GetBinLabel(1).find('IGH') != 0:
                print 'oops: %s %s' % (state.name, hist.GetXaxis().GetBinLabel(1))
                sys.exit()
            assert find_state_number(state.name) + 1 == find_state_number(hist.GetXaxis().GetBinLabel(1))
            # print 'only one transtition for %s' % state.name
            return

        print 'plotting %s' % state.name
        for ibin in range(1, hist.GetNbinsX()+1):
            hist.GetXaxis().SetBinLabel(ibin, simplify_state_name(hist.GetXaxis().GetBinLabel(ibin)) + '\nFOO' )

        hist.SetTitle(simplify_state_name(state.name))
        hist.Draw('hist')
        self.cvn.SaveAs(self.base_plotdir + '/' + gene_name + '/transitions/plots/' + state.name + '.svg')
        
    # ----------------------------------------------------------------------------------------
    def make_transition_hist(self, model):
        hframe = TH1F(model.name + '-frame', '', len(model.states) + 2, -1.5, len(model.states) + 1.5)
        hframe.SetNdivisions(202, 'y')
        hframe.Draw()
        ibin = 0
        for state in model.states:
            # if state.name != 'init':
            #     continue
            lines, texts = [], []
            total = 0.0
            for to_state, prob in state.transitions.iteritems():
                lines.append(TLine(-0.5 + ibin, total + prob, 0.5 + ibin, total + prob))
                # lines.append(TLine(-0.5 + ibin, prob + total, 100, prob + total))
                lines[-1].SetLineColor(kGreen+2)
                lines[-1].SetLineWidth(6)
                lines[-1].Draw()
                midpoint = 0.5*(prob + 2*total)
                texts.append(TPaveText(-0.5 + ibin, midpoint, 0.5 + ibin, midpoint + 0.05))
                texts[-1].AddText(-0.5 + ibin, midpoint, simplify_state_name(to_state))
                texts[-1].SetBorderSize(0)
                texts[-1].SetFillColor(0)
                texts[-1].SetFillStyle(0)
                texts[-1].Draw()
                total += prob
    
            ibin += 1
        self.cvn.SaveAs(os.getenv('www') + '/foo.png')

    # # ----------------------------------------------------------------------------------------
    # def plot_emissions(self, gene_name, state):
        

# plot transitions

# plot emissions

# plot pair emissions

if __name__ == '__main__':
    mplot = ModelPlotter(os.getenv('HOME') + '/work/partis/caches/fixed/hmm_parameters/hmms', os.getenv('www') + '/modelplots')
    # /IGHD1-14_star_01.yaml
