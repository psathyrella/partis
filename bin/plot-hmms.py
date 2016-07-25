#!/usr/bin/env python
import os
import argparse
from collections import OrderedDict
import glob
import operator
import yaml
import sys
from subprocess import check_call
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')
import plotting
import paramutils
import utils

# ----------------------------------------------------------------------------------------
def simplify_state_name(state_name):
    if state_name.find('IG') == 0:
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
    def __init__(self, args, base_plotdir):
        self.base_plotdir = base_plotdir

        self.eps_to_skip = 1e-3
        print 'skipping eps %f' % self.eps_to_skip

        plot_types = ('transitions', 'emissions')
        for ptype in plot_types:
            plotdir = self.base_plotdir + '/' + ptype
            utils.prep_dir(plotdir, wildlings=['*.png', '*.svg'])

        if args.hmmdir != None:
            filelist = glob.glob(args.hmmdir + '/*.yaml')
        else:
            filelist = utils.get_arg_list(args.infiles)
        if len(filelist) == 0:
            raise Exception('zero files passed to modelplotter')
        for infname in filelist:
            gene_name = os.path.basename(infname).replace('.yaml', '')  # the sanitized name, actually
            with open(infname) as infile:
                model = yaml.load(infile)
                self.make_transition_plot(gene_name, model)
                # self.make_emission_plot(gene_name, model)

    # ----------------------------------------------------------------------------------------
    def make_transition_plot(self, gene_name, model):
        fig, ax = plotting.mpl_init()
        fig.set_size_inches(plotting.plot_ratios[utils.get_region(gene_name)])

        ibin = 0
        drawn_name_texts, lines, texts = {}, {}, {}
        print gene_name
        legend_colors = set()  # add a color to this the first time you plot it
        for state in model.states:

            ax.text(-0.5 + ibin, -0.075, simplify_state_name(state.name), rotation='vertical', size=8)

            sorted_to_states = {}
            for name in state.transitions.keys():
                if name.find('IGH') == 0:
                    sorted_to_states[name] = int(simplify_state_name(name))
                else:
                    sorted_to_states[name] = name
            sorted_to_states = sorted(sorted_to_states.items(), key=operator.itemgetter(1))

            total = 0.0
            for to_state, simple_to_state in sorted_to_states:

                prob = state.transitions[to_state]

                alpha = 0.6
                style = '-'

                if 'insert' in str(simple_to_state):
                    label = 'insert'
                    color = '#3498db'  # blue
                    width = 5
                elif str(simple_to_state) == 'end':
                    label = 'end'
                    color = 'red'
                    width = 3
                else:  # regional/internal states
                    assert to_state.find('IG') == 0
                    label = 'internal'
                    color = 'green'
                    width = 3

                label_to_use = None
                if color not in legend_colors:
                    label_to_use = label
                    legend_colors.add(color)

                ax.plot([-0.5 + ibin, 0.5 + ibin], [total + prob, total + prob], color=color, linewidth=width, linestyle=style, alpha=alpha, label=label_to_use)

                midpoint = 0.5*(prob + 2*total)
                # ax.text(ibin, midpoint, simplify_state_name(to_state))

                total += prob
    
            ibin += 1

        ax.get_xaxis().set_visible(False)
        # plt.gcf().subplots_adjust(bottom=0.36)
        plotting.mpl_finish(ax, self.base_plotdir + '/transitions', gene_name, ybounds=(-0.01, 1.01), xbounds=(-3, len(model.states) + 3), leg_loc=(0.5, 0.1))

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
parser.add_argument('--outdir', default=os.getenv('www'))
args = parser.parse_args()

if __name__ == '__main__':
    if not os.path.exists(args.outdir):
        raise Exception('output directory %s does not exist' % args.outdir)
    mplot = ModelPlotter(args, os.getenv('www') + '/modelplots')
