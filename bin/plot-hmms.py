#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import argparse
from collections import OrderedDict
import glob
import operator
import yaml
import sys
from subprocess import check_call
import matplotlib as mpl
from io import open
mpl.use('Agg')
import matplotlib.pyplot as plt

import partis.plotting as plotting
import partis.paramutils as paramutils
import partis.utils as utils

# ----------------------------------------------------------------------------------------
class ModelPlotter(object):
    def __init__(self, args, base_plotdir):
        self.base_plotdir = base_plotdir

        self.eps_to_skip = 1e-3
        print('skipping eps %f' % self.eps_to_skip)

        plot_types = ('transitions', 'emissions')
        for ptype in plot_types:
            plotdir = self.base_plotdir + '/' + ptype
            utils.prep_dir(plotdir, wildlings=['*.png', '*.svg'])

        if args.hmmdir != None:
            self.filelist = glob.glob(args.hmmdir + '/*.yaml')
        else:
            self.filelist = utils.get_arg_list(args.infiles)
        if len(self.filelist) == 0:
            raise Exception('zero files passed to modelplotter')

    # ----------------------------------------------------------------------------------------
    def plot(self):
        for infname in self.filelist:
            gene_name = os.path.basename(infname).replace('.yaml', '')  # the sanitized name, actually
            with open(infname) as infile:
                model = yaml.load(infile, Loader=yaml.Loader)
                self.make_transition_plot(gene_name, model)
                self.make_emission_plot(gene_name, model)

    # ----------------------------------------------------------------------------------------
    def make_transition_plot(self, gene_name, model):
        """ NOTE shares a lot with make_mutefreq_plot() in python/paramutils.py """
        fig, ax = plotting.mpl_init()
        fig.set_size_inches(plotting.plot_ratios[utils.get_region(gene_name)])

        ibin = 0
        print(utils.color_gene(utils.unsanitize_name(gene_name)))
        legend_colors = set()  # add a color to this the first time you plot it
        for state in model.states:

            # bin label
            ax.text(-0.5 + ibin, -0.075, paramutils.simplify_state_name(state.name), rotation='vertical', size=8)

            sorted_to_states = {}
            for name in state.transitions.keys():
                if name.find('IG') == 0 or name.find('TR') == 0:
                    sorted_to_states[name] = int(paramutils.simplify_state_name(name))
                else:
                    sorted_to_states[name] = name
            sorted_to_states = sorted(list(sorted_to_states.items()), key=lambda x: str(x[1]))

            total = 0.0
            for to_state, simple_to_state in sorted_to_states:

                prob = state.transitions[to_state]

                alpha = 0.6
                width = 3

                if 'insert' in str(simple_to_state):
                    label = 'insert'
                    color = '#3498db'  # blue
                elif str(simple_to_state) == 'end':
                    label = 'end'
                    color = 'red'
                else:  # regional/internal states
                    assert to_state.find('IG') == 0 or to_state.find('TR') == 0
                    label = 'internal'
                    color = 'green'

                label_to_use = None
                if color not in legend_colors:
                    label_to_use = label
                    legend_colors.add(color)

                # horizontal line at height total+prob
                ax.plot([-0.5 + ibin, 0.5 + ibin], [total + prob, total + prob], color=color, linewidth=width, alpha=alpha, label=label_to_use)

                # vertical line from total to total + prob
                ax.plot([ibin, ibin], [total + 0.01, total + prob], color=color, alpha=alpha, linewidth=width)

                midpoint = 0.5*(prob + 2*total)
                # ax.text(ibin, midpoint, paramutils.simplify_state_name(to_state))  # nicely labels the midpoint of the chunk between lines, but there isn't really room for it

                total += prob
    
            ibin += 1

        ax.get_xaxis().set_visible(False)
        plotting.mpl_finish(ax, self.base_plotdir + '/transitions', gene_name, ybounds=(-0.01, 1.01), xbounds=(-3, len(model.states) + 3), leg_loc=(0.95, 0.1), adjust={'left' : 0.1, 'right' : 0.8}, leg_prop={'size' : 8})

    # ----------------------------------------------------------------------------------------
    def make_emission_plot(self, gene_name, model):
        plotting_info = []
        for state in model.states:
            if state.emissions is None:
                assert state.name == 'init'
                continue
            plotting_info.append({
                'name' : state.name,
                'nuke_freqs' : state.emissions['probs'],
                'gl_nuke' : state.extras['germline'] if 'germline' in state.extras else None
            })

        paramutils.make_mutefreq_plot(self.base_plotdir + '/emissions', gene_name, plotting_info, debug=True)

    # ----------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument('--hmmdir', help='directory with .yaml hmm model files, e.g. test/reference-results/test/parameters/simu/hmm/hmms')
parser.add_argument('--infiles', help='colon-separated list of .yaml hmm model files (either set this, or set --hmmdir)')
parser.add_argument('--outdir', required=True)
args = parser.parse_args()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

if args.hmmdir is None and args.infiles is None:
    raise Exception('have to specify either --hmmdir or --infiles')

if __name__ == '__main__':
    print('  %s the top line in the emission plots is usually yellow because the three non-germline bases are equally likely, and G comes last when sorted alphabetically' % utils.color('red', 'note'))
    if not os.path.exists(args.outdir):
        raise Exception('output directory %s does not exist' % args.outdir)
    mplot = ModelPlotter(args, args.outdir) # + '/modelplots')
    mplot.plot()
