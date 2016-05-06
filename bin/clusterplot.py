#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import argparse
import sys
import math
import numpy
import os
import glob
import csv
import seaborn as sns
import matplotlib.pyplot as plt
from subprocess import check_call
sns.set_style("ticks")

current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
if not os.path.exists(current_script_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % current_script_dir
sys.path.insert(1, current_script_dir)

from utils import get_arg_list, get_partition_from_str, correct_cluster_fractions
import plotting

"""
NOTE this is for plotting when you have multiple paths, i.e. when smc_particles > 1
"""

class ClusterPlot(object):
    def __init__(self, args, rebin=None, reco_info=None):
        self.args = args
        self.logprobs, self.n_clusters, self.adj_mis, self.logweights = {}, {}, {}, {}
        # self.clusters = {}
        # self.hists = {}
        # self.correct_cluster_fractions = {}
        self.final_logweights = []
        # self.tmp_n_true_clusters = None
        for fname in self.args.infnames:
            if self.args.debug:
                print fname
            with open(fname) as infile:
                for line in csv.DictReader(infile):
                    if not self.args.use_all_steps and int(line['n_procs']) != 1:  # skip any preliminary steps with more than one process
                        continue
                    ipath = int(line.get('path_index', 0))
                    if ipath != 0:
                        raise Exception('need to update some things for smc. see TODOs')
                    if ipath not in self.logprobs:
                        self.logprobs[ipath], self.n_clusters[ipath], self.adj_mis[ipath], self.logweights[ipath] = [], [], [], []
                        # self.clusters[ipath] = []
                        # self.correct_cluster_fractions[ipath] = []
                        # self.hists[ipath] = []
                    self.logprobs[ipath].append(float(line['logprob']))
                    self.n_clusters[ipath].append(int(line['n_clusters']))
                    self.adj_mis[ipath].append(float(line.get('adj_mi', -1)))
                    # partition = get_partition_from_str(line['partition'])  # TODO this is really wasteful to do this for every partition
                    # self.clusters[ipath].append(partition)
                    # self.correct_cluster_fractions[ipath].append(correct_cluster_fractions(partition, reco_info))
                    self.logweights[ipath].append(float(line.get('logweight', 1)))
                    # self.hists[ipath].append(plotting.get_cluster_size_hist(partition, rebin=rebin))
                    # if self.tmp_n_true_clusters is None:
                    #     self.tmp_n_true_clusters = int(line.get('n_true_clusters', -1))
    
        if len(self.logprobs) == 0:
            raise Exception('didn\'t read any lines from %s' % ' '.join(self.args.infnames))
        self.process_paths()
        self.plot()

    # ----------------------------------------------------------------------------------------
    def process_paths(self):
        max_logprob = {'logprobs' : [], 'adj_mis' : [], 'imaxes' : [], 'weights' : []}  # values at the point of maximum logprob
        max_adj_mi = {'logprobs' : [], 'adj_mis' : [], 'imaxes' : [], 'weights' : []}  # values at the point of maximum adj_mi
        self.imaxes = {'logprob' : None, 'adj_mi' : None}
        for ipath in self.logprobs.keys():
            imax_logprob = self.logprobs[ipath].index(max(self.logprobs[ipath]))
            max_logprob['logprobs'].append(self.logprobs[ipath][imax_logprob])
            max_logprob['adj_mis'].append(self.adj_mis[ipath][imax_logprob])
            max_logprob['imaxes'].append(imax_logprob)
            max_logprob['weights'].append(self.logweights[ipath][imax_logprob])
    
            imax_adj_mi = self.adj_mis[ipath].index(max(self.adj_mis[ipath]))
            max_adj_mi['logprobs'].append(self.logprobs[ipath][imax_adj_mi])
            max_adj_mi['adj_mis'].append(self.adj_mis[ipath][imax_adj_mi])
            max_adj_mi['imaxes'].append(imax_adj_mi)
            max_adj_mi['weights'].append(self.logweights[ipath][imax_adj_mi])
    
        def weighted_mean_rootvariance(vals, wgts):
            mean = numpy.average(vals, weights=wgts)
            variance = numpy.average((vals - mean)**2, weights=wgts)
            return (mean, math.sqrt(variance))
    
        self.imaxes['logprob'] = weighted_mean_rootvariance(max_logprob['imaxes'], max_logprob['weights'])
        self.imaxes['adj_mi'] = weighted_mean_rootvariance(max_adj_mi['imaxes'], max_adj_mi['weights'])
        delta_imaxes = [i_adj_mi - i_logprob for i_adj_mi, i_logprob in zip(max_adj_mi['imaxes'], max_logprob['imaxes'])]
        delta_wgts = [0.5*(adj_mi_wgt + logprob_wgt) for adj_mi_wgt, logprob_wgt in zip(max_adj_mi['weights'], max_logprob['weights'])]
        mean_delta_imax = weighted_mean_rootvariance(delta_imaxes, delta_wgts)
    
        max_logprob_means = {'logprobs' : weighted_mean_rootvariance(max_logprob['logprobs'], max_logprob['weights']),
                             'adj_mis' : weighted_mean_rootvariance(max_logprob['adj_mis'], max_logprob['weights'])}
        max_adj_mi_means = {'logprobs' : weighted_mean_rootvariance(max_adj_mi['logprobs'], max_adj_mi['weights']),
                            'adj_mis' : weighted_mean_rootvariance(max_adj_mi['adj_mis'], max_adj_mi['weights'])}

        if self.args.debug:
            print 'max logprob step %.1f:' % self.imaxes['logprob'][0]
            print '     logprob: %.1f +/- %.1e' % max_logprob_means['logprobs']
            print '      adj_mi: %.4f +/- %.1e' % max_logprob_means['adj_mis']
            print 'max adj_mi step %.1f:' % self.imaxes['adj_mi'][0]
            print '     logprob: %.1f +/- %.1e' % max_adj_mi_means['logprobs']
            print '      adj_mi: %.4f +/- %.1e' % max_adj_mi_means['adj_mis']
            print 'change: %.1f steps' % mean_delta_imax[0]
        logprob_at_adj_mi = max_adj_mi_means['logprobs'][0]
        logprob_at_logprob = max_logprob_means['logprobs'][0]
        if self.args.debug:
            print '     logprob: %.1f (delta %.3f)' % (logprob_at_adj_mi - logprob_at_logprob, (logprob_at_adj_mi - logprob_at_logprob) / logprob_at_logprob)
            print '      adj_mi: %.4f' % (max_adj_mi_means['adj_mis'][0] - max_logprob_means['adj_mis'][0])

        self.adj_mi_at_max_logprob = max_logprob_means['adj_mis'][0]
        # tmp_ipath = 0  # NOTE only for the first path a.t.m. TODO fix that
        # self.tmp_n_clusters = self.n_clusters[tmp_ipath][int(self.imaxes['logprob'][0])]  # [0] is because imaxes is (mean, err)
        # # self.tmp_cluster_size_hist = self.hists[tmp_ipath][int(self.imaxes['logprob'][0])]
        # if self.args.debug:
        #     print '   n_clusters at max logprob for zeroth path: %d' % self.tmp_n_clusters

    # ----------------------------------------------------------------------------------------
    def plot(self):

        max_length = -1
        for ipath in self.logprobs.keys():
            if len(self.logprobs[ipath]) > max_length:
                max_length = len(self.logprobs[ipath])
        
        min_logprob, max_logprob = None, None
        min_logprobs, max_logprobs = {}, {}  #[None for _ in range(len(logprobs.keys()))], [None for _ in range(len(logprobs.keys()))]
        if self.args.debug:
            print '%d paths' % len(self.logprobs.keys())
        for ipath in self.logprobs.keys():
            # while len(logprobs[ipath]) < max_length:
            #     logprobs[ipath].append(None)
            #     adj_mis[ipath].append(None)
            min_logprobs[ipath], max_logprobs[ipath] = None, None
            for il in range(len(self.logprobs[ipath])):
                # min/max for this path
                if min_logprobs[ipath] is None or self.logprobs[ipath][il] < min_logprobs[ipath]:
                    min_logprobs[ipath] = self.logprobs[ipath][il]
                if max_logprobs[ipath] is None or self.logprobs[ipath][il] > max_logprobs[ipath]:
                    max_logprobs[ipath] = self.logprobs[ipath][il]
                # global min/max
                if min_logprob is None or self.logprobs[ipath][il] < min_logprob:
                    min_logprob = self.logprobs[ipath][il]
                if max_logprob is None or self.logprobs[ipath][il] > max_logprob:
                    max_logprob = self.logprobs[ipath][il]
        
        if 'logprob' in self.args.normalize_axes:
            for ipath in self.logprobs.keys():
                # print ipath, min_logprobs[ipath], max_logprobs[ipath]
                for il in range(len(self.logprobs[ipath])):
                    if self.logprobs[ipath][il] is not None:
                        self.logprobs[ipath][il] = max_logprobs[ipath] / self.logprobs[ipath][il]
                        # print logprobs[ipath][il]
        
            min_logprob = 0.
            max_logprob = 1.
        
        # fig = plt.figure(1)
        # fig.clf()
        fig, ax = plt.subplots()
        # sns.despine(fig=fig,offset=10, trim=True)
        adj_mi_color = '#980000'
        logprob_color = '#1947A3'
        
        ax2 = ax.twinx()
        nxbins = 8
        nybins = 5
        
        yextrafactor = 1.
        xmin = 0
        min_adj_mi = 0.
        markersize = 30
        
        max_adj_mi = 1.
        if self.args.xbounds is None:
            self.args.xbounds = (0, 1 if 'step' in self.args.normalize_axes else max_length)
        if self.args.adjmi_bounds is not None:
            min_adj_mi = self.args.adjmi_bounds[0]
            max_adj_mi = self.args.adjmi_bounds[1]
        if self.args.logprob_bounds is not None:
            min_logprob = self.args.logprob_bounds[0]
            max_logprob = self.args.logprob_bounds[1]
        
        if __name__ == '__main__':
            ax.set_xlim(self.args.xbounds[0], self.args.xbounds[1])
            ax.set_ylim(min_adj_mi, yextrafactor * max_adj_mi)
            ax2.set_ylim(min_logprob, yextrafactor * max_logprob)
            ax.set_xlabel('agglomeration step', fontweight='bold')
            ax.set_ylabel('adjusted MI', color=adj_mi_color, fontweight='bold')
            ax2.set_ylabel('log prob', color=logprob_color, fontweight='bold')
            fig.tight_layout()
            plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)
            
            ax.locator_params(nbins=nxbins, axis='x')
            ax.locator_params(nbins=nybins, axis='y')
            
            for ipath in self.logprobs.keys():
                if 'step' in self.args.normalize_axes:
                    steps = [float(i) / len(self.logprobs[ipath]) for i in range(len(self.logprobs[ipath]))]
                else:
                    steps = [i for i in range(len(self.logprobs[ipath]))]
                sizes = [markersize for i in range(len(self.logprobs[ipath]))]
                fig_logprob = ax2.plot(steps, self.logprobs[ipath], color=logprob_color, alpha=.5, linewidth=1)
                if not self.args.is_data:
                    fig_adj_mi = ax.plot(steps, self.adj_mis[ipath], color=adj_mi_color, alpha=1, linewidth=1)
            
                fig_logprob_sc = ax2.scatter(steps, self.logprobs[ipath], color=logprob_color, alpha=.5, s=sizes)
                if not self.args.is_data:
                    fig_adj_mi_sc = ax.scatter(steps, self.adj_mis[ipath], color=adj_mi_color, alpha=1, s=sizes)
            
            
            xlinepos_logprob = [self.imaxes['logprob'][0] for _ in range(2)]
            if not self.args.is_data:
                xlinepos_adj_mi = [self.imaxes['adj_mi'][0] for _ in range(2)]
            if 'step' in self.args.normalize_axes:
                assert False
                # xlinepos_logprob = [ float(x) / len(logprobs xlinepos_logprob
            plt.plot(xlinepos_logprob, [min_logprob, max_logprob], color=logprob_color, linestyle='--')  #, color='k', linestyle='-', linewidth=2)
            if not self.args.is_data:
                plt.plot(xlinepos_adj_mi, [min_adj_mi, max_adj_mi], color=adj_mi_color, linestyle='--')  #, color='k', linestyle='-', linewidth=2)
            
            plotdir = os.getenv('www') + '/tmp'
            plotname = 'foo'
            plt.savefig(plotdir + '/' + plotname + '.svg')
            plt.close()  # wtf? doesn't eliminate warning about too many figures
            check_call(['permissify-www', plotdir])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', action='store_true')
    parser.add_argument('--infnames')
    parser.add_argument('--outfname')
    parser.add_argument('--normalize-axes', default=[])
    parser.add_argument('--use-all-steps', action='store_true')
    parser.add_argument('--is-data', action='store_true')
    parser.add_argument('--xbounds')
    parser.add_argument('--logprob-bounds')
    parser.add_argument('--adjmi-bounds')
    args = parser.parse_args()
    args.infnames = get_arg_list(args.infnames)
    if args.xbounds is not None:
        args.xbounds = get_arg_list(args.xbounds, floatify=True)
    if args.logprob_bounds is not None:
        args.logprob_bounds = get_arg_list(args.logprob_bounds, floatify=True)
    if args.adjmi_bounds is not None:
        args.adjmi_bounds = get_arg_list(args.adjmi_bounds, floatify=True)
    if len(args.normalize_axes) > 0:
        args.normalize_axes = get_arg_list(args.normalize_axes)
    
    fsize = 26
    mpl.rcParams.update({
        'font.size': 26,
        'axes.labelsize': 26,
        'xtick.labelsize':20,
        'ytick.labelsize':20,
        'font.family': 'Lato',
        'font.weight': 600,
        'axes.labelweight': 600,
        'legend.fontsize': fsize,
        'axes.titlesize': fsize
    })
    
    cplot = ClusterPlot(args)
