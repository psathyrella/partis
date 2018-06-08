#!/usr/bin/env python
import argparse
from collections import OrderedDict
import os
import glob
import sys
current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
if not os.path.exists(current_script_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % current_script_dir
sys.path.insert(1, current_script_dir)

import plotconfig
import plotting
import utils
import glutils
from hist import Hist

# ----------------------------------------------------------------------------------------
def get_hists_from_dir(dirname, histname, string_to_ignore=None):
    hists = {}
    for fname in glob.glob(dirname + '/*.csv'):
        varname = os.path.basename(fname).replace('.csv', '')
        if string_to_ignore is not None:
            varname = varname.replace(string_to_ignore, '')
        hists[varname] = Hist(fname=fname, title=histname)
    if len(hists) == 0:
        raise Exception('no csvs in %s' % dirname)
    return hists


# ----------------------------------------------------------------------------------------
def compare_directories(args, plotdirlist, outdir):
    utils.prep_dir(outdir, wildlings=['*.png', '*.svg', '*.csv'])

    # read hists from <plotdirlist>
    allhists = OrderedDict()
    allvars = set()  # all variables that appeared in any dir
    for idir in range(len(plotdirlist)):
        dirhists = get_hists_from_dir(plotdirlist[idir], args.names[idir])
        allvars |= set(dirhists.keys())
        allhists[args.names[idir]] = dirhists
    # then loop over all the <varname>s we found
    for varname in allvars:
        hlist = [allhists[dname].get(varname, Hist(1, 0, 1, title='null')) for dname in allhists]
        plot_single_variable(args, varname, hlist, outdir, pathnameclues=plotdirlist[0])

    plotting.make_html(outdir, n_columns=4)

# ----------------------------------------------------------------------------------------
def plot_single_variable(args, varname, hlist, outdir, pathnameclues):
    if varname in plotconfig.gene_usage_columns:
        hlist = plotting.add_bin_labels_not_in_all_hists(hlist)

    no_labels = False
    xline, bounds, figsize = None, None, None
    stats = args.extra_stats
    translegend = [0.0, -0.2]
    xtitle, ytitle = hlist[0].xtitle, hlist[0].ytitle
    if xtitle == '':  # arg, plotting.py thinks default should be None, hist.py thinks it's ''
        xtitle = None
    if '-mean-bins' in varname:
        raise Exception('darn, I was hoping I wasn\'t making these plots any more')
    plottitle = plotconfig.plot_titles[varname] if varname in plotconfig.plot_titles else varname

    ytitle = 'frequency' if args.normalize else 'counts'

    if 'mute-freqs/v' in pathnameclues or 'mute-freqs/d' in pathnameclues or 'mute-freqs/j' in pathnameclues:
        assert not args.normalize
        ytitle = 'mutation freq'

    if varname in plotconfig.gene_usage_columns:
        xtitle = 'allele'
        if hlist[0].n_bins == 2:
            stats = ' 0-bin'  # print the fraction of entries in the zero bin into the legend (i.e. the fraction correct)
    # elif hlist[0].bin_labels.count('') == hlist[0].n_bins + 2:
    #     xtitle = '???'

    line_width_override = None
    if args.performance_plots:
        if 'hamming_to_true_naive' in varname:
            xtitle = 'hamming distance'
            if '_normed' in varname:
                xtitle = 'fractional ' + xtitle
        elif '_vs_mute_freq' in varname:
            xtitle = 'mutation freq'
            ytitle = 'fraction correct'
            if varname[0] == 'v' or varname[0] == 'j':
                translegend = [-0.4, -0.4]
        elif varname.find('_gene') == 1:
            xtitle = ''
            ytitle = 'fraction correct'
        else:
            xtitle = 'inferred - true'
        bounds = plotconfig.true_vs_inferred_hard_bounds.setdefault(varname, None)
    else:
        bounds = plotconfig.default_hard_bounds.setdefault(varname, None)
        if bounds is None and 'insertion' in varname:
            bounds = plotconfig.default_hard_bounds.setdefault('all_insertions', None)
        if varname in plotconfig.gene_usage_columns:
            no_labels = True
            if 'j_' not in varname:
                figsize = (10, 5)
            line_width_override = 1
        elif 'per-gene-per-position/v' in pathnameclues:
            figsize = (20, 5)
            bounds = plotconfig.default_hard_bounds.setdefault(utils.unsanitize_name(varname), None)

    if 'IG' in varname or 'TR' in varname:
        if 'mute-freqs' in pathnameclues:
            gene = utils.unsanitize_name(varname)
            plottitle = gene  # + ' -- mutation frequency'
            xtitle = 'position'
            if utils.get_region(gene) == 'j':
                translegend = [0.1, 0.]  #(-0.35, -0.02)
            else:
                translegend = [0.15, -0.02]
            xline = None
            if args.glfo is not None:
                if utils.get_region(gene) in utils.conserved_codons[args.locus]:
                    xline = args.glfo[utils.conserved_codons[args.locus][utils.get_region(gene)] + '-positions'][gene]
        else:
            ilastdash = varname.rfind('-')
            gene = utils.unsanitize_name(varname[:ilastdash])
            base_varname = varname[ilastdash + 1 :]
            base_plottitle = plotconfig.plot_titles[base_varname] if base_varname in plotconfig.plot_titles else ''
            plottitle = gene + ' -- ' + base_plottitle

    if len(hlist) > 9:  # skootch it down so they (maybe) all fit
        translegend[1] -= 0.5
    if args.translegend is not None:  # override with the command line
        translegend = args.translegend
    if args.extra_stats == 'auto':  # kind of hackey
        if xtitle == 'inferred - true':
            stats = 'absmean'
        else:
            stats = 'mean'
    # draw that little #$*(!
    linewidths = [line_width_override, ] if line_width_override is not None else args.linewidths
    alphas = [0.6 for _ in range(len(hlist))]
    plotting.draw_no_root(hlist[0], plotname=varname, plotdir=outdir, more_hists=hlist[1:], write_csv=False, stats=stats, bounds=bounds,
                          shift_overflows=(os.path.basename(outdir) != 'gene-call'), plottitle=plottitle, colors=args.colors,
                          xtitle=xtitle, ytitle=ytitle, xline=xline, normalize=(args.normalize and '_vs_mute_freq' not in varname),
                          linewidths=linewidths, alphas=alphas, errors=True,
                          figsize=figsize, no_labels=no_labels, log=args.log, translegend=translegend)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--outdir', required=True)
parser.add_argument('--plotdirs', required=True)
parser.add_argument('--names', required=True)
parser.add_argument('--performance-plots', action='store_true')
parser.add_argument('--colors', default=':'.join(plotting.default_colors))
parser.add_argument('--linewidths', default=':'.join(plotting.default_linewidths))
parser.add_argument('--gldir', default='data/germlines/human')
parser.add_argument('--locus', default='igh')
parser.add_argument('--normalize', action='store_true')
parser.add_argument('--extra-stats')
parser.add_argument('--translegend')
parser.add_argument('--log', default='')

args = parser.parse_args()
args.plotdirs = utils.get_arg_list(args.plotdirs)
args.names = utils.get_arg_list(args.names)
args.colors = utils.get_arg_list(args.colors)
args.linewidths = utils.get_arg_list(args.linewidths)
args.translegend = utils.get_arg_list(args.translegend, floatify=True)
for iname in range(len(args.names)):
    args.names[iname] = args.names[iname].replace('@', ' ')

# if you just pass in one parent directory, we assume <args.names> contains the desired subdirs
if len(args.plotdirs) == 1:
    parentdir = args.plotdirs[0]
    args.plotdirs = [parentdir + '/' + n for n in args.names]

if len(args.plotdirs) != len(args.names):
    raise Exception('poorly formatted args:\n  %s\n  %s' % (' '.join(args.plotdirs), ' '.join(args.names)))

# if args.gldir is not 'none':
args.glfo = None
if os.path.exists(args.gldir):
    args.glfo = glutils.read_glfo(args.gldir, args.locus)

# figure out if there's subdirs we need to deal with
listof_plotdirlists, listof_outdirs = [], []
firstdir = args.plotdirs[0]
if len(glob.glob(firstdir + '/*.csv')) > 0:  # add the parent dirs if they've got csvs
    listof_plotdirlists.append(args.plotdirs)
    listof_outdirs.append(args.outdir)
for subdir in [d for d in os.listdir(firstdir) if os.path.isdir(firstdir + '/' + d)]:
    listof_plotdirlists.append([d + '/' + subdir for d in args.plotdirs])
    listof_outdirs.append(args.outdir + '/' + subdir)

for dlist, outdir in zip(listof_plotdirlists, listof_outdirs):
    compare_directories(args, dlist, outdir)
