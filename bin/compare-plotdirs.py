#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
from collections import OrderedDict
import os
import glob
import sys
import colored_traceback.always
import copy
import partis.plotconfig as plotconfig
import partis.plotting as plotting
import partis.utils as utils
import partis.glutils as glutils
from partis.hist import Hist
import partis.treeutils as treeutils

xtitledict = copy.deepcopy(plotting.legends)
xtitledict.update(plotconfig.xtitles)
xtitledict.update(treeutils.legtexts)

ptitledict = copy.deepcopy(plotting.legends)
ptitledict.update(plotconfig.plot_titles)
ptitledict.update(treeutils.legtexts)

# ----------------------------------------------------------------------------------------
def get_hists_from_dir(dirname, histname, string_to_ignore=None):
    hists = {}
    file_patterns = args.file_glob_str.split(',')
    for pattern in file_patterns:
        for fname in glob.glob('%s/%s' % (dirname, pattern)):
            varname = os.path.basename(fname)
            for rstr in args.file_replace_strs:
                varname = varname.replace(rstr, '')
            if string_to_ignore is not None:
                varname = varname.replace(string_to_ignore, '')
            hists[varname] = Hist(fname=fname, title=histname)
    if len(hists) == 0:
        print('    no csvs found%s in %s' % ('' if args.file_glob_str is None else ' with --file-glob-str \'%s\''%args.file_glob_str, dirname))
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

    plotting.make_html(outdir, n_columns=4, new_table_each_row=True)

# ----------------------------------------------------------------------------------------
def plot_single_variable(args, varname, hlist, outdir, pathnameclues):
    if varname in plotconfig.gene_usage_columns:
        hlist = plotting.add_bin_labels_not_in_all_hists(hlist)

    no_labels = False
    xline, bounds, figsize = None, None, None
    stats = args.extra_stats
    translegend = [0.0, -0.2]
    xtitle, ytitle = hlist[0].xtitle, hlist[0].ytitle
    bounds, xticks, xticklabels = args.xbounds, args.xticks, None
    if xtitle == '':  # arg, plotting.py thinks default should be None, hist.py thinks it's ''
        xtitle = None
    if '-mean-bins' in varname:
        raise Exception('darn, I was hoping I wasn\'t making these plots any more')
    plottitle = args.plottitle
    if plottitle is None:
        plottitle = ptitledict.get(varname, varname)

    ytitle = 'fraction of total' if args.normalize else 'counts'

    if 'mute-freqs/v' in pathnameclues or 'mute-freqs/d' in pathnameclues or 'mute-freqs/j' in pathnameclues:
        assert not args.normalize
        ytitle = 'mutation freq'

    if varname in plotconfig.gene_usage_columns:
        xtitle = 'allele'
        if hlist[0].n_bins == 2:
            stats = '0-bin'  # print the fraction of entries in the zero bin into the legend (i.e. the fraction correct)
            xtitle = None
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
        if bounds is None:
            bounds = plotconfig.default_hard_bounds.setdefault(varname, None)
        if varname in plotconfig.gene_usage_columns:
            # no_labels = True  # not sure why i wanted these labels turned off?
            if 'j_' not in varname:
                figsize = (10, 5)
            line_width_override = 1
        elif 'per-gene-per-position/v' in pathnameclues:
            figsize = (20, 5)
            if bounds is None:
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

    if varname == 'cluster-sizes':
        xtitle = 'cluster size'
        ytitle = 'fraction of clusters' if args.normalize else 'N clusters'
        plottitle = ''
        xticks, xticklabels = plotting.get_cluster_size_xticks(hlist=hlist)  # it would be better to use all the hists, but i think it'll just screw up the ticks
        import matplotlib.pyplot as plt
    if varname in ['func-per-drop', 'nonfunc-per-drop']:
        bounds = (-0.5, 15.5)
    if 'subtree-purity' in varname:
        if 'size' in varname:
            if args.log == '':
                args.log = 'xy'
            xticks = [1, 2, 3, 5, 10, 15, 20]
            xticklabels = ['1', '2', '3', '5', '10', '15', '20']
    if bounds is not None and any(h.xmin > bounds[1] or h.xmax < bounds[0] for h in hlist):  # if any hist is entirely outside of <bounds>, widen the <bounds>
        hist_bounds = [h.get_filled_bin_xbounds() for h in hlist]
        overall_min, overall_max = [f(vals+(b,)) for f, vals, b in zip([min, max], zip(*hist_bounds), bounds)]
        if overall_min < bounds[0] or overall_max > bounds[1]:
            print('    %s (%s) overriding bounds %s with filled hist bounds to get %s' % (utils.wrnstr(), varname, bounds, (overall_min, overall_max)))
        bounds = (overall_min, overall_max)

    if xtitle is None:
        xtitle = xtitledict.get(varname)

    if args.add_to_title is not None:
        plottitle += args.add_to_title

    if len(hlist) > 9:  # skootch it down so they (maybe) all fit
        translegend[1] -= 0.5
    if args.translegend is not None:  # override with the command line
        translegend = args.translegend
    if varname == 'paired-uids-per-uid':
        translegend = [translegend[0] + 0.15, translegend[1] - 0.3]
    if args.extra_stats == 'auto':  # kind of hackey
        if xtitle == 'inferred - true':
            stats = 'absmean'
        else:
            stats = 'mean'
    # draw that little #$*(!
    linewidths = [line_width_override, ] if line_width_override is not None else args.linewidths
    if args.alphas is None or len(args.alphas) != len(hlist):
        if args.alphas is not None and len(args.alphas) != len(hlist):
            print('  %s --alphas wrong length, using first entry for all' % utils.wrnstr())
        args.alphas = [0.6 if args.alphas is None else args.alphas[0] for _ in range(len(hlist))]
    shift_overflows = os.path.basename(outdir) != 'gene-call' and 'func-per-drop' not in varname
    plotting.draw_no_root(hlist[0], plotname=varname, plotdir=outdir, more_hists=hlist[1:], write_csv=False, stats=stats, bounds=bounds, ybounds=args.ybounds,
                          shift_overflows=shift_overflows, plottitle=plottitle, colors=args.colors,
                          xtitle=xtitle if args.xtitle is None else args.xtitle, ytitle=ytitle if args.ytitle is None else args.ytitle, xline=xline, normalize=(args.normalize and '_vs_mute_freq' not in varname),
                          linewidths=linewidths, linestyles=args.linestyles, markersizes=args.markersizes, alphas=args.alphas, errors=not args.no_errors, remove_empty_bins=True, #='y' in args.log,
                          figsize=figsize, no_labels=no_labels, log=args.log, translegend=translegend, xticks=xticks, xticklabels=xticklabels, square_bins=args.square_bins)

    if args.swarm_meta_key is not None:
        plotvals = {h.title : [h.get_bin_centers()[i] for i in h.ibiniter(True) for _ in range(int(h.bin_contents[i]))] for h in hlist}
        plotting.stack_meta_hists(varname, outdir, args.swarm_meta_key, plotvals, colors={h.title : c for h, c in zip(hlist, args.colors)}, xtitle=xtitle, swarm_plots=True, no_hist=True, xticks=xticks)

# ----------------------------------------------------------------------------------------
helpstr = """
Compare csv histogram plot files across multiple directories
    ./bin/compare-plotdirs.py --outdir _output/tmp-plots --plotdirs docs/example-plots/sw/mute-freqs/overall:docs/example-plots/hmm/mute-freqs/overall:docs/example-plots/multi-hmm/mute-freqs/overall --names sw:hmm:multi-hmm --normalize
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('--outdir', required=True, help='Output directory to which to write the resulting comparison plots. A summary .html file is also written to <outdir>.html')
parser.add_argument('--plotdirs', required=True, help='Colon-separated list of input plot directories, each of which must have identical structure. Looks for svgs first in each dir, but then also in the subdirs of each dir (so e.g. if each of them have a/, b/, and c/ subdirs, this script will make a separate comparison of a/, b/, and c/)')
parser.add_argument('--names', required=True, help='colon-separated list of names/labels corresponding to --plotdirs (use @ as space)')
parser.add_argument('--performance-plots', action='store_true', help='set to true if these are annotation performance plots, i.e. made with --plot-annotation-performance (this makes the axis labels more sensible)')
parser.add_argument('--colors', default=':'.join(plotting.default_colors), help='color-separated list of colors to cycle through for the plotdirs')
parser.add_argument('--alphas')
parser.add_argument('--linewidths', default=':'.join(plotting.default_linewidths), help='colon-separated list of linewidths to cycle through')
parser.add_argument('--markersizes', default=':'.join(plotting.default_markersizes), help='colon-separated list of linewidths to cycle through')
parser.add_argument('--linestyles', help='colon-separated list of linestyles to cycle through')
parser.add_argument('--gldirs', help='On plots showing mutation vs individual gene positions, if you\'d like a dashed veritcal line showing conserved codon positions, set this as a colon-separated list of germline info dirs corresponding to each plotdir') #, default=['data/germlines/human'])
parser.add_argument('--locus', default='igh')
parser.add_argument('--normalize', action='store_true', help='If set, the histograms from each plotdir are normalized (each bin contents divided by the integral) before making the comparison (e.g. for comparing different size samples).')
parser.add_argument('--extra-stats', help='if set, adds extra stat to legend, e.g. \'mean\', \'absmean\', \'auto\'')
parser.add_argument('--translegend', help='colon-separated list of x, y values with which to translate all the legends')
parser.add_argument('--log', default='', help='Display these axes on a log scale, set to either \'x\', \'y\', or \'xy\'')
parser.add_argument('--make-parent-html', action='store_true', help='after doing everything within subdirs, make a single html in the main/parent dir with all plots from subdirs')
parser.add_argument('--add-to-title', help='string to append to existing title (use @ as space)')
parser.add_argument('--file-glob-str', default='*.csv', help='Shell glob style regex for matching plot files. Separate multiple patterns with \',\' (and no curly braces {})')
parser.add_argument('--file-replace-strs', default='.csv', help='colon-separated list of strings to remove frome file base name to get variable name')
parser.add_argument('--xbounds')
parser.add_argument('--ybounds')
parser.add_argument('--xticks')
parser.add_argument('--plottitle')
parser.add_argument('--xtitle')
parser.add_argument('--ytitle')
parser.add_argument('--no-errors', action='store_true')
parser.add_argument('--single-plotdir', action='store_true')
parser.add_argument('--square-bins', action='store_true')
parser.add_argument('--swarm-meta-key', help='if set, also make swarm plots, pretending that each hist\'s title is the value for this "fake" meta info key, and treat each bin\'s entries as observations at the bin center\'s value')

args = parser.parse_args()
args.plotdirs = utils.get_arg_list(args.plotdirs)
args.names = utils.get_arg_list(args.names)
args.alphas = utils.get_arg_list(args.alphas, floatify=True)
args.colors = utils.get_arg_list(args.colors)
args.linewidths = utils.get_arg_list(args.linewidths, intify=True)
args.linestyles = utils.get_arg_list(args.linestyles)
args.markersizes = utils.get_arg_list(args.markersizes, intify=True)
args.gldirs = utils.get_arg_list(args.gldirs)
args.translegend = utils.get_arg_list(args.translegend, floatify=True)
args.xbounds = utils.get_arg_list(args.xbounds, floatify=True)
args.ybounds = utils.get_arg_list(args.ybounds, floatify=True)
args.xticks = utils.get_arg_list(args.xticks, floatify=True)
args.file_replace_strs = utils.get_arg_list(args.file_replace_strs)
for iname in range(len(args.names)):
    args.names[iname] = args.names[iname].replace('@', ' ')
if args.add_to_title is not None:
    args.add_to_title = args.add_to_title.replace('@', ' ')

if len(args.plotdirs) == 1 and not args.single_plotdir:
    print('  --plotdirs is length 1 (and --single-plotdir wasn\'t set), so assuming --names has the desired subdirs')
    parentdir = args.plotdirs[0]
    args.plotdirs = [parentdir + '/' + n for n in args.names]

if len(args.plotdirs) != len(args.names):
    raise Exception('poorly formatted args:\n  %s\n  %s' % (' '.join(args.plotdirs), ' '.join(args.names)))

# make a merged glfo from all the gldirs
args.glfo = None
if args.gldirs is not None:
    for gldir in [gd for gd in args.gldirs if os.path.exists(gd)]:
        tmpglfo = glutils.read_glfo(gldir, args.locus)
        if args.glfo is None:
            args.glfo = tmpglfo
        else:
            args.glfo, _ = glutils.get_merged_glfo(args.glfo, tmpglfo)  # discard name_mapping since we don't have annotations here

if any(not os.path.isdir(d) for d in args.plotdirs):
    print('   at least one of --plotdirs doesn\'t exist: %s' % ' '.join(d for d in args.plotdirs if not os.path.isdir(d)))
    sys.exit(0)

listof_plotdirlists, listof_outdirs = [], []
# first add the main/parent dir, if it has csvs
firstdir = args.plotdirs[0]
if len(glob.glob(firstdir + '/*.csv')) > 0:
    listof_plotdirlists.append(args.plotdirs)
    listof_outdirs.append(args.outdir)
else:
    print('    no csvs in main/parent dir %s' % firstdir)
# then figure out if there's subdirs we need to deal with
added_subds = []

for subdir in [d for d in os.listdir(firstdir) if os.path.isdir(firstdir + '/' + d)]:
    listof_plotdirlists.append([d + '/' + subdir for d in args.plotdirs])
    listof_outdirs.append(args.outdir + '/' + subdir)
    added_subds.append(subdir)
if len(added_subds) > 0:
    print('  added %d subdirs: %s' % (len(added_subds), ' '.join(added_subds)))

for dlist, outdir in zip(listof_plotdirlists, listof_outdirs):
    compare_directories(args, dlist, outdir)

if args.make_parent_html:  # didn't really test this very well
    fnoutstr, _ = utils.simplerun('find %s -type f -name *.svg' % args.outdir, return_out_err=True)
    plotting.make_html(args.outdir, fnames=[fnoutstr.strip().split('\n')], new_table_each_row=True)
