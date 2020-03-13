#!/usr/bin/env python
# has to be its own script, since ete3 requires its own god damn python version, installed in a separated directory
import time
import yaml
import itertools
import glob
import argparse
import copy
import random
import os
import tempfile
import subprocess
import sys
import colored_traceback.always
from collections import OrderedDict
import numpy
import math
try:
    import ete3
except ImportError:
    raise Exception('couldn\'t find the ete3 module. Either:\n          - it isn\'t installed (use instructions at http://etetoolkit.org/download/) or\n          - $PATH needs modifying (typically with the command \'% export PATH=~/anaconda_ete/bin:$PATH\')')

# ----------------------------------------------------------------------------------------
scolors = {
    'novel' : '#ffc300',  # 'Gold'
    'data' : 'LightSteelBlue',
    'pale-green' : '#85ad98',
    'pale-blue' : '#94a3d1',
    'tigger-default' : '#d77c7c', #'#c32222',  # red
    'igdiscover' : '#85ad98', #'#29a614',  # green
    'partis' : '#94a3d1', #'#2455ed',  # blue
    'lbi' : '#94a3d1',
}

# listcolors = [plotting.getgrey('medium') for _ in range(10)]
listfaces = [
    'red',
    'blue',
    'green',
]
used_colors, used_faces = {}, {}
simu_colors = OrderedDict((
    ('ok', 'DarkSeaGreen'),
    ('missing', '#d77c7c'),
    ('spurious', '#a44949'),
))
def get_scale_min(metric, vals):  # only make the color scale go down to here
    if metric == 'cons-dist-aa':
        return max(vals) - 10
    else:
        return min(vals)

# ----------------------------------------------------------------------------------------
def read_input(args):
    with open(args.treefname) as treefile:
        treestr = treefile.read().strip()
    treestr = treestr.replace('[&R] ', '').replace('\'', '')

    return {'treestr' : treestr}

# ----------------------------------------------------------------------------------------
min_size = 1.5
max_size = 10
opacity = 0.65
fsize = 7

# ----------------------------------------------------------------------------------------
def set_delta_affinities(etree, affyfo):  # set change in affinity from parent for each node, and return a list of all such affinity changes (for normalizing the cmap)
    delta_affyvals = []
    for node in etree.traverse():
        if node.up is None or any(n.name not in affyfo or affyfo[n.name] is None for n in (node, node.up)):
            node.add_feature('affinity_change', None)
            continue
        node.add_feature('affinity_change', affyfo[node.name] - affyfo[node.up.name])
        delta_affyvals.append(affyfo[node.name] - affyfo[node.up.name])

    return delta_affyvals

# ----------------------------------------------------------------------------------------
def get_size(vmin, vmax, val):
    return min_size + (val - vmin) * (max_size - min_size) / (vmax - vmin)

# ----------------------------------------------------------------------------------------
def add_legend(tstyle, varname, all_vals, smap, info, start_column, add_missing=False, add_sign=None, reverse_log=False, n_entries=5, fsize=4, no_opacity=False):  # NOTE very similar to add_smap_legend() in plot_2d_scatter() in python/lbplotting.py
    if len(all_vals) == 0:
        return
    assert add_sign in [None, '-', '+']
    tstyle.legend.add_face(ete3.TextFace('   %s ' % varname, fsize=fsize), column=start_column)
    min_val, max_val = get_scale_min(args.lb_metric, all_vals), max(all_vals)
    if min_val == max_val:
        return
    max_diff = (max_val - min_val) / float(n_entries - 1)
    val_list = list(numpy.arange(min_val, max_val + utils.eps, max_diff))  # first value is exactly <min_val>, last value is exactly <max_val> (eps is to keep it from missing the last one)
    # if add_sign is not None and add_sign == '-':  # for negative changes, we have the cmap using abs() and want to legend order to correspond
    #     val_list = reversed(val_list)  # arg, this breaks something deep in the legend maker, not sure what
    key_list = [None for _ in val_list]
    if add_missing:
        val_list += [None]
        key_list += ['missing!']  # doesn't matter what the last one is as long as it isn't in <affyfo>
    for val, key in zip(val_list, key_list):
        tstyle.legend.add_face(ete3.TextFace('', fsize=fsize), column=start_column)
        if smap is None:
            sz = get_size(min_val, max_val, val)
            rface = ete3.RectFace(sz, sz, bgcolor=plotting.getgrey(), fgcolor=None)
        else:
            rface = ete3.RectFace(6, 6, bgcolor=plotting.get_smap_color(smap, info, key=key, val=val), fgcolor=None)
        if not no_opacity:
            rface.opacity = opacity
        tstyle.legend.add_face(rface, column=start_column + 1)
        fstr = '%.1f' if args.lb_metric == 'cons-dist-aa' else '%.2f'
        tstyle.legend.add_face(ete3.TextFace((('  %s'+fstr) % (add_sign if add_sign is not None else '', math.exp(val) if reverse_log else val)) if key is None else '  missing', fsize=fsize), column=start_column + 2)

# ----------------------------------------------------------------------------------------
def set_meta_styles(args, etree, tstyle):
    lbfo = args.metafo[args.lb_metric]
    if args.lb_metric == 'lbr':  # remove zeroes
        lbfo = {u : (math.log(v) if args.log_lbr else v) for u, v in lbfo.items() if v > 0}
    lbvals = lbfo.values()
    if len(lbvals) == 0:
        return
    lb_smap = plotting.get_normalized_scalar_map(lbvals, 'viridis', hard_min=get_scale_min(args.lb_metric, lbvals) if args.lb_metric=='cons-dist-aa' else None)
    lb_min, lb_max = min(lbvals), max(lbvals)

    affyfo = None
    if args.affy_key in args.metafo and set(args.metafo[args.affy_key].values()) != set([None]):
        affyfo = args.metafo[args.affy_key]
        if args.lb_metric in affy_metrics:
            affyvals = affyfo.values()
            affy_smap = plotting.get_normalized_scalar_map([a for a in affyvals if a is not None], 'viridis')
        elif args.lb_metric in delta_affy_metrics:
            delta_affyvals = set_delta_affinities(etree, affyfo)
            delta_affy_increase_smap = plotting.get_normalized_scalar_map([v for v in delta_affyvals if v > 0], 'Reds', remove_top_end=True) if len(delta_affyvals) > 0 else None
            delta_affy_decrease_smap = plotting.get_normalized_scalar_map([abs(v) for v in delta_affyvals if v < 0], 'Blues', remove_top_end=True) if len(delta_affyvals) > 0 else None
        else:
            assert False

    for node in etree.traverse():
        node.img_style['size'] = 0
        rfsize = 0
        bgcolor = plotting.getgrey()
        if args.lb_metric in affy_metrics:
            if node.name not in lbfo:  # really shouldn't happen
                print '  %s missing lb info for node \'%s\'' % (utils.color('red', 'warning'), node.name)
                continue
            if affyfo is not None:
                rfsize = get_size(lb_min, lb_max, lbfo[node.name])
                if node.name in affyfo:
                    bgcolor = plotting.get_smap_color(affy_smap, affyfo, key=node.name)
            else:
                rfsize = 5
                bgcolor = plotting.get_smap_color(lb_smap, lbfo, key=node.name)
        elif args.lb_metric in delta_affy_metrics:
            node.img_style['vt_line_color'] = plotting.getgrey()  # if they're black, it's too hard to see the large changes in affinity, since they're very dark (at least with current color schemes)
            # rfsize = get_size(lb_min, lb_max, lbfo[node.name]) if node.name in lbfo else 1.5
            rfsize = 5 if node.name in lbfo else 1.5
            bgcolor = plotting.get_smap_color(lb_smap, lbfo, key=node.name)
            if affyfo is not None and delta_affy_increase_smap is not None and node.affinity_change is not None:
                # tface = ete3.TextFace(('%+.4f' % node.affinity_change) if node.affinity_change != 0 else '0.', fsize=3)
                # node.add_face(tface, column=0)
                if node.affinity_change > 0:  # increase
                    node.img_style['hz_line_color'] = plotting.get_smap_color(delta_affy_increase_smap, None, val=node.affinity_change)
                    node.img_style['hz_line_width'] = 1.2
                elif node.affinity_change < 0:  # decrease
                    node.img_style['hz_line_color'] = plotting.get_smap_color(delta_affy_decrease_smap, None, val=abs(node.affinity_change))
                    node.img_style['hz_line_width'] = 1.2
                else:
                    node.img_style['hz_line_color'] = plotting.getgrey()
        if args.queries_to_include is not None and node.name in args.queries_to_include:
            tface = ete3.TextFace(node.name, fsize=3, fgcolor='red')
            node.add_face(tface, column=0)
        rface = ete3.RectFace(width=rfsize, height=rfsize, bgcolor=bgcolor, fgcolor=None)
        rface.opacity = opacity
        node.add_face(rface, column=0)

    affy_label = args.affy_key.replace('_', ' ')
    if args.lb_metric in affy_metrics:
        if affyfo is None:
            add_legend(tstyle, args.lb_metric, lbvals, lb_smap, lbfo, 0, n_entries=4)
        else:
            add_legend(tstyle, args.lb_metric, lbvals, None, lbfo, 0, n_entries=4)
            add_legend(tstyle, affy_label, [a for a in affyvals if a is not None], affy_smap, affyfo, 3)
    elif args.lb_metric in delta_affy_metrics:
        add_legend(tstyle, args.lb_metric, lbvals, lb_smap, lbfo, 0, reverse_log=args.log_lbr)
        if affyfo is not None:
            add_legend(tstyle, '%s decrease' % affy_label, [abs(v) for v in delta_affyvals if v < 0], delta_affy_decrease_smap, affyfo, 3, add_sign='-', no_opacity=True)
            add_legend(tstyle, '%s increase' % affy_label, [v for v in delta_affyvals if v > 0], delta_affy_increase_smap, affyfo, 6, add_sign='+', no_opacity=True)

# ----------------------------------------------------------------------------------------
def plot_trees(args):
    treefo = read_input(args)

    treestr = treefo['treestr']
    if len(treestr.split()) == 2 and treestr.split()[0] in ['[&U]', '[&R']:  # dumbest #$!#$#ing format in the goddamn world (ete barfs on other programs' rooting information)
        treestr = treefo['treestr'].split()[1]
    etree = ete3.Tree(treestr, format=1)  # , quoted_node_names=True)

    tstyle = ete3.TreeStyle()
    tstyle.mode = args.tree_style[0]
    # tstyle.show_scale = False
    tstyle.scale_length = 1. / treeutils.typical_bcr_seq_len
    # tstyle.show_branch_length = True
    # tstyle.complete_branch_lines_when_necessary = True

    if args.metafo is not None:
        set_meta_styles(args, etree, tstyle)

    # print '      %s' % args.outfname
    tstyle.show_leaf_name = False
    etree.render(args.outfname, tree_style=tstyle)

# ----------------------------------------------------------------------------------------
affy_metrics = ['lbi', 'cons-dist-aa', 'cons-dist-nuc', 'shm']  # it would be nice to instead use the info at the top of treeutils/lbplotting
delta_affy_metrics = ['lbr']
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--treefname', required=True)
parser.add_argument('--outfname', required=True)
parser.add_argument('--lb-metric', default='lbi', choices=affy_metrics+delta_affy_metrics)
parser.add_argument('--affy-key', default='affinity', choices=['affinity', 'relative_affinity'])
# parser.add_argument('--lb-tau', required=True, type=float)
parser.add_argument('--metafname')
parser.add_argument('--queries-to-include')
parser.add_argument('--tree-style', default='rectangular', choices=['rectangular', 'circular'])
parser.add_argument('--partis-dir', default=utils.get_partis_dir(), help='path to main partis install dir')
parser.add_argument('--log-lbr', action='store_true')
args = parser.parse_args()

sys.path.insert(1, args.partis_dir + '/python')
try:
    import utils
    import treeutils
    import glutils
    import plotting
except ImportError as e:
    print e
    raise Exception('couldn\'t import from main partis dir \'%s\' (set with --partis-dir)' % args.partis_dir)

args.queries_to_include = utils.get_arg_list(args.queries_to_include)
args.metafo = None
if args.metafname is not None:
    with open(args.metafname) as metafile:
        args.metafo = yaml.load(metafile, Loader=yaml.Loader)

plot_trees(args)
