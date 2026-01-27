#!/usr/bin/env python3
# has to be its own script, since ete3 requires its own god damn python version, installed in a separated directory
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import time
import yaml
import itertools
import glob
import argparse
import copy
import random
import os
import subprocess
import sys
import colored_traceback.always
from collections import OrderedDict
import numpy
import math
import re
from io import open
import ete3
from pathlib import Path

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
opacity = 0.65
node_fsize = 7

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
    if vmin == vmax:
        return args.min_face_size
    if args.use_node_area:
        val = math.sqrt(val)
    rfsize = args.min_face_size + (val - vmin) * (args.max_face_size - args.min_face_size) / float(vmax - vmin)
    return rfsize

# ----------------------------------------------------------------------------------------
def add_legend(tstyle, varname, all_vals, smap, info, start_column, add_missing=False, add_sign=None, reverse_log=False, n_entries=5, fsize=4, no_opacity=False):  # NOTE very similar to add_smap_legend() in plot_2d_scatter() in python/lbplotting.py
    if len(all_vals) == 0:
        all_vals = [-1, 1]  # um, maybe this is ok?
        # return  # NOTE you *cannot* return here, since if we don't actually add the stuff then it later (when rendering) crashes with a key error deep within ete due to the <start_column> being wrong/inconsistent
    assert add_sign in [None, '-', '+']
    tstyle.legend.add_face(ete3.TextFace('   %s ' % varname, fsize=fsize), column=start_column)
    min_val, max_val = get_scale_min(varname, all_vals), max(all_vals)
    if min_val == max_val:
        min_val, max_val = plotting.expand_bounds([min_val, max_val], only_down=True)  # <only_down> is for affinity increase scale: expand downward if only one value so the one value shows up as dark red (rather than super light red)
    val_list = plotting.get_leg_entries(n_entries=n_entries, min_val=min_val, max_val=max_val)
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
        if reverse_log:
            val = math.exp(val)
        def vstr():
            if varname == 'cons-dist-aa': return '%.1f' % val
            elif 'affinity' in varname: return '%s' % utils.round_to_n_digits(val, 2)
            else: return '%.2f' % val
        if key is None:
            tfstr = '  %s%s' % (utils.non_none([add_sign, '']), vstr())
        else:
            ftstr = '  missing'
        tstyle.legend.add_face(ete3.TextFace(tfstr, fsize=fsize), column=start_column + 2)

# ----------------------------------------------------------------------------------------
def label_node(lnode, root_node):
    # ----------------------------------------------------------------------------------------
    def meta_emph(tname):
        if args.meta_info_to_emphasize is not None:
            key, val = list(args.meta_info_to_emphasize.items())[0]
            if key in args.metafo and tname in args.metafo[key] and utils.meta_info_equal(key, val, args.metafo[key][tname], formats=args.meta_emph_formats):
                return True
        return False
    # ----------------------------------------------------------------------------------------
    def use_node_name(tname, tnode=None):
        if args.label_all_nodes:
            return True
        if tnode is not None and args.label_leaf_nodes and tnode.is_leaf():
            return True
        if args.queries_to_include is not None and tname in args.queries_to_include:
            return True
        if tnode is not None and args.label_root_node and tnode is root_node:
            return True
        if meta_emph(tname):
            return True
        if args.uid_translations is not None and tname in args.uid_translations:
            return True
        if args.node_label_regex is not None and len(re.findall(args.node_label_regex, tname)) > 0:
            return True
        return False
    # ----------------------------------------------------------------------------------------
    def split_line(lstr):  # split into two rows if more than 3 entries
        if lstr.count(',') < 3:
            return lstr
        blist = lstr.split(', ')
        return '%s\n%s' % (', '.join(blist[:len(blist)//2]), ', '.join(blist[len(blist)//2:]))
    # ----------------------------------------------------------------------------------------
    def get_nlabel(tname):
        if not use_node_name(tname):
            return None
        nlabel = tname
        if args.uid_translations is not None and nlabel in args.uid_translations:
            nlabel = args.uid_translations[nlabel]
        if args.node_label_regex is not None:
            mstrs = re.findall(args.node_label_regex, nlabel)
            if len(mstrs) > 0:
                nlabel = '+'.join(mstrs)
        return nlabel
    # ----------------------------------------------------------------------------------------
    def get_ncolor(tname):
        return 'red' if meta_emph(tname) or args.meta_info_to_emphasize is None else 'black'
    # ----------------------------------------------------------------------------------------
    use_name = use_node_name(lnode.name, tnode=lnode)
    if 'duplicates' in args.metafo and lnode.name in args.metafo['duplicates']:
        use_name |= any(use_node_name(n) for n in args.metafo['duplicates'][lnode.name])
    if use_name:
        nlabels, ncolors = [[tfn(lnode.name)] for tfn in (get_nlabel, get_ncolor)]
        if 'duplicates' in args.metafo and lnode.name in args.metafo['duplicates']:
            nlabels += [get_nlabel(n) for n in args.metafo['duplicates'][lnode.name]]
            ncolors += [get_ncolor(n) for n in args.metafo['duplicates'][lnode.name]]
        nlabel = ', '.join(sorted(set(l for l in nlabels if l is not None)))
        ncolor = 'red' if 'red' in ncolors else 'black'
    else:
        if 'labels' in args.metafo:
            nlabel, ncolor = '', 'red'
        else:
            return
    if 'labels' in args.metafo:
        mlabel = args.metafo['labels'].get(lnode.name, '')
        blabels, tlabels, tcolors, bcolors = ['', ''], ['', ''], ['black' for _ in range(2)], ['black' for _ in range(2)]
        label_list = mlabel.split('\n')
        if 'h:' in mlabel or 'l:' in mlabel:
            for lstr in label_list:
                if 'nuc' in lstr and 'aa' in lstr:
                    assert lstr.count(',') == 1  # e.g. '3 nuc, 1 aa'
                    tlabels[0], blabels[0] = lstr.split(',')
                elif lstr.find('h:') == 0:
                    tlabels[1], tcolors[1] = ' '+lstr, 'blue'
                elif lstr.find('l:') == 0:
                    blabels[1], bcolors[1] = ' '+lstr, 'green'
                else:
                    raise Exception('couldn\'t parse \'%s\'' % mlabel)
        elif '\n' in mlabel:
            tlabels[1], blabels[1] = label_list[0], '\n'.join(label_list[1:])
            blabels[1] = split_line(blabels[1])
        else:
            tlabels[0] = mlabel
        for il, (blab, bcol) in enumerate(zip(blabels, bcolors)):  # blabels are branch bottom labels
            lnode.add_face(ete3.TextFace(blab, fsize=node_fsize, fgcolor=bcol), column=il, position='branch-bottom')
        for il, (tlab, tcol) in enumerate(zip(tlabels, tcolors)):  # branch top labels
            lnode.add_face(ete3.TextFace(tlab, fsize=node_fsize, fgcolor=tcol), column=il, position='branch-top')
    if nlabel != '':
        tface = ete3.TextFace(' '+nlabel, fsize=node_fsize, fgcolor=ncolor)  # <nlabel> is usually the uid/sequence name
        lnode.add_face(tface, column=1) # position='branch-bottom')

# ----------------------------------------------------------------------------------------
def set_lb_styles(args, etree, tstyle):
    # ----------------------------------------------------------------------------------------
    lbfo = args.metafo[args.lb_metric]
    if 'lbr' in args.lb_metric or 'lbf' in args.lb_metric:  # remove zeros + maybe apply log()
        lbfo = {u : (math.log(v) if args.log_lbr else v) for u, v in lbfo.items() if v > 0}
    lbvals = list(lbfo.values())
    if len(lbvals) == 0:
        return
    lb_smap = plotting.get_normalized_scalar_map(lbvals, 'viridis', hard_min=get_scale_min(args.lb_metric, lbvals) if args.lb_metric=='cons-dist-aa' else None)
    lb_min, lb_max = min(lbvals), max(lbvals)

    affyfo = None
    if args.affy_key in args.metafo and set(args.metafo[args.affy_key].values()) != set([None]):
        affyfo = args.metafo[args.affy_key]
        if args.lb_metric in treeutils.affy_metrics:
            affyvals = list(affyfo.values())
            affy_smap = plotting.get_normalized_scalar_map([a for a in affyvals if a is not None], 'viridis')
        elif args.lb_metric in treeutils.daffy_metrics:
            delta_affyvals = set_delta_affinities(etree, affyfo)
            affy_increases = [v for v in delta_affyvals if v > 0]
            if len(set(affy_increases)) == 1:  # if there's only one affinity increase, expand downward so color is dark red for actual observed value
                affy_increases = plotting.expand_bounds([affy_increases[0] for _ in range(2)], only_down=True)
            delta_affy_increase_smap = plotting.get_normalized_scalar_map(affy_increases, 'Reds', remove_top_end=True) if len(delta_affyvals) > 0 else None
            delta_affy_decrease_smap = plotting.get_normalized_scalar_map([abs(v) for v in delta_affyvals if v < 0], 'Blues', remove_top_end=True) if len(delta_affyvals) > 0 else None
        else:
            assert False

    for node in etree.traverse():
        node.img_style['size'] = 0
        rfsize = 0
        bgcolor = plotting.getgrey()
        if args.lb_metric in treeutils.affy_metrics:
            if node.name not in lbfo:  # really shouldn't happen
                print('  %s missing lb info for node \'%s\'' % (utils.color('red', 'warning'), node.name))
                continue
            if affyfo is not None:
                rfsize = get_size(lb_min, lb_max, lbfo[node.name])
                if node.name in affyfo:
                    bgcolor = plotting.get_smap_color(affy_smap, affyfo, key=node.name)
            else:
                rfsize = 5
                bgcolor = plotting.get_smap_color(lb_smap, lbfo, key=node.name)
        elif args.lb_metric in treeutils.daffy_metrics:
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
        label_node(node, etree.get_tree_root())
        rface = ete3.RectFace(width=rfsize, height=rfsize, bgcolor=bgcolor, fgcolor=None)
        rface.opacity = opacity
        node.add_face(rface, column=0)

    affy_label = args.affy_key.replace('_', ' ')
    mleg = lbplotting.mtitlestr('per-seq', args.lb_metric)
    if args.lb_metric in treeutils.affy_metrics:
        if affyfo is None:
            add_legend(tstyle, mleg, lbvals, lb_smap, lbfo, 0, n_entries=4)
        else:
            add_legend(tstyle, mleg, lbvals, None, lbfo, 0, n_entries=4)
            add_legend(tstyle, affy_label, [a for a in affyvals if a is not None], affy_smap, affyfo, 3)
    elif args.lb_metric in treeutils.daffy_metrics:
        add_legend(tstyle, mleg, lbvals, lb_smap, lbfo, 0, reverse_log=args.log_lbr)
        if affyfo is not None:
            add_legend(tstyle, '%s decrease' % affy_label, [abs(v) for v in delta_affyvals if v < 0], delta_affy_decrease_smap, affyfo, 3, add_sign='-', no_opacity=True)
            add_legend(tstyle, '%s increase' % affy_label, [v for v in delta_affyvals if v > 0], delta_affy_increase_smap, affyfo, 6, add_sign='+', no_opacity=True)

# ----------------------------------------------------------------------------------------
def set_meta_styles(args, etree, tstyle):
    all_emph_vals, emph_colors = None, None
    if args.meta_info_key_to_color is not None:
        mvals = args.metafo.get(args.meta_info_key_to_color, {})
        all_emph_vals, emph_colors = plotting.meta_emph_init(args.meta_info_key_to_color, formats=args.meta_emph_formats, all_emph_vals=set(mvals.values()))
        mcolors = {v : c for v, c in emph_colors}
        plotting.make_meta_info_legend(os.path.dirname(args.outfname), utils.getprefix(os.path.basename(args.outfname)), args.meta_info_key_to_color, emph_colors, all_emph_vals, meta_emph_formats=args.meta_emph_formats, alpha=0.6)
    if args.node_size_key is not None:
        nsvals = set(args.metafo[args.node_size_key].values()) - set([None])
        min_nsval, max_nsval = [mfcn(nsvals) for mfcn in [min, max]]
        if args.use_node_area:
            min_nsval, max_nsval = [math.sqrt(v) for v in [min_nsval, max_nsval]]
    if args.branch_color_key is not None:
        bcvals = args.metafo[args.branch_color_key]
        smvals = [v for v in bcvals.values() if v is not None]
        if 'vrc01' in args.branch_color_key:
            smvals = [v for v in smvals if v > 0]
        bc_smap = plotting.get_normalized_scalar_map(smvals, 'Reds', remove_top_end=True)
        add_legend(tstyle, plotting.legends.get(args.branch_color_key, args.branch_color_key), smvals, bc_smap, bcvals, 3, no_opacity=True)
    for node in etree.traverse():
        node.img_style['size'] = 0
        rfsize = 5
        if args.node_size_key is not None:
            rfsize = get_size(min_nsval, max_nsval, args.metafo[args.node_size_key].get(node.name, min_nsval))
        bgcolor = plotting.getgrey()
        if args.meta_info_key_to_color is not None and node.name in mvals:
            bgcolor = mcolors.get(mvals[node.name], bgcolor)

        label_node(node, etree.get_tree_root())
        ftypes = {'rect' : ete3.RectFace, 'circle' : ete3.CircleFace}
        if args.face_type == 'rect':
            rface = ete3.RectFace(width=rfsize, height=rfsize, bgcolor=bgcolor, fgcolor=None)
        elif args.face_type == 'circle':
            rface = ete3.CircleFace(radius=rfsize, color=bgcolor)
        else:
            assert False
        rface.opacity = opacity
        node.add_face(rface, column=0)

        if args.branch_color_key is not None:
            bval = bcvals.get(node.name)
            bcol = plotting.getgrey() if 'vrc01' in args.branch_color_key and bval==0 else plotting.get_smap_color(bc_smap, bcvals, key=node.name)
            node.img_style['hz_line_color'] = bcol
            node.img_style['hz_line_width'] = 1.2

# ----------------------------------------------------------------------------------------
def plot_trees(args):
    treefo = read_input(args)

    treestr = treefo['treestr']
    if len(treestr.split()) == 2 and treestr.split()[0] in ['[&U]', '[&R]']:  # dumbest #$!#$#ing format in the goddamn world (ete barfs on other programs' rooting information)
        treestr = treefo['treestr'].split()[1]
    etree = ete3.Tree(treestr, format=1)  # , quoted_node_names=True)

    tstyle = ete3.TreeStyle()
    tstyle.mode = args.tree_style[0]
    # tstyle.show_scale = False
    if getattr(tstyle, 'scale_length', None) is not None:
        tstyle.scale_length = 1. / treeutils.typical_bcr_seq_len
    # tstyle.show_branch_length = True
    # tstyle.complete_branch_lines_when_necessary = True

    if args.metafo is not None:
        if args.lb_metric is not None:
            set_lb_styles(args, etree, tstyle)
        if args.meta_info_key_to_color is not None or args.meta_info_to_emphasize:
            set_meta_styles(args, etree, tstyle)
    else:
        print('  %s --metafo is not set, so no node formats (e.g. labels) will be set)' % utils.wrnstr())

    # print '      %s' % args.outfname
    tstyle.show_leaf_name = False
    etree.render(args.outfname, tree_style=tstyle)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--treefname', required=True)
parser.add_argument('--outfname', required=True)
parser.add_argument('--lb-metric') #, default='lbi') #, choices=treeutils.affy_metrics+treeutils.daffy_metrics)
parser.add_argument('--affy-key', default='affinity', choices=['affinity', 'relative_affinity'])
# parser.add_argument('--lb-tau', required=True, type=float)
parser.add_argument('--metafname')
parser.add_argument('--queries-to-include')
parser.add_argument('--label-all-nodes', action='store_true')
parser.add_argument('--label-leaf-nodes', action='store_true')
parser.add_argument('--label-root-node', action='store_true')
parser.add_argument('--node-label-regex', help='portion of node label to keep (rest is discarded if regex is found, if no regex label is left unchanged). E.g. \'ig.\' reduces them all to the locus')
parser.add_argument('--tree-style', default='rectangular', choices=['rectangular', 'circular'])
parser.add_argument('--partis-dir', default=str(Path(__file__).parent.parent), help='path to main partis install dir')
parser.add_argument('--log-lbr', action='store_true')
parser.add_argument('--seq-len', type=int)
parser.add_argument('--uid-translations', help='colon-separated list of comma-separated pairs of uid:translated-id pairs')
parser.add_argument('--meta-info-to-emphasize', help='see partis help')
parser.add_argument('--meta-info-key-to-color', help='see partis help')
parser.add_argument('--meta-emph-formats', help='see partis help')
parser.add_argument('--node-size-key', help='annotation key with which to scale the node size')
parser.add_argument('--use-node-area', action='store_true', help='for --node-size-key scale by area rather than edge/radius')
parser.add_argument('--branch-color-key', help='annotation key with which to scale the branch length color')
parser.add_argument('--face-type', choices=['rect', 'circle'], default='rect', help='shape of symbol for each node')
parser.add_argument('--min-face-size', type=float, default=1.5, help='min size for node symbol')
parser.add_argument('--max-face-size', type=float, default=15, help='min size for node symbol')
args = parser.parse_args()
if args.meta_info_key_to_color is None and not args.meta_info_to_emphasize:  # it'd be nice to move this stuff, but whatevs
    print('  note: if neither --meta-info-key-to-color or --meta-info-to-emphasize are set, other style attributes may not be set')

import partis.utils as utils
import partis.treeutils as treeutils
import partis.glutils as glutils
import partis.plotting as plotting
import partis.lbplotting as lbplotting
args.meta_info_to_emphasize = utils.get_arg_list(args.meta_info_to_emphasize, key_val_pairs=True)
args.meta_emph_formats = utils.get_arg_list(args.meta_emph_formats, key_val_pairs=True)
utils.meta_emph_arg_process(args)
args.uid_translations = utils.get_arg_list(args.uid_translations, key_val_pairs=True)

if args.node_label_regex is not None and not args.label_all_nodes and not args.label_leaf_nodes:
    print('  note: --node-label-regex set, but neither --label-all-nodes nor --label-leaf-nodes were set, so they may not actually get labeled')
    # print('  --node-label-regex: turning on --label-all-nodes')
    # args.label_all_nodes = True

args.queries_to_include = utils.get_arg_list(args.queries_to_include)
args.metafo = None
if args.metafname is not None:
    with open(args.metafname) as metafile:
        args.metafo = yaml.load(metafile, Loader=yaml.CLoader)

plot_trees(args)
