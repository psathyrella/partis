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
try:
    import ete3
except ImportError:
    raise Exception('couldn\'t find the ete3 module. Either:\n          - it isn\'t installed (use instructions at http://etetoolkit.org/download/) or\n          - $PATH needs modifying (typically with the command \'% export PATH=~/anaconda_ete/bin:$PATH\')')

# ----------------------------------------------------------------------------------------
def getgrey(gtype='medium'):
    if gtype == 'medium':
        return '#929292'
    elif gtype == 'light-medium':
        return '#cdcdcd'
    elif gtype == 'light':
        return '#d3d3d3'
    elif gtype == 'white':
        return '#ffffff'
    else:
        assert False

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

listcolors = [getgrey('medium') for _ in range(10)]
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


# ----------------------------------------------------------------------------------------
def read_input(args):
    with open(args.treefname) as treefile:
        treestr = treefile.read().strip()
    treestr = treestr.replace('[&R] ', '').replace('\'', '')

    return {'treestr' : treestr}

# ----------------------------------------------------------------------------------------
def get_color(smap, info, key, val=None):
    if val is None:
        assert key is not None
        if key not in info:
            return getgrey()
        val = info[key]
    rgb_code = smap.to_rgba(val)[:3]
    return plotting.rgb_to_hex(rgb_code)

# ----------------------------------------------------------------------------------------
min_size = 1.5
max_size = 10
opacity = 0.65
lb_metric = 'lbi'
fsize = 7

# ----------------------------------------------------------------------------------------
def get_size(vmin, vmax, val):
    return min_size + (val - vmin) * (max_size - min_size) / (vmax - vmin)

# ----------------------------------------------------------------------------------------
def draw_tree(args, treefo):
    etree = ete3.Tree(treefo['treestr'], format=1)

    tstyle = ete3.TreeStyle()
    # tstyle.show_scale = False
    # tstyle.scale_length = args.lb_tau
    # tstyle.show_branch_length = True
    # tstyle.complete_branch_lines_when_necessary = True

    if args.metafo is not None:
        lbfo = args.metafo[lb_metric]
        lbvals = lbfo.values()
        lb_min, lb_max = min(lbvals), max(lbvals)
        # lb_smap = plotting.get_normalized_scalar_map(lbvals, None)
        
        if 'affinity' in args.metafo:
            affyfo = args.metafo['affinity']
            affyvals = affyfo.values()
            affy_min, affy_max = min(affyvals), max(affyvals)
            affy_smap = plotting.get_normalized_scalar_map(affyvals, 'viridis')

        for node in etree.traverse():
            node.img_style['size'] = 0
            if node.name in lbfo:
                rfsize = get_size(lb_min, lb_max, lbfo[node.name])
                if affyfo is not None:
                    if node.name in affyfo:
                        bgcolor = get_color(affy_smap, affyfo, node.name)
                    else:
                        bgcolor = getgrey()
                rface = ete3.RectFace(width=rfsize, height=rfsize, bgcolor=bgcolor, fgcolor=None)
                rface.opacity = opacity
                node.add_face(rface, column=0)

        # lb legend
        tstyle.legend.add_face(ete3.TextFace(lb_metric + '   ', fsize=fsize), column=0)
        tstyle.legend.add_face(ete3.TextFace('', fsize=fsize), column=0)
        tstyle.legend.add_face(ete3.TextFace('', fsize=fsize), column=0)
        middle_val = affy_min + (affy_max - affy_min) / 2.
        tstyle.legend.add_face(ete3.RectFace(min_size, min_size, bgcolor=get_color(affy_smap, affyfo, key=None, val=middle_val), fgcolor=None), column=1)
        tstyle.legend.add_face(ete3.RectFace(max_size, max_size, bgcolor=get_color(affy_smap, affyfo, key=None, val=middle_val), fgcolor=None), column=1)
        tstyle.legend.add_face(ete3.TextFace('  %.4f' % lb_min, fsize=fsize), column=2)
        tstyle.legend.add_face(ete3.TextFace('  %.4f' % lb_max, fsize=fsize), column=2)

        # affy legend
        tstyle.legend.add_face(ete3.TextFace('   affinity ', fsize=fsize), column=3)
        delta_affy = (affy_max - affy_min) / 5.
        affy_val_list = list(numpy.arange(affy_min, affy_max + utils.eps, delta_affy))  # first value is exactly <affy_min>, last value is exactly <affy_max>
        affy_key_list = [None for _ in affy_val_list]
        affy_val_list += [None]
        affy_key_list += ['missing!']  # doesn't matter what the last one is as long as it isn't in <affyfo>
        for aval, akey in zip(affy_val_list, affy_key_list):
            tstyle.legend.add_face(ete3.TextFace('', fsize=fsize), column=3)
            middle_size = min_size + (max_size - min_size) / 2.
            rface = ete3.RectFace(middle_size, middle_size, bgcolor=get_color(affy_smap, affyfo, key=akey, val=aval), fgcolor=None)
            rface.opacity = opacity
            tstyle.legend.add_face(rface, column=4)
            tstyle.legend.add_face(ete3.TextFace(('  %.4f' % aval) if akey is None else '  missing', fsize=fsize), column=5)

    suffix = '.svg'
    imagefname = args.plotdir + '/' + args.plotname + suffix
    print '      %s' % imagefname
    tstyle.show_leaf_name = False
    etree.render(imagefname, tree_style=tstyle)

# ----------------------------------------------------------------------------------------
def plot_trees(args):
    treefo = read_input(args)
    draw_tree(args, treefo)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--treefname', required=True)
parser.add_argument('--plotdir', required=True)
parser.add_argument('--lb-tau', required=True, type=float)
parser.add_argument('--plotname', default='test')
parser.add_argument('--metafname')
parser.add_argument('--partis-dir', default=os.getcwd(), help='path to main partis install dir')
args = parser.parse_args()

sys.path.insert(1, args.partis_dir + '/python')
try:
    import utils
    import glutils
    import plotting
except ImportError as e:
    print e
    raise Exception('couldn\'t import from main partis dir \'%s\' (set with --partis-dir)' % args.partis_dir)

args.metafo = None
if args.metafname is not None:
    with open(args.metafname) as metafile:
        args.metafo = yaml.load(metafile)

plot_trees(args)
