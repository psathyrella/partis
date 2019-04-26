#!/usr/bin/env python
import os
import sys
import argparse
import colored_traceback.always
import yaml

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import treeutils

class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter)
parser.add_argument('infname', help='input newick tree file')
parser.add_argument('outfname', help='output yaml file')
parser.add_argument('--input-metafname', help='optional input yaml meta file with multiplicity information. Example: yaml.dump({n.taxon.label : {\'multiplicity\' : 1} for n in dtree.postorder_node_iter()}, tfile, width=150)')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()


if not os.path.exists(args.infname):
    raise Exception('infname %s does not exist' % args.infname)
with open(args.infname) as infile:
    dtree = treeutils.get_dendro_tree(treefname=args.infname, debug=args.debug)

input_metafo = None
if args.input_metafname is not None:
    if not os.path.exists(args.input_metafname):
        raise Exception('--input-metafname %s does not exist' % args.input_metafname)
    with open(args.input_metafname) as imfile:
        input_metafo = yaml.load(imfile, Loader=yaml.Loader)

lbvals = treeutils.calculate_lb_values(dtree, treeutils.default_lb_tau, treeutils.default_lbr_tau_factor * treeutils.default_lb_tau, input_metafo=input_metafo, use_multiplicities=(input_metafo is not None), debug=args.debug)

with open(args.outfname, 'w') as outfile:
    yaml.dump(lbvals, outfile, width=200)
