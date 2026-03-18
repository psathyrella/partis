#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import csv
import sys
import os
import copy
import collections
import colored_traceback.always

import partis.utils as utils
import partis.glutils as glutils

parser = argparse.ArgumentParser(
    description='Compare two germline sets. Paths can be germline-sets dirs or '
    'parameter dirs (auto-detected by presence of germline-sets/ subdir and '
    'v_gene-probs.csv). When parameter dirs are given, also computes '
    'sequence-aware usage distance (EMD).'
)
parser.add_argument('dir1', help='germline-sets dir or parameter dir')
parser.add_argument('dir2', help='germline-sets dir or parameter dir')
parser.add_argument('--names', default='a:b', help='colon-separated list of length 2 with labels for dir1 and dir2')
parser.add_argument('--colors', default='blue:red')
parser.add_argument('--locus', default='igh')
parser.add_argument('--debug', type=int, default=1, help='debug level for usage comparison (0=silent, 1=table, 2=+alignment)')
args = parser.parse_args()
args.names = utils.get_arg_list(args.names)
args.colors = utils.get_arg_list(args.colors)

# ----------------------------------------------------------------------------------------
def read_gene_probs_csv(csv_path):
    """Read a partis {region}_gene-probs.csv file (columns: {region}_gene, count)."""
    usage = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            # column name is like v_gene, d_gene, j_gene
            gene = [v for k, v in row.items() if k.endswith('_gene')][0]
            count = float(row['count'])
            if count > 0:
                usage[gene] = count
    return usage

# ----------------------------------------------------------------------------------------
def resolve_dir(path):
    """Return (gldir, is_param_dir). Auto-detects parameter dirs."""
    gs_subdir = os.path.join(path, 'germline-sets')
    probs_file = os.path.join(path, 'v_gene-probs.csv')
    if os.path.isdir(gs_subdir) and os.path.isfile(probs_file):
        return gs_subdir, True
    return path, False

# ----------------------------------------------------------------------------------------
gldirs, param_dirs = [], []
for name, d, color in zip(args.names, [args.dir1, args.dir2], args.colors):
    gldir, is_param_dir = resolve_dir(d)
    gldirs.append(gldir)
    param_dirs.append(d if is_param_dir else None)
    label = '%s (parameter dir)' % d if is_param_dir else gldir
    print('%s: %s' % (utils.color(color, name), label))

glfos = [glutils.read_glfo(gldir, args.locus, debug=True) for gldir in gldirs]

# glutils.compare_glfos(glfos, args.names, args.locus)

if all(param_dirs):
    print('\n%s' % utils.color('green', 'sequence-aware usage comparison (EMD)'))
    for region in [r for r in utils.regions if r in glfos[0]['seqs']]:
        probs_file = '%s_gene-probs.csv' % region
        region_usage = []
        for pdir in param_dirs:
            csv_path = os.path.join(pdir, probs_file)
            if os.path.isfile(csv_path):
                region_usage.append(read_gene_probs_csv(csv_path))
            else:
                region_usage.append({})
        if not any(region_usage):
            continue
        glutils.compare_germline_usage(glfos, region_usage, region, names=args.names, debug=args.debug)
