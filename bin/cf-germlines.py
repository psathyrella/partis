#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import sys
import os
import copy
import collections
import colored_traceback.always

import partis.utils as utils
import partis.glutils as glutils

parser = argparse.ArgumentParser()
parser.add_argument('gldir1')
parser.add_argument('gldir2')
parser.add_argument('--names', default='+gl-1:+gl-2', help='colon-separated list of length 2 with labels for gldir1 and gldir2, which will be appended to each gene name in the ascii output')
parser.add_argument('--locus', default='igh')
args = parser.parse_args()
args.names = utils.get_arg_list(args.names)

# ----------------------------------------------------------------------------------------
def clrname(name):
    return utils.color('blue', name)

# ----------------------------------------------------------------------------------------
glfos = []
for name, gldir in zip(args.names, [args.gldir1, args.gldir2]):
    print('%s:' % clrname(name))
    glfos.append(glutils.read_glfo(gldir, args.locus, debug=True))

for region in [r for r in utils.regions if r in glfos[0]['seqs']]:
    aseqs, bseqs = [{s : n for n, s in g['seqs'][region].items()} for g in glfos]  # dict of names keyed by seqs
    a_only_seqs, b_only_seqs = set(aseqs) - set(bseqs), set(bseqs) - set(aseqs)

    print('%s' % utils.color('green', region))

    common_seqs = set(aseqs) & set(bseqs)
    common_name_seqs = [aseqs[s] for s in common_seqs if aseqs[s]==bseqs[s]]
    print('    %3d seqs in common with same name: %s' % (len(common_name_seqs), utils.color_genes(sorted(common_name_seqs))))
    dnamed_seqs = [(aseqs[s], bseqs[s]) for s in common_seqs if aseqs[s] != bseqs[s]]
    if len(dnamed_seqs) > 0:
        print('      %s %d common seq%s with different names: %s' % (utils.wrnstr(), len(dnamed_seqs), utils.plural(len(dnamed_seqs)), ',  '.join(utils.color_genes([an,bn]) for an, bn in dnamed_seqs)))
    print('    only in:\n      %12s: %3d  %s\n      %12s: %3d  %s' % (clrname(args.names[0]), len(a_only_seqs), utils.color_genes(sorted(aseqs[s] for s in a_only_seqs)),
                                                                      clrname(args.names[1]), len(b_only_seqs), utils.color_genes(sorted(bseqs[s] for s in b_only_seqs))))

    tmpfo = glutils.get_empty_glfo(args.locus)  # make a new glfo that will only have non-shared genes
    for gname, oname, only_seqs, allseqs, ogfo in zip(args.names, reversed(args.names), [a_only_seqs, b_only_seqs], [aseqs, bseqs], reversed(glfos)):  # <gset> is the genes that're only in <gname>
        print('  finding nearest seq in %s for %d seqs only in %s' % (clrname(oname), len(only_seqs), clrname(gname)))
        for oseq in only_seqs:
            glutils.find_nearest_gene_in_glfo(ogfo, oseq, new_name=allseqs[oseq], region=region, debug=True)
