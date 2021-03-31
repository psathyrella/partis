#!/usr/bin/env python
import argparse
import sys
import os
import copy
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')
import utils
import glutils
import collections
import colored_traceback.always

parser = argparse.ArgumentParser()
parser.add_argument('gldir1')
parser.add_argument('gldir2')
parser.add_argument('--names', default='+gl-1:+gl-2', help='colon-separated list of length 2 with labels for gldir1 and gldir2, which will be appended to each gene name in the ascii output')
parser.add_argument('--locus', default='igh')
args = parser.parse_args()
args.names = utils.get_arg_list(args.names)

glfos = []
for name, gldir in zip(args.names, [args.gldir1, args.gldir2]):
    print '%s:' % utils.color('yellow', name)
    glfos.append(glutils.read_glfo(gldir, args.locus, debug=True))

for region in [r for r in utils.regions if r in glfos[0]['seqs']]:
    aset, bset = [set(g['seqs'][region]) for g in glfos]

    tmpfo = glutils.get_empty_glfo(args.locus)  # make a new glfo that will only have non-shared genes
    for glabel, gset, gfo in zip(args.names, [aset - bset, bset - aset], glfos):  # <gset> is the genes that're only in <glabel>
        for ogene in gset:
            glutils.add_new_allele(tmpfo, {'gene' : '+'.join([ogene, glabel]), 'seq' : gfo['seqs'][region][ogene], 'cpos' : utils.cdn_pos(gfo, region, ogene)}, use_template_for_codon_info=False)

    # eh, maybe this doesn't really add anything?
    # # add the nearest genes that they both have for comparison NOTE this gives one comparison gene for *each* gene, so usually you get a bunch of comparison/'both' genes in each block in the ascii output
    # for bgene in aset & bset:
    #     _, nearest_gene, _ = glutils.find_nearest_gene_with_same_cpos(glfos[0], glfos[0]['seqs'][region][bgene], new_cpos=utils.cdn_pos(glfos[0], region, bgene))  # i think it doesn't matter which glfo we get it from, so arbitrarily choose the first one
    #     glutils.add_new_allele(tmpfo, {'gene' : '+'.join([nearest_gene, 'both']), 'seq' : glfos[0]['seqs'][region][nearest_gene], 'cpos' : utils.cdn_pos(glfos[0], region, bgene)}, use_template_for_codon_info=False)

    print '%s: only in:\n      %12s: %2d  %s\n      %12s: %2d  %s' % (utils.color('green', region), args.names[0], len(aset - bset), utils.color_genes(sorted(aset - bset)), args.names[1], len(bset - aset), utils.color_genes(sorted(bset - aset)))
    if len(tmpfo['seqs'][region]) > 0:
        print ' comparing to nearest genes that were in both (labeled \'both\'):'
        glutils.print_glfo(tmpfo, only_region=region)
