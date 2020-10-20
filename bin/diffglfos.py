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
    tmp_glfo = copy.deepcopy(glfos[0])  # make a new glfo that will only have non-shared genes

    # first remove any that are in both
    for gene, seq in tmp_glfo['seqs'][region].items():
        if gene not in glfos[1]['seqs'][region]:  # keep it, under a new name, if it's not in <glfos[1]>
            glutils.add_new_allele(tmp_glfo, {'gene' : '+'.join([gene, args.names[0]]), 'seq' : seq, 'template-gene' : gene})  # add an extra str to the name so we know which one it came from
        glutils.remove_gene(tmp_glfo, gene)

    # then add any that are only in the second one
    for gene, seq in glfos[1]['seqs'][region].items():
        if gene not in glfos[0]['seqs'][region]:
            cpos = glfos[1][utils.cdn(glfos[1], region)+'-positions'][gene] if utils.cdn(glfos[1], region) is not None else None
            glutils.add_new_allele(tmp_glfo, {'gene' : '+'.join([gene, args.names[1]]), 'seq' : seq, 'cpos' : cpos}, use_template_for_codon_info=False)  # can't use template cause we might've deleted it in the first loop

    # then add the nearest genes that they both have for comparison
    for gene, seq in tmp_glfo['seqs'][region].items():
        _, nearest_gene, _ = glutils.find_nearest_gene_with_same_cpos(glfos[0], seq, new_cpos=utils.cdn_pos(tmp_glfo, region, gene))  # i think it doesn't matter which glfo we get it from, so arbitrarily choose the first one
        glutils.add_new_allele(tmp_glfo, {'gene' : '+'.join([nearest_gene, 'both']), 'seq' : glfos[0]['seqs'][region][nearest_gene], 'template-gene' : gene})

    if len(tmp_glfo['seqs'][region]) > 0:
        print ' comparing to nearest genes that were in both (labeled \'both\'):'
        glutils.print_glfo(tmp_glfo, only_region=region)
