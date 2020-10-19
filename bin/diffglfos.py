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

glfo_1, glfo_2 = [glutils.read_glfo(gd, args.locus, debug=True) for gd in [args.gldir1, args.gldir2]]
# glutils.print_glfo(glfo_1)
# glutils.print_glfo(glfo_2)

for region in [r for r in utils.regions if r in glfo_1['seqs']]:
    tmp_glfo = copy.deepcopy(glfo_1)  # make a new glfo that will only have non-shared genes
    for gene, seq in tmp_glfo['seqs'][region].items():
        if gene not in glfo_2['seqs'][region]:  # keep it, under a new name, if it's not in <glfo_2>
            glutils.add_new_allele(tmp_glfo, {'gene' : gene+args.names[0], 'seq' : seq, 'template-gene' : gene})  # add an extra str to the name so we know which one it came from
        glutils.remove_gene(tmp_glfo, gene)
    for gene, seq in glfo_2['seqs'][region].items():
        if gene not in tmp_glfo['seqs'][region]:
            cpos = glfo_2[utils.cdn(glfo_2, region)+'-positions'][gene] if utils.cdn(glfo_2, region) is not None else None
            glutils.add_new_allele(tmp_glfo, {'gene' : gene+args.names[1], 'seq' : seq, 'cpos' : cpos}, use_template_for_codon_info=False)  # can't use template cause we might've deleted it in the first loop
    glutils.print_glfo(tmp_glfo, only_region=region)
