#!/usr/bin/env python
from distutils.version import StrictVersion
import argparse
import copy
import time
import random
import sys
import subprocess
import multiprocessing
import numpy
import scipy
if StrictVersion(scipy.__version__) < StrictVersion('0.17.0'):
    raise RuntimeError("scipy version 0.17 or later is required (found version %s)." % scipy.__version__)
import colored_traceback.always
import os
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils

# ----------------------------------------------------------------------------------------
def trim_and_remove_genes(region, gene, seq, glfo, template_glfo, debug=False):
    nearest_template_gene = glutils.find_nearest_gene_using_names(template_glfo, gene)
    nearest_template_seq = template_glfo['seqs'][region][nearest_template_gene]
    # extra_bases = glfo['cyst-positions'][gene] - template_glfo['cyst-positions'][nearest_template_gene]  # not right if there's some internal gaps in the alignment
    aligned_nearest_template_seq, aligned_seq = utils.align_seqs(nearest_template_seq, seq)

    if debug:
        print '    %s' % utils.color_gene(gene)
        utils.color_mutants(aligned_nearest_template_seq, aligned_seq, print_result=True, ref_label='template ', extra_str='       ')

    if aligned_seq[0] not in utils.gap_chars and aligned_nearest_template_seq[0] not in utils.gap_chars:
        if debug:
            print '      ok'
    elif aligned_seq[0] in utils.gap_chars:
        if debug:
            print '      %s, removing' % utils.color('red', 'too small')
        glutils.remove_gene(glfo, gene)
    else:
        if debug:
            print '        extra bases %s' % utils.color_gene(gene)
        extra_bases = len(aligned_nearest_template_seq) - len(aligned_nearest_template_seq.lstrip('-'))
        seq = seq[extra_bases:]
        glfo['seqs'][region][gene] = seq
        glfo['cyst-positions'][gene] -= extra_bases
        if debug:
            print '          removed %d bases' % extra_bases
        # utils.color_mutants(nearest_template_seq, seq, print_result=True, ref_label='template ', align=True, extra_str='            ')
        assert utils.codon_unmutated('cyst', glfo['seqs'][region][gene], glfo['cyst-positions'][gene], debug=True)

# ----------------------------------------------------------------------------------------
def read_ramesh_file(fname, outdir, debug=False):
    seqfos = utils.read_fastx(fname)
    glseqs = {l : {r : {} for r in utils.loci[l]} for l in utils.loci if 'ig' in l}
    for sfo in seqfos:
        if os.path.basename(fname) == 'coding.fa':
            meta = [x.strip('[]').split('=') for x in sfo['infostrs']]
            mdict = {m[0] : m[1] for m in meta if len(m) == 2}
            if 'gene' not in mdict:
                print 'no gene for %s' % sfo['infostrs']
                continue
            gene = mdict['gene']
        else:
            mdict = {}
            gene = sfo['name']
        if debug:
            print gene
        if utils.is_constant_gene(gene):
            if debug:
                print '  constant'
            continue
        region = utils.get_region(gene)
        utils.split_gene(gene)
        # if 'partial' in mdict:
        #     gene += '_partial_%s' % mdict['partial'].replace('\'', '').replace(',', '')
        if sfo['seq'] in glseqs[utils.get_locus(gene)][region].values():
            if debug:
                print '  duplicate'
            continue
        glseqs[utils.get_locus(gene)][region][gene] = sfo['seq']

    return glseqs

# ----------------------------------------------------------------------------------------
def parse_ramesh_seqs(glseqs, outdir, debug=False):
    for locus in glseqs:
        glutils.remove_glfo_files(outdir, locus)
        # write to a glfo dir without extra info
        for region in glseqs[locus]:
            fn = glutils.get_fname(outdir, locus, region)
            if not os.path.exists(os.path.dirname(fn)):
                os.makedirs(os.path.dirname(fn))
            with open(fn, 'w') as ofile:
                for gene, seq in glseqs[locus][region].items():
                    ofile.write('>%s\n%s\n' % (gene, seq))

        # figure out extra info
        template_glfo = glutils.read_glfo('data/germlines/macaque', locus)
        glfo = glutils.read_glfo(outdir, locus, template_glfo=template_glfo, remove_bad_genes=True, debug=True)

        # trim non-coding stuff upstream of v (and remove non-full-length ones)
        gene_groups = {}
        for region in ['v']:
            group_labels = sorted(set([utils.gene_family(g) for g in glfo['seqs'][region]]))
            gene_groups[region] = [(glabel, {g : glfo['seqs'][region][g] for g in glfo['seqs'][region] if utils.gene_family(g) == glabel}) for glabel in group_labels]
        for region in [r for r in utils.regions if r in gene_groups]:
            if debug:
                print '%s' % utils.color('reverse_video', utils.color('green', region))
            for group_label, group_seqs in gene_groups[region]:  # ok, this isn't really doing anything any more
                if debug:
                    print '  %s' % utils.color('blue', group_label)
                for gene, seq in group_seqs.items():
                    trim_and_remove_genes(region, gene, seq, glfo, template_glfo, debug=debug)

        # remove any seqs with ambiguous bases
        for region in [r for r in utils.regions if r in glfo['seqs']]:
            for gene, seq in glfo['seqs'][region].items():
                if utils.ambig_frac(seq) > 0.:
                    if debug:
                        print '   %d ambiguous bases: %s' % (len(seq) * utils.ambig_frac(seq), utils.color_gene(gene))
                    glutils.remove_gene(glfo, gene)

        glutils.print_glfo(glfo)

        # write final result
        glutils.write_glfo(outdir, glfo, debug=True)

fname = 'macaque/ramesh-v1/coding.fa'
outdir = 'macaque/ramesh-cleaned'
# parse_ramesh_seqs(read_ramesh_file(fname, outdir), outdir, debug=True)
# sys.exit()

# ----------------------------------------------------------------------------------------
for locus in ['igh', 'igk', 'igl']:
    # glfo = glutils.read_glfo('macaque/imgt-downloaded', locus, debug=True)
    # glutils.write_glfo('macaque/imgt-parsed', glfo, debug=True)
    # glfo = glutils.read_glfo('data/germlines/macaque', locus, debug=True)
    glfo = glutils.read_glfo(outdir, locus, debug=True)
    glutils.print_glfo(glfo)
sys.exit()

# ----------------------------------------------------------------------------------------
for locus in ['igh']:  # ['igh', 'igk', 'igl']:
    glfo = glutils.read_glfo('macaque/ramesh', locus, debug=True)
    # glutils.write_glfo('macaque/imgt-parsed', glfo, debug=True)
    # glfo = glutils.read_glfo('data/germlines/macaque', locus, debug=True)
    # glutils.print_glfo(glfo)
