#!/usr/bin/env python

import os
import utils
from hmmwriter import HmmWriter

n_max_versions = 0  # only look at the first n gene versions (speeds things up for testing)
only_genes = 'IGHV3-64*04:IGHV1-18*01:IGHV3-23*04:IGHV3-72*01:IGHV5-51*01:IGHD4-23*01:IGHD3-10*01:IGHD4-17*01:IGHD6-19*01:IGHD3-22*01:IGHJ4*02_F:IGHJ5*02_F:IGHJ6*02_F:IGHJ3*02_F:IGHJ2*01_F'.split(':')
germline_seqs = utils.read_germlines('/home/dralph/Dropbox/work/recombinator')
human = 'A'
naivety = 'N'

for region in utils.regions:
    outdir = 'bcell/' + region
    if os.path.exists(outdir):
        for hmmfile in os.listdir(outdir):
            if hmmfile.endswith(".hmm"):
                os.remove(outdir + "/" + hmmfile)
    else:
        os.makedirs(outdir)
    
    igene = 0
    for gene_name in germline_seqs[region]:
        if n_max_versions != 0 and igene >= n_max_versions:
            print 'breaking after %d gene versions' % n_max_versions
            break
        if only_genes != '' and gene_name not in only_genes:
            continue
        print '  %d / %d (%s)' % (igene, len(germline_seqs[region]), gene_name)
        igene += 1
        writer = HmmWriter('/home/dralph/Dropbox/work/recombinator/data/human-beings/' + human + '/' + naivety, 'bcell', region, gene_name, naivety, germline_seqs[region][gene_name])
        writer.write()
