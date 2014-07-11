#!/usr/bin/env python

import os
import sys
import utils
from hmmwriter import HmmWriter

only_genes = ''
try:
    ionly = sys.argv.index('--only_genes')
    only_genes = sys.argv[ionly+1]
except:
    pass

outfname = ''
try:
    ioutfname = sys.argv.index('--outfname')
    outfname = sys.argv[ioutfname+1]
except:
    pass

# only_genes = 'IGHV3-64*04:IGHV1-18*01:IGHV3-23*04:IGHV3-72*01:IGHV5-51*01:IGHD4-23*01:IGHD3-10*01:IGHD4-17*01:IGHD6-19*01:IGHD3-22*01:IGHJ4*02_F:IGHJ5*02_F:IGHJ6*02_F:IGHJ3*02_F:IGHJ2*01_F'.split(':')
germline_seqs = utils.read_germlines('/home/dralph/Dropbox/work/recombinator')
human = 'A'
naivety = 'M'

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
        if only_genes != '' and gene_name not in only_genes:
            continue
        print '  %d / %d (%s)' % (igene, len(germline_seqs[region]), gene_name)
        igene += 1
        writer = HmmWriter('/home/dralph/Dropbox/work/recombinator/data/human-beings/' + human + '/' + naivety, 'bcell', gene_name, naivety, germline_seqs[region][gene_name])
        writer.write(outfname)
