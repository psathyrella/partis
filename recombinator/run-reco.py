#!/usr/bin/env python
import os
import sys
from recombinator import Recombinator
import utils

#----------------------------------------------------------------------------------------
#assert len(sys.argv) == 2
#human = sys.argv[1]  #'C'
naivety = 'N'
only_genes=''  #'IGHV3-64*04:IGHV1-18*01:IGHV3-23*04:IGHV3-72*01:IGHV5-51*01:IGHD4-23*01:IGHD3-10*01:IGHD4-17*01:IGHD6-19*01:IGHD3-22*01:IGHJ4*02_F:IGHJ5*02_F:IGHJ6*02_F:IGHJ3*02_F:IGHJ2*01_F',
for human in utils.humans:
    reco = Recombinator('../parameters', '../data', human, naivety, only_genes=only_genes, total_length_from_right=130, apply_shm=(naivety=='M'))
    for _ in range(100):
        outdir = 'output/' + human + '/' + naivety
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outfname = outdir + '/simu.csv'
        success = False
        while not success:  # returns False on failure, so keep trying (failure usually means we chose inconsistent cdr3 length and gene choices, or something similar)
            success = reco.combine(outfname, 'append')
    # sys.exit()
