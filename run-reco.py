#!/usr/bin/env python
from Recombinator import Recombinator

#----------------------------------------------------------------------------------------
reco = Recombinator('data/human-beings/C/N/probs.csv.bz2')
#reco = Recombinator('data/hack-probs.csv.bz2')  # hacked this file by hand to use one v, j, and cdr3 length but vary the d gene.
for _ in range(100):
    print 'combining'
#    raw_reco_seq = reco.combine('out.csv')
    raw_reco_seq = reco.combine()
