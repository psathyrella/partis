#!/usr/bin/env python
from recombinator import Recombinator

#----------------------------------------------------------------------------------------
reco = Recombinator('data/human-beings','C','M')
#reco = Recombinator('data/hack-probs.csv.bz2')  # hacked this file by hand to use one v, j, and cdr3 length but vary the d gene.
for _ in range(1):
    print 'combining'
    raw_reco_seq = reco.combine('simulated-seqs.csv')
