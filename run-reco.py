#!/usr/bin/env python
from recombinator import Recombinator

#----------------------------------------------------------------------------------------
reco = Recombinator('data/human-beings/C/N/probs.csv.bz2')
for _ in range(100):
    print 'combining'
#    raw_reco_seq = reco.combine('out.csv')
    raw_reco_seq = reco.combine()
