#!/usr/bin/env python
from recombinator import Recombinator

#----------------------------------------------------------------------------------------
reco = Recombinator()
for _ in range(100):
    print 'combining'
    raw_reco_seq = reco.combine('out.csv')
#    raw_reco_seq = reco.combine()
