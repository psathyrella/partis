#!/usr/bin/env python

from Recombinator import Recombinator

#----------------------------------------------------------------------------------------
reco = Recombinator()
for _ in range(1):
    print 'combining'
    raw_reco_seq = reco.combine('out.csv')
#    print 'final %s' % raw_reco_seq
