#!/usr/bin/env python
from Recombinator import Recombinator
#----------------------------------------------------------------------------------------
#
reco = Recombinator()
for _ in range(20):
    print 'combining...'
    raw_reco_seq = reco.combine()
    print 'final %s' % raw_reco_seq
