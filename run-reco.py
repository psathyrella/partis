#!/usr/bin/env python
import os
import sys
from recombinator import Recombinator
import utils

#----------------------------------------------------------------------------------------
#assert len(sys.argv) == 2
#human = sys.argv[1]  #'C'
naivety = 'M'
for human in utils.humans:
    reco = Recombinator('data/human-beings', human, naivety)
    for _ in range(1):
        outdir = 'output/' + human + '/' + naivety
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outfname = outdir + '/simu.csv'
        success = False
        while not success:
            success = reco.combine()  #outfname, 'append')
        sys.exit()
