#!/usr/bin/env python

import argparse
import recombinator


def make_events(args, n_events, iproc, random_ints):
    # NOTE all the different seeds! this sucks but is necessary
    reco = Recombinator(args, seed=args.seed+iproc, sublabel=str(iproc))
    for ievt in range(n_events):
        # print ievt,
        # sys.stdout.flush()
        reco.combine(random_ints[ievt])
