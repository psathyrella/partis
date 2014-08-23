#!/usr/bin/python

import utils
from versioncounter import integrate_version_freqs

naivety = 'M'
for human in utils.humans:
    integrate_version_freqs('data/human-beings/' + human+ '/' + naivety+ '/probs.csv.bz2')
