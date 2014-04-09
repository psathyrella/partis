#!/usr/bin/python

from versioncounter import VersionCounter

vc = VersionCounter('data/human-beings/01-C-N_filtered.vdjcdr3.csv.bz2', 'data/human-beings/A/N')
vc.count()

