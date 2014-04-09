#!/usr/bin/python

from versioncounter import VersionCounter

person = 'C'
stage = 'N'
vc = VersionCounter('data/human-beings/01-' + person + '-' + stage + '_filtered.vdjcdr3.csv.bz2', 'data/human-beings/' + person + '/' + stage)
vc.count()
