#!/usr/bin/python
""" Count the number of occurrences of each v/d/j/cdr3_length combination """
# TEST
from versioncounter import VersionCounter

cryptic_strings = (('04-A-M',),
                   ('02-C-M',),
                   ('01-C-N',),
                   ('06-B-M',),
                   ('03-A-N1', '03-A-N2'),
                   ('03-A-N2',),
                   ('05-B-N',))

for crypt in cryptic_strings:
    crypt_split = crypt[0].split('-')
    person = crypt_split[1]
    stage = crypt_split[2].rstrip('1').rstrip('2')
    fnames = ('data/human-beings/' + subcrypt + '_filtered.vdjcdr3.csv.bz2' for subcrypt in crypt)
    vc = VersionCounter(fnames, 'data/human-beings/' + person + '/' + stage)
    vc.count()
