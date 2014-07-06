#!/usr/bin/python
""" Count the number of occurrences of each v/d/j/cdr3_length combination """

import sys
from versioncounter import VersionCounter
import utils

cryptic_strings = (('04-A-M',),
                   ('02-C-M',),
                   ('01-C-N',),
                   ('06-B-M',),
                   ('03-A-N1', '03-A-N2'),
                   ('03-A-N2',),
                   ('05-B-N',))

indir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw'

for crypt in cryptic_strings:
    crypt_split = crypt[0].split('-')
    human = crypt_split[1]
    stage = crypt_split[2].rstrip('1').rstrip('2')
    fnames = list(indir + '/' + human + '/' + subcrypt + '_filtered_train.vdjcdr3.csv.bz2' for subcrypt in crypt)
    for column_dependency in utils.column_dependencies:
        vc = VersionCounter(fnames, 'data/human-beings/' + human + '/' + stage, min_counts=3, index_columns=column_dependency)
        vc.count()
