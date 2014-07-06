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

# correlations are taken from figure in the bcellap repo
# first entry is the column of interest, and it depends upon the following entries
column_lists = []
column_lists.append(('v_gene',))  # TODO v choice actually depends on everything... but not super strongly, so a.t.m. I ignore it
column_lists.append(('v_3p_del',      'v_gene'))
column_lists.append(('d_gene',        'd_5p_del', 'd_3p_del'))
column_lists.append(('d_5p_del',      'd_3p_del', 'd_gene'))
column_lists.append(('d_3p_del',      'd_5p_del', 'd_gene'))
column_lists.append(('j_5p_del',))  # strange but seemingly true: does not depend on j choice
column_lists.append(('vd_insertion',))
column_lists.append(('dj_insertion',  'j_gene'))
column_lists.append(utils.index_columns)

for crypt in cryptic_strings:
    crypt_split = crypt[0].split('-')
    human = crypt_split[1]
    stage = crypt_split[2].rstrip('1').rstrip('2')
    fnames = list(indir + '/' + human + '/' + subcrypt + '_filtered_train.vdjcdr3.csv.bz2' for subcrypt in crypt)
    for column_list in column_lists:
        vc = VersionCounter(fnames, 'data/human-beings/' + human + '/' + stage, min_counts=10, index_columns=column_list)
        vc.count()
    sys.exit()
