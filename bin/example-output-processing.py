#!/usr/bin/env python
import csv
import sys

partis_path = '.'  # edit this if you're not running from the main partis dir
sys.path.insert(1, partis_path + '/python')
import utils
import glutils
from clusterpath import ClusterPath

# read default germline info
glfo = glutils.read_glfo(partis_path + '/data/germlines/human', locus='igh')

print 'first parse an annotation csv file:'
with open(partis_path + '/test/reference-results/annotate-new-simu.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for line in reader:
        if line['v_gene'] == '':  # failed (i.e. couldn't find an annotation)
            continue
        utils.process_input_line(line)  # converts strings in the csv file to floats/ints/dicts/etc.
        utils.add_implicit_info(glfo, line)  # add stuff to <line> that's useful, isn't written to the csv since it's redundant
        utils.print_reco_event(line)  # print ascii-art representation of the rearrangement event
        print '\navailable annotation info for each line (see manual for descriptions):'
        for key, val in line.items():
            print '%20s %s' % (key, val)
        break

print '\n\nthen parse a partition csv file:'
cp = ClusterPath()
cp.readfile(partis_path + '/test/reference-results/seed-partition-new-simu.csv')
cp.print_partitions(abbreviate=True)
