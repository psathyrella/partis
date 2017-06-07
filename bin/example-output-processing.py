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
        utils.process_input_line(line)
        utils.add_implicit_info(glfo, line)
        utils.print_reco_event(line)
        break

print 'then parse a partition csv file:'
cp = ClusterPath()
cp.readfile(partis_path + '/test/reference-results/seed-partition-new-simu.csv')
cp.print_partitions(abbreviate=True)
