#!/usr/bin/env python
import csv
import sys

partis_path = '.'  # edit this if you're not running from the main partis dir
sys.path.insert(1, partis_path + '/python')
import utils
import glutils
from clusterpath import ClusterPath

# read default germline info
glfo = glutils.read_glfo(partis_path + '/data/germlines/human', chain='h')

print 'first parse an annotation csv file:'
with open(partis_path + '/test/reference-results/annotate-new-simu.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for line in reader:
        utils.process_input_line(line)
        utils.add_implicit_info(glfo, line)
        utils.print_reco_event(line)
        cdr3_bounds = (line['codon_positions']['v'], line['codon_positions']['j'] + 3)
        print ''
        print '  should match the above:'
        print '    %s naive cdr3' % line['naive_seq'][cdr3_bounds[0] : cdr3_bounds[1]]
        print '    %s mature' % line['indel_reversed_seqs'][0][cdr3_bounds[0] : cdr3_bounds[1]]
        print ''
        break

print 'then parse a partition csv file:'
cp = ClusterPath()
cp.readfile(partis_path + '/test/reference-results/seed-partition-new-simu.csv')
cp.print_partitions(abbreviate=True)
