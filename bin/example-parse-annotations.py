#!/usr/bin/env python
import csv
import os
import sys
import argparse

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils

parser = argparse.ArgumentParser()
parser.add_argument('--annotation-file', default=partis_dir + '/test/reference-results/annotate-ref-simu.yaml')
args = parser.parse_args()

glfo, annotations = utils.read_yaml_annotations(args.annotation_file)
print 'print one annotation, then break:'
for line in annotations.values():
    utils.print_reco_event(line)  # print ascii-art representation of the rearrangement event
    print '\n\navailable keys:'
    for key, val in line.items():
        print '%20s %s' % (key, val)
    break

# for old csv files:
# parser.add_argument('--annotation-file', default=partis_dir + '/test/reference-results/annotate-ref-simu.csv')
# parser.add_argument('--germline-info-dir', default=partis_dir + '/test/reference-results/test/parameters/simu/hmm/germline-sets')
# parser.add_argument('--locus', default='igh')
# glfo = glutils.read_glfo(args.germline_info_dir, locus=args.locus)  # read germline info (it'll be nice to switch to yaml output files so we can store this in the same file as the rest of the output)
# with open(args.annotation_file) as csvfile:
#     reader = csv.DictReader(csvfile)
#     for line in reader:  # one line for each annotation
#         if line['v_gene'] == '':  # failed (i.e. couldn't find an annotation)
#             continue
#         utils.process_input_line(line)  # converts strings in the csv file to floats/ints/dicts/etc.
#         utils.add_implicit_info(glfo, line)  # add stuff to <line> that's useful, but isn't written to the csv since it's redundant
#         print 'print one annotation, then break:'
#         utils.print_reco_event(line)  # print ascii-art representation of the rearrangement event
#         print '\n\navailable keys:'
#         for key, val in line.items():
#             print '%20s %s' % (key, val)
#         break
