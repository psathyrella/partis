#!/usr/bin/env python
import csv
import argparse
import sys
sys.path.insert(1, './python')

import seqfileopener
from glomerator import Glomerator
import utils

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')
parser.add_argument('--label')
parser.add_argument('--indir')
parser.add_argument('--simfname')
args = parser.parse_args()

input_info, reco_info = seqfileopener.get_seqfile_info(args.simfname, is_data=False)

infname = args.indir + '/' + args.label.replace('-', '_') + '_db-pass_parse-select_clone-pass.tab'

id_clusters = {}  # map from cluster id to list of seq ids
with open(infname) as chfile:
    reader = csv.DictReader(chfile, delimiter='\t')
    for line in reader:
        clid = line['CLONE']
        uid = line['SEQUENCE_ID']
        if clid not in id_clusters:
            id_clusters[clid] = []
        id_clusters[clid].append(uid)

partition = [ids for ids in id_clusters.values()]
glom = Glomerator(reco_info)
adj_mi = glom.mutual_information(partition, debug=True)
