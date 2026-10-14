#!/usr/bin/env python3
"""Write ogrdb-functional.csv (label, functional, sequence) for the OGRDB macaque germline sets G00091 (IGH), G00092 (IGK), and G00093 (IGL).

The sets themselves were staged in this dir as ungapped fastas (see data/germlines/README.md); here we get OGRDB's 'functional' field, which is only in
the AIRR-format download. 'sequence' is the ungapped coding sequence, i.e. the same sequence as in the staged fastas.
"""
import csv
import json
import os
import urllib.request

set_ids = ['G00091', 'G00092', 'G00093']
outfname = os.path.dirname(os.path.realpath(__file__)) + '/ogrdb-functional.csv'
rows = []
for set_id in set_ids:
    url = 'https://ogrdb.airr-community.org/api/germline/set/%s/published/airr' % set_id
    with urllib.request.urlopen(url) as response:
        gset = json.load(response)['GermlineSet'][0]
    print('  %s %s release %s (%s): %d alleles' % (set_id, gset['germline_set_name'], gset['release_version'], gset['release_date'], len(gset['allele_descriptions'])))
    for afo in gset['allele_descriptions']:
        rows.append({'label' : afo['label'], 'functional' : afo['functional'], 'sequence' : afo['coding_sequence'].upper().replace('.', ''), 'set' : '%s-release-%s' % (set_id, gset['release_version'])})
with open(outfname, 'w') as ofile:
    writer = csv.DictWriter(ofile, ['label', 'functional', 'sequence', 'set'])
    writer.writeheader()
    writer.writerows(rows)
print('  wrote %d alleles to %s' % (len(rows), outfname))
