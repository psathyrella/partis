#!/usr/bin/env python
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always

# parse fasta files downloaded from ogrdb https://ogrdb.airr-community.org/germline_sets/Mouse
# (run from main partis dir)

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/data/germlines', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils

debug = False
for locus in ['igh', 'igk', 'igl']:
    print utils.color('blue', locus)
    merged_glfo = None  # glfo merged from original (imgt) mouse/ plus c57bl and balbc
    for mouse_type in ['c57bl', 'balbc']:
        print '  %s' % utils.color('blue_bkg', mouse_type)
        old_glfo = glutils.read_glfo('data/germlines/mouse', locus, debug=debug)
        ogdir = 'data/germlines/ogrdb-download/%s' % mouse_type
        new_glfo = glutils.read_glfo(ogdir, locus, template_glfo=old_glfo, add_dummy_name_components=True, debug=debug)
        glutils.write_glfo('data/germlines/%s' % mouse_type, new_glfo, debug=debug)
        if merged_glfo is None:
            merged_glfo = old_glfo
        merged_glfo = glutils.get_merged_glfo(merged_glfo, new_glfo, debug=True)
    glutils.write_glfo('data/germlines/merged-mouse', merged_glfo, debug=debug)
