import csv
import sys
import os
import re

import python.utils as utils

# ----------------------------------------------------------------------------------------
def get_gcid(prn, mouse, node, gcn):
    # return 'pr-%s-m-%s-gc-%s' % (prn, mouse, gcn)
    return 'PR%s-%s-%s-%s-GC' % (prn, mouse, node, gcn)
# PR%s*-%s-%s-%s-GC
# btt-PR-1-8-1-LP-4-GC

# ----------------------------------------------------------------------------------------
def reverse_gcid(gcid):
    if 'btt' in gcid:  # UGH
        gcid = fix_btt_id(gcid)
    prstr, mouse, node, gcn, gcstr = gcid.split('-')
    prn = prstr.replace('PR', '')
    assert gcstr == 'GC'
    return prn, mouse, node, gcn

# ----------------------------------------------------------------------------------------
def fix_btt_id(gcid):
    mstr = utils.get_single_entry(re.findall('btt-PR-.-.', gcid))
    btstr, prstr, prn1, prn2  = mstr.split('-')
    assert btstr == 'btt' and prstr == 'PR'
    return gcid.replace(mstr, 'PR%d.%02d' % (int(prn1), int(prn2)))

# ----------------------------------------------------------------------------------------
def read_gcreplay_metadata(gcreplay_dir):
    # ----------------------------------------------------------------------------------------
    def readcfn(prn):
        with open('%s/metadata.PR%s.csv' % (gcreplay_dir, prn)) as cfile:
            reader = csv.DictReader(cfile)
            for line in reader:
                if prn == '1':
                    line['time'] = line['imm_duration']
                    del line['imm_duration']
                elif prn == '2':
                    line['time'] = 'd15'
                    line['strain'] = 'wt'
                else:
                    assert False
                gcid = get_gcid(line['PR'], line['mouse'], line['node'], line['gc'])
                assert gcid not in rpmeta
                rpmeta[gcid] = line
    # ----------------------------------------------------------------------------------------
    print('  reading gcreplay meta info from %s' % gcreplay_dir)
    rpmeta = {}
    for prn in ['1', '2']:
        readcfn(prn)
    return rpmeta
