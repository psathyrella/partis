import csv
import sys
import os
import re

import partis.utils as utils

replay_hl_naive_seqs = [
    'GAGGTGCAGCTTCAGGAGTCAGGACCTAGCCTCGTGAAACCTTCTCAGACTCTGTCCCTCACCTGTTCTGTCACTGGCGACTCCATCACCAGTGGTTACTGGAACTGGATCCGGAAATTCCCAGGGAATAAACTTGAGTACATGGGGTACATAAGCTACAGTGGTAGCACTTACTACAATCCATCTCTCAAAAGTCGAATCTCCATCACTCGAGACACATCCAAGAACCAGTACTACCTGCAGTTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCAAGGGACTTCGATGTCTGGGGCGCAGGGACCACGGTCACCGTCTCCTCA',
    'GACATTGTGATGACTCAGTCTCAAAAATTCATGTCCACATCAGTAGGAGACAGGGTCAGCGTCACCTGCAAGGCCAGTCAGAATGTGGGTACTAATGTAGCCTGGTATCAACAGAAACCAGGGCAATCTCCTAAAGCACTGATTTACTCGGCATCCTACAGGTACAGTGGAGTCCCTGATCGCTTCACAGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAATGTGCAGTCTGAAGACTTGGCAGAGTATTTCTGTCAGCAATATAACAGCTATCCTCTCACGTTCGGCTCGGGGACTAAGCTAGAAATAAAA'
]
replay_naive_seq = ''.join(replay_hl_naive_seqs)

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
    if 'btt-' not in gcid:
        return gcid
    mstr = utils.get_single_entry(re.findall('btt-PR-.-[0-9][0-9]*', gcid))
    btstr, prstr, prn1, prn2  = mstr.split('-')
    assert btstr == 'btt' and prstr == 'PR'
    return gcid.replace(mstr, 'PR%d.%02d' % (int(prn1), int(prn2)))

# ----------------------------------------------------------------------------------------
def read_gcreplay_metadata(gcreplay_dir, old_style=False):
    # ----------------------------------------------------------------------------------------
    def readcfn(prn=None):
        mfn = ('%s/gc_metadata.csv'%gcreplay_dir) if prn is None else '%s/metadata.PR%s.csv' % (gcreplay_dir, prn)
        with open(mfn) as cfile:
            reader = csv.DictReader(cfile)
            for line in reader:
                if prn in ['1', None]:
                    line['time'] = line['imm_duration']
                    del line['imm_duration']
                elif prn == '2':
                    line['time'] = 'd15'
                    line['strain'] = 'wt'
                else:
                    assert False
                if prn is None:  # new style file
                    gcid = line['uid']
                else:
                    gcid = get_gcid(line['PR'], line['mouse'], line['node'], line['gc'])
                assert gcid not in rpmeta
                rpmeta[gcid] = line
    # ----------------------------------------------------------------------------------------
    print('  reading gcreplay meta info from %s' % gcreplay_dir)
    rpmeta = {}
    if old_style:  # UGH
        for prn in ['1', '2']:
            readcfn(prn=prn)
    else:
        readcfn()
    return rpmeta

