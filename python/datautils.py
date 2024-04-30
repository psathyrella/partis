import csv
import sys
import os

# ----------------------------------------------------------------------------------------
def get_gcid(prn, mouse, node, gcn):
    # print(prn, mouse, gcn)
    # sys.exit()
    # return 'pr-%s-m-%s-gc-%s' % (prn, mouse, gcn)
    return 'PR%s-%s-%s-%s-GC' % (prn, mouse, node, gcn)
# PR%s*-%s-%s-%s-GC
# btt-PR-1-8-1-LP-4-GC

# ----------------------------------------------------------------------------------------
def reverse_gcid(gcid):
    # return 'pr-%s-m-%s-gc-%s' % (prn, mouse, gcn)
    return prn, mouse, gcn

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

