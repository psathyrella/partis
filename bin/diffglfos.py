#!/usr/bin/env python
import argparse
import sys
import os
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')
import utils
import glutils

parser = argparse.ArgumentParser()
parser.add_argument('file1')
parser.add_argument('file2')
args = parser.parse_args()

genes = {fn : set() for fn in (args.file1, args.file2)}

def readfile(fname):
    with open(fname) as fastafile:
        n_skipped_pseudogenes = 0
        for line in fastafile:
            line = line.strip()
            if line == '':
                continue
            if line[0] != '>':
                continue
            linefo = [p.replace('>', '').strip() for p in line.split('|')]
            gene = None
            for piece in linefo:
                if piece[:2] == 'IG':
                    gene = piece
            if gene is None:
                raise Exception('couldn\'t fine gene in %s' % line)

            if len(linefo) > 1:
                functionality = linefo[glutils.imgt_info_indices.index('functionality')]
                if functionality not in glutils.functionalities:
                    raise Exception('unexpected functionality %s' % functionality)
                if functionality == 'P':
                    n_skipped_pseudogenes += 1
                    continue

            genes[fname].add(gene)
        if n_skipped_pseudogenes > 0:
            print '    skipped %d pseudogenes' % n_skipped_pseudogenes

readfile(args.file1)
readfile(args.file2)

print 'file1: %d' % len(genes[args.file1])
print 'file2: %d' % len(genes[args.file2])
print 'both: %d' % len(genes[args.file1] & genes[args.file2])
only_file1 = genes[args.file1] - genes[args.file2]
print 'only file1: %d  (%s)' % (len(only_file1), ' '.join([utils.color_gene(g) for g in only_file1]))
only_file2 = genes[args.file2] - genes[args.file1]
print 'only file2: %d  (%s)' % (len(only_file2), ' '.join([utils.color_gene(g) for g in only_file2]))
