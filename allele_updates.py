#!/usr/bin/env python

from __future__ import division
import collections
import csv
import json
import itertools
import operator

from vdjalign.util import readfq

NUM_COLS = ['A', 'C', 'G', 'T', 'n_reads', 'alignment_position']

def load_csv(fp):
    rows = csv.DictReader(fp)
    for row in rows:
        for c in NUM_COLS:
            row[c] = int(row[c])
        yield row


def list_changes(group):
    r = []
    od = collections.OrderedDict
    for i in group:
        n = float(i['n_reads'])
        for k in 'ACGT':
            if k != i['ref_base'] and i[k] / n > 0.1:
                r.append(od([('position', i['alignment_position']),
                             ('from', i['ref_base']),
                             ('to', k),
                             ('n', i[k]),
                             ('t', n),
                             ('p', i[k]/n)]))
    return r


def introduce_mutations(sequence, changes):
    sequence = list(sequence)
    for change in changes:
        pos = change['position']
        assert(sequence[pos] == change['from'])
        sequence[pos] = change['to']
    return ''.join(sequence)



def main():
    with open('python/vdjalign/imgt/data/ighv_aligned.fasta') as fp:
        sequences = {i[0].split('|')[1].strip(): i[:2] for i in readfq(fp)}

    with open('poss_alleles.csv') as fp:
        rows = list(load_csv(fp))

    grouped = itertools.groupby(rows, operator.itemgetter('subject', 'reference'))
    new_alleles = collections.Counter()
    changes = collections.OrderedDict()
    seen = set()
    for (_, r), g in grouped:
        base = r.rpartition('*')[0]
        new_allele = new_alleles[base] + 91
        new_ref = '{0}*{1:02d}'.format(base, new_allele)
        c = list_changes(g)
        cset = (base, frozenset((i['position'], i['to']) for i in c))
        if cset in seen:
            print 'skipping', cset
            continue
        seen.add(cset)

        changes[new_ref] = {'base': r, 'changes': c}
        new_alleles[base] += 1

    with open('allele_updates.json', 'w') as ofp:
        json.dump(changes, ofp, indent=2)

    with open('allele_updates.fasta', 'w') as ofp:
        for k, v in changes.iteritems():
            orig_id, orig_sequence = sequences[v['base']]
            new_id = orig_id.replace(v['base'], k).split('|')
            new_id[0] += '_update'
            new_id = '|'.join(new_id)
            sequence = introduce_mutations(orig_sequence, v['changes'])
            ofp.write('>{0}\n{1}\n'.format(new_id, sequence))




if __name__ == '__main__':
    main()
