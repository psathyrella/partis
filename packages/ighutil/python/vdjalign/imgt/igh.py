"""
Access to IMGT germline abs
"""
import csv
import contextlib
import itertools
import os.path
import shutil
import subprocess

from pkg_resources import resource_stream

from lxml import etree

from .. import util

_PKG = __package__ + '.data'

# 0-based index of cysteine in NT alignment (position 104 in AA)
_CYSTEINE_POSITION = 3 * (104 - 1)

_GAP = '.-'

def _remove_allele(s):
    return s.rpartition('*')[0]

def _position_lookup(sequence):
    result = {}
    seq_idx = itertools.count().next

    for i, base in enumerate(sequence):
        if base not in _GAP:
            result[i] = seq_idx()

    return result

def ighv_position_lookup(fname='ighv_aligned.fasta'):
    result = {}
    with resource_stream(_PKG, fname) as fp:
        sequences = ((name.split('|')[1], seq)
                     for name, seq, _ in util.readfq(fp))

        for name, s in sequences:
            p = [i for i, b in enumerate(s) if b not in _GAP]
            assert len(p) == len(s.translate(None, _GAP))
            result[name] = p

    return result

def cysteine_map(fname='ighv_aligned.fasta'):
    """
    Calculate cysteine position in each sequence of the ighv alignment
    """
    result = {}
    with resource_stream(_PKG, fname) as fp:
        sequences = ((name.split('|')[1], seq)
                     for name, seq, _ in util.readfq(fp))

        for name, s in sequences:
            name = _remove_allele(name)
            p = _position_lookup(s).get(_CYSTEINE_POSITION)
            if p is None:
                continue
            result[name] = max(p, result.get(name, -1))

    return result

def cdr_gff3(out_fp, fname='ighv_aligned.fasta', gffname='ighv_aligned.gff3'):
    """
    Calculate cysteine position in each sequence of the ighv alignment
    """

    from .. import gff3

    with resource_stream(_PKG, fname) as fp, resource_stream(_PKG, gffname) as gff_fp:
        sequences = ((name.split('|')[1], seq)
                     for name, seq, _ in util.readfq(fp))

        records = list(gff3.parse(gff_fp))

        out_fp.write('##gff-version\t3\n')
        w = csv.writer(out_fp, delimiter='\t', quoting=csv.QUOTE_NONE, lineterminator='\n')
        for name, s in sequences:
            pl = _position_lookup(s)

            for feature in records:
                if feature.start0 not in pl and feature.end - 1 not in pl:
                    continue

                start0 = pl.get(feature.start0, min(q for r, q in pl.items() if r >= feature.start0))
                end = pl.get(feature.end - 1, max(q for r, q in pl.items() if r <= feature.end - 1))

                f = feature._replace(seqid=name,
                                     start=start0 + 1,
                                     end=end + 1)
                w.writerow(f.update_attributes(ID=f.attribute_dict()['Name'] + '_' + name))

def load_ighj():
    with resource_stream(_PKG, 'ighj_adaptive.xml') as fp:
        doc = etree.parse(fp)

        for elem in doc.xpath('/region/germline'):
            record = dict(elem.attrib)
            for k in ('length', 'aminoacidindex', 'cdr3index'):
                record[k] = int(record[k])

            record['trp_index'] = record['length'] + record['cdr3index'] - 3
            record['cdr3end'] = record['trp_index'] + 3

            yield record


def consensus_by_allele(file_name):
    with resource_stream(_PKG, file_name) as fp:
        sequences = ((name.split('|')[1], seq)
                     for name, seq, _ in util.readfq(fp))
        sequences = sorted(sequences)
        grouped = itertools.groupby(sequences,
                                    lambda (name, _): _remove_allele(name))
        for g, seqs in grouped:
            seqs = list(seqs)
            drop_indices = frozenset([i for i in xrange(len(seqs[0][1]))
                                      if all(s[i] in _GAP for _, s in seqs)])

            def remove_dropped(s):
                return ''.join(c for i, c in enumerate(s)
                               if i not in drop_indices)

            if len(seqs) == 1:
                yield g, remove_dropped(seqs[0][1])
                continue

            # Build consensus
            cmd = ['consambig', '-name', g, '-filter']
            p = subprocess.Popen(cmd,
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
            with p.stdin as ifp:
                for name, seq in seqs:
                    ifp.write('>{0}\n{1}\n'.format(name, seq.replace('.', '-')))
            with p.stdout as fp:
                name, seq, _ = next(util.readfq(fp))
            returncode = p.wait()
            assert returncode == 0, returncode
            yield name, remove_dropped(seq).upper()

