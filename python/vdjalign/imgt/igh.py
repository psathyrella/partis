"""
Access to IMGT germline abs
"""
import itertools
import subprocess

from pkg_resources import resource_stream

from lxml import etree

from .. import util

_PKG = 'vdjalign.imgt.data'

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


def tryptophan_map():
    """
    Calculate T position in each sequence of the ighj alignment
    """
    result = {}
    with resource_stream(_PKG, 'ighj_orig.fasta') as fp:
        sequences = ((name, seq.upper()) for name, seq, _ in util.readfq(fp))

        for desc, s in sequences:
            splt = desc.split('|')
            name = splt[1]
            codon_start = int(splt[7])

            for i in xrange(codon_start - 1, len(s), 3):
                if s[i:i + 3] == 'TGG':
                    result[name] = i
                    break
            else:
                raise ValueError(s)
    return result


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
