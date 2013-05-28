"""
Access to IMGT germline abs
"""
import itertools
from pkg_resources import resource_stream

from .. import util

_PKG = 'vdjalign.imgt.data'

# 0-based index of cysteine in NT alignment (position 104 in AA)
_CYSTEINE_POSITION = 3 * 104 - 1

_GAP = '.-'

def _position_lookup(sequence):
    result = {}
    seq_idx = itertools.count().next

    for i, base in enumerate(sequence):
        if base not in _GAP:
            result[i] = seq_idx()

    return result

def cysteine_map():
    """
    Calculate cysteine position in each sequence of the ighv alignment
    """
    result = {}
    with resource_stream(_PKG, 'ighv_aligned.fasta') as fp:
        sequences = ((name.split('|')[1], seq) for name, seq, _ in util.readfq(fp))

        for name, s in sequences:
            result[name] = _position_lookup(s).get(104 * 3 - 1)
    return result

def phe_trp_map():
    """
    Calculate T/P position in each sequence of the ighj alignment
    """
    result = {}
    with resource_stream(_PKG, 'ighj_orig.fasta') as fp:
        sequences = ((name, seq.upper()) for name, seq, _ in util.readfq(fp))

        for desc, s in sequences:
            splt = desc.split('|')
            name = splt[1]
            codon_start = int(splt[7])

            for i in xrange(codon_start - 1, len(s), 3):
                if s[i:i+3] in ('TTT', 'TTC', 'TGG'):
                    result[name] = i
                    break
            else:
                raise ValueError(s)
    return result


