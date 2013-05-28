"""
Access to IMGT germline abs
"""
import itertools
from pkg_resources import resource_stream

from .. import util

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
    Calculate cysteine position in each sequence
    """
    result = {}
    with resource_stream('vdjalign.imgt.data', 'ighv_aligned.fasta') as fp:
        sequences = ((name.split('|')[1], seq) for name, seq, _ in util.readfq(fp))

        for name, s in sequences:
            result[name] = _position_lookup(s).get(104 * 3 - 1)
    return result
