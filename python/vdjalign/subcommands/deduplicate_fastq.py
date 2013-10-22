"""
Collapse identical reads, adding a count to the first record
"""
import logging
import operator
import sys

from vdjalign.util import opener, readfq

log = logging.getLogger(__name__)

class ReadWithAbundance(object):
    __slots__ = ['name', 'sequence', 'count']
    def __init__(self, name, sequence, count=1):
        self.name = name
        self.sequence = sequence
        self.count = count

    def __repr__(self):
        return '<{0} [name={1}; sequence={2}; count={3}]>'.format(
            self.__class__.__name__, self.name, self.sequence, self.count)


def build_parser(p):
    p.add_argument('infile', type=opener('r'), help="""Input FAST[AQ]""")
    p.add_argument('-o', '--outfile', default=sys.stdout, type=opener('w'))
    p.set_defaults(func=action)

def action(a):
    with a.infile as ifp:
        seqs = readfq(ifp)
        result = {}
        for name, seq, _ in seqs:
            try:
                result[seq].count += 1
            except KeyError:
                result[seq] = ReadWithAbundance(name, seq)

    dedup = sorted(result.values(), key=operator.attrgetter('count'),
                   reverse=True)

    with a.outfile as ofp:
        for read in dedup:
            ofp.write('>{0} {1}\n{2}\n'.format(read.name, read.count, read.sequence))
