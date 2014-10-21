"""
Extract the best scoring V, D, J
"""
import contextlib
import csv
import itertools
import logging
import operator
import sys

import pysam

from .. import util

log = logging.getLogger('extract-best-hits')


def build_parser(p):
    p.add_argument('bam_file', type=pysam.Samfile)
    p.add_argument('-o', '--outfile', type=util.opener('w'),
                   default=sys.stdout)

    p.set_defaults(func=action)


def action(a):
    with contextlib.closing(a.bam_file) as bam, a.outfile as ofp:
        grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))

        references = bam.references

        segments = [i[:4] for i in references]

        def get_segment(read):
            return segments[read.tid]

        w = csv.DictWriter(ofp, ['seq_id', 'IGHV', 'IGHD', 'IGHJ'],
                            lineterminator='\n')
        w.writeheader()
        for qname, records in grouped:
            records = sorted(records, key=get_segment)
            by_segment = itertools.groupby(records, get_segment)
            d = {'seq_id': qname}
            d.update({k: references[next(v).tid] for k, v in by_segment})
            w.writerow(d)
