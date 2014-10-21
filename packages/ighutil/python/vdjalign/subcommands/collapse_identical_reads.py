"""
Collapse identical reads.
"""
import csv
import itertools
import logging
import operator
import sys

from vdjalign.util import opener


KEEP = ['nucleotide', 'copy', 'cdr3Length', 'jGeneName', 'jTies', 'vIndex',
        'jIndex', 'source', 'sequenceStatus']

OUT = ['name'] + KEEP + ['sources', 'n_sources']


log = logging.getLogger(__name__)


def load_csv_rows(fname, conversions=None, **kwargs):
    log.info("Loading %s", fname)
    with opener('rU')(fname) as fp:
        lines = (i for i in fp if not i.startswith('#'))
        r = csv.DictReader(lines, **kwargs)
        for row in r:
            row['source'] = fname
            if conversions:
                for k, f in conversions:
                    row[k] = f(row[k])
            yield row


def subset_dict(d, keys):
    return {k: d[k] for k in keys}


def build_parser(p):
    p.add_argument('infiles', metavar='infile', nargs='+')
    p.add_argument('--prefix', default='')
    p.add_argument('-o', '--outfile', default=sys.stdout, type=opener('w'))
    p.set_defaults(func=action)


def action(a):
    rows = (row for f in a.infiles for row in load_csv_rows(f, [('copy', int)],
                                                            delimiter='\t'))
    rows = (subset_dict(d, KEEP) for d in rows)
    sort_key = operator.itemgetter('nucleotide', 'copy', 'source')
    rows = sorted(rows, key=sort_key, reverse=True)
    logging.info("All data loaded and sorted.")

    with a.outfile as ofp:
        writer = csv.DictWriter(ofp, OUT, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        grouped = itertools.groupby(rows, operator.itemgetter('nucleotide'))
        grouped = (list(v) for _, v in grouped)
        grouped = sorted(grouped,
                         key=lambda v: (len(v), sum(i['copy'] for i in v)),
                         reverse=True)

        for i, values in enumerate(grouped):
            representative = values[0]
            representative['n_sources'] = len(values)
            representative['sources'] = ','.join(i['source'] for i in values)
            representative['copy'] += sum(i['copy'] for i in values)
            representative['name'] = '{0}{1:07d}'.format(a.prefix, i)
            writer.writerow(representative)
