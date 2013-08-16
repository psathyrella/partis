"""
Graph of uncertainty
"""
import collections
import csv
import logging
import functools
import operator
import sys

import networkx

from .. import util

log = logging.getLogger('adaptive-vj-uncertainty')

Combination = collections.namedtuple('Combination', ['vs', 'js', 'cdr3_length', 'count'])

def load_csv_rows(fname, conversions=None, **kwargs):
    log.info("Loading %s", fname)
    with util.opener('rU')(fname) as fp:
        lines = (i for i in fp if not i.startswith('#'))
        r = csv.DictReader(lines, **kwargs)
        for row in r:
            row['source'] = fname
            if conversions:
                for k, f in conversions:
                    row[k] = f(row[k])
            yield row

def or_none(f):
    @functools.wraps(f)
    def inner(s):
        if s is not None and s != '':
            return f(s)
    return inner

def int_set(s):
    return frozenset(int(i) for i in s.split(','))

def distinct_by(iterable, key=None):
    seen = set()
    for i in iterable:
        k = key(i) if key else i
        if k not in seen:
            seen.add(k)
            yield i

def row_to_combination(row):
    if row['vGeneName'] != -1:
        vs = frozenset([row['vGeneName']])
    else:
        vs = frozenset(row['vTies'])
    if row['jGeneName'] != -1:
        js = frozenset([row['jGeneName']])
    else:
        js = frozenset(row['jTies'])
    return Combination(vs, js, row['cdr3Length'], 1)

def build_parser(p):
    p.add_argument('infiles', metavar='infile', nargs='+')
    p.add_argument('--prefix', default='')
    p.add_argument('-o', '--outfile', default=sys.stdout, type=util.opener('w'))
    p.set_defaults(func=action)

def resolutions(comb):
    for v in comb.vs:
        for j in comb.js:
            yield v, j, comb.cdr3_length

def freqs_to_graph(freqs):
    resolution_nodes = collections.defaultdict(set)
    g = networkx.Graph()
    for n in freqs:
        g.add_node(n)
        for resolution in resolutions(n):
            neighbors = resolution_nodes[resolution]
            for neighbor in neighbors:
                g.add_edge(n, neighbor)
            neighbors.add(n)
    return g

def action(a):
    rows = (row for f in a.infiles
            for row in load_csv_rows(f, [('vGeneName', or_none(int)),
                                         ('jGeneName', or_none(int)),
                                         ('vTies', or_none(int_set)),
                                         ('jTies', or_none(int_set)),
                                         ('cdr3Length', or_none(int))],
                                      delimiter='\t'))
    rows = distinct_by(rows, operator.itemgetter('nucleotide'))
    rows = (row_to_combination(row) for row in rows)
    freqs = (k._replace(count=v) for k, v in collections.Counter(rows).iteritems())
    logging.info('creating graph')
    g = freqs_to_graph(freqs)
    logging.info('%d nodes, %d edges', g.number_of_nodes(), g.number_of_edges())

    with a.outfile as fp:
        w = csv.writer(fp, lineterminator='\n')
        w.writerow(['vs', 'js', 'cdr3_length', 'n_nodes', 'total_count'])
        for cc in networkx.connected_components(g):
            vs = '|'.join(set(str(j) for i in cc for j in i.vs))
            js = '|'.join(set(str(j) for i in cc for j in i.js))
            cdr3s = '|'.join(set(str(i.cdr3_length) for i in cc))
            w.writerow([vs, js, cdr3s, len(cc), sum(i.count for i in cc)])
