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

def resolutions(comb):
    for v in comb.vs:
        for j in comb.js:
            yield v, j, comb.cdr3_length

def combination_str(c):
    js = '|'.join(set(str(j) for j in c.js))
    vs = '|'.join(set(str(v) for v in c.vs))
    return '{0}-{1}-{2}'.format(vs, c.cdr3_length, js)

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

def build_parser(p):
    p.add_argument('infiles', metavar='infile', nargs='+')
    p.add_argument('-o', '--outfile', default=sys.stdout, type=util.opener('w'))
    p.add_argument('-c', '--frequency-cutoffs', default=[0], type=int_set)
    p.set_defaults(func=action)

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
    log.info('creating graph')
    g = freqs_to_graph(freqs)
    log.info('%d nodes, %d edges', g.number_of_nodes(), g.number_of_edges())

    with a.outfile as fp:
        w = csv.writer(fp, lineterminator='\n')
        w.writerow(['cutoff', 'vs', 'js', 'cdr3_length', 'n_nodes', 'total_count'])
        for cutoff in sorted(a.frequency_cutoffs):
            g2 = g.copy()
            g2.remove_nodes_from(i for i in g.nodes() if i.count < cutoff)
            connected = networkx.connected_components(g2)
            log.info('cutoff=%d: %d connected components  %d nodes  %d edges',
                     cutoff,
                     len(connected),
                     g2.number_of_nodes(),
                     g2.number_of_edges())
            for cc in connected:
                vs = '|'.join(set(str(j) for i in cc for j in i.vs))
                js = '|'.join(set(str(j) for i in cc for j in i.js))
                cdr3s = '|'.join(set(str(i.cdr3_length) for i in cc))
                w.writerow([cutoff, vs, js, cdr3s, len(cc), sum(i.count for i in cc)])
            #if cutoff == 1000:
                #cc = max((i for i in connected if len(i) < 1000), key=len)
                #g2.remove_nodes_from(i for i in g2.nodes() if i not in set(cc))
                #networkx.write_dot(g2, 'hairball.dot')
