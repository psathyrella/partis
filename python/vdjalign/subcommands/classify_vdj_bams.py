"""
Annotate sequences based on VDJ alignments
"""
import collections
import csv
import itertools
import logging
import sys

from pkg_resources import resource_stream
import pysam

from .. import gff3, util
from ..classify_vdj import classify_record

log = logging.getLogger('classify-vdj-bams')

def _load_annotations_by_seqid(name):
    result = collections.defaultdict(dict)
    with resource_stream('vdjalign.imgt.data', name) as fp:
        records = gff3.parse(fp)

        for record in records:
            d = result[record.seqid]
            attr = record.attribute_dict()
            if attr['Name'] in d:
                raise ValueError('{0} already in {1}'.format(attr['Name'], d))

            d[attr['Name']] = record

    return result


def load_v_annotations():
    return _load_annotations_by_seqid('ighv.gff3')

def load_j_annotations():
    return _load_annotations_by_seqid('ighj_adaptive.gff3')

def build_parser(p):
    p.add_argument('v_bam', type=pysam.Samfile)
    p.add_argument('d_bam', type=pysam.Samfile)
    p.add_argument('j_bam', type=pysam.Samfile)
    p.add_argument('-o', '--outfile', default=sys.stdout, type=util.opener('w'))
    p.set_defaults(func=action)


def action(a):
    v_annot = load_v_annotations()
    all_loci = frozenset(i for v in v_annot.values() for i in v)
    j_annot = load_j_annotations()

    it = itertools.izip_longest(a.v_bam, a.d_bam, a.j_bam)

    v_targets = a.v_bam.references
    d_targets = a.d_bam.references
    j_targets = a.j_bam.references

    with a.outfile as ofp:
        w = None
        for v, d, j in it:
            if not all([v, d, j]):
                raise ValueError('Number of records does not match: {0}'.format([v, d, j]))
            if not d.qname == v.qname and j.qname == v.qname:
                raise ValueError('Names do not match: {0.qname} {1.qname} {2.qname}'.format(v, d, j))


            # Everything checks out
            classification = classify_record(v, d, j,
                                             v_annot=v_annot[v_targets[v.tid]],
                                             j_annot=j_annot[j_targets[j.tid]])
            classification['qname'] = v.qname
            classification['v'] = v_targets[v.tid]
            classification['d'] = d_targets[d.tid]
            classification['j'] = j_targets[j.tid]

            for i in all_loci:
                for k in ['begin', 'end', 'complete']:
                    key = '_'.join((i, k))
                    if key not in classification:
                        classification[key] = None

            if w is None:
                w = csv.DictWriter(ofp, classification.keys(), lineterminator='\n')
                w.writeheader()

            w.writerow(classification)
