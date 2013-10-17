"""
Align V and J starting from an adaptive CSV file
"""
import collections
import contextlib
import csv
import functools
import itertools
import logging
import operator
import tempfile


import pysam

from .. import util
from . import align_fastq

log = logging.getLogger('vdjalign')
TAG_COUNT = 'XC'
TAG_CDR3_LENGTH = 'XL'
TAG_STATUS = 'XS'

def add_tags(reads, rows):
    """Tags should be the length of reads, in same order"""
    reads = itertools.groupby(reads, operator.attrgetter('qname'))
    for r, (g, v) in itertools.izip_longest(rows, reads):
        assert r.name == g, (g, r.name)
        for read in v:
            t = read.tags
            t.extend(r.tags.iteritems())
            read.tags = t
            yield read

def or_none(fn):
    @functools.wraps(fn)
    def apply_or_none(s):
        if s:
            return fn(s)
    return apply_or_none


def build_parser(p):
    p.add_argument('csv_file', type=util.opener('rU'))
    p.add_argument('-d', '--delimiter', default='\t',
                   help="""Delimiter [default: tab]""")
    p.add_argument('-c', '--count-column', default='n_sources')
    p.add_argument('-j', '--threads', default=1, type=int, help="""Number of
                   threads [default: %(default)d]""")
    p.add_argument('out_bamfile', help="""Path for V alignments""")
    p.add_argument('-r', '--read-group')
    agrp = p.add_argument_group('Alignment options')
    agrp.add_argument('-k', '--keep', help="""Number of alignments to keep per
                      gene [default: 15V, 5J, 10D]""", type=int)
    agrp.add_argument('-m', '--match', default=1, type=int, help="""Match score
                      [default: %(default)d]""")
    agrp.add_argument('-u', '--mismatch', default=1, type=int, help="""Match score
                      [default: %(default)d]""")
    agrp.add_argument('-o', '--gap-open', default=7, type=int, help="""Gap
                      opening penalty [default: %(default)d]""")
    agrp.add_argument('-e', '--gap-extend', default=1, type=int, help="""Gap
                      extension penalty [default: %(default)d]""")
    agrp.add_argument('--max-drop', default=0, type=int, help="""Maximum
                      alignment score drop [default: %(default)d]""")
    p.set_defaults(func=action)

def action(a):
    with a.csv_file as ifp:
        csv_lines = (i for i in ifp if not i.startswith('#'))
        r = csv.DictReader(csv_lines, delimiter=a.delimiter)
        int_or_none = or_none(int)
        Row = collections.namedtuple('Row', ['name', 'sequence', 'v_index',
                                             'j_index', 'tags'])
        log.info('Loading sequences.')
        rows= [Row(name=row.get('name', str(i)),
                   sequence=row['nucleotide'],
                   v_index=int_or_none(row.get('vIndex')),
                   j_index=int_or_none(row.get('jIndex')),
                   tags={TAG_COUNT: int_or_none(row[a.count_column]),
                         TAG_CDR3_LENGTH: int_or_none(row['cdr3Length']),
                         TAG_STATUS: int(row['sequenceStatus'])})
               for i, row in enumerate(r)]
        log.info('Done.')

        sequences = [(i.name, i.sequence) for i in rows]

        align = functools.partial(align_fastq.sw_to_bam, n_threads=a.threads,
                                  read_group=a.read_group,
                                  match=a.match,
                                  mismatch=a.mismatch,
                                  gap_open=a.gap_open,
                                  gap_extend=a.gap_extend,
                                  max_drop=a.max_drop)

        closing = contextlib.closing
        with align_fastq.resource_fasta('ighv_functional.fasta') as vf, \
                align_fastq.resource_fasta('ighj_adaptive.fasta') as jf, \
                align_fastq.resource_fasta('ighd.fasta') as df, \
                tempfile.NamedTemporaryFile(suffix='.bam') as tbam, \
                tempfile.NamedTemporaryFile(suffix='.fasta') as tf:
            for name, sequence in sequences:
                tf.write('>{0}\n{1}\n'.format(name, sequence))
            tf.flush()

            align(vf, tf.name, tbam.name, extra_ref_paths=[df, jf])
            with closing(pysam.Samfile(tbam.name, 'rb')) as in_sam, \
                    closing(pysam.Samfile(a.out_bamfile, 'wb', template=in_sam)) as out_sam:
                for read in add_tags(in_sam, rows):
                    out_sam.write(read)
