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

from .. import imgt, util
from . import align_fastq

log = logging.getLogger('vdjalign')
TAG_COUNT = 'XC'
TAG_CDR3_LENGTH = 'XL'
TAG_STATUS = 'XS'

def add_tags(reads, rows):
    """Tags should be the length of reads, in same order"""
    reads = itertools.groupby(reads, operator.attrgetter('qname'))
    for (name, tags), (g, v) in itertools.izip_longest(rows, reads):
        assert name == g, (g, name)
        for read in v:
            t = read.tags
            t.extend(tags.iteritems())
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
    p.add_argument('-c', '--count-column', default='n_sources', help="""Column
                   containing (integer) counts [default: %(default)s]""")

    align_fastq.fill_targets_alignment_options(p)

    p.set_defaults(func=action)

def action(a):
    with a.csv_file as ifp:
        csv_lines = (i for i in ifp if not i.startswith('#'))
        r = csv.DictReader(csv_lines, delimiter=a.delimiter)
        int_or_none = or_none(int)
        Row = collections.namedtuple('Row', ['name', 'sequence', 'v_index',
                                             'j_index', 'tags'])
        log.info('Loading sequences.')
        rows = (Row(name=row.get('name', str(i)),
                    sequence=row['nucleotide'],
                    v_index=int_or_none(row.get('vIndex')),
                    j_index=int_or_none(row.get('jIndex')),
                    tags={TAG_COUNT: int_or_none(row[a.count_column]),
                          TAG_CDR3_LENGTH: int_or_none(row['cdr3Length']),
                          TAG_STATUS: int(row['sequenceStatus'])})
               for i, row in enumerate(r))

        sequences = []
        tags = []
        for i in rows:
            sequences.append((i.name, i.sequence))
            tags.append((i.name, i.tags))
        log.info('Done: read %d', len(sequences))

        align = functools.partial(align_fastq.sw_to_bam, n_threads=a.threads,
                                  read_group=a.read_group,
                                  match=a.match,
                                  mismatch=a.mismatch,
                                  gap_open=a.gap_open,
                                  gap_extend=a.gap_extend,
                                  min_score=a.min_score,
                                  bandwidth=a.bandwidth,
                                  max_drop=a.max_drop)

        closing = contextlib.closing

        with imgt.temp_fasta(a.locus, 'v', a.v_subset) as vf, \
                imgt.temp_fasta(a.locus, 'j', a.j_subset) as jf, \
                util.with_if(a.locus == 'IGH', imgt.temp_fasta, a.locus, 'd',
                        collection=a.d_subset) as df, \
                tempfile.NamedTemporaryFile(suffix='.bam') as tbam, \
                tempfile.NamedTemporaryFile(suffix='.fasta') as tf:
            for name, sequence in sequences:
                tf.write('>{0}\n{1}\n'.format(name, sequence))
            tf.flush()
            ex_refs = [i for i in [df, jf] if i is not None]
            align(vf, tf.name, tbam.name, extra_ref_paths=ex_refs)

            with closing(pysam.Samfile(tbam.name, 'rb')) as in_sam, \
                    closing(pysam.Samfile(a.out_bamfile, 'wb', template=in_sam)) as out_sam:
                for read in add_tags(in_sam, tags):
                    out_sam.write(read)
