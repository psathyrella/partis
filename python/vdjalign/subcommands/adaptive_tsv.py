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
import os
import shutil
import subprocess
import tempfile

from pkg_resources import resource_stream

import pysam

from .. import util, sw

log = logging.getLogger('vdjalign')
TAG_COUNT = 'XC'
TAG_CDR3_LENGTH = 'XL'
TAG_STATUS = 'XS'

@contextlib.contextmanager
def resource_fasta(names):
    if isinstance(names, basestring):
        names = [names]
    if not names:
        raise ValueError('Missing sequence files')
    with util.tempdir(prefix=names[0]) as j, \
            open(j(names[0]), 'w') as fp:
        for name in names:
            with resource_stream('vdjalign.imgt.data', name) as in_fasta:
                shutil.copyfileobj(in_fasta, fp)
        fp.close()
        yield fp.name

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

@contextlib.contextmanager
def tmpfifo(**kwargs):
    fifo_name = kwargs.pop('name', 'fifo')
    with util.tempdir(**kwargs) as j:
        f = j(fifo_name)
        os.mkfifo(f)
        yield f

def sw_to_bam(ref_path, sequence_iter, bam_dest, n_threads,
              read_group=None, n_keep=-1, max_drop=100,
              match=1, mismatch=1, gap_open=7, gap_extend=1):
    with tmpfifo(prefix='pw-to-bam', name='samtools-view-fifo') as fifo_path, \
            tempfile.NamedTemporaryFile(suffix='.fasta', prefix='pw_to_bam') as tf:
        for name, sequence in sequence_iter:
            tf.write('>{0}\n{1}\n'.format(name, sequence))
        tf.flush()
        cmd1 = ['samtools', 'view', '-o', bam_dest, '-Sb', fifo_path]
        log.info(' '.join(cmd1))
        p = subprocess.Popen(cmd1)
        sw.align(ref_path, tf.name, fifo_path, n_threads=n_threads,
                 read_group=read_group, n_keep=n_keep, max_drop=max_drop,
                 gap_open=gap_open, mismatch=mismatch, match=match,
                 gap_extend=gap_extend)
        returncode = p.wait()
        if returncode:
            raise subprocess.CalledProcessError(returncode, cmd1)

def build_parser(p):
    p.add_argument('csv_file', type=util.opener('rU'))
    p.add_argument('-d', '--delimiter', default='\t',
                   help="""Delimiter [default: tab]""")
    p.add_argument('-c', '--count-column', default='n_sources')
    p.add_argument('-j', '--threads', default=1, type=int, help="""Number of
                   threads [default: %(default)d]""")
    p.add_argument('v_bamfile', help="""Path for V alignments""")
    p.add_argument('d_bamfile', help="""Path for D alignmnets""")
    p.add_argument('j_bamfile', help="""Path for J alignments""")
    p.add_argument('-r', '--read-group')
    agrp = p.add_argument_group('Alignment options')
    agrp.add_argument('-k', '--keep', help="""Number of alignments to keep per
                      gene [default: 15V, 5J, 10D]""", type=int)
    agrp.add_argument('--max-drop', help="""Discard alignments with scores
                      VALUE below max score.""", default=5, type=int)
    agrp.add_argument('-m', '--match', default=1, type=int, help="""Match score
                      [default: %(default)d]""")
    agrp.add_argument('-u', '--mismatch', default=1, type=int, help="""Match score
                      [default: %(default)d]""")
    agrp.add_argument('-o', '--gap-open', default=7, type=int, help="""Gap
                      opening penalty [default: %(default)d]""")
    agrp.add_argument('-e', '--gap-extend', default=1, type=int, help="""Gap
                      extension penalty [default: %(default)d]""")
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
                   v_index=int_or_none(row['vIndex']),
                   j_index=int_or_none(row['jIndex']),
                   tags={TAG_COUNT: int_or_none(row[a.count_column]),
                         TAG_CDR3_LENGTH: int_or_none(row['cdr3Length']),
                         TAG_STATUS: int(row['sequenceStatus'])})
               for i, row in enumerate(r)]
        log.info('Done.')

        sequences = [(i.name, i.sequence) for i in rows]

        align = functools.partial(sw_to_bam, n_threads=a.threads,
                                  read_group=a.read_group,
                                  max_drop=a.max_drop,
                                  match=a.match,
                                  mismatch=a.mismatch,
                                  gap_open=a.gap_open,
                                  gap_extend=a.gap_extend)

        closing = contextlib.closing
        def align_gene(resource_name, outfile, keep=a.keep):
            with resource_fasta(resource_name) as fasta, \
                    tempfile.NamedTemporaryFile(suffix='.bam') as tf:
                align(fasta, sequences, tf.name, n_keep=keep)
                with closing(pysam.Samfile(tf.name, 'rb')) as in_sam, \
                        closing(pysam.Samfile(outfile, 'wb', template=in_sam)) as out_sam:
                    for read in add_tags(in_sam, rows):
                        out_sam.write(read)

        log.info('aligning V')
        align_gene('ighv_functional.fasta', a.v_bamfile, keep=a.keep or 15)
        log.info('aligning D')
        align_gene('ighd.fasta', a.d_bamfile, keep=a.keep or 10)
        log.info('aligning J')
        align_gene('ighj_adaptive.fasta', a.j_bamfile, keep=a.keep or 5)
