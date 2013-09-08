"""
Align V and J starting from an adaptive CSV file
"""
import collections
import contextlib
import csv
import functools
import logging
import os
import shutil
import subprocess
import tempfile

from pkg_resources import resource_stream

log = logging.getLogger('vdjalign')

from .. import util, sw


TAG_FRAME = 'XF'
TAG_CDR3_START = 'XB'
TAG_CDR3_END = 'XE'
TAG_COUNT = 'XC'
TAG_CDR3_LENGTH = 'XL'
TAG_JGENE = 'XJ'

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
              read_group=None, n_keep=-1):
    with tmpfifo(prefix='pw-to-bam', name='samtools-view-fifo') as fifo_path, \
            tempfile.NamedTemporaryFile(suffix='.fasta', prefix='pw_to_bam') as tf:
        for name, sequence in sequence_iter:
            tf.write('>{0}\n{1}\n'.format(name, sequence))
            tf.flush()
        cmd1 = ['samtools', 'view', '-o', bam_dest, '-Sb', fifo_path]
        log.info(' '.join(cmd1))
        p = subprocess.Popen(cmd1)
        sw.align(ref_path, tf.name, fifo_path, n_threads=n_threads,
                    read_group=read_group, n_keep=n_keep,
                    gap_open=10, mismatch=4, match=1)
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
    p.add_argument('v_bamfile')
    p.add_argument('d_bamfile')
    p.add_argument('j_bamfile')
    p.add_argument('-r', '--read-group')
    p.add_argument('--default-qual')
    p.set_defaults(func=action)

def action(a):
    with a.csv_file as ifp:
        csv_lines = (i for i in ifp if not i.startswith('#'))
        r = csv.DictReader(csv_lines, delimiter=a.delimiter)
        int_or_none = or_none(int)
        Row = collections.namedtuple('Row', ['name', 'sequence', 'v_index',
                                             'j_index', 'tags'])
        log.info('Loading sequences.')
        sequences = [Row(name=row['name'],
                         sequence=row['nucleotide'],
                         v_index=int_or_none(row['vIndex']),
                         j_index=int_or_none(row['jIndex']),
                         tags={TAG_COUNT: int_or_none(row[a.count_column]),
                               TAG_CDR3_LENGTH: int_or_none(row['cdr3Length']),
                               TAG_JGENE: row['jGeneName'] if row['jGeneName'] != '-1' else row['jTies']})
                     for i, row in enumerate(r)]
        log.info('Done.')

        #tags = {i.name: i.tags for i in sequences}
        sequences = ((i.name, i.sequence) for i in sequences)

        log.info('aligning V')
        with resource_fasta('ighv_functional.fasta') as fasta:
            sw_to_bam(fasta, sequences, a.v_bamfile, n_threads=a.threads,
                  read_group=a.read_group)
        log.info('aligning D')
        with resource_fasta('ighd.fasta') as fasta:
            sw_to_bam(fasta, sequences, a.d_bamfile, n_threads=a.threads,
                  read_group=a.read_group)
        log.info('aligning J')
        with resource_fasta('ighj_adaptive.fasta') as fasta:
            sw_to_bam(fasta, sequences, a.j_bamfile, n_threads=a.threads,
                  read_group=a.read_group)
