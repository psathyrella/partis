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
import pysam

log = logging.getLogger('vdjalign')

from .. import util, sw


TAG_FRAME = 'XF'
TAG_CDR3_START = 'XB'
TAG_CDR3_END = 'XE'
TAG_COUNT = 'XC'
TAG_CDR3_LENGTH = 'XL'
TAG_JGENE = 'XJ'

@contextlib.contextmanager
def resource_fasta(name):
    with resource_stream('vdjalign.imgt.data', name) as in_fasta, \
            tempfile.NamedTemporaryFile(prefix=name) as tf:
        shutil.copyfileobj(in_fasta, tf)
        in_fasta.close()
        tf.flush()
        yield tf.name

def annotate_alignments(input_bam, tags):
    for read in input_bam:
        if read.is_reverse:
            continue
        read.tags = ([(k, v) for k, v in read.tags if k and k != 'SA'] +
                     tags[read.qname].items())
        yield read


def or_none(fn):
    @functools.wraps(fn)
    def apply_or_none(s):
        if s:
            return fn(s)
    return apply_or_none

@contextlib.contextmanager
def tmpfifo(**kwargs):
    with util.tempdir(**kwargs) as j:
        f = j('fifo')
        os.mkfifo(f)
        yield f

def sw_to_bam(ref_path, sequence_iter, bam_dest, n_threads):
    with tmpfifo(prefix='pw-to-bam') as fifo_path, \
            tempfile.NamedTemporaryFile(suffix='.fasta', prefix='pw_to_bam') as tf:
        for name, sequence in sequence_iter:
            tf.write('>{0}\n{1}\n'.format(name, sequence))
            tf.flush()
        cmd1 = ['samtools', 'view', '-o', bam_dest, '-Sb', fifo_path]
        log.info(' '.join(cmd1))
        p = subprocess.Popen(cmd1)
        sw.align(ref_path, tf.name, fifo_path, n_threads=n_threads)
        returncode = p.wait()
        if returncode:
            raise subprocess.CalledProcessError(cmd1, returncode)

def build_parser(p):
    p.add_argument('csv_file', type=util.opener('rU'))
    p.add_argument('-d', '--delimiter', default='\t',
                   help="""Delimiter [default: tab]""")
    p.add_argument('-c', '--count-column', default='n_sources')
    p.add_argument('-j', '--threads', default=1, type=int, help="""Number of
            threads [default: %(default)d]""")
    p.add_argument('v_bamfile')
    p.add_argument('j_bamfile')
    p.add_argument('-r', '--read-group')
    p.add_argument('--default-qual')
    p.set_defaults(func=action)

def action(a):
    ntf = functools.partial(tempfile.NamedTemporaryFile, suffix='.bam')
    closing = contextlib.closing

    with resource_fasta('ighv.fasta') as v_fasta, \
            resource_fasta('ighj_adaptive.fasta') as j_fasta, \
            a.csv_file as ifp, \
            ntf(prefix='v_alignments-') as v_tf, \
            ntf(prefix='j_alignments-') as j_tf:

        csv_lines = (i for i in ifp if not i.startswith('#'))
        r = csv.DictReader(csv_lines, delimiter=a.delimiter)
        int_or_none = or_none(int)
        Row = collections.namedtuple('Row', ['name', 'sequence', 'v_index',
                                             'j_index', 'tags'])
        sequences = [Row(name=row['name'],
                         sequence=row['nucleotide'],
                         v_index=int_or_none(row['vIndex']),
                         j_index=int_or_none(row['jIndex']),
                         tags={TAG_COUNT: int_or_none(row[a.count_column]),
                               TAG_CDR3_LENGTH: int_or_none(row['cdr3Length']),
                               TAG_JGENE: row['jGeneName'] if row['jGeneName'] != '-1' else row['jTies']})
                     for i, row in enumerate(r)]

        tags = {i.name: i.tags for i in sequences}
        v_sequences = ((i.name, i.sequence[:i.v_index])
                       for i in sequences if i.v_index is not None)
        j_sequences = ((i.name, i.sequence[i.j_index:])
                       for i in sequences if i.j_index is not None)

        # TODO: read group
        sw_to_bam(v_fasta, v_sequences, v_tf.name, n_threads=a.threads)

        with pysam.Samfile(v_tf.name, 'rb') as v_tmp_bam, \
                pysam.Samfile(a.v_bamfile, 'wb', template=v_tmp_bam) as v_bam:
            log.info('Identifying frame and CDR3 start')
            reads = annotate_alignments(v_tmp_bam, tags)
            for read in reads:
                if not (read.is_unmapped or read.is_secondary or read.is_reverse):
                    if a.default_qual and read.rlen:
                        read.qual = a.default_qual * read.rlen
                v_bam.write(read)

        log.info('Aligning J-region')
        sw_to_bam(j_fasta, j_sequences, j_tf.name, n_threads=a.threads)
        logging.info('Annotating J-region')
        with closing(pysam.Samfile(j_tf.name, 'rb')) as j_tmp_bam, \
                closing(pysam.Samfile(a.j_bamfile, 'wb',
                                      template=j_tmp_bam)) as j_bam:
            reads = annotate_alignments(j_tmp_bam, tags)
            for read in reads:
                if not (read.is_unmapped or read.is_secondary or read.is_reverse):
                    if a.default_qual and read.rlen:
                        read.qual = a.default_qual * read.rlen
                j_bam.write(read)
