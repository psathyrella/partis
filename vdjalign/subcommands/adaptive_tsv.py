"""
Align V and J starting from an adaptive CSV file
"""
import collections
import contextlib
import csv
import logging
import functools
import shutil
import subprocess
import tempfile

from pkg_resources import resource_stream
import pysam

log = logging.getLogger('vdjalign')

from .. import util

BWA_OPTS = ['-k', '6', '-O', '10', '-L', '0', '-v', '2', '-T', '10', '-M',
            '-a']

TAG_FRAME = 'XF'
TAG_CDR3_START = 'XB'
TAG_CDR3_END = 'XE'
TAG_COUNT = 'XC'
TAG_CDR3_LENGTH = 'XL'
TAG_JGENE = 'XJ'


@contextlib.contextmanager
def bwa_index_of_package_imgt_fasta(file_name):
    with util.tempdir(prefix='bwa-') as td, \
            resource_stream('vdjalign.imgt.data', file_name) as in_fasta, \
            tempfile.TemporaryFile(prefix='stderr', dir=td()) as dn:
        with open(td(file_name), 'w') as fp:
            shutil.copyfileobj(in_fasta, fp)
        cmd = ['bwa', 'index', file_name]
        log.info('Running "%s" in %s', ' '.join(cmd), td())
        subprocess.check_call(cmd, cwd=td(), stderr=dn)
        yield td(file_name)


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


def bwa_to_bam(index_path, sequence_iter, bam_dest, bwa_opts=None, sam_opts=None, alg='mem'):
    cmd1 = ['bwa', alg, index_path, '-'] + (bwa_opts or [])
    cmd2 = ['samtools', 'view', '-o', bam_dest, '-Sb', '-'] + (sam_opts or [])
    log.info('BWA: %s | %s', ' '.join(cmd1), ' '.join(cmd2))
    p1 = subprocess.Popen(cmd1, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout)

    with p1.stdin as fp:
        for name, seq in sequence_iter:
            fp.write('>{0}\n{1}\n'.format(name, seq))

    exit_code2 = p2.wait()
    if exit_code2:
        raise subprocess.CalledProcessError(exit_code2, cmd2)

    exit_code1 = p1.wait()
    if exit_code1:
        raise subprocess.CalledProcessError(exit_code1, cmd1)


def align_bwa(index_path, sequence_iter, bam_dest, threads=1, rg=None):
    opts = BWA_OPTS[:]
    if threads > 1:
        opts.extend(('-t', str(threads)))
    if rg is not None:
        opts.extend(('-R', rg))
    return bwa_to_bam(index_path, sequence_iter, bam_dest,
                      bwa_opts=opts)

def build_parser(p):
    p.add_argument('csv_file', type=util.opener('rU'))
    p.add_argument('-d', '--delimiter', default='\t',
                   help="""Delimiter [default: tab]""")
    p.add_argument('-c', '--count-column', default='n_sources')
    p.add_argument('-j', '--threads', default=1, type=int, help="""Number of
            threads for BWA [default: %(default)d]""")
    p.add_argument('v_bamfile')
    p.add_argument('j_bamfile')
    p.add_argument('-r', '--read-group')
    p.add_argument('--default-qual')
    p.set_defaults(func=action)


def action(a):
    indexed = bwa_index_of_package_imgt_fasta
    ntf = functools.partial(tempfile.NamedTemporaryFile, suffix='.bam')
    closing = contextlib.closing

    with indexed('ighv.fasta') as v_index, \
         indexed('ighj_adaptive.fasta') as j_index, \
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

        align_bwa(v_index, v_sequences, v_tf.name, threads=a.threads,
                  rg=a.read_group)
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
        align_bwa(j_index, j_sequences, j_tf.name, threads=a.threads,
                    rg=a.read_group)
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
