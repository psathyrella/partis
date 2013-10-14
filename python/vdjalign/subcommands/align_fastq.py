"""
Align V D, J starting from a FASTQ file
"""
import contextlib
import functools
import logging
import shutil
import subprocess

from pkg_resources import resource_stream

from .. import util, sw

log = logging.getLogger('vdjalign')

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

def sw_to_bam(ref_path, sequence_path, bam_dest, n_threads,
              read_group=None, extra_ref_paths=[],
              match=1, mismatch=1, gap_open=7, gap_extend=1):
    with util.tmpfifo(prefix='pw-to-bam', name='samtools-view-fifo') as fifo_path:
        cmd1 = ['samtools', 'view', '-@', n_threads, '-o', bam_dest, '-Sb', fifo_path]
        log.info(' '.join(cmd1))
        p = subprocess.Popen(cmd1)
        sw.ig_align(ref_path, sequence_path, fifo_path, n_threads=n_threads,
                    read_group=read_group, extra_ref_paths=extra_ref_paths,
                    gap_open=gap_open, mismatch=mismatch, match=match,
                    gap_extend=gap_extend)
        returncode = p.wait()
        if returncode:
            raise subprocess.CalledProcessError(returncode, cmd1)

def build_parser(p):
    p.add_argument('fastq')
    p.add_argument('-j', '--threads', default=1, type=int, help="""Number of
                   threads [default: %(default)d]""")
    p.add_argument('out_bamfile', help="""Path for alignments""")
    p.add_argument('-r', '--read-group')
    agrp = p.add_argument_group('Alignment options')
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
    align = functools.partial(sw_to_bam, n_threads=a.threads,
                              read_group=a.read_group,
                              match=a.match,
                              mismatch=a.mismatch,
                              gap_open=a.gap_open,
                              gap_extend=a.gap_extend)

    log.info('aligning')
    with resource_fasta('ighv_functional.fasta') as vf, \
            resource_fasta('ighj_adaptive.fasta') as jf, \
            resource_fasta('ighd.fasta') as df:
        align(vf, a.fastq, a.out_bamfile, extra_ref_paths=[df, jf])
