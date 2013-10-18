"""
Align V D, J starting from a FASTQ file
"""
import contextlib
import functools
import logging
import subprocess

from .. import imgt, util, sw

log = logging.getLogger('vdjalign')


def sw_to_bam(ref_path, sequence_path, bam_dest, n_threads,
              read_group=None, extra_ref_paths=[],
              match=1, mismatch=1, gap_open=7, gap_extend=1,
              max_drop=0):
    with util.tmpfifo(prefix='pw-to-bam', name='samtools-view-fifo') as fifo_path:
        cmd1 = ['samtools', 'view', '-@', str(n_threads), '-o', bam_dest, '-Sb', fifo_path]
        log.info(' '.join(cmd1))
        p = subprocess.Popen(cmd1)
        sw.ig_align(ref_path, sequence_path, fifo_path, n_threads=n_threads,
                    read_group=read_group, extra_ref_paths=extra_ref_paths,
                    gap_open=gap_open, mismatch=mismatch, match=match,
                    gap_extend=gap_extend, max_drop=max_drop)
        returncode = p.wait()
        if returncode:
            raise subprocess.CalledProcessError(returncode, cmd1)

def fill_targets_alignment_options(p):
    p.add_argument('-j', '--threads', default=1, type=int, help="""Number of
                   threads [default: %(default)d]""")
    p.add_argument('out_bamfile', help="""Path for alignments""")
    p.add_argument('-r', '--read-group', help="""read group header line such as
                   '@RG\tID:foo\tSM:bar'""")

    tgrp = p.add_argument_group('Target options')
    tgrp.add_argument('-l', '--locus', default='IGH', choices=('IGH', 'IGK', 'IGL'),
             help="""Locus to align reads against [default: %(default)s]""")
    tgrp.add_argument('--v-subset', help="""Subset of V reads to use (e.g.,
                      functional) [default: %(default)s]""")
    tgrp.add_argument('--d-subset', help="""Subset of D reads to use
                      (currently: No options) [default: %(default)s]""")
    tgrp.add_argument('--j-subset', help="""Subset of J reads to use (e.g.,
                      adaptive) [default: %(default)s]""")

    agrp = p.add_argument_group('Alignment options')
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

def build_parser(p):
    p.add_argument('fastq')
    fill_targets_alignment_options(p)
    p.set_defaults(func=action)

def action(a):
    align = functools.partial(sw_to_bam, n_threads=a.threads,
                              read_group=a.read_group,
                              match=a.match,
                              mismatch=a.mismatch,
                              gap_open=a.gap_open,
                              gap_extend=a.gap_extend,
                              max_drop=a.max_drop)

    log.info('aligning')
    with imgt.temp_fasta(a.locus, 'v', a.v_subset) as vf, \
            imgt.temp_fasta(a.locus, 'j', a.j_subset) as jf, \
            util.with_if(a.locus == 'IGH', imgt.temp_fasta, a.locus, 'd',
                    collection=a.d_subset) as df:
        ex_refs = [i for i in [df, jf] if i is not None]
        align(vf, a.fastq, a.out_bamfile, extra_ref_paths=ex_refs)
