"""
Align V D, J starting from a FASTQ file
"""
import functools
import logging
import subprocess

from .. import imgt, util, sw

log = logging.getLogger('vdjalign')


def sw_to_bam(ref_path, sequence_path, bam_dest, n_threads,
              read_group=None, extra_ref_paths=[],
              match=1, mismatch=1, gap_open=7, gap_extend=1,
              max_drop=0, min_score=0, bandwidth=150, samtools_dir=None):
    with util.tmpfifo(prefix='pw-to-bam', name='samtools-view-fifo') as fifo_path:
        samtool_binary = 'samtools'
        if samtools_dir is not None:
            samtool_binary = samtools_dir + '/' + samtool_binary
        cmd1 = [samtool_binary, 'view', '-@', str(n_threads), '-o', bam_dest, '-Sb', fifo_path]
        log.info(' '.join(cmd1))
        p = subprocess.Popen(cmd1)
        sw.ig_align(ref_path, sequence_path, fifo_path, n_threads=n_threads,
                    read_group=read_group, extra_ref_paths=extra_ref_paths,
                    gap_open=gap_open, mismatch=mismatch, match=match,
                    gap_extend=gap_extend, max_drop=max_drop,
                    bandwidth=bandwidth, min_score=min_score)
        returncode = p.wait()
        if returncode:
            raise subprocess.CalledProcessError(returncode, cmd1)

def fill_targets_alignment_options(p):
    p.add_argument('-j', '--threads', default=1, type=int, help="""Number of
                   threads [default: %(default)d]""")
    p.add_argument('out_bamfile', help="""Destination path""")
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
    tgrp.add_argument('--vdj-dir', help="""Directory from which to read
                      germline genes) [default: %(default)s]""")
    tgrp.add_argument('--samtools-dir', help="""Path to samtools binary.""")

    agrp = p.add_argument_group('Alignment options')
    agrp.add_argument('-m', '--match', default=1, type=int, help="""Match score
                      [default: %(default)d]""")
    agrp.add_argument('-u', '--mismatch', default=1, type=int, help="""Match score
                      [default: %(default)d]""")
    agrp.add_argument('-o', '--gap-open', default=7, type=int, help="""Gap
                      opening penalty [default: %(default)d]""")
    agrp.add_argument('-e', '--gap-extend', default=1, type=int, help="""Gap
                      extension penalty [default: %(default)d]""")
    agrp.add_argument('--max-drop', default=0, type=int, help="""Only output
                      alignments with scores within MAX_DROP of the best
                      scoring alignment for each segment [default:
                      %(default)d]""")
    agrp.add_argument('--min-score', default=0, type=int, help="""Minimum (V)
                      alignment score [default: %(default)d]""")
    agrp.add_argument('--bandwidth', default=150, type=int, help="""Bandwidth
                      for global alignment [default: %(default)d]""")

def build_parser(p):
    p.add_argument('fastq')
    fill_targets_alignment_options(p)
    p.set_defaults(func=action)

def action(a):
    align = functools.partial(sw_to_bam, n_threads=a.threads,
                              read_group=a.read_group,
                              match=a.match,
                              mismatch=a.mismatch,
                              min_score=a.min_score,
                              gap_open=a.gap_open,
                              gap_extend=a.gap_extend,
                              bandwidth=a.bandwidth,
                              max_drop=a.max_drop,
                              samtools_dir=a.samtools_dir)

    log.info('aligning')
    with imgt.temp_fasta(a.locus, 'v', a.v_subset, vdj_dir=a.vdj_dir) as vf, \
            imgt.temp_fasta(a.locus, 'j', a.j_subset, vdj_dir=a.vdj_dir) as jf, \
            util.with_if(a.locus == 'IGH', imgt.temp_fasta, a.locus, 'd',
                    collection=a.d_subset, vdj_dir=a.vdj_dir) as df:
        ex_refs = [i for i in [df, jf] if i is not None]
        align(vf, a.fastq, a.out_bamfile, extra_ref_paths=ex_refs)
