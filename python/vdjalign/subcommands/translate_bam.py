"""
Translate a BAM file
"""
import contextlib
import logging
import sys

import pysam

from .. import codons, util
from .adaptive_tsv import TAG_CDR3_START

log = logging.getLogger(__name__)

def build_parser(p):
    p.add_argument('bam_files', nargs='+')
    p.add_argument('-o', '--outfile', default=sys.stdout, type=util.opener('w'))
    p.set_defaults(func=action)

def action(a):
    with a.outfile as ofp:
        for bam_path in a.bam_files:
            with contextlib.closing(pysam.Samfile(bam_path)) as bam:
                for read in bam:
                    if read.is_secondary or read.is_unmapped:
                        continue
                    try:
                        frame = read.opt('XF')
                    except KeyError:
                        continue

                    r = read.seq[frame - 1:]
                    ofp.write('>{0}\n{1}\n'.format(read.qname,
                                                   codons.translate(r)))
                    cysteine_pos = (read.opt(TAG_CDR3_START) - (frame - 1)) // 3
                    ofp.write('>{0}_ann\n{1}\n'.format(read.qname, ('-' * cysteine_pos) + '*'))
                    ofp.write('>{0}_nt\n{1}\n'.format(read.qname, r))
