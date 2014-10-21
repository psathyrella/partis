
"""
Align V D, J starting from a FASTQ file
"""
import logging
import shutil

from .. import imgt, util


log = logging.getLogger('vdjalign')

def fill_targets_alignment_options(p):
    p.add_argument('fasta_out', help="""Path for alignments""")
    p.add_argument('gff3_out', help="""Path for alignments""", nargs='?')

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

def comma_separated(s):
    return [i.trim() for i in s.split(',')]

def build_parser(p):
    p.add_argument('fasta_out', help="""Path for alignments""", type=util.opener('w'))
    p.add_argument('gff3_out', help="""Path for alignments""", nargs='?', type=util.opener('w'))

    tgrp = p.add_argument_group('Target options')
    tgrp.add_argument('-l', '--locus', default='IGH', choices=('IGH', 'IGK', 'IGL'),
             help="""Locus to align reads against [default: %(default)s]""")
    tgrp.add_argument('--v-subset', help="""Subset of V reads to use (e.g.,
                      functional) [default: %(default)s]""")
    tgrp.add_argument('--d-subset', help="""Subset of D reads to use
                      (currently: No options) [default: %(default)s]""")
    tgrp.add_argument('--j-subset', help="""Subset of J reads to use (e.g.,
                      adaptive) [default: %(default)s]""")
    tgrp.add_argument('-s', '--segments', default=('v', 'd', 'j'), type=list,
                      help="""Segments to output, no delimiters [default:
                      vdj]""")

    p.set_defaults(func=action)

def action(a):
    pieces = [('v', a.v_subset),
              ('d', a.d_subset),
              ('j', a.j_subset)]
    if a.locus != 'IGH':
        pieces.pop(1) # No d
    pieces = [(s, c) for s, c in pieces if s in a.segments]

    with a.fasta_out as ofp:
        for segment, collection in pieces:
            with imgt.fasta_handle(a.locus, segment, collection=collection) as ifp:
                shutil.copyfileobj(ifp, ofp)
    if a.gff3_out:
        with a.gff3_out as ofp:
            for segment, collection in pieces:
                if imgt.has_gff3(a.locus, segment, collection=collection):
                    with imgt.gff3_handle(a.locus, segment,
                                          collection=collection) as ifp:
                        shutil.copyfileobj(ifp, ofp)
