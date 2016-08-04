"""
Collapse identical reads, adding a count to the first record
"""

from vdjalign.dedup_fq import deduplicate

def build_parser(p):
    p.add_argument('infile', help="""Input FAST[AQ]""")
    p.add_argument('outfile', help="""Output [will be gzipped]""")
    p.set_defaults(func=action)

def action(a):
    deduplicate(a.infile, a.outfile)
