"""
Align V and J starting from an adaptive XsV file
"""
import contextlib
import collections
import csv
import logging
import os
import os.path
import shutil
import subprocess
import sys
import tempfile

from pkg_resources import resource_stream
import pysam

log = logging.getLogger(__name__)

from .. import util
from .adaptive_tsv import TAG_COUNT

@contextlib.contextmanager
def indexed_bam(bam_path):
    idx_path = bam_path + '.bai'
    if os.path.exists(idx_path):
        yield
    else:
        try:
            cmd = ['samtools', 'index', bam_path]
            log.info('Running "%s"', ' '.join(cmd))
            subprocess.check_call(cmd)
            yield
        finally:
            if os.path.exists(idx_path):
                os.remove(idx_path)

@contextlib.contextmanager
def indexed_fasta_of_resource(file_name):
    with util.tempdir(prefix='faidx-') as td, \
         resource_stream('vdjalign.imgt.data', file_name) as in_fasta, \
         tempfile.TemporaryFile(prefix='stderr', dir=td()) as dn:
        with open(td(file_name), 'w') as fp:
            shutil.copyfileobj(in_fasta, fp)
        cmd = ['samtools', 'faidx', file_name]
        log.info('Running "%s" in %s', ' '.join(cmd), td())
        subprocess.check_call(cmd, cwd=td(), stderr=dn)
        with contextlib.closing(pysam.Fastafile(td(file_name))) as fasta:
            yield fasta

FILE_MAP = {'v': 'ighv.fasta', 'j': 'ighj.fasta'}

def build_parser(p):
    p.add_argument('gene', choices=('v', 'j'))
    p.add_argument('bam_files', nargs='+')
    p.add_argument('-o', '--outfile', default=sys.stdout, type=util.opener('w'))
    p.set_defaults(func=action)

def action(a):
    closing = contextlib.closing

    with indexed_fasta_of_resource(FILE_MAP[a.gene]) as fasta:
        w = csv.DictWriter(a.outfile,
                ['gene', 'read_group', 'pos', 'coverage', 'wt', 'A', 'C', 'G',
                 'T', 'N', 'Ins', 'Del'],
                lineterminator='\n')
        w.writeheader()
        for bam_path in a.bam_files:
            with indexed_bam(bam_path), \
                 closing(pysam.Samfile(bam_path)) as bam:
                pileups = bam.pileup(stepper='all', fastafile=fasta)
                for pileup in pileups:
                    rows = collections.defaultdict(collections.Counter)
                    pos = pileup.pos
                    ref = bam.getrname(pileup.tid)
                    rbase = fasta.fetch(ref, pos, pos+1)
                    assert len(rbase) == 1
                    for r in pileup.pileups:
                        if r.alignment.is_unmapped or r.alignment.is_secondary:
                            continue
                        rg = r.alignment.opt('RG')
                        if r.indel > 0:
                            base = 'Ins'
                        elif r.indel < 0:
                            base = 'Del'
                        else:
                            base = r.alignment.seq[r.qpos]

                        try:
                            c = r.alignment.opt(TAG_COUNT)
                        except:
                            raise KeyError("{0} missing count".format(r.alignment.qname))

                        rows[rg][base] += c

                    for rg, counts in rows.iteritems():
                        row = {'gene': ref, 'read_group': rg, 'pos': pos, 'wt': rbase,
                                'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'Ins': 0, 'Del': 0,
                               'coverage': pileup.n}
                        row.update(counts)
                        w.writerow(row)
