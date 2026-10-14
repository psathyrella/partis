#!/usr/bin/env python3
""" Merge the newly-downloaded germline set with VDJbase alleles downloaded for the vanwinkle-170 cohort,
    remove duplicate sequences, and write the combined set to germlines/processed/. Also compares the two input sets.
    The idea is to get every gene that could possibly be in the changeo annotations for vanwinkle (somehow we still
    end up missing some).

    This built /fh/fast/matsen_e/data/vanwinkle-170/germlines/processed/, which merge-germline-set.py then merged into the default human set
    (issue #396). Copied from psathyrella/datascripts meta/vanwinkle-170/process-germlines.py (b05d89c), with paths made into arguments.
    NOTE it reads imgt with skip_orfs=False, i.e. it keeps imgt ORFs, and it doesn't remove genes with bad conserved codons (see issue #421).
    To reproduce processed/ (fasta files byte-for-byte, extras.csv up to row order), run split-germlines.py in this dir, and pass as
    --template-dir the default human set from before 2026, e.g.:
      mkdir -p /tmp/$USER/gl-pre2026 && git archive d0e4f57b2^ data/germlines/human | tar -x -C /tmp/$USER/gl-pre2026
      ./data/germlines/vdjbase-download/human/process-germlines.py --template-dir /tmp/$USER/gl-pre2026/data/germlines/human --outdir <outdir>
    NOTE read_glfo(..., remove_duplicates=True) appends any new groups of duplicate names to data/germlines/duplicate-names.csv.
"""
import argparse
import os
import sys

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/data/germlines/vdjbase-download/human', '')
if partis_dir not in sys.path:
    sys.path.insert(1, partis_dir)
import partis.utils as utils
import partis.glutils as glutils
import colored_traceback.always

parser = argparse.ArgumentParser()
parser.add_argument('--template-dir', required=True, help='germline set used as the template for codon positions and duplicate names (originally the then-default data/germlines/human)')
parser.add_argument('--imgt-dir', default='data/germlines/imgt-download/human')
parser.add_argument('--vdjbase-dir', default=os.path.dirname(os.path.realpath(__file__)), help='dir with <locus>/<locus>{v,d,j}.fasta, as written by split-germlines.py')
parser.add_argument('--outdir', required=True)
args = parser.parse_args()

for locus in utils.sub_loci('ig'):
    print('%s' % utils.color('red', locus))
    processed_ref_glfo = glutils.read_glfo(args.template_dir, locus)
    new_imgt_glfo = glutils.read_glfo(args.imgt_dir, locus, template_glfo=processed_ref_glfo, remove_duplicates=True, skip_orfs=False, debug=True)
    print(new_imgt_glfo['seqs']['d'].keys())
    vdj_glfo = glutils.read_glfo(args.vdjbase_dir, locus, template_glfo=new_imgt_glfo, remove_duplicates=True, debug=True)
    print(utils.color('blue', 'merge'))
    merged_glfo, _ = glutils.get_merged_glfo(new_imgt_glfo, vdj_glfo, debug=True)
    print(utils.color('blue', 'write'))
    glutils.write_glfo(args.outdir, merged_glfo, debug=True)
    print(utils.color('blue', 'read'))
    merged_glfo = glutils.read_glfo(args.outdir, locus, debug=True)
    print(utils.color('blue', 'comparing %s' % locus))
    glutils.compare_glfos([new_imgt_glfo, vdj_glfo], ['default', 'vdjbase'], locus)
