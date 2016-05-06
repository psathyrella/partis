#!/usr/bin/env python
import os
import argparse
import csv
import re
import shutil
import glob
import sys
from Bio import SeqIO
from subprocess import Popen, PIPE

current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
if not os.path.exists(current_script_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % current_script_dir
sys.path.insert(1, current_script_dir)

import utils

# ----------------------------------------------------------------------------------------
def clean_dir():
    if os.path.exists(args.dirname):
        for fname in files_to_copy + [unaligned_fname, aligned_fname, cyst_fname, error_fname]:
            if os.path.exists(args.dirname + '/' + fname):
                os.remove(args.dirname + '/' + fname)
        remaining_files = glob.glob(args.dirname + '/*')
        if len(remaining_files) > 0:
            raise Exception('unexpected files in %s: %s' % (args.dirname, ' '.join(remaining_files)))
    else:
        os.makedirs(args.dirname)

# ----------------------------------------------------------------------------------------
def write_new_alignments(aligned_germlines, all_new_genes):
    """ write an alignment for each sequence in <all_new_genes>, taking the aligned sequences from <aligned_germlines> """
    with open(args.dirname + '/' + aligned_fname, 'w') as tmpfile:
        for gene, seq in aligned_germlines.items():
            if gene not in all_new_genes:
                continue
            tmpfile.write('>%s\n%s\n' % (gene, seq.replace('-', '.').upper()))

# ----------------------------------------------------------------------------------------
def get_alignments(old_aligned_genes, genes_without_alignments, all_new_genes):
    if len(genes_without_alignments) == 0:
        print 'already have alignments for all genes'
        write_new_alignments(old_aligned_genes, all_new_genes)
        return

    print 'missing alignments for %d genes' % len(genes_without_alignments)
    old_aligned_fname = args.dirname + '/old-aligned.fasta'
    missing_fname = args.dirname + '/missing-alignments.fasta'
    msa_table_fname = args.dirname + '/msa-table.txt'
    all_fname = args.dirname + '/all.fa'
    with open(old_aligned_fname, 'w') as tmpfile:
        for gene, seq in old_aligned_genes.items():
            tmpfile.write('>%s\n%s\n' % (gene, seq.replace('.', '-')))
    with open(missing_fname, 'w') as tmpfile:
        for gene, seq in genes_without_alignments.items():
            tmpfile.write('>%s\n%s\n' % (gene, seq.replace('.', '-')))

    def run(cmd):
        print 'RUN %s' % cmd
        proc = Popen(cmd, shell=True, stderr=PIPE)
        out, err = proc.communicate()  # the commands all redirect stdout to a file
        err = err.replace('\r', '\n')
        printstrs = []
        for errstr in err.split('\n'):  # remove the stupid progress bar things
            matches = re.findall('[0-9][0-9]* / [0-9][0-9]*', errstr)
            if len(matches) == 1 and errstr.strip() == matches[0]:
                continue
            printstrs.append(errstr)
        print '        ' + '\n        '.join(printstrs)

    cmd = 'ruby bin/makemergetable.rb ' + old_aligned_fname + ' 1>' + msa_table_fname
    run(cmd)
    cmd = 'cat ' + old_aligned_fname + ' ' + missing_fname + ' >' + all_fname
    run(cmd)
    cmd = 'mafft --merge ' + msa_table_fname + ' ' + all_fname + ' >' + args.dirname + '/' + aligned_fname  # options=  # "--localpair --maxiterate 1000"
    run(cmd)

    # then rewrite aligned file with only new genes, converting to upper case and dots for gaps
    all_aligned_germlines = utils.read_germline_seqs(args.dirname, only_region='v', aligned=True)['v']
    write_new_alignments(all_aligned_germlines, all_new_genes)

    os.remove(old_aligned_fname)
    os.remove(missing_fname)
    os.remove(msa_table_fname)
    os.remove(all_fname)

# ----------------------------------------------------------------------------------------
def get_n_gaps_up_to_cpos(aligned_seq, cpos):
    """ return number of gapped positions in <aligned_seq> before <cpos> """
    ipos = 0  # position in unaligned sequence
    n_gaps_passed = 0  # number of gapped positions in the aligned sequence that we pass before getting to cpos (i.e. while ipos < cpos)
    while ipos < cpos:
        if aligned_seq[ipos + n_gaps_passed] in utils.gap_chars:
            n_gaps_passed += 1
        else:
            ipos += 1

    return n_gaps_passed

# ----------------------------------------------------------------------------------------
def get_cpos_in_alignment(aligned_seq, seq, cpos):
    """ given <cpos> in <seq>, find the cysteine's position in <aligned_seq> """
    utils.check_conserved_cysteine(seq, cpos)
    cpos_in_alignment = cpos + get_n_gaps_up_to_cpos(aligned_seq, cpos)
    utils.check_conserved_cysteine(aligned_seq, cpos_in_alignment)
    return cpos_in_alignment

# ----------------------------------------------------------------------------------------
def write_cyst_file(known_cyst_positions):
    unaligned_genes = utils.read_germline_seqs(args.dirname, only_region='v')['v']
    aligned_genes = utils.read_germline_seqs(args.dirname, only_region='v', aligned=True)['v']

    known_gene = None  # we need to find at least one gene that's in the old and the new sets, so we know how to convert cyst positions
    for gene, info in known_cyst_positions.items():
        if gene in aligned_genes:
            known_gene = gene
            break
    if known_gene is None:
        raise Exception('couldn\'t find any genes in common between %s and %s, so can\'t write new cyst position file' % (args.reference, args.dirname + '/' + aligned_fname))

    known_cpos = known_cyst_positions[known_gene]
    cpos_in_alignment = get_cpos_in_alignment(aligned_genes[known_gene], unaligned_genes[known_gene], known_cpos)
    cyst_positions = {}
    errors = []
    for gene, seq in unaligned_genes.items():
        unaligned_cpos = cpos_in_alignment - utils.count_gaps(aligned_genes[gene], istop=cpos_in_alignment)
        try:
            utils.check_conserved_cysteine(seq, unaligned_cpos, debug=True)
        except:
            print '  %s cysteine not found in %s, skipping' % (utils.color('red', 'warning'), gene)
            # print gene, unaligned_cpos
            # print seq
            # print aligned_genes[gene]
            errors.append(gene)
            continue
        cyst_positions[gene] = unaligned_cpos

    with open(args.dirname + '/' + cyst_fname, 'w') as cystfile:
        writer = csv.DictWriter(cystfile, ('gene', 'istart'))
        writer.writeheader()
        for gene, cpos in cyst_positions.items():
            writer.writerow({'gene' : gene, 'istart' : cpos})
    with open(args.dirname + '/' + error_fname, 'w') as errorfile:
        for gene in errors:
            errorfile.write('%s\n' % gene)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('ighv_fname', help='input germline v set (presumably a new one), in fasta')
parser.add_argument('--dirname', help='directory name for output (if not specified, we use <infname> with its suffix removed)')
parser.add_argument('--reference-dir', default='data/imgt', help='directory with reference/old germline sets')
args = parser.parse_args()
if args.dirname is None:
    args.dirname = os.path.os.path.splitext(args.ighv_fname)[0]

files_to_copy = ['ighd.fasta', 'ighj.fasta', 'tryp-positions.csv']
unaligned_fname = 'ighv.fasta'
aligned_fname = 'ighv-aligned.fasta'
cyst_fname = 'cyst-positions.csv'
error_fname = 'bad-genes.txt'

# ----------------------------------------------------------------------------------------
clean_dir()
shutil.copyfile(args.ighv_fname, args.dirname + '/' + unaligned_fname)

# figure out which v genes we need to align
old_aligned_genes = utils.read_germline_seqs(args.reference_dir, only_region='v', aligned=True)
all_new_genes = utils.read_germline_seqs(args.dirname, only_region='v')  # all genes in ighv_fname, not just the new ones
genes_without_alignments = {}
for gene in all_new_genes['v']:
    if gene in old_aligned_genes['v']:
        old_seq = old_aligned_genes['v'][gene]
        for gc in utils.gap_chars:
            old_seq = old_seq.replace(gc, '')
        if old_seq != all_new_genes['v'][gene]:
            raise Exception('new gene named %s has a different sequence to the gene of the same name in %s/' % (gene, args.reference_dir))
    else:
        genes_without_alignments[gene] = all_new_genes['v'][gene]

get_alignments(old_aligned_genes['v'], genes_without_alignments, all_new_genes['v'])

print 'cp %s --> %s' % (args.reference_dir, args.dirname)
for fname in files_to_copy:
    print '    %s' % fname
    shutil.copyfile(args.reference_dir + '/' + fname, args.dirname + '/' + fname)

known_cyst_positions = utils.read_codon_positions(args.reference_dir + '/cyst-positions.csv')
write_cyst_file(known_cyst_positions)
