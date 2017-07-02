#!/usr/bin/env python
import time
import colored_traceback.always
import argparse
import subprocess
import sys
import os

sys.path.insert(1, './python')
import utils
import glutils

# ----------------------------------------------------------------------------------------
def get_glfname(region, aligned):  # igblast uses unaligned ones
    if aligned:
        return args.igbdir + '/gl-fastas/human_igh_' + region + '.fasta'
    else:
        return args.igbdir + '/human_gl_' + region.upper() + '.fa'

# ----------------------------------------------------------------------------------------
def run_igblast(infname, outfname):
    if utils.output_exists(args, outfname, offset=8):
        return

    cmd = './igblastn'
    cmd += ' -germline_db_V human_gl_V -germline_db_D human_gl_V -germline_db_J human_gl_J'
    cmd += ' -auxiliary_data optional_file/human_gl.aux'
    cmd += ' -domain_system imgt -ig_seqtype Ig -organism human -outfmt \'7 std qseq sseq btop\''
    cmd += ' -num_threads %d' % args.n_procs
    cmd += ' -query ' + infname + ' -out ' + outfname
    
    cmd = 'cd %s; %s' % (args.igbdir, cmd)
    # utils.simplerun(cmd, shell=True, print_time='igblast')

# ----------------------------------------------------------------------------------------
def run_changeo(infname, igblast_outfname, outfname):
    if utils.output_exists(args, outfname, offset=8):
        return

    glfnames = [get_glfname(r, aligned=True) for r in utils.regions]
    cmd = args.changeo_path + '/bin/MakeDb.py igblast'
    cmd += ' -i %s -s %s -r %s --regions --scores' % (igblast_outfname, infname, ' '.join(glfnames))
    utils.simplerun(cmd, print_time='changeo')

# ----------------------------------------------------------------------------------------
def run_partis(infname, outfname):
    if utils.output_exists(args, outfname, offset=8):
        return

    aligned_gl_seqs = {}  # keyed by seq so it's easy to check for duplicates
    for r in utils.regions:  # deduplicate before passing to partis
        for seqfo in utils.read_fastx(get_glfname(r, aligned=True)):
            if seqfo['seq'] in aligned_gl_seqs:
                continue
            aligned_gl_seqs[seqfo['seq']] = '|'.join(seqfo['infostrs'])
    aligned_germline_fname = args.workdir + '/all-aligned-gl-seqs.fa'
    with open(aligned_germline_fname, 'w') as merged_file:
        for seq, gene in aligned_gl_seqs.items():
            merged_file.write('>%s\n%s\n' % (gene, seq))

    cmd = './bin/partis cache-parameters'
    cmd += ' --infname ' + infname
    cmd += ' --dont-remove-unlikely-alleles --dont-allele-cluster --dont-find-new-alleles'
    cmd += ' --presto-output --only-smith-waterman'
    cmd += ' --outfname ' + outfname
    cmd += ' --aligned-germline-fname ' + aligned_germline_fname
    cmd += ' --n-procs ' + str(args.n_procs)

    utils.simplerun(cmd, print_time='partis annotation')

    os.remove(aligned_germline_fname)

# ----------------------------------------------------------------------------------------
def fiddle_with_gene_info(gfo):
    if '_' in gfo['name']:
        print '      ', gfo['name']
        gfo['name'] = gfo['name'].replace('_', '+')
        print '      ', gfo['name']
        _, _, _ = utils.split_gene(gfo['name'])
        try:
            _, _, _ = utils.split_gene(gfo['name'])
        except:
            print '    couldn\'t parse \'%s\'' % gfo['name']
        print '    new tigger gene: %s' % gfo['name']

# ----------------------------------------------------------------------------------------
def run_tigger(infname, outfname):
    if utils.output_exists(args, outfname, offset=8):
        return

    rcmds = ['library(tigger)', 'library(dplyr)']
    # rcmds += ['data(sample_db, germline_ighv)']

    db_name = 'annotations'
    gls_name = 'gls'
    rcmds += ['%s = read.csv("%s", sep="\t")' % (db_name, infname)]
    rcmds += ['%s = readIgFasta("%s")' % (gls_name, get_glfname('v', aligned=True))]

    rcmds += ['novel_df = findNovelAlleles(%s, %s, germline_min=2, nproc=%d)' % (db_name, gls_name, args.n_procs)]  #
    rcmds += ['geno = inferGenotype(%s, find_unmutated = FALSE, germline_db = %s, novel_df = novel_df)' % (db_name, gls_name)]
    rcmds += ['genotype_seqs = genotypeFasta(geno, %s, novel_df)' % (gls_name)]
    rcmds += ['writeFasta(genotype_seqs, "%s")' % outfname.replace('.fasta', '-tigger.fasta')]
    cmdfname = args.workdir + '/tigger-in.cmd'
    with open(cmdfname, 'w') as cmdfile:
        cmdfile.write('\n'.join(rcmds) + '\n')
    cmdstr = 'R --slave -f ' + cmdfname
    utils.simplerun(cmdstr, shell=True, print_time='tigger')

    # post-process tigger .fa (mostly to make sure we're handling the gene names properly)
    tgfo = utils.read_fastx(outfname.replace('.fasta', '-tigger.fasta'))
    with open(outfname, 'w') as outfile:
        for gfo in tgfo:
            fiddle_with_gene_info(gfo)
            outfile.write('>%s\n%s\n' % (gfo['name'], gfo['seq']))
    os.remove(cmdfname)

# ----------------------------------------------------------------------------------------
def run_alignment(args, outdir):
    infbase = utils.getprefix(os.path.basename(args.infname))
    if args.aligner == 'igblast':
        igblast_outfname = outdir + '/' + infbase + '-igblast.fmt7'
        changeo_outfname = outdir + '/' + infbase + '-igblast_db-pass.tab'
        run_igblast(args.infname, igblast_outfname)
        run_changeo(args.infname, igblast_outfname, changeo_outfname)
        return changeo_outfname
    elif args.aligner == 'partis':
        outfname = outdir + '/' + infbase + '-partis-sw-annotations.tsv'
        run_partis(args.infname, outfname)  # can't change it in changeo (I think), so may as well use the same name here
        return outfname
    else:
        assert False

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--gls-gen', action='store_true')
parser.add_argument('--infname', required=True)
parser.add_argument('--outfname', required=True)
parser.add_argument('--workdir', required=True)
parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--aligner', choices=['igblast', 'partis'], default='partis')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--igbdir', default='./packages/ncbi-igblast-1.6.1/bin')
parser.add_argument('--glfo-dir')
parser.add_argument('--changeo-path', default=os.getenv('HOME') + '/.local')
args = parser.parse_args()

if not args.gls_gen:
    print '%s can\'t really run without --gls-gen, since you\'d need to work out how to change j parameters' % utils.color('red', 'warning')
if args.glfo_dir is not None:
    print '%s not using --glfo-dir a.t.m. (would need to rebuild igblast db)' % utils.color('red', 'warning')

# ----------------------------------------------------------------------------------------
outdir = os.path.dirname(args.outfname)  # kind of annoying having <args.workdir> and <outdir>, but the former is for stuff we don't want to keep (not much...  maybe just .cmd file), and the latter is for stuff we do
utils.prep_dir(args.workdir, wildlings=['*.cmd', '*.fa']) #'*.fmt7'])
utils.prep_dir(outdir, allow_other_files=True)

outfname = run_alignment(args, outdir)
run_tigger(outfname, args.outfname)
