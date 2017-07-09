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

    if args.glfo_dir is not None:
        print '%s --glfo-dir isn\'t getting plugged in to igblast/changeo (would need to rebuild igblast db)' % utils.color('red', 'warning')

    cmd = './igblastn'
    cmd += ' -germline_db_V human_gl_V -germline_db_D human_gl_V -germline_db_J human_gl_J'
    cmd += ' -auxiliary_data optional_file/human_gl.aux'
    cmd += ' -domain_system imgt -ig_seqtype Ig -organism human -outfmt \'7 std qseq sseq btop\''
    cmd += ' -num_threads %d' % args.n_procs
    cmd += ' -query ' + infname + ' -out ' + outfname
    
    cmd = 'cd %s; %s' % (args.igbdir, cmd)
    utils.simplerun(cmd, shell=True, print_time='igblast')

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
    cmd += ' --leave-default-germline'
    cmd += ' --presto-output --only-smith-waterman'
    cmd += ' --outfname ' + outfname
    if args.glfo_dir is not None:
        cmd += ' --initial-germline-dir ' + args.glfo_dir
    cmd += ' --aligned-germline-fname ' + aligned_germline_fname
    cmd += ' --n-procs ' + str(args.n_procs)

    utils.simplerun(cmd, print_time='partis annotation')

    os.remove(aligned_germline_fname)

# ----------------------------------------------------------------------------------------
def run_tigger(infname, outfname, outdir):
    if utils.output_exists(args, outfname, offset=8):
        return

    rcmds = ['library(tigger)', 'library(dplyr)']
    # rcmds += ['data(sample_db, germline_ighv)']

    db_name = 'annotations'
    gls_name = 'gls'
    rcmds += ['%s = read.csv("%s", sep="\t")' % (db_name, infname)]
    rcmds += ['%s = readIgFasta("%s")' % (gls_name, get_glfname('v', aligned=True))]

    tigger_outfname = outdir + '/tigger.fasta'
    rcmds += ['novel_df = findNovelAlleles(%s, %s, germline_min=2, nproc=%d)' % (db_name, gls_name, args.n_procs)]  #
    rcmds += ['geno = inferGenotype(%s, find_unmutated = FALSE, germline_db = %s, novel_df = novel_df)' % (db_name, gls_name)]
    rcmds += ['genotype_seqs = genotypeFasta(geno, %s, novel_df)' % (gls_name)]
    rcmds += ['writeFasta(genotype_seqs, "%s")' % tigger_outfname]
    cmdfname = args.workdir + '/tigger-in.cmd'
    with open(cmdfname, 'w') as cmdfile:
        cmdfile.write('\n'.join(rcmds) + '\n')
    cmdstr = 'R --slave -f ' + cmdfname
    utils.simplerun(cmdstr, shell=True, print_time='tigger')

    # post-process tigger .fa
    gldir = args.glfo_dir if args.glfo_dir is not None else 'data/germlines/human'
    glfo = glutils.read_glfo(gldir, args.locus)
    tigger_alleles = set()
    for seqfo in utils.read_fastx(tigger_outfname):
        seq = seqfo['seq'].replace(utils.gap_chars[0], '')  # it should be just dots...
        tigger_alleles.add(seqfo['name'])
        if seqfo['name'] not in glfo['seqs'][args.region]:
            newfo = {'gene' : seqfo['name'], 'seq' : seq}
            use_template_for_codon_info = False
            if '+' in newfo['gene']:
                newfo['template-gene'] = newfo['gene'].split('+')[0]
                use_template_for_codon_info = True
            glutils.add_new_allele(glfo, newfo, use_template_for_codon_info=use_template_for_codon_info, debug=True)
        elif glfo['seqs'][args.region][seqfo['name']] != seq:
            print '%s different sequences in glfo and tigger output for %s:\n    %s\n    %s' % (utils.color('red', 'error'), seqfo['name'], glfo['seqs'][args.region][seqfo['name']], seqfo['seq'])
    for gene in glfo['seqs'][args.region]:  # remove them afterwards so we can use existing ones to get codon info
        if gene not in tigger_alleles:
            glutils.remove_gene(glfo, gene)

    out_gldir = os.path.dirname(outfname).rstrip('/' + args.locus)
    assert glutils.get_fname(out_gldir, args.locus, args.region) == outfname
    glutils.write_glfo(out_gldir, glfo)

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
parser.add_argument('--locus', default='igh')
parser.add_argument('--region', default='v')
parser.add_argument('--changeo-path', default=os.getenv('HOME') + '/.local')
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
outdir = os.path.dirname(args.outfname)  # kind of annoying having <args.workdir> and <outdir>, but the former is for stuff we don't want to keep (not much...  maybe just .cmd file), and the latter is for stuff we do
assert outdir.split('/')[-1] == args.locus
outdir = outdir.rstrip('/' + args.locus)

utils.prep_dir(args.workdir, wildlings=['*.cmd', '*.fa']) #'*.fmt7'])
utils.prep_dir(outdir, allow_other_files=True)

outfname = run_alignment(args, outdir)
run_tigger(outfname, args.outfname, outdir)
