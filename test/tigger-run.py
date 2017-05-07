#!/usr/bin/env python
import time
import colored_traceback.always
import argparse
import subprocess
import sys
import os

sys.path.insert(1, './python')
import utils

# ----------------------------------------------------------------------------------------
def get_glfname(region, aligned):  # igblast uses unaligned ones
    if aligned:
        return args.igbdir + '/gl-fastas/human_igh_' + region + '.fasta'
    else:
        return args.igbdir + '/human_gl_' + region.upper() + '.fa'

# ----------------------------------------------------------------------------------------
def run_igblast(infname, outfname):
    if utils.output_exists(args, outfname):
        return

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
    if utils.output_exists(args, outfname):
        return

    changeo_path = os.getenv('HOME') + '/.local/bin'
    glfnames = [get_glfname(r, aligned=True) for r in utils.regions]
    cmd = changeo_path + '/MakeDb.py igblast'
    cmd += ' -i %s -s %s -r %s --regions --scores' % (igblast_outfname, infname, ' '.join(glfnames))
    utils.simplerun(cmd, print_time='changeo')

# ----------------------------------------------------------------------------------------
def run_tigger(infname, outfname):
    if utils.output_exists(args, outfname):
        return

    rcmds = ['library(tigger)', 'library(dplyr)']
    # rcmds += ['data(sample_db, germline_ighv)']

    db_name = 'annotations'
    gls_name = 'gls'
    rcmds += ['%s = read.csv("%s", sep="\t")' % (db_name, infname)]
    rcmds += ['%s = readIgFasta("%s")' % (gls_name, get_glfname('v', aligned=True))]

    rcmds += ['novel_df = findNovelAlleles(%s, %s, nproc=%d)' % (db_name, gls_name, args.n_procs)]  # , germline_min=2
    rcmds += ['geno = inferGenotype(%s, find_unmutated = FALSE, germline_db = %s, novel_df = novel_df)' % (db_name, gls_name)]
    rcmds += ['genotype_seqs = genotypeFasta(geno, %s, novel_df)' % (gls_name)]
    rcmds += ['writeFasta(genotype_seqs, "%s")' % outfname]
    cmdfname = args.workdir + '/tigger-in.cmd'
    with open(cmdfname, 'w') as cmdfile:
        cmdfile.write('\n'.join(rcmds) + '\n')
    # subprocess.check_call(['cat', cmdfname])
    cmdstr = 'R --slave -f ' + cmdfname
    utils.simplerun(cmdstr, shell=True, print_time='tigger')
    os.remove(cmdfname)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--gls-gen', action='store_true')
parser.add_argument('--infname', required=True)
parser.add_argument('--outfname', required=True)
parser.add_argument('--workdir', required=True)
parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--igbdir', default='./packages/ncbi-igblast-1.6.1/bin')
args = parser.parse_args()
if not args.gls_gen:
    print '%s can\'t really run without --gls-gen, since you\'d need to work out how to change j parameters' % utils.color('red', 'warning')

# ----------------------------------------------------------------------------------------
outdir = os.path.dirname(args.outfname)  # kind of annoying having <args.workdir> and <outdir>, but the former is for stuff we don't want to keep (not much...  maybe just .cmd file), and the latter is for stuff we do
infbase = utils.getprefix(os.path.basename(args.infname))
igblast_outfname = outdir + '/' + infbase + '-igblast.fmt7'
changeo_outfname = outdir + '/' + infbase + '-igblast_db-pass.tab'
utils.prep_dir(args.workdir, wildlings=['*.cmd']) #'*.fmt7'])
utils.prep_dir(outdir, allow_other_files=True)

run_igblast(args.infname, igblast_outfname)
run_changeo(args.infname, igblast_outfname, changeo_outfname)
run_tigger(changeo_outfname, args.outfname)
