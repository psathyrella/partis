#!/usr/bin/env python
import subprocess
import sys
import os

sys.path.insert(1, './python')
import utils

infname = os.getcwd() + '/tmp.fa'
tigger_outfname = 'ighv.fasta'
igbdir = './packages/ncbi-igblast-1.6.1/bin'
n_procs = 4
workdir = '/tmp/dkralph/tigger'
utils.prep_dir(workdir, wildlings='*.cmd')

def changeo_outfname(infname):
    assert infname.count('.') == 1
    return infname.split('.')[0] + '_db-pass.tab'

def get_glfname(region):
    return igbdir + '/human_gl_' + region.upper() + '.fa'

def run_igblast(infname, outfname):
    cmd = './igblastn'
    cmd += ' -germline_db_V human_gl_V -germline_db_D human_gl_V -germline_db_J human_gl_J'
    cmd += ' -auxiliary_data optional_file/human_gl.aux'
    cmd += ' -domain_system imgt -ig_seqtype Ig -organism human -outfmt \'7 std qseq sseq btop\''
    cmd += ' -num_threads %d' % n_procs
    cmd += ' -query ' + infname + ' -out ' + outfname
    
    cmd = 'cd %s; %s' % (igbdir, cmd)
    utils.simplerun(cmd, shell=True)

def run_changeo(infname, igblast_outfname):
    changeo_path = os.getenv('HOME') + '/.local/bin'
    glfnames = [get_glfname(r) for r in utils.regions]
    cmd = changeo_path + '/MakeDb.py igblast'
    cmd += ' -i %s -s %s -r %s --regions --scores' % (igblast_outfname, infname, ' '.join(glfnames))
    utils.simplerun(cmd)

def run_tigger(infname, outfname):
    rcmds = ['library(tigger)', 'library(dplyr)']
    # rcmds += ['data(sample_db, germline_ighv)']

    db_name = 'annotations'
    gls_name = 'gls'
    rcmds += ['%s = read.csv("%s", sep="\t")' % (db_name, changeo_outfname(infname))]
    rcmds += ['%s = readIgFasta("%s")' % (gls_name, get_glfname('v'))]

    rcmds += ['novel_df = findNovelAlleles(%s, %s, nproc=%d, germline_min=2)' % (db_name, gls_name, n_procs)]
    rcmds += ['geno = inferGenotype(%s, find_unmutated = FALSE, germline_db = %s, novel_df = novel_df)' % (db_name, gls_name)]
    rcmds += ['genotype_seqs = genotypeFasta(geno, %s, novel_df)' % (gls_name)]
    rcmds += ['writeFasta(genotype_seqs, "%s")' % outfname]
    cmdfname = workdir + '/tigger-in.cmd'
    with open(cmdfname, 'w') as cmdfile:
        cmdfile.write('\n'.join(rcmds) + '\n')
    # subprocess.check_call(['cat', cmdfname])
    cmdstr = 'R --slave -f ' + cmdfname
    utils.simplerun(cmdstr, shell=True)
    os.remove(cmdfname)

igblast_outfname = os.getcwd() + '/%s.fmt7' % os.path.basename(os.path.splitext(infname)[0])
run_igblast(infname, igblast_outfname)
run_changeo(infname, igblast_outfname)
run_tigger(infname, tigger_outfname)
