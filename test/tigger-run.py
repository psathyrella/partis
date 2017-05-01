#!/usr/bin/env python
import subprocess
import sys
import os

sys.path.insert(1, './python')
import utils

igbdir = './packages/ncbi-igblast-1.6.1/bin'
n_procs = 4

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
    cmd += ' -query ' + infname + ' -out ' + outfname
    
    cmd = 'cd %s; %s' % (igbdir, cmd)
    utils.simplerun(cmd, shell=True)

def run_changeo(infname, igblast_outfname):
    changeo_path = os.getenv('HOME') + '/.local/bin'
    glfnames = [get_glfname(r) for r in utils.regions]
    cmd = changeo_path + '/MakeDb.py igblast'
    cmd += ' -i %s -s %s -r %s --regions --scores' % (igblast_outfname, infname, ' '.join(glfnames))
    utils.simplerun(cmd)

infname = os.getcwd() + '/tmp.fa'
igblast_outfname = os.getcwd() + '/%s.fmt7' % os.path.basename(os.path.splitext(infname)[0])
# run_igblast(infnmae, igblast_outfname)
# run_changeo(infname, igblast_outfname)

    
cmd = ['library(tigger)', 'library(dplyr)']
# 'data(sample_db, germline_ighv)'
# cmd += 'novel_df = findNovelAlleles(sample_db, germline_ighv, nproc=1)'
cmd += 'annotations = read.csv("%s")' % changeo_outfname(infname)
cmd += 'germlines = readIgFasta("%s")' % get_glfname('v')
cmd += 'novel_df = findNovelAlleles(annotations, germlines, nproc=%d)' % n_procs
cmd += 'novel = selectNovel(novel_df)'
# cmd += 'geno = inferGenotype(sample_db, find_unmutated = TRUE, germline_db = germline_ighv, novel_df = novel_df)'
# 'genotype_seqs = genotypeFasta(geno, germline_ighv, novel_df)'
# 'writeFasta(genotype_seqs, "tmp.fa")'

# # ----------------------------------------------------------------------------------------
# require(devtools, quietly=TRUE)
# require(ggplot2, quietly=TRUE)
# require(dplyr, quietly=TRUE)
# load_all("alakazam")
# load_all("shm")
# load_all("tigger")
# annotations = read.csv('/home/dralph/work/partis-dev/_tmp/tigger/run-viterbi.csv')
# germlines = readIgFasta('/home/dralph/work/partis-dev/_tmp/tigger/germlines/ighv-aligned.fasta')
# # germlines = readIgFasta('/home/dralph/work/partis-dev/data/imgt/ighv-aligned.fasta')
# novel_df = findNovelAlleles(annotations, germlines, nproc=1, germline_min=2)
# novel = selectNovel(novel_df)
# glimpse(novel)
# # ----------------------------------------------------------------------------------------
# # gl = load('gl.rda')
# #data(head, gl)
# # save(headdata, file='head.rda')
# # headdata = load('head.rda')
# # ----------------------------------------------------------------------------------------
