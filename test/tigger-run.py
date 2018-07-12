#!/usr/bin/env python
import numpy
import random
import glob
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
    return '%s/%s/%s%s%s.fasta' % (args.igbdir, args.locus, args.locus, region, '' if aligned else '-unaligned')

# ----------------------------------------------------------------------------------------
def initialize_germline_info(outdir):
    if len(glob.glob('%s/%s/%s*-unaligned.fasta.nsq' % (args.igbdir, args.locus, args.locus))) == len(utils.regions):
        return

    cmds = ['#!/bin/bash']
    cmds += ['export PATH=%s:$PATH' % args.condapath]
    # cmds += ['igblastn -help']
    cmds += ['cd %s/%s' % (args.igbdir, args.locus)]
    for tmpreg in utils.regions:  # need a dummy d for igk/igl
        cmds += ['perl ../edit_imgt_file.pl %s%s.fasta > %s%s-unaligned.fasta' % (args.locus, tmpreg, args.locus, tmpreg)]
        cmds += ['makeblastdb -parse_seqids -dbtype nucl -in %s%s-unaligned.fasta' % (args.locus, tmpreg)]
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname=outdir + '/run.sh')

# ----------------------------------------------------------------------------------------
def run_igblast(infname, outfname):
    if utils.output_exists(args, outfname, offset=8):
        return

    if args.glfo_dir is not None:
        print '%s --glfo-dir isn\'t getting plugged in to igblast/changeo (would need to rebuild igblast db)' % utils.color('red', 'warning')

    if args.n_random_queries is not None:
        sub_infname = os.path.dirname(outfname) + '/' + os.path.basename(infname.replace(utils.getsuffix(infname), '-n-random-queries-%d%s' % (args.n_random_queries, utils.getsuffix(infname))))
        if os.path.exists(sub_infname):
            print '    --n-random-queries: leaving existing fasta for igblast (hopefully it has %d queries)' % args.n_random_queries
        else:
            print '    --n-random-queries: writing new fasta for igblast (%d queries)' % args.n_random_queries
            seqfos = utils.read_fastx(infname, n_random_queries=args.n_random_queries)
            with open(sub_infname, 'w') as sub_infile:
                for seqfo in seqfos:
                    sub_infile.write('>%s\n%s\n' % (seqfo['name'], seqfo['seq']))
        infname = sub_infname

    cmds = ['#!/bin/bash']
    cmds += ['cd %s/%s' % (args.igbdir, args.locus)]
    cmds += ['export PATH=%s:$PATH' % args.condapath]
    cmds += ['igblastn']
    for tmpreg in utils.regions:
        cmds[-1] += ' -germline_db_%s %s%s-unaligned.fasta' % (tmpreg.upper(), args.locus, tmpreg)
    cmds[-1] += ' -auxiliary_data optional_file/%s_gl.aux' % args.species
    cmds[-1] += ' -domain_system imgt -ig_seqtype Ig -organism %s -outfmt \'7 std qseq sseq btop\'' % args.species
    cmds[-1] += ' -num_threads %d' % utils.auto_n_procs()
    cmds[-1] += ' -query ' + infname + ' -out ' + outfname
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname=args.workdir + '/run.sh')

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
    cmd += ' --dont-write-parameters'
    cmd += ' --parameter-dir ' + os.path.dirname(outfname) + '/partis-parameters'  # sw cache file gets written to the parameter dir even with --dont-write-parameters
    if args.glfo_dir is not None:
        cmd += ' --initial-germline-dir ' + args.glfo_dir
    cmd += ' --aligned-germline-fname ' + aligned_germline_fname
    cmd += ' --n-procs ' + str(args.n_procs)
    if args.n_random_queries is not None:
        cmd += ' --n-random-queries %d' % args.n_random_queries
    if args.slurm:
        cmd += ' --batch-system slurm'
    cmd += ' --locus ' + args.locus

    utils.simplerun(cmd, print_time='partis annotation')

    os.remove(aligned_germline_fname)

# ----------------------------------------------------------------------------------------
def run_tigger(infname, outfname, outdir):
    if utils.output_exists(args, outfname, offset=8):
        return

    rcmds = ['library(ggplot2)', 'library(tigger, warn.conflicts=FALSE)', 'library(dplyr, warn.conflicts=FALSE)']
    # rcmds += ['data(sample_db, germline_ighv)']

    db_name = 'annotations'
    gls_name = 'gls'
    rcmds += ['%s = read.csv("%s", sep="\t")' % (db_name, infname)]
    rcmds += ['%s = readIgFasta("%s")' % (gls_name, get_glfname('v', aligned=True))]

    tigger_outfname = outdir + '/tigger.fasta'
    find_novel_argstr = '%s, %s, nproc=%d' % (db_name, gls_name, utils.auto_n_procs())
    if args.tuned_tigger_params:
        germline_min = 5 # only analyze genes which correspond to at least this many V calls (default 200)
        min_seqs = 5  # minimum number of total sequences
        j_max = 0.95  # of sequences which align perfectly (i.e. zero mutation?) to a new allele, no more than this fraction can correspond to each junction length + j gene combination (default 0.15)
        find_novel_argstr += ', germline_min=%d, min_seqs=%d, j_max=%f' % (germline_min, min_seqs, j_max)
    rcmds += ['novel_df = findNovelAlleles(%s)' % find_novel_argstr]
    # rcmds += ['sessionInfo()']
    rcmds += ['print(novel_df)']
    rcmds += ['geno = inferGenotype(%s, find_unmutated = TRUE, germline_db = %s, novel_df = novel_df)' % (db_name, gls_name)]
    rcmds += ['genotype_seqs = genotypeFasta(geno, %s, novel_df)' % (gls_name)]
    rcmds += ['writeFasta(genotype_seqs, "%s")' % tigger_outfname]
    cmdfname = args.workdir + '/tigger-in.cmd'
    with open(cmdfname, 'w') as cmdfile:
        cmdfile.write('\n'.join(rcmds) + '\n')
    cmdstr = 'R --slave -f ' + cmdfname

    cmdfo = {'cmd_str' : cmdstr, 'logdir' : args.workdir, 'env' : os.environ}
    # TODO maybe switch to utils.run_r()?
    proc = utils.run_cmd(cmdfo)
    while proc.poll() is None:
        time.sleep(0.01)
    if proc.returncode != 0:  # damn thing crashes if it thinks the sample size is small
        with open(args.workdir + '/err') as ferr:
            errstr = ''.join(ferr.readlines())
        if 'Not enough sample sequences were assigned to any germline' in errstr:
            with open(tigger_outfname, 'w') as dummy_outfasta:
                dummy_outfasta.write('')
        else:
            subprocess.check_call(['cat', args.workdir + '/out'])
            subprocess.check_call(['cat', args.workdir + '/err'])
            sys.exit(proc.returncode)

    for oe in ['err', 'out']:
        with open(args.workdir + '/' + oe) as oefile:
            print ''.join(oefile.readlines())
        os.remove(args.workdir + '/' + oe)

    # post-process tigger .fa
    template_gldir = args.glfo_dir if args.glfo_dir is not None else 'data/germlines/human'
    glfo = glutils.create_glfo_from_fasta(tigger_outfname, args.locus, args.region, template_gldir, simulation_germline_dir=args.simulation_germline_dir)
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
def install():
    rcmds = ['install.packages("tigger", repos="http://cran.rstudio.com/")']
    workdir = '/tmp/%s/%d' % (os.getenv('USER'), random.randint(0,999999))
    os.makedirs(workdir)
    utils.run_r(rcmds, workdir)
    os.rmdir(workdir)

# install()
# sys.exit()

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--gls-gen', action='store_true')
parser.add_argument('--tuned-tigger-params', action='store_true')
parser.add_argument('--infname', required=True)
parser.add_argument('--outfname', required=True)
parser.add_argument('--workdir', required=True)
parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--n-random-queries', type=int)
parser.add_argument('--aligner', choices=['igblast', 'partis'], default='partis')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--igbdir', default='./packages/ncbi-igblast-1.6.1/bin')
parser.add_argument('--glfo-dir')
parser.add_argument('--simulation-germline-dir')
parser.add_argument('--locus', default='igh')
parser.add_argument('--region', default='v')
parser.add_argument('--species', default='human')
parser.add_argument('--slurm', action='store_true')
parser.add_argument('--changeo-path', default=os.getenv('HOME') + '/.local')
parser.add_argument('--condapath', default=os.getenv('HOME') + '/miniconda3/bin')
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
outdir = os.path.dirname(args.outfname)  # kind of annoying having <args.workdir> and <outdir>, but the former is for stuff we don't want to keep (not much...  maybe just .cmd file), and the latter is for stuff we do
assert outdir.split('/')[-1] == args.locus
outdir = outdir.rstrip('/' + args.locus)

utils.prep_dir(args.workdir, wildlings=['*.cmd', '*.fa', '*.sh']) #'*.fmt7'])
utils.prep_dir(outdir, allow_other_files=True)

initialize_germline_info(outdir)  # deal with igblast germline crap
outfname = run_alignment(args, outdir)  # get the alignments, either with igblast or partis
run_tigger(outfname, args.outfname, outdir)

os.rmdir(args.workdir)
