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
def run_tigger(infname, outfname):
    if utils.output_exists(args, outfname):
        return

    cmd = 'export PATH=%s:$PATH' % args.condapath
    cmd += 'igdiscover init --db %s/db --reads %s discovertest' % (args.workdir, args.infname)
    'cp igdiscover-testdata-*/igdiscover.yaml discovertest/'
    'cd discovertest && igdiscover run'

    cmd = './igblastn'
    cmd += ' -germline_db_V human_gl_V -germline_db_D human_gl_V -germline_db_J human_gl_J'
    cmd += ' -auxiliary_data optional_file/human_gl.aux'
    cmd += ' -domain_system imgt -ig_seqtype Ig -organism human -outfmt \'7 std qseq sseq btop\''
    cmd += ' -num_threads %d' % args.n_procs
    cmd += ' -query ' + infname + ' -out ' + outfname
    
    cmd = 'cd %s; %s' % (args.igbdir, cmd)
    utils.simplerun(cmd, shell=True, print_time='igblast')

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--infname', required=True)
parser.add_argument('--outfname', required=True)
parser.add_argument('--workdir', required=True)
# parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--condapath', default=os.getenv('HOME') + '/miniconda3/bin')
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
outdir = os.path.dirname(args.outfname)  # kind of annoying having <args.workdir> and <outdir>, but the former is for stuff we don't want to keep (not much...  maybe just .cmd file), and the latter is for stuff we do
# '/tmp/igdiscover-testdata-b4d969a8f58a/'
# infbase = utils.getprefix(os.path.basename(args.infname))
# igblast_outfname = outdir + '/' + infbase + '-igblast.fmt7'
# changeo_outfname = outdir + '/' + infbase + '-igblast_db-pass.tab'
# utils.prep_dir(args.workdir, wildlings=['*.cmd']) #'*.fmt7'])
# utils.prep_dir(outdir, allow_other_files=True)

# run_igblast(args.infname, igblast_outfname)
# run_changeo(args.infname, igblast_outfname, changeo_outfname)
# run_tigger(changeo_outfname, args.outfname)
