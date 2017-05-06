#!/usr/bin/env python
import yaml
import time
import colored_traceback.always
import argparse
import subprocess
import sys
import os

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/test', '')
sys.path.insert(1, './python')
import utils
import glutils

# ----------------------------------------------------------------------------------------
def get_igd_glsfname(outdir, region):
    return outdir + '/db/' + region.upper() + '.fasta'

# ----------------------------------------------------------------------------------------
def prepare_igdiscover_outdir(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if os.path.exists(outdir + '/db'):
        for fn in [get_igd_glsfname(outdir, r) for r in utils.regions]:
            print fn
            print os.path.exists(fn)
            if os.path.exists(fn):
                os.remove(fn)
    else:
        os.makedirs(outdir + '/db')
    for region in utils.regions:
        subprocess.check_call(['ln', '-s', glutils.get_fname(args.glfo_dir, args.locus, region), get_igd_glsfname(outdir, region)])

    cfgfname = outdir + '/' + os.path.basename(args.yamlfname)  # this is the .yaml in igdiscover/ (but *not* in igdiscover/work/) have to write it in the parent workdir, then cp to work/, because... meh, who cares why, just do it like this so shit works
    if os.path.exists(cfgfname):
        os.remove(cfgfname)
    with open(args.yamlfname) as cfgfile:  # whereas this is the template .yaml in partis/test/
        cfgdata = yaml.load(cfgfile)
    if not args.gls_gen:
        for filtername in ['pre_germline_filter', 'germline_filter']:
            for cfgvar in ['unique_js', 'unique_cdr3s']:
                cfgdata[filtername][cfgvar] = 0
    with open(cfgfname, 'w') as cfgfile:
        yaml.dump(cfgdata, cfgfile, width=200)

    if os.path.exists(outdir + '/work'):  # sigh, it spams out too much different output, can't get away without a -r
        subprocess.check_call(['rm', '-rv', outdir + '/work'])

# ----------------------------------------------------------------------------------------
def run_igdiscover(infname, outfname, outdir):
    if utils.output_exists(args, outfname):
        return

    cmds = ['#!/bin/bash']
    cmds += ['export PATH=%s:$PATH' % args.condapath]
    cmds += ['export PYTHONNOUSERSITE=True']  # otherwise it finds the pip-installed packages in .local and breaks (see https://github.com/conda/conda/issues/448)
    cmds += ['cd %s' % outdir]
    cmds += ['igdiscover init --db db --single-reads %s work' % args.infname]  # prepares to run, putting files into <outdir>
    cmds += ['cp %s work/' % os.path.basename(args.yamlfname)]
    cmds += ['cd work']
    cmds += ['igdiscover run']
    cmdfname = outdir + '/run.sh'
    with open(cmdfname, 'w') as cmdfile:
        for cmd in cmds:
            cmdfile.write(cmd + '\n')
    subprocess.check_call(['chmod', '+x', cmdfname])
    cmdfos = [{'cmd_str' : cmdfname,
               'workdir' : outdir,
               'outfname' : outdir + '/work/final/%s_usage.tab' % 'v'.upper()}]
    # utils.simplerun(cmdfname, shell=True, print_time='igdiscover')
    utils.prepare_cmds(cmdfos)
    start = time.time()
    utils.run_cmds(cmdfos, ignore_stderr=True)
    print '      igdiscover time: %.1f' % (time.time()-start)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--gls-gen', action='store_true')
parser.add_argument('--infname', required=True)
parser.add_argument('--outfname', required=True)
parser.add_argument('--glfo-dir', required=True)
parser.add_argument('--yamlfname', default=partis_dir + '/test/igdiscover.yaml')
parser.add_argument('--locus', default='igh')
parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--condapath', default=os.getenv('HOME') + '/miniconda3/bin')
args = parser.parse_args()

print '%s not using --n-procs for igdiscover yet' % utils.color('red', 'note')

# ----------------------------------------------------------------------------------------
outdir = os.path.dirname(args.outfname)
assert outdir.split('/')[-1] == args.locus  # otherwise will need to update things
outdir = outdir.replace('/igh', '')
prepare_igdiscover_outdir(outdir)
run_igdiscover(args.infname, args.outfname, outdir)
