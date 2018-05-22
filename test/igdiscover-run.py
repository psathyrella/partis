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
            if os.path.exists(fn):
                os.remove(fn)
    else:
        os.makedirs(outdir + '/db')
    for region in utils.regions:
        targetname = glutils.get_fname(args.glfo_dir, args.locus, region)
        linkname = get_igd_glsfname(outdir, region)
        if region in utils.getregions(args.locus):
            if not os.path.exists(targetname):
                raise Exception('gl file %s d.n.e.' % targetname)
            if not os.path.islink(linkname):
                subprocess.check_call(['ln', '-s', targetname, linkname])
        else:
            with open(linkname, 'w') as dummy_d_file:
                dummy_d_file.write('>%sDx-x*x\n%s\n' % (args.locus.upper(), 'aa'))

    cfgfname = outdir + '/' + os.path.basename(args.yamlfname)  # this is the .yaml in igdiscover/ (but *not* in igdiscover/work/) have to write it in the parent workdir, then cp to work/, because... meh, who cares why, just do it like this so shit works
    if os.path.exists(cfgfname):
        os.remove(cfgfname)
    with open(args.yamlfname) as cfgfile:  # whereas this is the template .yaml in partis/test/
        cfgdata = yaml.load(cfgfile)
    if True: #not args.gls_gen:
        for filtername in ['pre_germline_filter', 'germline_filter']:
            for cfgvar in ['unique_js', 'unique_cdr3s']:
                cfgdata[filtername][cfgvar] = 0
    if args.species != 'human':
        if args.species == 'macaque':
            cfgdata['species'] = 'rhesus_monkey'
        else:
            assert False
    with open(cfgfname, 'w') as cfgfile:
        yaml.dump(cfgdata, cfgfile, width=200)

    if os.path.exists(outdir + '/work'):  # sigh, it spams out too much different output, can't get away without a '-r'
        subprocess.check_call(['rm', '-r', outdir + '/work'])

# ----------------------------------------------------------------------------------------
def getpathcmd(env=None):
    cmds = ['#!/bin/bash']
    cmds += ['. %s/etc/profile.d/conda.sh' % args.condapath]  # NOTE have to update conda (using the old version in the next two lines) in order to get this to work
    cmds += ['export PATH=%s/bin:$PATH' % args.condapath]
    # cmds += ['export PYTHONNOUSERSITE=True']  # otherwise it finds the pip-installed packages in .local and breaks (see https://github.com/conda/conda/issues/448)
    return cmds

# ----------------------------------------------------------------------------------------
def update_igdiscover():
    cmds = getpathcmd('test')

    # args.env_label = 'igdiscover'
    # # install:
    # cmds += ['conda config --add channels defaults']
    # cmds += ['conda config --add channels conda-forge']
    # cmds += ['conda config --add channels bioconda']
    # cmds += ['conda create -n %s igdiscover' % args.env_label]
    # cmds += ['conda activate %s' % args.env_label]
    # # update:
    # cmds += ['conda activate %s' % args.env_label]
    # cmds += ['igdiscover --version']
    # cmds += ['conda update igdiscover']
    # cmds += ['igdiscover --version']

    # args.env_label = 'igdiscover-dev'
    # install_dir = partis_dir + '/packages'
    # # install dev version:
    # if not os.path.exists(install_dir):
    #     os.makedirs(install_dir)
    # cmds += ['cd %s' % install_dir]
    # cmds += ['git clone https://github.com/NBISweden/IgDiscover.git']
    # cmds += ['cd IgDiscover']
    # cmds += ['conda env create -n %s -f environment.yml' % args.env_label]
    # cmds += ['source activate %s' % args.env_label]
    # cmds += ['python3 -m pip install -e .']
    # # update dev version:
    # cmds += ['cd IgDiscover']
    # cmds += ['git pull']

    utils.simplerun('\n'.join(cmds) + '\n', cmdfname='/tmp/tmprun.sh', debug=True)

# ----------------------------------------------------------------------------------------
def run_igdiscover(infname, outfname, outdir):
    if utils.output_exists(args, outfname):
        return

    prepare_igdiscover_outdir(outdir)

    if args.n_random_queries is not None:
        sub_infname = outdir + '/' + os.path.basename(infname.replace(utils.getsuffix(infname), '-n-random-queries-%d%s' % (args.n_random_queries, utils.getsuffix(infname))))
        if os.path.exists(sub_infname):
            print '    --n-random-queries: leaving existing fasta for igdiscover (hopefully it has %d queries)' % args.n_random_queries
        else:
            print '    --n-random-queries: writing new fasta for igdiscover (%d queries)' % args.n_random_queries
            seqfos = utils.read_fastx(infname, n_random_queries=args.n_random_queries)
            with open(sub_infname, 'w') as sub_infile:
                for seqfo in seqfos:
                    sub_infile.write('>%s\n%s\n' % (seqfo['name'], seqfo['seq']))
        infname = sub_infname

    igdiscover_outfname = outdir + '/work/final/database/%s.fasta' % args.region.upper()

    cmds = getpathcmd()
    cmds += ['conda activate %s' % args.env_label]
    cmds += ['cd %s' % outdir]
    cmds += ['igdiscover init --db db --single-reads %s work' % infname]  # prepares to run, putting files into <outdir>
    cmds += ['cp %s work/' % os.path.basename(args.yamlfname)]
    cmds += ['cd work']
    cmds += ['igdiscover run']
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname=outdir + '/run.sh', print_time='igdiscover', debug=True)

    template_gldir = args.glfo_dir  # if args.glfo_dir is not None else 'data/germlines/ XXX human'  # can probably delete this now that --glfo-dir is required (but leaving for now, to show how it used to be in case it comes up)
    glfo = glutils.create_glfo_from_fasta(igdiscover_outfname, args.locus, args.region, template_gldir, simulation_germline_dir=args.simulation_germline_dir)
    out_gldir = os.path.dirname(outfname).rstrip('/' + args.locus)
    assert glutils.get_fname(out_gldir, args.locus, args.region) == outfname
    glutils.write_glfo(out_gldir, glfo, debug=True)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--update-igdiscover', action='store_true', help='update to the latest igdiscover version')
parser.add_argument('--gls-gen', action='store_true')
parser.add_argument('--slurm', action='store_true', help='doesn\'t do shit, it\'s just for compatibility with tigger-run')
parser.add_argument('--infname', required=True)
parser.add_argument('--outfname', required=True)
parser.add_argument('--glfo-dir', required=True)
parser.add_argument('--simulation-germline-dir')
parser.add_argument('--yamlfname', default=partis_dir + '/test/igdiscover.yaml')
parser.add_argument('--region', default='v')
parser.add_argument('--locus', default='igh')
parser.add_argument('--species', default='human')
parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--n-random-queries', type=int)
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--condapath', default=os.getenv('HOME') + '/miniconda3')
parser.add_argument('--env-label', default='igdiscover-dev')
args = parser.parse_args()

if args.update_igdiscover:
    update_igdiscover()
    sys.exit()

print '%s not using --n-procs for igdiscover yet (it uses all the available threads by default)' % utils.color('red', 'note')
if not args.gls_gen:
    print '%s can\'t really run without --gls-gen, since you\'d need to work out how to change j parameters' % utils.color('red', 'warning')

if args.slurm:
    print '    note: --slurm doesn\'t do anything'

# ----------------------------------------------------------------------------------------
outdir = os.path.dirname(args.outfname)
assert outdir.split('/')[-1] == args.locus  # otherwise will need to update things
outdir = outdir.replace('/' + args.locus, '')
run_igdiscover(args.infname, args.outfname, outdir)
