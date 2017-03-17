#!/usr/bin/env python
import os
import random
import argparse
import sys
from subprocess import check_call
sys.path.insert(1, './python')

import utils
import glutils

base_cmd = './bin/test-allele-finding.py'
fsdir = '/fh/fast/matsen_e'
baseoutdir = fsdir + '/dralph/partis/allele-finder'
locus = 'igh'
region = 'v'

# ----------------------------------------------------------------------------------------
sys.path.insert(1, './datascripts')
import heads
for study in ['kate-qrs-2016-09-09', 'laura-mb-2016-12-22', 'chaim-donor-45-2016-08-04', 'adaptive-billion-read-2016-04-07', 'vollmers-2016-04-08']:
# for study in ['adaptive-billion-read-2016-04-07']:
    names, dirs = [], []
    metafo = heads.read_metadata(study)
    print study
    for dset in metafo:
        bdir = fsdir + '/processed-data/partis/' + study + '/vz/' + dset
        if metafo[dset]['timepoint'] == 'merged':
            continue
        if dset == 'Hs-LN3-5RACE-IgG':  # bad one
            continue
        if not os.path.exists(bdir):
            print '    %s missing' % dset
            continue
        names.append(metafo[dset]['shorthand'])
        dirs.append(bdir + '/plots/sw/mute-freqs/overall')
    check_call(['./bin/compare-plotdirs.py', '--outdir', fsdir + '/dralph/partis/tmp/lots-of-mfreqs/' + study, '--plotdirs', ':'.join(dirs), '--names', ':'.join(names)])
# check_call(['./bin/compare-plotdirs.py', '--outdir', fsdir + '/dralph/partis/tmp/lots-of-mfreqs/zeroths', '--plotdirs', ':'.join(dirs), '--names', ':'.join(names)])
sys.exit()

# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print '%s %s' % (utils.color('red', 'run'), cmd_str)
    sys.stdout.flush()
    check_call(cmd_str.split())

# ----------------------------------------------------------------------------------------
def get_performance(outdir, debug=False):
    sglfo = glutils.read_glfo(outdir + '/germlines/simulation', locus=locus)
    iglfo = glutils.read_glfo(outdir + '/simu-test/sw/germline-sets', locus=locus)
    missing_alleles = set(sglfo['seqs'][region]) - set(iglfo['seqs'][region])
    spurious_alleles = set(iglfo['seqs'][region]) - set(sglfo['seqs'][region])
    if debug:
        if len(missing_alleles) > 0:
            print '    %2d  missing %s' % (len(missing_alleles), ' '.join([utils.color_gene(g) for g in missing_alleles]))
        if len(spurious_alleles) > 0:
            print '    %2d spurious %s' % (len(spurious_alleles), ' '.join([utils.color_gene(g) for g in spurious_alleles]))
        if len(missing_alleles) == 0 and len(spurious_alleles) == 0:
            print '    none missing'
    return len(missing_alleles), len(spurious_alleles)

# ----------------------------------------------------------------------------------------
def cf_nsnps(args, original_glfo):
    if args.plot:
        for nsnp in args.nsnp_list:
            missing, spurious = zip(*[get_performance(baseoutdir + '/nsnp-' + str(nsnp) + '/' + str(iproc) + '/') for iproc in range(args.n_tests)])
            assert len(set(missing + spurious) - set([0, 1])) == 0  # should only be zeroes and/or ones
            print '  missing:  %2d / %-2d = %.2f' % (missing.count(1), len(missing), float(missing.count(1)) / len(missing))
            print '  spurious: %2d / %-2d = %.2f' % (spurious.count(1), len(spurious), float(spurious.count(1)) / len(spurious))
        return
    v_gene = args.v_genes[0]
    for nsnp in args.nsnp_list:
        cmd = base_cmd + ' --n-procs 5 --n-tests ' + str(args.n_tests) + ' --n-sim-events 1000' # --slurm'
        cmd += ' --sim-v-genes ' + v_gene
        cmd += ' --inf-v-genes ' + v_gene
        cmd += ' --nsnp-list ' + str(nsnp)
        cmd += ' --outdir ' + baseoutdir + '/' + 'nsnp-' + str(nsnp)
        run(cmd)

# need to add sample size as a subdir for everybody
# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['nsnp'])
parser.add_argument('--nsnp-list', default='1')
parser.add_argument('--n-tests', type=int, default=5)
parser.add_argument('--plot', action='store_true')
parser.add_argument('--v-genes', default='IGHV4-39*01')
args = parser.parse_args()

args.nsnp_list = utils.get_arg_list(args.nsnp_list, intify=True)
args.v_genes = utils.get_arg_list(args.v_genes)

original_glfo = glutils.read_glfo('data/germlines/human', locus=locus)

# ----------------------------------------------------------------------------------------
if args.action == 'nsnp':
    cf_nsnps(args, original_glfo)
