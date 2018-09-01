#!/usr/bin/env python
import os
import sys

current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
sys.path.insert(1, current_script_dir)
import utils

base_outdir = '_output'
extrastr = 'xxx'
ete_path = '/home/' + os.getenv('USER') + '/anaconda_ete/bin'
bcr_phylo_path = os.getenv('PWD') + '/packages/bcr-phylo-benchmark'

# ----------------------------------------------------------------------------------------
def simdir(stype):
  return '%s/bcr-phylo-%s' % (base_outdir, stype)
def plotdir(stype):
  return '%s/partis/bcr-phylo-%s' % (os.getenv('fs'), stype)
def simfname(stype):
  return '%s/%s.fasta' % (simdir(stype), extrastr)

# ----------------------------------------------------------------------------------------
def simulate(stype):
  os.environ['TMPDIR'] = '/tmp'
  cmd = 'export PATH=%s:$PATH && xvfb-run -a python %s/bin/simulator.py' % (ete_path, bcr_phylo_path)
  cmd += ' --mutability %s/S5F/Mutability.csv --substitution %s/S5F/Substitution.csv' % (bcr_phylo_path, bcr_phylo_path)

  if stype == 'neutral':
    cmd += ' --lambda %f --lambda0 %f --N %d' % (1.5, 0.365, 1000)
  elif stype == 'selection':
    cmd += ' --selection'
    cmd += ' --lambda %f --lambda0 %f' % (2., 0.365)
    # cmd += ' --n %d' % 60
    cmd += ' --T %d --target_dist 5 --target_count 100' % 100
    cmd += ' --carry_cap 1000 --skip_update 100'
  else:
    assert False

  cmd += ' --outbase %s/%s' % (simdir(stype), extrastr)
  cmd += ' --random_seq %s/sequence_data/AbPair_naive_seqs.fa' % bcr_phylo_path
  if not os.path.exists(simdir(stype)):
    os.makedirs(simdir(stype))

  utils.simplerun(cmd, shell=True, debug=True) #, dryrun=True)

  if stype == 'selection':  # extract kd values from pickle file (since it requires ete/anaconda to read)
    cmd = 'export PATH=%s:$PATH && xvfb-run -a python ./bin/view-trees.py %s/%s_lineage_tree.p %s/kd-vals.csv' % (ete_path, simdir(stype), extrastr, simdir(stype))
    utils.simplerun(cmd, shell=True)

# ----------------------------------------------------------------------------------------
def partition(stype):
    # cmd = './bin/partis cache-parameters --infname %s --parameter-dir %s/params --n-procs 8' % (simfname(stype), simdir(stype))
    # cmd = './bin/partis partition --infname %s --parameter-dir %s/params --plotdir %s --n-procs 8 --outfname %s/partition.yaml' % (simfname(stype), simdir(stype), plotdir(stype), simdir(stype))
    # cmd = './bin/partis plot-partitions --plotdir %s --outfname %s/partition.yaml' % (plotdir(stype), simdir(stype))
    cmd = './bin/partis view-output --outfname %s/partition.yaml --abb' % simdir(stype)
    utils.simplerun(cmd, debug=True) #, dryrun=True)

for stype in ['selection', 'neutral']:
    simulate(stype)
    # partition(stype)
    break
