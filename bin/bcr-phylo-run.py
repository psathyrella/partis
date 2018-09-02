#!/usr/bin/env python
import os
import sys

current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
sys.path.insert(1, current_script_dir)
import utils

base_outdir = '%s/partis/bcr-phylo' % os.getenv('fs')
extrastr = 'simu'
label = 'test'
n_sim_seqs = 10
selection_time = 5  # 35
carry_cap = 1000  # 1000
n_mutations = 5

ete_path = '/home/' + os.getenv('USER') + '/anaconda_ete/bin'
bcr_phylo_path = os.getenv('PWD') + '/packages/bcr-phylo-benchmark'

# ----------------------------------------------------------------------------------------
def simdir(stype):
  return '%s/%s/%s/simu' % (base_outdir, stype, label)
def infdir(stype):
  return '%s/%s/%s/partis' % (base_outdir, stype, label)
def simfname(stype):
  return '%s/%s.fasta' % (simdir(stype), extrastr)

# ----------------------------------------------------------------------------------------
def simulate(stype):
  os.environ['TMPDIR'] = '/tmp'
  prof_cmds = '' #'-m cProfile -s tottime -o prof.out'
  cmd = 'export PATH=%s:$PATH && xvfb-run -a python %s %s/bin/simulator.py' % (ete_path, prof_cmds, bcr_phylo_path)
  cmd += ' --mutability %s/S5F/Mutability.csv --substitution %s/S5F/Substitution.csv' % (bcr_phylo_path, bcr_phylo_path)
  # from maual: --T 35 --n 60 --selection --target_dist 5 --target_count 100 --carry_cap 1000 --skip_update 100

  if stype == 'neutral':
    cmd += ' --lambda %f --lambda0 %f --N %d' % (1.5, 0.365, n_sim_seqs)
  elif stype == 'selection':
    cmd += ' --selection'
    cmd += ' --lambda %f --lambda0 %f' % (2., 0.365)
    # cmd += ' --n %d' % 60  # cells downsampled (default None)
    cmd += ' --T %d' % selection_time
    cmd += ' --target_dist %d --target_count %d' % (n_mutations, n_sim_seqs)  # um, I think target count is maybe the number of sequences?
    cmd += ' --carry_cap %d --skip_update 100' % (carry_cap)
  else:
    assert False

  cmd += ' --verbose'
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
    cmd = './bin/partis cache-parameters --infname %s --parameter-dir %s/params --n-procs 8' % (simfname(stype), infdir(stype))
    # cmd = './bin/partis partition --infname %s --parameter-dir %s/params --plotdir %s/plots --n-procs 8 --outfname %s/partition.yaml' % (simfname(stype), infdir(stype), infdir(stype), infdir(stype))
    # cmd = './bin/partis plot-partitions --plotdir %s/plots --outfname %s/partition.yaml' % (infdir(stype), infdir(stype))
    # cmd = './bin/partis view-output --outfname %s/partition.yaml --abb' % infdir(stype)
    utils.simplerun(cmd, debug=True) #, dryrun=True)

for stype in ['selection', 'neutral']:
    # simulate(stype)
    partition(stype)
    break
