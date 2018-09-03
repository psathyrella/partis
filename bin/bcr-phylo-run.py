#!/usr/bin/env python
import csv
import colored_traceback.always
import collections
import copy
import os
import sys

current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
sys.path.insert(1, current_script_dir)
import utils
import indelutils

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

  cmd = './bin/partis simulate --simulate-from-scratch --mutation-multiplier 0.0001 --debug 1 --n-leaves 1 --constant-number-of-leaves --outfname %s/naive-simu.yaml' % simdir(stype)
  utils.simplerun(cmd, debug=True)
  glfo, annotation_list, cpath = utils.read_output('%s/naive-simu.yaml' % simdir(stype))
  assert len(annotation_list) == 1  # see below
  naive_line = annotation_list[0]
  # # don't need to write them a.t.m., since we're only running bcr-phylo one family at a time
  # with open('%s/naive-simu.fa' % simdir(stype), 'w') as nfile:
  #   for line in annotation_list:
  #     if len(line['unique_ids']) > 1:
  #       raise Exception('should just be simulating naive sequences with partis above')
  #     nfile.write('>%s\n%s\n' % (line['unique_ids'][0], line['naive_seq']))  # kind of weird to write the naive seq instead of the input seqs, but it's what we want, just don't forget that this is how you're arranging things

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
  # cmd += ' --random_seq %s/sequence_data/AbPair_naive_seqs.fa' % bcr_phylo_path
  # cmd += ' --random_seq %s/naive-simu.fa' % simdir(stype)  # see note above
  cmd += ' --sequence %s' % naive_line['naive_seq']
  if not os.path.exists(simdir(stype)):
    os.makedirs(simdir(stype))

  utils.simplerun(cmd, shell=True, debug=True) #, dryrun=True)

  seqfos = utils.read_fastx('%s/%s.fasta' % (simdir(stype), extrastr))  # output mutated sequences from bcr-phylo
  seqfos = [sfo for sfo in seqfos if sfo['name'] != 'naive']  # don't have a kd value for it, anyway

  assert len(annotation_list) == 1  # see notes above
  assert len(naive_line['unique_ids']) == 1
  assert not indelutils.has_indels(naive_line['indelfos'][0])  # would have to handle this below
  utils.print_reco_event(naive_line)
  reco_info = collections.OrderedDict()
  for sfo in seqfos:
    mline = copy.deepcopy(naive_line)
    utils.remove_all_implicit_info(mline)
    mline['unique_ids'] = [sfo['name']]
    mline['seqs'] = [sfo['seq']]  # it's really important to set both the seqs (since they're both already in there from the naive line)
    mline['input_seqs'] = [sfo['seq']]  # it's really important to set both the seqs (since they're both already in there from the naive line)
    reco_info[sfo['name']] = mline
    utils.add_implicit_info(glfo, mline)
  final_line = utils.synthesize_multi_seq_line_from_reco_info([sfo['name'] for sfo in seqfos], reco_info)
  utils.print_reco_event(final_line)

  if stype == 'selection':  # extract kd values from pickle file (since it requires ete/anaconda to read)
    cmd = 'export PATH=%s:$PATH && xvfb-run -a python ./bin/view-trees.py %s/%s_lineage_tree.p %s/kd-vals.csv' % (ete_path, simdir(stype), extrastr, simdir(stype))
    utils.simplerun(cmd, shell=True)
    kdvals = {}
    with open('%s/kd-vals.csv' % simdir(stype)) as kdfile:
      reader = csv.DictReader(kdfile)
      for line in reader:
        kdvals[line['uid']] = float(line['kd'])
    if set(kdvals) != set(final_line['unique_ids']):
      print '        missing from final_line: %s' % ' '.join(set(kdvals) - set(final_line['unique_ids']))
      print '        missing from kdvals: %s' % ' '.join(set(final_line['unique_ids']) - set(kdvals))
      print '      %s kd file uids don\'t match final line (see above, it\'s maybe just internal nodes?)' % utils.color('red', 'note:')
    final_line['affinities'] = [kdvals[u] for u in final_line['unique_ids']]

  assert len(cpath.partitions) == 0
  utils.write_annotations('%s/mutated-simu.yaml' % simdir(stype), glfo, [final_line], utils.simulation_headers)

# ----------------------------------------------------------------------------------------
def partition(stype):
    # cmd = './bin/partis cache-parameters --infname %s --parameter-dir %s/params --n-procs 8' % (simfname(stype), infdir(stype))
    cmd = './bin/partis partition --calculate-tree-metrics --infname %s --parameter-dir %s/params --n-procs 1 --outfname %s/partition.yaml' % (simfname(stype), infdir(stype), infdir(stype))
    # cmd = './bin/partis plot-partitions --plotdir %s/plots --outfname %s/partition.yaml' % (infdir(stype), infdir(stype))
    # cmd = './bin/partis view-output --outfname %s/partition.yaml --abb' % infdir(stype)
    utils.simplerun(cmd, debug=True) #, dryrun=True)

for stype in ['selection', 'neutral']:
    simulate(stype)
    # partition(stype)
    break
