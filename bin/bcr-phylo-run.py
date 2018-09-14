#!/usr/bin/env python
import argparse
import csv
import colored_traceback.always
import collections
import copy
import os
import sys

current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
sys.path.insert(1, current_script_dir)
sys.path.insert(1, current_script_dir + '/../packages/baltic')
import utils
import indelutils
import treeutils
from event import RecombinationEvent

base_outdir = '%s/partis/bcr-phylo' % os.getenv('fs')
label = 'test'

ete_path = '/home/' + os.getenv('USER') + '/anaconda_ete/bin'
bcr_phylo_path = os.getenv('PWD') + '/packages/bcr-phylo-benchmark'

# ----------------------------------------------------------------------------------------
def simdir(stype):
    return '%s/%s/%s/simu' % (base_outdir, stype, label)
def infdir(stype):
    return '%s/%s/%s/partis' % (base_outdir, stype, label)
def simfname(stype):
    return '%s/mutated-simu.yaml' % simdir(stype)

# ----------------------------------------------------------------------------------------
def rearrange():
    cmd = './bin/partis simulate --simulate-from-scratch --mutation-multiplier 0.0001 --seed %d --debug 1 --n-leaves 1 --constant-number-of-leaves --outfname %s/naive-simu.yaml' % (args.seed, simdir(args.stype))
    utils.simplerun(cmd, debug=True)

# ----------------------------------------------------------------------------------------
def read_naive_simulation():
    glfo, annotation_list, cpath = utils.read_output('%s/naive-simu.yaml' % simdir(args.stype))
    assert len(annotation_list) == 1  # would need to change some things
    assert len(cpath.partitions) == 0
    naive_line = annotation_list[0]
    return glfo, naive_line

# ----------------------------------------------------------------------------------------
def run_bcr_phylo():
    _, naive_line = read_naive_simulation()
    os.environ['TMPDIR'] = '/tmp'
    prof_cmds = ''  # '-m cProfile -s tottime -o prof.out'
    cmd = 'export PATH=%s:$PATH && xvfb-run -a python %s %s/bin/simulator.py' % (ete_path, prof_cmds, bcr_phylo_path)

    if args.run_help:
        cmd += ' --help'
    elif args.stype == 'neutral':
        assert False  # needs updating (well, maybe not, but I'm not thinking about it when I move the selection parameters to command line args)
        cmd += ' --lambda %f --lambda0 %f' % (1.5, 0.365)
        cmd += ' --n_final_seqs %d' % args.n_sim_seqs
    elif args.stype == 'selection':
        cmd += ' --selection'
        cmd += ' --lambda %f' % args.branching_parameter
        cmd += ' --lambda0 %f' % args.base_mutation_rate
        cmd += ' --obs_times %d' % args.obs_time
        cmd += ' --target_dist %d' % args.target_distance
        cmd += ' --carry_cap %d' % args.carry_cap
        cmd += ' --n_to_downsample %d' % args.n_sim_seqs  # Number of cells kept (not discarded) during final downsampling step (default: None)
        # cmd += ' --stop_dist %d'  % xxx  # Stop when any simulated sequence is closer than this (hamming distance) to any of the target sequences.
        # can't set both N and T, but need to set T for selection:
        # cmd += ' --n_final_seqs %d' % args.n_sim_seqs  # Desired number of final sequences (actual number may be less due to removal of nonsense sequences)
    else:
        assert False

    cmd += ' --no_context'
    cmd += ' --verbose'
    cmd += ' --outbase %s/%s' % (simdir(args.stype), args.extrastr)
    # cmd += ' --random_seq %s/sequence_data/AbPair_naive_seqs.fa' % bcr_phylo_path
    cmd += ' --sequence %s' % naive_line['naive_seq']
    cmd += ' --random_seed %d' % args.seed
    if not os.path.exists(simdir(args.stype)):
        os.makedirs(simdir(args.stype))

    utils.simplerun(cmd, shell=True, debug=True) #, dryrun=True)

# ----------------------------------------------------------------------------------------
def parse_bcr_phylo_output():
    glfo, naive_line = read_naive_simulation()

    outdir = simdir(args.stype)
    seqfos = utils.read_fastx('%s/%s.fasta' % (outdir, args.extrastr))  # output mutated sequences from bcr-phylo
    seqfos = [sfo for sfo in seqfos if sfo['name'] != 'naive']  # don't have a kd value for the naive sequence, so may as well throw it out

    assert len(naive_line['unique_ids']) == 1  # enforces that we ran naive-only, 1-leaf partis simulation above
    assert not indelutils.has_indels(naive_line['indelfos'][0])  # would have to handle this below
    if args.debug:
        utils.print_reco_event(naive_line)
    reco_info = collections.OrderedDict()
    for sfo in seqfos:
        mline = copy.deepcopy(naive_line)
        utils.remove_all_implicit_info(mline)
        del mline['tree']
        mline['unique_ids'] = [sfo['name']]
        mline['seqs'] = [sfo['seq']]  # it's really important to set both the seqs (since they're both already in there from the naive line)
        mline['input_seqs'] = [sfo['seq']]  # it's really important to set both the seqs (since they're both already in there from the naive line)
        reco_info[sfo['name']] = mline
        utils.add_implicit_info(glfo, mline)
    final_line = utils.synthesize_multi_seq_line_from_reco_info([sfo['name'] for sfo in seqfos], reco_info)
    tmp_event = RecombinationEvent(glfo)  # I don't want to move the function out of event.py right now
    tmp_event.set_reco_id(final_line)  # I _think_ I don't need to set <irandom>
    if args.debug:
        utils.print_reco_event(final_line)

    # extract kd values from pickle file (use a separate script since it requires ete/anaconda to read)
    if args.stype == 'selection':
        cmd = 'export PATH=%s:$PATH && xvfb-run -a python ./bin/view-trees.py --pickle-tree-file %s/%s_lineage_tree.p --kdfile %s/kd-vals.csv --newick-tree-file %s/simu.nwk' % (ete_path, outdir, args.extrastr, outdir, outdir)
        utils.simplerun(cmd, shell=True)
        kdvals = {}
        with open('%s/kd-vals.csv' % outdir) as kdfile:
            reader = csv.DictReader(kdfile)
            for line in reader:
                kdvals[line['uid']] = float(line['kd'])
        # if len(set(kdvals) - set(final_line['unique_ids'])) > 0:  # uids in the kd file but not the <line> (i.e. not in the newick/fasta files) are probably just bcr-phylo discarding internal nodes
        #     print '        missing from final_line (probably just internal nodes that bcr-phylo wrote to the tree without names): %s' % ' '.join(set(kdvals) - set(final_line['unique_ids']))
        if len(set(final_line['unique_ids']) - set(kdvals)) > 0:
            print '        missing from kdvals: %s' % ' '.join(set(final_line['unique_ids']) - set(kdvals))
        final_line['affinities'] = [kdvals[u] for u in final_line['unique_ids']]
        tree = treeutils.get_dendro_tree(treefname='%s/simu.nwk' % outdir , ignore_existing_internal_node_labels=True)  # bcr-phylo sets all the internal node labels to '1', so we have to relabel them so dendropy doesn't later barf
        if args.debug:
            print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=tree), padwidth=12)
        final_line['tree'] = tree.as_string(schema='newick')

    print '  writing annotations to %s' % simfname(args.stype)
    utils.write_annotations(simfname(args.stype), glfo, [final_line], utils.simulation_headers)

# ----------------------------------------------------------------------------------------
def simulate():
    # rearrange()
    # run_bcr_phylo()
    parse_bcr_phylo_output()

# ----------------------------------------------------------------------------------------
def partition():
    n_procs = 1
    cmd = './bin/partis cache-parameters --infname %s --parameter-dir %s/params --n-procs %d --seed %d' % (simfname(args.stype), infdir(args.stype), n_procs, args.seed)
    utils.simplerun(cmd, debug=True) #, dryrun=True)
    cmd = './bin/partis partition --debug 0 --n-final-clusters 1 --write-additional-cluster-annotations 0:5 --is-simu --calculate-tree-metrics --infname %s --parameter-dir %s/params --n-procs %d --outfname %s/partition.yaml --seed %d' % (simfname(args.stype), infdir(args.stype), n_procs, infdir(args.stype), args.seed)
    utils.simplerun(cmd, debug=True) #, dryrun=True)

    # cmd = './bin/partis view-output --outfname %s/partition.yaml --calculate-tree-metrics' % infdir(args.stype)  # TODO not sure this works
    # utils.simplerun(cmd, debug=True) #, dryrun=True)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--stype', default='selection', choices=('selection', 'neutral'))
parser.add_argument('--debug', action='store_true')
parser.add_argument('--run-help', action='store_true')
parser.add_argument('--seed', type=int, default=1, help='random seed (note that bcr-phylo doesn\'t seem to support setting its random seed)')
parser.add_argument('--extrastr', default='simu', help='doesn\'t really do anything, but it\'s required by bcr-phylo')
parser.add_argument('--n-sim-seqs', type=int, default=35, help='desired number of final sequences (typically, we downsample to get to this number)')
parser.add_argument('--obs-time', type=int, default=35, help='number of rounds of reproduction')
parser.add_argument('--carry-cap', type=int, default=1000, help='carrying capacity of germinal center')
parser.add_argument('--target-distance', type=int, default=25, help='Desired distance (number of non-synonymous mutations) between the naive sequence and the target sequences. If larger than 50 (ish) it seems to have trouble finding target sequencess.')
parser.add_argument('--branching-parameter', type=float, default=2., help='')
parser.add_argument('--base-mutation-rate', type=float, default=0.365, help='')
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
simulate()
partition()
