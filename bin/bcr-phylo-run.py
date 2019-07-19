#!/usr/bin/env python
import argparse
import csv
import colored_traceback.always
import collections
import copy
import os
import sys
import numpy

current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
sys.path.insert(1, current_script_dir)
import utils
import indelutils
import treeutils
from event import RecombinationEvent

ete_path = '/home/' + os.getenv('USER') + '/anaconda_ete/bin'
bcr_phylo_path = os.getenv('PWD') + '/packages/bcr-phylo-benchmark'

# ----------------------------------------------------------------------------------------
def simdir():
    return '%s/%s/simu' % (args.base_outdir, args.stype)
def infdir():
    return '%s/%s/partis' % (args.base_outdir, args.stype)
def naive_fname():
    return '%s/naive-simu.yaml' % simdir()
def bcr_phylo_fasta_fname(outdir):
    return '%s/%s.fasta' % (outdir, args.extrastr)
def simfname():
    return '%s/mutated-simu.yaml' % simdir()
def param_dir():
    return '%s/params' % infdir()
def partition_fname():
    return '%s/partition.yaml' % infdir()

# ----------------------------------------------------------------------------------------
def rearrange():
    if utils.output_exists(args, naive_fname(), outlabel='naive simu', offset=4):
        return
    cmd = './bin/partis simulate --simulate-from-scratch --mutation-multiplier 0.0001 --n-leaves 1 --constant-number-of-leaves'  # tends to get in infinite loop if you actually pass 0. (yes, I should fix this)
    cmd += ' --debug %d --seed %d --outfname %s --n-sim-events %d' % (int(args.debug), args.seed, naive_fname(), args.n_sim_events)
    utils.simplerun(cmd, debug=True)

# ----------------------------------------------------------------------------------------
def run_bcr_phylo(naive_line, outdir, ievent, n_total_events):
    if utils.output_exists(args, bcr_phylo_fasta_fname(outdir), outlabel='bcr-phylo', offset=4):
        return

    cmd = '%s/bin/simulator.py' % bcr_phylo_path
    if args.run_help:
        cmd += ' --help'
    elif args.stype == 'neutral':
        assert False  # needs updating (well, maybe not, but I'm not thinking about it when I move the selection parameters to command line args)
        cmd += ' --lambda %f --lambda0 %f' % (1.5, 0.365)
        cmd += ' --n_final_seqs %d' % args.n_sim_seqs_per_generation
    elif args.stype == 'selection':
        cmd += ' --selection'
        cmd += ' --lambda %f' % args.branching_parameter
        cmd += ' --lambda0 %f' % args.base_mutation_rate
        cmd += ' --obs_times %s' % ' '.join(['%d' % t for t in args.obs_times])
        cmd += ' --n_to_sample %s' % ' '.join('%d' % n for n in args.n_sim_seqs_per_generation)
        cmd += ' --metric_for_target_dist %s' % args.metric_for_target_distance
        cmd += ' --target_dist %d' % args.target_distance
        cmd += ' --target_count %d' % args.target_count
        cmd += ' --carry_cap %d' % args.carry_cap
        cmd += ' --observe_common_ancestors'

        # cmd += ' --n_target_clusters 1'
        # cmd += ' --target_cluster_distance 1'

        # cmd += ' --observe_based_on_affinity'  # implementation in bcr-phylo needs some work
    else:
        assert False

    cmd += ' --debug %d' % args.debug
    cmd += ' --n_tries 30'
    cmd += ' --no_context'
    cmd += ' --no_plot'
    cmd += ' --outbase %s/%s' % (outdir, args.extrastr)
    cmd += ' --random_seed %d' % (args.seed + ievent)
    if n_total_events > 1:  # if the final sample's going to contain many trees, it's worth making the uids longer so there's fewer collisions/duplicates
        cmd += ' --uid_str_len 7'
    cmd += ' --naive_seq %s' % naive_line['naive_seq']

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    utils.run_ete_script(cmd, ete_path)  # NOTE kind of hard to add a --dry-run option, since we have to loop over the events we made in rearrange()

# ----------------------------------------------------------------------------------------
def parse_bcr_phylo_output(glfo, naive_line, outdir, ievent):
    seqfos = utils.read_fastx(bcr_phylo_fasta_fname(outdir))  # output mutated sequences from bcr-phylo

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
        mline['duplicates'] = [[]]
        reco_info[sfo['name']] = mline
        utils.add_implicit_info(glfo, mline)
    final_line = utils.synthesize_multi_seq_line_from_reco_info([sfo['name'] for sfo in seqfos], reco_info)
    if args.debug:
        utils.print_reco_event(final_line)

    # extract kd values from pickle file (use a separate script since it requires ete/anaconda to read)
    if args.stype == 'selection':
        cmd = './bin/read-bcr-phylo-trees.py --pickle-tree-file %s/%s_lineage_tree.p --kdfile %s/kd-vals.csv --newick-tree-file %s/simu.nwk' % (outdir, args.extrastr, outdir, outdir)
        utils.run_ete_script(cmd, ete_path)
        nodefo = {}
        with open('%s/kd-vals.csv' % outdir) as kdfile:
            reader = csv.DictReader(kdfile)
            for line in reader:
                nodefo[line['uid']] = {
                    'kd' : float(line['kd']),
                    'relative_kd' : float(line['relative_kd']),
                    'lambda' : line.get('lambda', None),
                    'target_index' : int(line['target_index']),
                }
        if len(set(nodefo) - set(final_line['unique_ids'])) > 0:  # uids in the kd file but not the <line> (i.e. not in the newick/fasta files) are probably just bcr-phylo discarding internal nodes
            print '        in kd file, but missing from final_line (probably just internal nodes that bcr-phylo wrote to the tree without names): %s' % (set(nodefo) - set(final_line['unique_ids']))
        if len(set(final_line['unique_ids']) - set(nodefo)) > 0:
            print '        in final_line, but missing from kdvals: %s' % ' '.join(set(final_line['unique_ids']) - set(nodefo))
        final_line['affinities'] = [1. / nodefo[u]['kd'] for u in final_line['unique_ids']]
        final_line['relative_affinities'] = [1. / nodefo[u]['relative_kd'] for u in final_line['unique_ids']]
        final_line['lambdas'] = [nodefo[u]['lambda'] for u in final_line['unique_ids']]
        final_line['nearest_target_indices'] = [nodefo[u]['target_index'] for u in final_line['unique_ids']]
        tree = treeutils.get_dendro_tree(treefname='%s/simu.nwk' % outdir)
        tree.scale_edges(1. / numpy.mean([len(s) for s in final_line['seqs']]))
        if args.debug:
            print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=tree), padwidth=12)
        final_line['tree'] = tree.as_string(schema='newick')
    tmp_event = RecombinationEvent(glfo)  # I don't want to move the function out of event.py right now
    tmp_event.set_reco_id(final_line, irandom=ievent)  # not sure that setting <irandom> here actually does anything

    # get target sequences
    target_seqfos = utils.read_fastx('%s/%s_targets.fa' % (outdir, args.extrastr))
    final_line['target_seqs'] = [tfo['seq'] for tfo in target_seqfos]

    return final_line

# ----------------------------------------------------------------------------------------
def simulate():

    rearrange()

    glfo, naive_event_list, cpath = utils.read_output(naive_fname())
    assert len(naive_event_list) == args.n_sim_events

    outdirs = ['%s/event-%d' % (simdir(), i) for i in range(len(naive_event_list))]

    for ievent, (naive_line, outdir) in enumerate(zip(naive_event_list, outdirs)):
        if args.n_sim_events > 1:
            print '  %s %d' % (utils.color('blue', 'ievent'), ievent)
        run_bcr_phylo(naive_line, outdir, ievent, len(naive_event_list))

    if utils.output_exists(args, simfname(), outlabel='mutated simu', offset=4):  # i guess if it crashes during the plotting just below, this'll get confused
        return

    mutated_events = []
    for ievent, (naive_line, outdir) in enumerate(zip(naive_event_list, outdirs)):
        mutated_events.append(parse_bcr_phylo_output(glfo, naive_line, outdir, ievent))

    print '  writing annotations to %s' % simfname()
    utils.write_annotations(simfname(), glfo, mutated_events, utils.simulation_headers)

    if not args.only_csv_plots:
        import lbplotting
        for outdir, event in zip(outdirs, mutated_events):
            lbplotting.plot_bcr_phylo_simulation(outdir, event, args.extrastr, args.metric_for_target_distance)
        # utils.simplerun('cp -v %s/simu_collapsed_runstat_color_tree.svg %s/plots/' % (outdir, outdir))

# ----------------------------------------------------------------------------------------
def cache_parameters():
    if utils.output_exists(args, param_dir() + '/hmm/hmms', outlabel='parameters', offset=4):
        return
    cmd = './bin/partis cache-parameters --infname %s --parameter-dir %s --n-procs %d --seed %d' % (simfname(), param_dir(), args.n_procs, args.seed)
    utils.simplerun(cmd, debug=True) #, dryrun=True)

# ----------------------------------------------------------------------------------------
def partition():
    if utils.output_exists(args, partition_fname(), outlabel='partition', offset=4):
        return
    cmd = './bin/partis partition --simultaneous-true-clonal-seqs --is-simu --infname %s --parameter-dir %s --n-procs %d --outfname %s --seed %d' % (simfname(), param_dir(), args.n_procs, partition_fname(), args.seed)
    #  --write-additional-cluster-annotations 0:5  # I don't think there was really a good reason for having this
    if not args.dont_get_tree_metrics:
        cmd += ' --get-tree-metrics --plotdir %s' % (infdir() + '/plots')
    if args.lb_tau is not None:
        cmd += ' --lb-tau %f' % args.lb_tau
    utils.simplerun(cmd, debug=True) #, dryrun=True)
    # cmd = './bin/partis get-tree-metrics --outfname %s/partition.yaml' % infdir()
    # utils.simplerun(cmd, debug=True) #, dryrun=True)

# ----------------------------------------------------------------------------------------
all_actions = ('simu', 'cache-parameters', 'partition')
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter)
parser.add_argument('--stype', default='selection', choices=('selection', 'neutral'))
parser.add_argument('--actions', default=':'.join(all_actions))
parser.add_argument('--base-outdir', default='%s/partis/bcr-phylo/test' % os.getenv('fs', default=os.getenv('HOME')))
parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
parser.add_argument('--run-help', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--only-csv-plots', action='store_true')
parser.add_argument('--dont-get-tree-metrics', action='store_true', help='Partition without getting tree metrics, presumably because you want to run them yourself later')
parser.add_argument('--seed', type=int, default=1, help='random seed (note that bcr-phylo doesn\'t seem to support setting its random seed)')
parser.add_argument('--n-procs', type=int, default=1)
parser.add_argument('--extrastr', default='simu', help='doesn\'t really do anything, but it\'s required by bcr-phylo')
parser.add_argument('--n-sim-seqs-per-generation', default='100', help='Number of sequences to sample at each time in --obs-times.')
parser.add_argument('--n-sim-events', type=int, default=1, help='number of simulated rearrangement events')
parser.add_argument('--obs-times', default='100:120', help='Times (reproductive rounds) at which to selection sequences for observation.')
parser.add_argument('--carry-cap', type=int, default=1000, help='carrying capacity of germinal center')
parser.add_argument('--target-distance', type=int, default=15, help='Desired distance (number of non-synonymous mutations) between the naive sequence and the target sequences.')
parser.add_argument('--metric-for-target-distance', default='aa', choices=['aa', 'nuc', 'aa-sim'], help='see bcr-phylo docs')
parser.add_argument('--target-count', type=int, default=1, help='Number of target sequences to generate.')
parser.add_argument('--branching-parameter', type=float, default=2., help='see bcr-phylo docs')
parser.add_argument('--base-mutation-rate', type=float, default=0.365, help='see bcr-phylo docs')
parser.add_argument('--lb-tau', type=float, help='')

args = parser.parse_args()

args.obs_times = utils.get_arg_list(args.obs_times, intify=True)
args.n_sim_seqs_per_generation = utils.get_arg_list(args.n_sim_seqs_per_generation, intify=True)
args.actions = utils.get_arg_list(args.actions, choices=all_actions)

# ----------------------------------------------------------------------------------------
if 'simu' in args.actions:
    simulate()
if 'cache-parameters' in args.actions:
    cache_parameters()
if 'partition' in args.actions:
    partition()
