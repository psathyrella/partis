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
    cmd = './bin/partis simulate --simulate-from-scratch --mutation-multiplier 0.0001 --n-leaves 1 --constant-number-of-leaves'  # tends to get in infinite loop if you actually pass 0. (yes, I should fix this)
    cmd += ' --debug %d --seed %d --outfname %s/naive-simu.yaml --n-sim-events %d' % (int(args.debug), args.seed, simdir(args.stype), args.n_sim_events)
    utils.simplerun(cmd, debug=True)

# ----------------------------------------------------------------------------------------
def run_bcr_phylo(naive_line, outdir, ievent):
    tmpdir = utils.choose_random_subdir('/tmp/%s' % os.getenv('USER'))  # this is I think just for xvfb-run
    os.makedirs(tmpdir)
    prof_cmds = ''  # '-m cProfile -s tottime -o prof.out'
    cmd = 'export TMPDIR=%s && export PATH=%s:$PATH && xvfb-run -a python %s %s/bin/simulator.py' % (tmpdir, ete_path, prof_cmds, bcr_phylo_path)

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
        cmd += ' --n_to_downsample %d' % args.n_sim_seqs_per_generation
        cmd += ' --target_dist %d' % args.target_distance
        cmd += ' --target_count %d' % args.target_count
        cmd += ' --carry_cap %d' % args.carry_cap
        cmd += ' --observe_common_ancestors'

        # cmd += ' --observe_based_on_affinity'  # implementation in bcr-phylo needs some work
        # cmd += ' --kd_fuzz_fraction 0.03'  # this makes the plots nicer, but also (of course, in retrospect) reduces the force of selection
    else:
        assert False


    # cmd += ' --verbose'
    cmd += ' --no_context'
    cmd += ' --outbase %s/%s' % (outdir, args.extrastr)
    cmd += ' --naive_seq %s' % naive_line['naive_seq']
    cmd += ' --random_seed %d' % (args.seed + ievent)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    utils.simplerun(cmd, shell=True, extra_str='        ', debug=True) #, dryrun=True)
    os.rmdir(tmpdir)

# ----------------------------------------------------------------------------------------
def parse_bcr_phylo_output(glfo, naive_line, outdir, ievent):
    seqfos = utils.read_fastx('%s/%s.fasta' % (outdir, args.extrastr))  # output mutated sequences from bcr-phylo

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
        if len(set(kdvals) - set(final_line['unique_ids'])) > 0:  # uids in the kd file but not the <line> (i.e. not in the newick/fasta files) are probably just bcr-phylo discarding internal nodes
            print '        in kd file, but missing from final_line (probably just internal nodes that bcr-phylo wrote to the tree without names): %s' % (set(kdvals) - set(final_line['unique_ids']))
        if len(set(final_line['unique_ids']) - set(kdvals)) > 0:
            print '        in final_line, but missing from kdvals: %s' % ' '.join(set(final_line['unique_ids']) - set(kdvals))
        final_line['affinities'] = [1. / kdvals[u] for u in final_line['unique_ids']]
        tree = treeutils.get_dendro_tree(treefname='%s/simu.nwk' % outdir)
        if args.debug:
            print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=tree), padwidth=12)
        final_line['tree'] = tree.as_string(schema='newick')
    tmp_event = RecombinationEvent(glfo)  # I don't want to move the function out of event.py right now
    tmp_event.set_reco_id(final_line, irandom=ievent)  # not sure that setting <irandom> here actually does anything

    # get target sequences
    target_seqfos = utils.read_fastx('%s/%s_targets.fa' % (outdir, args.extrastr))
    final_line['target_seqs'] = [tfo['seq'] for tfo in target_seqfos]
    from Bio.Seq import Seq
    final_line['nearest_target_indices'] = []
    aa_targets = [Seq(seq).translate() for seq in final_line['target_seqs']]
    for mseq in final_line['input_seqs']:
        aa_mseq = Seq(mseq).translate()
        aa_hdists = [utils.hamming_distance(aa_t, aa_mseq, amino_acid=True) for aa_t in aa_targets]
        imin = aa_hdists.index(min(aa_hdists))  # NOTE doesn't do anything differently if there's more than one min
        final_line['nearest_target_indices'].append(imin)

    return final_line

# ----------------------------------------------------------------------------------------
def simulate():

    rearrange()

    glfo, naive_event_list, cpath = utils.read_output('%s/naive-simu.yaml' % simdir(args.stype))
    assert len(naive_event_list) == args.n_sim_events

    outdirs = ['%s/event-%d' % (simdir(args.stype), i) for i in range(len(naive_event_list))]

    print '    running bcr-phylo for %d naive rearrangements' % len(naive_event_list)
    for ievent, (naive_line, outdir) in enumerate(zip(naive_event_list, outdirs)):
        run_bcr_phylo(naive_line, outdir, ievent)

    mutated_events = []
    for ievent, (naive_line, outdir) in enumerate(zip(naive_event_list, outdirs)):
        mutated_events.append(parse_bcr_phylo_output(glfo, naive_line, outdir, ievent))

    print '  writing annotations to %s' % simfname(args.stype)
    utils.write_annotations(simfname(args.stype), glfo, mutated_events, utils.simulation_headers)

    import plotting
    for outdir, event in zip(outdirs, mutated_events):
        plotting.plot_bcr_phylo_simulation(outdir, event, args.extrastr)
    # utils.simplerun('cp -v %s/simu_collapsed_runstat_color_tree.svg %s/plots/' % (outdir, outdir))

# ----------------------------------------------------------------------------------------
def partition():
    n_procs = 1
    cmd = './bin/partis cache-parameters --infname %s --parameter-dir %s/params --n-procs %d --seed %d' % (simfname(args.stype), infdir(args.stype), n_procs, args.seed)
    utils.simplerun(cmd, debug=True) #, dryrun=True)
    cmd = './bin/partis partition --n-final-clusters 1 --write-additional-cluster-annotations 0:5 --lb-tau %f --is-simu --get-tree-metrics --infname %s --parameter-dir %s/params --plotdir %s --n-procs %d --outfname %s/partition.yaml --seed %d' % (args.lb_tau, simfname(args.stype), infdir(args.stype), infdir(args.stype) + '/plots', n_procs, infdir(args.stype), args.seed)
    utils.simplerun(cmd, debug=True) #, dryrun=True)
    # cmd = './bin/partis get-tree-metrics --outfname %s/partition.yaml' % infdir(args.stype)
    # utils.simplerun(cmd, debug=True) #, dryrun=True)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--stype', default='selection', choices=('selection', 'neutral'))
parser.add_argument('--debug', action='store_true')
parser.add_argument('--run-help', action='store_true')
parser.add_argument('--seed', type=int, default=1, help='random seed (note that bcr-phylo doesn\'t seem to support setting its random seed)')
parser.add_argument('--extrastr', default='simu', help='doesn\'t really do anything, but it\'s required by bcr-phylo')
parser.add_argument('--n-sim-seqs-per-generation', type=int, default=100, help='Number of sequences to sample at each time in --obs-times.')
parser.add_argument('--n-sim-events', type=int, default=1, help='number of simulated rearrangement events')
parser.add_argument('--obs-times', default='100:120', help='Times (reproductive rounds) at which to selection sequences for observation.')
parser.add_argument('--carry-cap', type=int, default=1000, help='carrying capacity of germinal center')
parser.add_argument('--target-distance', type=int, default=15, help='Desired distance (number of non-synonymous mutations) between the naive sequence and the target sequences.')
parser.add_argument('--target-count', type=int, default=10, help='Number of target sequences to generate.')
parser.add_argument('--branching-parameter', type=float, default=2., help='')
parser.add_argument('--base-mutation-rate', type=float, default=0.365, help='')

parser.add_argument('--lb-tau', type=float, default=0.4, help='')
args = parser.parse_args()

args.obs_times = utils.get_arg_list(args.obs_times, intify=True)

# ----------------------------------------------------------------------------------------
simulate()
partition()
