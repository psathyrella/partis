#!/usr/bin/env python
import argparse
import time
import random
import sys
from multiprocessing import Process, active_children
import os
sys.path.insert(1, './python')

import utils
# merged data: /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_merged.tsv.bz2

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true', help='Passed on to ROOT when plotting')
sys.argv.append('-b')
parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
parser.add_argument('--no-clean', action='store_true', help='Don\'t remove the various temp files')

# basic actions
parser.add_argument('--action', choices=('cache-parameters', 'run-viterbi', 'run-forward', 'partition', 'simulate', 'build-hmms', 'generate-trees'), help='What do you want to do?')

# finer action control
# parser.add_argument('--pair', action='store_true', help='Run on every pair of sequences in the input')
parser.add_argument('--n-sets', type=int, default=1, help='Run on sets of sequences of size <n> (i.e. \"k-hmm\")')
parser.add_argument('--all-combinations', action='store_true', help='Run algorithm on *all* possible combinations of the input queries of length <n-sets> (otherwise we run on sequential sets of <n-sets> in the input file).')
parser.add_argument('--is-data', action='store_true', help='True if not simulation')
parser.add_argument('--skip-unproductive', action='store_true', help='Skip sequences which Smith-Waterman determines to be unproductive (they have stop codons, are out of frame, etc.)')
parser.add_argument('--plot-performance', action='store_true', help='Write out plots comparing true and inferred distributions')
parser.add_argument('--truncate-n-sets', action='store_true', help='If running on <n-sets> sequences, truncate such that they all have the same length to the left and right of the conserved cysteine')
parser.add_argument('--pad-sequences', action='store_true', help='Pad (with Ns) all input sequences out to, at minimum, the same length on either side of the conserved cysteine. Also pad out far enough so as to eliminate all v_5p and j_3p deletions.')
parser.add_argument('--naivety', default='M', choices=['N', 'M'], help='Naive or mature sequences?')
parser.add_argument('--seed', type=int, default=int(time.time()), help='Random seed for use (mostly) by recombinator (to allow reproducibility)')
parser.add_argument('--branch-length-multiplier', type=float, help='Multiply observed branch lengths by some factor when simulating, e.g. if in data it was 0.05, but you want ten percent in your simulation, set this to 2')
# parser.add_argument('--plot-all-best-events', action='store_true', help='Plot all of the <n-best-events>, i.e. sample from the posterior')
parser.add_argument('--plot-parameters', action='store_true', help='Plot inferred parameters?')
parser.add_argument('--dont-mimic-data-read-length', action='store_true', help='Simulate events with the entire v, d, and j regions? (Otherwise we mimic the read length observed in data)')
parser.add_argument('--no-plot', action='store_true', help='Don\'t write any plots (we write a *lot* of plots for debugging, which can be slow).')
parser.add_argument('--vollmers-clustering', action='store_true', help='Perform annotation-based clustering from Vollmers paper')
parser.add_argument('--force-dont-randomize-input-order', action='store_true', help='For scons test we want to be able to overide randomization of sequence order.')
parser.add_argument('--rescale-emissions', action='store_true')
# parser.add_argument('--use_mean_at_boundaries', action='store_true')

# input and output locations
parser.add_argument('--seqfile', help='input sequence file')
parser.add_argument('--parameter-dir', required=True, help='Directory to/from which to write/read sample-specific parameters')
parser.add_argument('--datadir', default='data/imgt', help='Directory from which to read non-sample-specific information (e.g. germline genes)')
parser.add_argument('--outfname')
parser.add_argument('--plotdir', default='_plots')
parser.add_argument('--ighutil-dir', default=os.getenv('HOME') + '/.local', help='Path to vdjalign executable. The default is where \'pip install --user\' typically puts things')
parser.add_argument('--workdir', default='/tmp/' + os.path.basename(os.getenv('HOME')) + '/hmms/' + str(random.randint(0, 99999)), help='Temporary working directory (see also <no-clean>)')
parser.add_argument('--initial-cachefname')

# run/batch control
parser.add_argument('--n-procs', default='1', help='Max/initial number of processes over which to parallelize (Can be colon-separated list: first number is procs for hmm, second (should be smaller) is procs for smith-waterman, hamming, etc.)')
parser.add_argument('--slurm', action='store_true', help='Run multiple processes with slurm, otherwise just runs them on local machine. NOTE make sure to set <workdir> to something visible on all batch nodes.')
parser.add_argument('--queries', help='Colon-separated list of query names to which we restrict ourselves')
parser.add_argument('--reco-ids', help='Colon-separated list of rearrangement-event IDs to which we restrict ourselves')  # or recombination events
parser.add_argument('--n-max-queries', type=int, default=-1, help='Maximum number of query sequences on which to run (except for simulator, where it\'s the number of rearrangement events)')
parser.add_argument('--only-genes', help='Colon-separated list of genes to which to restrict the analysis')
parser.add_argument('--n-best-events', type=int, default=3, help='Number of best events to print (i.e. n-best viterbi paths)')
parser.add_argument('--max_clusters_per_proc', type=int, default=50)

# tree generation (see also branch-length-fname)
# NOTE see also branch-length-multiplier, although that comes into play after the trees are generated
parser.add_argument('--n-trees', type=int, default=500, help='Number of trees to generate')
parser.add_argument('--n-leaves', type=int, default=5, help='Number of leaves per tree')
parser.add_argument('--random-number-of-leaves', action='store_true', help='For each tree choose a random number of leaves based on <n-leaves> (a.t.m. from a hacktified exponential). Otherwise give all trees <n-leaves> leaves')

# numerical inputs
parser.add_argument('--hamming-cluster-cutoff', type=float, default=0.2, help='Threshold for hamming distance single-linkage preclustering')  # See plots in this (https://github.com/psathyrella/partis-dev/issues/70) issue for justification. TODO set threshold dynamically (for each cluster pair) based on uncertainty derived from n-best viterbi paths
parser.add_argument('--pair-hmm-cluster-cutoff', type=float, default=0.0, help='Threshold for pair hmm single-linkage preclustering')
parser.add_argument('--min_observations_to_write', type=int, default=20, help='For hmmwriter.py, if we see a gene version fewer times than this, we sum over other alleles, or other versions, etc. (see hmmwriter)')
parser.add_argument('--n-max-per-region', default='3:5:2', help='Number of best smith-waterman matches (per region, in the format v:d:j) to pass on to the hmm')
parser.add_argument('--default-v-fuzz', type=int, default=5, help='Size of the k space region over which to sum in the v direction')
parser.add_argument('--default-d-fuzz', type=int, default=2, help='Size of the k space region over which to sum in the d direction')
parser.add_argument('--smc-particles', type=int, default=1, help='Number of particles (clustering paths) to simulate with SMC')
parser.add_argument('--gap-open-penalty', type=int, default=1000, help='Penalty for indel creation in Smith-Waterman step (reduce this if you want to look for indels).')

# temporary arguments (i.e. will be removed as soon as they're not needed)
# parser.add_argument('--tree-parameter-file', default='/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz', help='File from which to read inferred tree parameters (from mebcell analysis)')
parser.add_argument('--gtrfname', default='data/recombinator/gtr.txt', help='File with list of GTR parameters. Fed into bppseqgen along with the chosen tree')
parser.add_argument('--branch-length-fname', default='data/recombinator/branch-lengths.txt', help='Branch lengths from Connor\'s mebcell stuff')
# NOTE command to generate gtr parameter file: [stoat] partis/ > zcat /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz | jq .independentParameters | grep -v '[{}]' | sed 's/["\:,]//g' | sed 's/^[ ][ ]*//' | sed 's/ /,/' | sort >data/gtr.txt

# uncommon arguments
parser.add_argument('--apply-choice_probs_in_sw', action='store_true', help='Apply gene choice probs in Smith-Waterman step. Probably not a good idea (see comments in waterer.py).')
parser.add_argument('--insertion-base-content', default=True, action='store_true',help='Account for non-uniform base content in insertions. Slows us down by a factor around five and gives no performance benefit.')
parser.add_argument('--allow-unphysical-insertions', action='store_true', help='allow insertions on left side of v and right side of j. NOTE this is very slow.')
# parser.add_argument('--allow_external_deletions', action='store_true')     # ( " ) deletions (               "                     )
# parser.add_argument('--total-length-from-right', type=int, default=-1, help='Total read length you want for simulated sequences')
parser.add_argument('--joint-emission', action='store_true', help='Use information about both sequences when writing pair emission probabilities?')

args = parser.parse_args()
args.only_genes = utils.get_arg_list(args.only_genes)

args.n_procs = utils.get_arg_list(args.n_procs, intify=True)
if len(args.n_procs) == 1:
    args.n_fewer_procs = args.n_procs[0]
else:
    args.n_fewer_procs = args.n_procs[1]
args.n_procs = args.n_procs[0]

print '\n\nmove bad cysteine-skipping code somewhere better? or at least rationalize it a bit\n\n'
print '\n\nclean up raw_best stuff\n\n'
if args.slurm and '/tmp' in args.workdir:
    raise Exception('ERROR it appears that <workdir> isn\'t set to something visible to all slurm nodes')

if args.pad_sequences and not args.allow_unphysical_insertions:
    raise Exception('pad-sequences set without allowing v_5p and j_3p deletions... which isn\'t going to work')

print '\nsetting rescale emissions to true\n'
args.rescale_emissions = True

assert not args.truncate_n_sets  # disabled and deprecated (I'm breaking it to make N padding easier to implement)
if args.plot_performance:
    assert not args.is_data
    # assert args.algorithm == 'viterbi'
# if args.plot_all_best_events:
#     assert args.n_max_queries == 1  # at least for now

# ----------------------------------------------------------------------------------------
def run_simulation(args):
    print 'simulating'
    assert args.parameter_dir != None and args.outfname != None
    assert args.n_max_queries > 0
    random.seed(args.seed)
    n_per_proc = int(float(args.n_max_queries) / args.n_procs)
    all_random_ints = []
    for iproc in range(args.n_procs):  # have to generate these all at once, 'cause each of the subprocesses is going to reset its seed and god knows what happens to our seed at that point
        all_random_ints.append([random.randint(0, sys.maxint) for i in range(n_per_proc)])
    for iproc in range(args.n_procs):
        proc = Process(target=make_events, args=(args, n_per_proc, iproc, all_random_ints[iproc]))
        proc.start()
    while len(active_children()) > 0:
        # print ' wait %s' % len(active_children()),
        sys.stdout.flush()
        time.sleep(1)
    utils.merge_csvs(args.outfname, [args.workdir + '/recombinator-' + str(iproc) + '/' + os.path.basename(args.outfname) for iproc in range(args.n_procs)], cleanup=(not args.no_clean))

# ----------------------------------------------------------------------------------------
def make_events(args, n_events, iproc, random_ints):
    # NOTE all the different seeds! this sucks but is necessary
    reco = Recombinator(args, seed=args.seed+iproc, sublabel=str(iproc))  #, total_length_from_right=args.total_length_from_right)
    for ievt in range(n_events):
        # print ievt,
        # sys.stdout.flush()
        reco.combine(random_ints[ievt])

# ----------------------------------------------------------------------------------------
if args.action == 'simulate' or args.action == 'generate-trees':
    if args.action == 'generate-trees':
        from treegenerator import TreeGenerator, Hist
        treegen = TreeGenerator(args, args.parameter_dir + '/mean-mute-freqs.csv')
        treegen.generate_trees(self.args.outfname)
        sys.exit(0)
    # if not args.no_clean:
    #     os.rmdir(reco.workdir)
    from recombinator import Recombinator
    run_simulation(args)
else:
    from partitiondriver import PartitionDriver
    random.seed(args.seed)

    args.queries = utils.get_arg_list(args.queries, intify=False)
    args.reco_ids = utils.get_arg_list(args.reco_ids, intify=False)
    args.n_max_per_region = utils.get_arg_list(args.n_max_per_region)
    if len(args.n_max_per_region) != 3:
        raise Exception('ERROR n-max-per-region should be of the form \'x:y:z\', but I got' + str(args.n_max_per_region))

    parter = PartitionDriver(args)

    if args.action == 'build-hmms':  # just build hmms without doing anything else -- you wouldn't normally do this
        parter.write_hmms(args.parameter_dir)
    elif args.action == 'cache-parameters':
        parter.cache_parameters()
    elif 'run-' in args.action:
        parter.run_algorithm(args.action.replace('run-', ''))
    elif args.action == 'partition':
        parter.partition()
    else:
        raise Exception('ERROR bad action ' + args.action)

    parter.clean()
