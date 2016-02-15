#!/usr/bin/env python
import argparse
import time
import random
import sys
from multiprocessing import Process, active_children
from subprocess import check_output
import os
sys.path.insert(1, './python')

import utils
# merged data: /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_merged.tsv.bz2
parser = argparse.ArgumentParser()
parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
parser.add_argument('--sw-debug', type=int, default=0, choices=[0, 1, 2], help='debug level for smith-waterman')
parser.add_argument('--no-clean', action='store_true', help='Don\'t remove the various temp files')

# basic actions
parser.add_argument('--action', choices=('cache-parameters', 'run-viterbi', 'run-forward', 'partition', 'simulate', 'generate-trees'), help='What do you want to do?')

# finer action control
parser.add_argument('--n-sets', type=int, default=1, help='Run on sets of sequences of size <n> (i.e. \"k-hmm\")')
parser.add_argument('--all-combinations', action='store_true', help='Run algorithm on *all* possible combinations of the input queries of length <n-sets> (otherwise we run on sequential sets of <n-sets> in the input file).')
parser.add_argument('--is-data', action='store_true', help='True if not simulation')
parser.add_argument('--skip-unproductive', action='store_true', help='Skip sequences which Smith-Waterman determines to be unproductive (they have stop codons, are out of frame, etc.)')
parser.add_argument('--plot-performance', action='store_true', help='Write out plots comparing true and inferred distributions')
parser.add_argument('--seed', type=int, default=int(time.time()), help='Random seed for use (mostly) by recombinator (to allow reproducibility)')
parser.add_argument('--mutation-multiplier', type=float, help='Multiply observed branch lengths by some factor when simulating, e.g. if in data it was 0.05, but you want closer to ten percent in your simulation, set this to 2')
parser.add_argument('--mimic-data-read-length', action='store_true', help='trim V 5\' and D 3\' to mimic read lengths seen in data')
parser.add_argument('--annotation-clustering', help='Perform annotation-based clustering from Vollmers paper')
parser.add_argument('--rescale-emissions', action='store_true', default=True)
parser.add_argument('--print-partitions', action='store_true', help='Print partition info in <outfname> and then exit.')
# parser.add_argument('--use_mean_at_boundaries', action='store_true', help='see note in hmmwriter')
parser.add_argument('--annotation-clustering-thresholds', default='0.9', help='colon-separated list of thresholds for annotation-based (e.g. vollmers) clustering')
parser.add_argument('--naive-hamming-bounds')
parser.add_argument('--logprob-ratio-threshold', type=float, default=18., help='reaches a min value of <this> minus five for large clusters.')
parser.add_argument('--naive-vsearch', action='store_true')
parser.add_argument('--naive-swarm', action='store_true')
parser.add_argument('--synthetic-distance-based-partition', action='store_true')
parser.add_argument('--no-indels', action='store_true', help='Skip smith-waterman matches that include indels. Note that this just *skips* them, you probably also want to increase the gap-open penalty to prevent vdjalign from finding them in the first place.')
parser.add_argument('--n-partition-steps', type=int, default=99999, help='Instead of proceeding until we reach 1 process, stop after <n> partitioning steps.')
parser.add_argument('--no-random-divvy', action='store_true', help='Don\'t shuffle the order of the input sequences before passing on to ham')  # it's imperative to shuffle if you're partitioning on simulation, or if you're partitioning with more than one process. But it may also be kinda slow.
parser.add_argument('--naive-hamming', action='store_true', help='agglomerate purely with naive hamming distance, i.e. set the low and high preclustering bounds to the same value')
parser.add_argument('--naivety', default='M', choices=['M', 'N'])
parser.add_argument('--print-cluster-annotations', action='store_true', help='print annotation for each final cluster')
parser.add_argument('--presto-output', action='store_true', help='write output file in presto format')
parser.add_argument('--only-csv-plots', action='store_true', help='only write csv plots')

if 'ralph' in os.getenv('USER'):
    print 'TODO make sure all mpl figures are getting closed'
    print 'clean up print_reco_event, and get it working with indels'
    print 'TODO holy crap stop recursively calling add_implicit_info and its ilk from itselvothers'
    print 'TODO remove all the columns and all_columns checks (and seq/seqs assertions in add_qr_seqs)'
    print 'TODO stop rewriting germlines at every sw step'
    print 'TODO make sure reset_effective_erosions_and_effective_insertions is actually modifying everything it should be'
    print 'TODO remove hist under/overflow warnings when simulating'

# input and output locations
parser.add_argument('--seqfile', help='input sequence file')
parser.add_argument('--parameter-dir', required=True, help='Directory to/from which to write/read sample-specific parameters')
parser.add_argument('--datadir', default=os.getcwd() + '/data/imgt', help='Directory from which to read non-sample-specific information (e.g. germline genes)')
parser.add_argument('--outfname')
parser.add_argument('--plotdir', help='Base directory to which to write plots (no plots are written if this isn\'t set)')
parser.add_argument('--ighutil-dir', default=os.getenv('HOME') + '/.local', help='Path to vdjalign executable. The default is where \'pip install --user\' typically puts things')
parser.add_argument('--workdir', help='Temporary working directory (see also <no-clean>)')
parser.add_argument('--persistent-cachefname')
parser.add_argument('--cache-naive-hfracs', action='store_true')

# run/batch control
parser.add_argument('--n-procs', default='1', help='Max/initial number of processes over which to parallelize (Can be colon-separated list: first number is procs for hmm, second (should be smaller) is procs for smith-waterman, hamming, etc.)')
parser.add_argument('--n-max-procs', default=500, help='never allow more processes than this')
parser.add_argument('--slurm', action='store_true', help='Run multiple processes with slurm, otherwise just runs them on local machine. NOTE make sure to set <workdir> to something visible on all batch nodes.')
parser.add_argument('--queries', help='Colon-separated list of query names to which we restrict ourselves')
parser.add_argument('--reco-ids', help='Colon-separated list of rearrangement-event IDs to which we restrict ourselves')  # or recombination events
parser.add_argument('--seed-unique-id', help='only look for sequences that are clonally related to this unique id')
parser.add_argument('--n-max-queries', type=int, default=-1, help='Maximum number of query sequences on which to run (except for simulator, where it\'s the number of rearrangement events)')
parser.add_argument('--only-genes', help='Colon-separated list of genes to which to restrict the analysis')
parser.add_argument('--n-best-events', default=None, help='Number of best events to print (i.e. n-best viterbi paths). Default is set in bcrham.')

# simulation (see also gtr-fname)
# NOTE see also mutation-multiplier, although that comes into play after the trees are generated
parser.add_argument('--n-sim-events', type=int, default=1, help='Number of rearrangement events to simulate')
parser.add_argument('--n-trees', type=int, default=500, help='Number of trees to generate')
parser.add_argument('--n-leaves', type=int, default=5, help='Number of leaves per tree (used as the mean when drawing from a distribution)')
parser.add_argument('--constant-number-of-leaves', action='store_true', help='Give all trees the same number of leaves (default is to choose each tree\'s number of leaves from a geometric with mean <n_leaves>)')
parser.add_argument('--n-leaf-distribution', default='geometric', choices=['geometric', 'box'], help='distribution from which to draw the number of leaves for each tree')
parser.add_argument('--indel-frequency', default=0., type=float, help='fraction of simulated sequences with indels')
# disabled for now, but if you want multiple indels per sequence you can use this (you'd also need to uncomment a line in recombinator):
# parser.add_argument('--mean-n-indels', default=1, type=int, help='mean number of indels in each sequence which we\'ve already decided has indels (geometric distribution)')
parser.add_argument('--mean-indel-length', default=5, help='mean length of each indel (geometric distribution)')
parser.add_argument('--indel-location', choices=[None, 'v', 'cdr3'], help='where to put the indels')

# numerical inputs
parser.add_argument('--min_observations_to_write', type=int, default=20, help='For hmmwriter.py, if we see a gene version fewer times than this, we sum over other alleles, or other versions, etc. (see hmmwriter)')
parser.add_argument('--n-max-per-region', default='3:5:2', help='Number of best smith-waterman matches (per region, in the format v:d:j) to pass on to the hmm')
parser.add_argument('--default-v-fuzz', type=int, default=5, help='Size of the k space region over which to sum in the v direction')
parser.add_argument('--default-d-fuzz', type=int, default=2, help='Size of the k space region over which to sum in the d direction')
parser.add_argument('--smc-particles', type=int, default=1, help='Number of particles (clustering paths) to simulate with SMC')
parser.add_argument('--gap-open-penalty', type=int, default=30, help='Penalty for indel creation in Smith-Waterman step.')
parser.add_argument('--match-mismatch', default='5:1', help='match:mismatch scores for smith-waterman.')
parser.add_argument('--max-logprob-drop', type=float, default=5., help='stop glomerating when the total logprob has dropped by this much')
parser.add_argument('--n-partitions-to-write', type=int, default=100, help='')

parser.add_argument('--version', action='version', help='print version and exit', version='partis %s' % check_output(['git', 'tag']))

# temporary arguments (i.e. will be removed as soon as they're not needed)
parser.add_argument('--gtrfname', default='data/recombinator/gtr.txt', help='File with list of GTR parameters. Fed into bppseqgen along with the chosen tree')
# NOTE command to generate gtr parameter file: [stoat] partis/ > zcat /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz | jq .independentParameters | grep -v '[{}]' | sed 's/["\:,]//g' | sed 's/^[ ][ ]*//' | sed 's/ /,/' | sort >data/gtr.txt

# uncommon arguments
parser.add_argument('--apply-choice_probs_in_sw', action='store_true', help='Apply gene choice probs in Smith-Waterman step. Probably not a good idea (see comments in waterer.py).')
parser.add_argument('--joint-emission', action='store_true', help='Use information about both sequences when writing pair emission probabilities?')

args = parser.parse_args()
args.only_genes = utils.get_arg_list(args.only_genes)
args.n_procs = utils.get_arg_list(args.n_procs, intify=True)
args.n_fewer_procs = args.n_procs[0] if len(args.n_procs) == 1 else args.n_procs[1]
args.n_procs = args.n_procs[0]

if args.slurm and '/tmp' in args.workdir:
    raise Exception('it appears that <workdir> isn\'t set to something visible to all slurm nodes')

if args.smc_particles != 1:
    raise Exception('sequential monte carlo is not supported at this juncture.')

if args.no_indels and args.gap_open_penalty < 1000:
    print 'forcing gap open to 1000 to prevent indels'
    args.gap_open_penalty = 1000

if args.workdir is None:  # set default here so we know whether it was set by hand or not
    args.workdir = '/tmp/' + os.path.basename(os.getenv('HOME')) + '/hmms/' + str(random.randint(0, 999999))
if os.path.exists(args.workdir):
    raise Exception('workdir %s already exists' % args.workdir)

if args.plot_performance:
    if args.plotdir is None:
        raise Exception('can\'t plot performance unless --plotdir is specified')

# ----------------------------------------------------------------------------------------
def run_simulation(args):
    if args.outfname is None:
        raise Exception('have to specify --outfname for simulation')
    if not os.path.exists(args.parameter_dir):
        raise Exception('parameter dir %s d.n.e.' % args.parameter_dir)
    if not args.n_sim_events > 0:
        raise Exception('--n-sim-events has to be a positivie number')
    if args.slurm:
        raise Exception('simulator parallelization does not handle slurm')

    def make_events(n_events, iproc, random_ints):
        assert n_events > 0
        # NOTE all the different seeds! this sucks but is necessary
        reco = Recombinator(args, seed=args.seed+iproc, sublabel=str(iproc))
        for ievt in range(n_events):
            # print ievt,
            # sys.stdout.flush()
            failed = True
            itry = 0
            while failed:
                if itry > 0:
                    print 'try again: %d' % itry
                failed = not reco.combine(random_ints[ievt] + itry)
                itry += 1

    print 'simulating'
    random.seed(args.seed)
    n_per_proc = int(float(args.n_sim_events) / args.n_procs)
    all_random_ints = []
    for iproc in range(args.n_procs):  # have to generate these all at once, 'cause each of the subprocesses is going to reset its seed and god knows what happens to our seed at that point
        all_random_ints.append([random.randint(0, sys.maxint) for i in range(n_per_proc)])
    for iproc in range(args.n_procs):
        proc = Process(target=make_events, args=(n_per_proc, iproc, all_random_ints[iproc]))
        proc.start()
    while len(active_children()) > 0:
        # print ' wait %s' % len(active_children()),
        sys.stdout.flush()
        time.sleep(1)
    utils.merge_csvs(args.outfname, [args.workdir + '/recombinator-' + str(iproc) + '/' + os.path.basename(args.outfname) for iproc in range(args.n_procs)], cleanup=(not args.no_clean))

if args.action == 'simulate' or args.action == 'generate-trees':
    if args.action == 'generate-trees':
        from treegenerator import TreeGenerator, Hist
        treegen = TreeGenerator(args, args.parameter_dir + '/mean-mute-freqs.csv')
        treegen.generate_trees(args.outfname)
        sys.exit(0)
    # if not args.no_clean:
    #     os.rmdir(reco.workdir)
    from recombinator import Recombinator
    run_simulation(args)
else:
    start = time.time()
    from partitiondriver import PartitionDriver
    random.seed(args.seed)
    args.queries = utils.get_arg_list(args.queries)
    args.reco_ids = utils.get_arg_list(args.reco_ids)
    args.n_max_per_region = utils.get_arg_list(args.n_max_per_region, intify=True)
    args.match_mismatch = utils.get_arg_list(args.match_mismatch, intify=True)
    args.annotation_clustering_thresholds = utils.get_arg_list(args.annotation_clustering_thresholds, floatify=True)
    args.naive_hamming_bounds = utils.get_arg_list(args.naive_hamming_bounds, floatify=True)
    if len(args.n_max_per_region) != 3:
        raise Exception('n-max-per-region should be of the form \'x:y:z\', but I got ' + str(args.n_max_per_region))
    if len(args.match_mismatch) != 2:
        raise Exception('match-mismatch should be of the form \'match:mismatch\', but I got ' + str(args.n_max_per_region))

    parter = PartitionDriver(args)

    if args.action == 'cache-parameters':
        parter.cache_parameters()
    elif 'run-' in args.action:
        parter.run_algorithm(args.action.replace('run-', ''))
    elif args.action == 'partition':
        parter.partition()
    else:
        raise Exception('ERROR bad action ' + args.action)

    parter.clean()
    print '      total time: %.3f' % (time.time()-start)
