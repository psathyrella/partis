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

# basic actions:
parser.add_argument('--cache-parameters', action='store_true', help='cache parameter counts and hmm files in dir specified by <parameter_dir>')
parser.add_argument('--run-algorithm', choices=['viterbi', 'forward'], help='Run the specified algorithm once')
parser.add_argument('--partition', action='store_true', help='Find the best partition for the given sequences')
parser.add_argument('--simulate', action='store_true', help='Create simulated rearrangement events')
parser.add_argument('--build-hmms', action='store_true', help='just build hmms (and write \'em out) from existing parameter csvs')
parser.add_argument('--generate-trees', action='store_true', help='generate trees to pass to bppseqgen for simulation')

# finer action control
# parser.add_argument('--pair', action='store_true', help='Run on every pair of sequences in the input')
parser.add_argument('--n-sets', type=int, default=1, help='Run on sets of sequences of size <n> (i.e. \"k-hmm\")')
parser.add_argument('--all-combinations', action='store_true', help='Run algorithm on *all* possible combinations of the input queries of length <n-sets>')
parser.add_argument('--is-data', action='store_true', help='True if not simulation')
parser.add_argument('--skip-unproductive', action='store_true', help='Skip sequences which Smith-Waterman determines to be unproductive (they have stop codons, are out of frame, etc.)')
parser.add_argument('--plot-performance', action='store_true', help='Write out plots comparing true and inferred distributions')
parser.add_argument('--truncate-pairs', action='store_true', help='If pairing two sequences (for hamming distance or hmm pair scores) of different length, truncate the left side of the longer one.')
parser.add_argument('--naivety', default='M', choices=['N', 'M'], help='Naive or mature sequences?')
parser.add_argument('--seed', type=int, default=int(time.time()), help='Random seed for use by recombinator (to allow reproducibility)')
parser.add_argument('--branch-length-multiplier', type=float, help='Multiply observed branch lengths by some factor when simulting, e.g. if in data it was 0.05, but you want ten percent in your simulation, set this to 2')
parser.add_argument('--plot-all-best-events', action='store_true', help='Plot all of the <n-best-events>, i.e. sample from the posterior')
parser.add_argument('--plot-parameters', action='store_true', help='Plot inferred parameters?')
parser.add_argument('--mimic-data-read-length', action='store_true', help='Simulate events with the same read length as ovserved in data? (Otherwise use the entire v and j genes)')
parser.add_argument('--baum-welch-iterations', type=int, default=1, help='Number of Baum-Welch-like iterations.')

# input and output locations
parser.add_argument('--seqfile', help='input sequence file')
parser.add_argument('--parameter-dir', required=True, help='Directory to write sample-specific parameters to and/or read \'em from (e.g. mutation freqs)')
parser.add_argument('--datadir', default='data/imgt', help='Directory from which to read non-sample-specific information (e.g. germline genes)')
parser.add_argument('--outfname')
parser.add_argument('--partitionfname')
parser.add_argument('--plotdir', default='/tmp/partis/plots')
parser.add_argument('--ighutil-dir', default=os.getenv('HOME') + '/.local', help='Path to vdjalign executable. The default is where \'pip install --user\' typically puts things')
parser.add_argument('--workdir', default='/tmp/' + os.path.basename(os.getenv('HOME')) + '/hmms/' + str(os.getpid()), help='Temporary working directory (see also <no-clean>)')

# run/batch control
parser.add_argument('--n-procs', type=int, default=1, help='Max number of processes over which to parallelize')
parser.add_argument('--n-fewer-procs', type=int, help='Number of processes for Smith-Waterman and hamming distance')
parser.add_argument('--slurm', action='store_true', help='Run multiple processes with slurm, otherwise just runs them on local machine. NOTE make sure to set <workdir> to something visible on all batch nodes.')
parser.add_argument('--queries', help='Colon-separated list of query names to which we restrict ourselves')
parser.add_argument('--reco-ids', help='Colon-separated list of rearrangement-event IDs to which we restrict ourselves')  # or recombination events
parser.add_argument('--n-max-queries', type=int, default=-1, help='Maximum number of query sequences on which to run (except for simulator, where it\'s the number of rearrangement events)')
parser.add_argument('--only-genes', help='Colon-separated list of genes to which to restrict the analysis')
parser.add_argument('--n-best-events', type=int, default=3, help='Number of best events to print (i.e. n-best viterbi paths)')

# tree generation (see also branch-length-fname)
# NOTE see also branch-length-multiplier, although that comes into play after the trees are generated
parser.add_argument('--n-trees', type=int, default=100, help='Number of trees to generate')
parser.add_argument('--n-leaves', type=int, default=5, help='Number of leaves per tree')
parser.add_argument('--random-number-of-leaves', action='store_true', help='For each tree choose a random number of leaves in [2, <n-leaves>] (inclusive!). Otherwise give all trees <n-leaves> leaves')

# numerical inputs
parser.add_argument('--hamming-cluster-cutoff', type=float, default=0.5, help='Threshold for hamming distance single-linkage preclustering')
parser.add_argument('--pair-hmm-cluster-cutoff', type=float, default=0.0, help='Threshold for pair hmm single-linkage preclustering')
parser.add_argument('--min_observations_to_write', type=int, default=20, help='For hmmwriter.py, if we see a gene version fewer times than this, we sum over other alleles, or other versions, etc. (see hmmwriter)')
parser.add_argument('--n-max-per-region', default='3:5:2', help='Number of best smith-waterman matches (per region, in the format v:d:j) to pass on to the hmm')
parser.add_argument('--default-v-fuzz', type=int, default=5, help='Size of the k space region over which to sum in the v direction')
parser.add_argument('--default-d-fuzz', type=int, default=2, help='Size of the k space region over which to sum in the d direction')

# temporary arguments (i.e. will be removed as soon as they're not needed)
# parser.add_argument('--tree-parameter-file', default='/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz', help='File from which to read inferred tree parameters (from mebcell analysis)')
parser.add_argument('--gtrfname', default='data/recombinator/gtr.txt', help='File with list of GTR parameters. Fed into bppseqgen along with the chosen tree')
parser.add_argument('--branch-length-fname', default='data/recombinator/branch-lengths.txt', help='Branch lengths from Connor\'s mebcell stuff')
# NOTE command to generate gtr parameter file: [stoat] partis/ > zcat /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz | jq .independentParameters | grep -v '[{}]' | sed 's/["\:,]//g' | sed 's/^[ ][ ]*//' | sed 's/ /,/' | sort >data/gtr.txt

# uncommon arguments
parser.add_argument('--apply-choice_probs_in_sw', action='store_true', help='Apply gene choice probs in Smith-Waterman step. Probably not a good idea (see comments in waterer.py).')
parser.add_argument('--insertion-base-content', default=True, action='store_true',help='Account for non-uniform base content in insertions. Slows us down by a factor around five and gives no performance benefit.')
parser.add_argument('--allow_unphysical_insertions', action='store_true', help='allow insertions on left side of v and right side of j. NOTE this is very slow.')
# parser.add_argument('--allow_external_deletions', action='store_true')     # ( " ) deletions (               "                     )
# parser.add_argument('--total-length-from-right', type=int, default=-1, help='Total read length you want for simulated sequences')
parser.add_argument('--joint-emission', action='store_true', help='Use information about both sequences when writing pair emission probabilities?')

args = parser.parse_args()
args.only_genes = utils.get_arg_list(args.only_genes)
if args.n_fewer_procs == None:
    args.n_fewer_procs = args.n_procs
if args.slurm and '/tmp' in args.workdir:
    print 'ERROR it appears that <workdir> isn\'t set to something visible to all slurm nodes'
    sys.exit()

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
if args.simulate or args.generate_trees:
    if args.generate_trees:
        from treegenerator import TreeGenerator, Hist
        treegen = TreeGenerator(args, args.parameter_dir + '/mean-mute-freqs.csv')
        treegen.generate_trees(self.args.outfname)
        sys.exit(0)
    # if not args.no_clean:
    #     os.rmdir(reco.workdir)
    from recombinator import Recombinator
    run_simulation(args)
else:
    # assert args.cache_parameters or args.point_estimate or args.partition
    from partitiondriver import PartitionDriver

    args.queries = utils.get_arg_list(args.queries, intify=True)
    args.reco_ids = utils.get_arg_list(args.reco_ids, intify=True)
    args.n_max_per_region = utils.get_arg_list(args.n_max_per_region)
    if len(args.n_max_per_region) != 3:
        print 'ERROR n-max-per-region should be form \'x:y:z\', but I got', args.n_max_per_region
        sys.exit()

    utils.prep_dir(args.workdir)
    parter = PartitionDriver(args)

    if args.build_hmms:        
        parter.write_hmms(args.parameter_dir, None)
        sys.exit()

    if args.cache_parameters:
        parter.cache_parameters()
    elif args.run_algorithm != None:
        parter.run_algorithm()
    else:
        parter.partition()
