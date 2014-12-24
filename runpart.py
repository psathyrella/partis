#!/usr/bin/env python
import argparse
import time
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
parser.add_argument('--point-estimate', action='store_true', help='Find the <n-best-events>-best viterbi paths')
parser.add_argument('--partition', action='store_true', help='Find the best partition for the given sequences')
parser.add_argument('--simulate', action='store_true', help='Create simulated rearrangement events')
parser.add_argument('--build-hmms', action='store_true', help='just build hmms (and write \'em out) from existing parameter csvs')
parser.add_argument('--generate-trees', action='store_true', help='generate trees to pass to bppseqgen for simulation')

# finer action control
parser.add_argument('--pair', action='store_true', help='Run on every pair of sequences in the input')
parser.add_argument('--is-data', action='store_true', help='True if not simulation')
parser.add_argument('--skip-unproductive', action='store_true', help='Skip sequences which Smith-Waterman determines to be unproductive (they have stop codons, are out of frame, etc.)')
parser.add_argument('--plot-performance', action='store_true', help='Write out plots comparing true and inferred distributions')
parser.add_argument('--truncate-pairs', action='store_true', help='If pairing two sequences (for hamming distance or hmm pair scores) of different length, truncate the left side of the longer one.')
parser.add_argument('--naivety', default='M', choices=['N', 'M'], help='Naive or mature sequences?')
parser.add_argument('--seed', type=int, help='Random seed for use by recombinator (to allow reproducibility)')

# input and output locations
parser.add_argument('--seqfile', help='input sequence file')
parser.add_argument('--parameter-dir', required=True, help='Directory to write sample-specific parameters to and/or read \'em from (e.g. mutation freqs)')
parser.add_argument('--datadir', required=True, help='Directory from which to read non-sample-specific information (e.g. germline genes)')
parser.add_argument('--outfname')
parser.add_argument('--plotdir', default='/tmp/partis/plots')
parser.add_argument('--ighutil-dir', default=os.getenv('HOME') + '/.local', help='Path to vdjalign executable. The default is where \'pip install --user\' typically puts things')
parser.add_argument('--workdir', default='/tmp/' + os.path.basename(os.getenv('HOME')) + '/hmms/' + str(os.getpid()), help='Temporary working directory (see also <no-clean>)')

# run control
parser.add_argument('--n-procs', type=int, default=1, help='number of processes over which to parallelize')
parser.add_argument('--slurm', action='store_true', help='Run multiple processes with slurm, otherwise just runs them on local machine. NOTE make sure to set <workdir> to something visible on all batch nodes.')
parser.add_argument('--queries', help='Colon-separated list of query names to which we restrict ourselves')
parser.add_argument('--reco-ids', help='Colon-separated list of rearrangement-event IDs to which we restrict ourselves')  # or recombination events
parser.add_argument('--n-max-queries', type=int, default=-1, help='Maximum number of query sequences on which to run (except for simulator, where it\'s the number of rearrangement events)')
parser.add_argument('--only-genes', help='Colon-separated list of genes to which to restrict the analysis')
parser.add_argument('--n-best-events', type=int, default=3, help='Number of best events to print (i.e. n-best viterbi paths)')

# numerical inputs
parser.add_argument('--hamming-cluster-cutoff', type=int, default=0.5, help='Threshold for hamming distance single-linkage preclustering')
parser.add_argument('--pair-hmm-cluster-cutoff', type=int, default=0.0, help='Threshold for pair hmm single-linkage preclustering')
parser.add_argument('--min_observations_to_write', type=int, default=20, help='For hmmwriter.py, if we see a gene version fewer times than this, we sum over other alleles, or other versions, etc. (see hmmwriter)')
parser.add_argument('--n-max-per-region', default='3:5:2', help='Number of best smith-waterman matches (per region, in the format v:d:j) to pass on to the hmm')
parser.add_argument('--default-v-fuzz', type=int, default=5, help='Size of the k space region over which to sum in the v direction')
parser.add_argument('--default-d-fuzz', type=int, default=2, help='Size of the k space region over which to sum in the d direction')

# temporary arguments (i.e. will be removed as soon as they're not needed)
parser.add_argument('--tree-parameter-file', default='/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz', help='File from which to read inferred tree parameters (from mebcell analysis)')
parser.add_argument('--treefname', default='data/recombinator/trees.tre', help='File with list of newick-formated trees. Each rearrangement event chooses one at random')
parser.add_argument('--gtrfname', default='data/recombinator/gtr.txt', help='File with list of GTR parameters. Fed into bppseqgen along with the chosen tree')
# NOTE command to generate gtr parameter file: [stoat] partis/ > zcat /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz | jq .independentParameters | grep -v '[{}]' | sed 's/["\:,]//g' | sed 's/^[ ][ ]*//' | sed 's/ /,/' | sort >data/gtr.txt

# uncommon arguments
parser.add_argument('--apply-choice_probs_in_sw', action='store_true', help='Apply gene choice probs in Smith-Waterman step. Probably not a good idea (see comments in waterer.py).')
parser.add_argument('--insertion-base-content', default=True, action='store_true',help='Account for non-uniform base content in insertions. Slows us down by a factor around five and gives no performance benefit.')
parser.add_argument('--allow_unphysical_insertions', action='store_true', help='allow insertions on left side of v and right side of j. NOTE this is very slow.')
# parser.add_argument('--allow_external_deletions', action='store_true')     # ( " ) deletions (               "                     )
parser.add_argument('--total-length-from-right', type=int, default=-1, help='Total read length you want for simulated sequences')
parser.add_argument('--joint-emission', action='store_true', help='Use information about both sequences when writing pair emission probabilities?')

args = parser.parse_args()
args.only_genes = utils.get_arg_list(args.only_genes)

def run_simulation(args, iproc):
    reco = Recombinator(args, iprocess=iproc, total_length_from_right=args.total_length_from_right)
    for ievt in range(args.n_max_queries):
        print ievt,
        sys.stdout.flush()
        reco.combine()

# ----------------------------------------------------------------------------------------
if args.simulate:
    if args.generate_trees:
        from treegenerator import TreeGenerator, Hist
        treegen = TreeGenerator(args)
        sys.exit(0)
    from recombinator import Recombinator
    assert args.parameter_dir != None and args.outfname != None
    for iproc in range(args.n_procs):
        proc = Process(target=run_simulation, args=(args, iproc))
        proc.start()
    while len(active_children()) > 0:
        # print ' wait %s' % len(active_children()),
        sys.stdout.flush()
        time.sleep(1)
    utils.merge_csvs(args.outfname, [args.workdir + '/recombinator-' + str(iproc) + '/' + os.path.basename(args.outfname) for iproc in range(args.n_procs)], cleanup=(not args.no_clean))
    # if not args.no_clean:
    #     os.rmdir(reco.workdir)
        
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
    elif args.point_estimate:
        parter.point_estimate()
    else:
        parter.partition()
