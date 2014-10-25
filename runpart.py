#!/usr/bin/env python
import argparse
import sys
import os
sys.path.insert(1, './python')

import utils

# ----------------------------------------------------------------------------------------
def get_arg_list(arg):  # make lists from args that are passed as strings of colon-separated values
    if arg == None:
        return arg
    else:
        return arg.strip().split(':')  # to allow ids with minus signs, need to add a space, which you then have to strip() off

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
sys.argv.append('-b')
parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
parser.add_argument('--no_clean', action='store_true')  # don't remove the various temp files

# basic actions:
parser.add_argument('--simulate', action='store_true')  # make simulated recombination events
parser.add_argument('--generate_trees', action='store_true')  # generate trees to pass to bppseqgen for simulation
parser.add_argument('--cache_parameters', action='store_true')  # cache parameter counts and hmm files in dir specified by <parameter_dir>
parser.add_argument('--point_estimate', action='store_true')
parser.add_argument('--partition', action='store_true')

# action flags
parser.add_argument('--pair', action='store_true')
parser.add_argument('--is_data', action='store_true')
parser.add_argument('--skip_unproductive', action='store_true')  # skip unproductive rearrangements
parser.add_argument('--apply_choice_probs_in_sw', action='store_true')
parser.add_argument('--plot_performance', action='store_true')

parser.add_argument('--seqfile')  # input file
parser.add_argument('--parameter_dir', required=True)  # sample-specific parameters (mutation rates, gene version freqs, ...)
parser.add_argument('--datadir', default='./data')  # non-sample-specific information (e.g. germline gene versions)
parser.add_argument('--outdir')
parser.add_argument('--plotdir')

parser.add_argument('--n_bases_skip', type=int, default=0)  # number of bases to skip on the left side of the sequence
parser.add_argument('--n_procs', type=int, default=1)  # number of processes over which to parallelize
parser.add_argument('--naivety', default='M', choices=['N', 'M'])
parser.add_argument('--queries')  # restrict to certain query seqs
parser.add_argument('--reco_ids')  # or recombination events
parser.add_argument('--n_max_queries', type=int, default=-1)  # stop after this many queries
parser.add_argument('--only_genes')  # skip all gene matches except for these when parsing the sw output  #'IGHV3-64*04:IGHV1-18*01:IGHV3-23*04:IGHV3-72*01:IGHV5-51*01:IGHD4-23*01:IGHD3-10*01:IGHD4-17*01:IGHD6-19*01:IGHD3-22*01:IGHJ4*02_F:IGHJ5*02_F:IGHJ6*02_F:IGHJ3*02_F:IGHJ2*01_F',

parser.add_argument('--n_max_per_region', type=int, default=5)  # number of best smith-waterman matches (per region) to keep and pass on to the hmm
parser.add_argument('--n_best_events', type=int, default=3)
parser.add_argument('--default_v_fuzz', type=int, default=2)  # TODO play around with these default fuzzes
parser.add_argument('--default_d_fuzz', type=int, default=2)
parser.add_argument('--ighutil_dir', default=os.getenv('HOME') + '/.local')  # this is where '% pip install --user' puts things by default
parser.add_argument('--workdir', default='/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid()))

# temporary arguments (i.e. will be removed as soon as they're not needed)
parser.add_argument('--hackey_extra_data_dir', default='data/recombinator')  # dir for tree parameters that I'm not yet inferring. TODO fix that, obviously
parser.add_argument('--tree_parameter_file', default='/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz')
# NOTE command to generate gtr parameter file: [stoat] partis/ > zcat /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz | jq .independentParameters | grep -v '[{}]' | sed 's/["\:,]//g' | sed 's/^[ ][ ]*//' | sed 's/ /,/' | sort >data/gtr.txt

args = parser.parse_args()
args.only_genes = get_arg_list(args.only_genes)

# ----------------------------------------------------------------------------------------
if args.simulate:
    if args.generate_trees:
        from treegenerator import TreeGenerator, Hist
        treegen = TreeGenerator(args)
        sys.exit(0)
    from recombinator import Recombinator
    assert args.parameter_dir != None and args.outdir != None
    reco = Recombinator(args, total_length_from_right=130)
    for ievt in range(args.n_max_queries):
        print ievt,
        sys.stdout.flush()
        reco.combine()
    os.rmdir(reco.workdir)
        
else:
    # assert args.cache_parameters or args.point_estimate or args.partition
    from partitiondriver import PartitionDriver

    args.queries = get_arg_list(args.queries)
    args.reco_ids = get_arg_list(args.reco_ids)

    utils.prep_dir(args.workdir)
    parter = PartitionDriver(args)

    # parter.write_hmms(args.parameter_dir, None)
    # sys.exit()

    if args.cache_parameters:
        parter.cache_parameters()
    elif args.point_estimate:
        parter.point_estimate()
    else:
        parter.partition()
