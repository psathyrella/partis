#!/usr/bin/env python
import argparse
import sys
import os
sys.path.insert(1, './python')

import utils
# merged data: /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_merged.tsv.bz2

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
sys.argv.append('-b')
parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
parser.add_argument('--no-clean', action='store_true')  # don't remove the various temp files

# basic actions:
parser.add_argument('--simulate', action='store_true')  # make simulated recombination events
parser.add_argument('--build-hmms', action='store_true')  # just build hmms (and write 'em out) from existing parameter csvs
parser.add_argument('--generate-trees', action='store_true')  # generate trees to pass to bppseqgen for simulation
parser.add_argument('--cache-parameters', action='store_true')  # cache parameter counts and hmm files in dir specified by <parameter_dir>
parser.add_argument('--point-estimate', action='store_true')
parser.add_argument('--partition', action='store_true')

# action flags
parser.add_argument('--pair', action='store_true')
parser.add_argument('--is-data', action='store_true')
parser.add_argument('--skip-unproductive', action='store_true')  # skip unproductive rearrangements
parser.add_argument('--apply-choice_probs_in_sw', action='store_true')
parser.add_argument('--plot-performance', action='store_true')
parser.add_argument('--insertion-base-content', default=True, action='store_true')
# TODO tell waterer about these allowances
parser.add_argument('--allow_unphysical_insertions', action='store_true')  # allow insertions on left side of v and right side of j
# parser.add_argument('--allow_external_deletions', action='store_true')     # ( " ) deletions (               "                     )

parser.add_argument('--seqfile')  # input file
parser.add_argument('--parameter-dir', required=True)  # sample-specific parameters (mutation rates, gene version freqs, ...)
parser.add_argument('--datadir', required=True)  # non-sample-specific information (e.g. germline gene versions)
parser.add_argument('--outfname')  # for recombinator, write simulated events to this file. For waterer and partitiondriver, write output to this file
parser.add_argument('--plotdir')

# parser.add_argument('--simu_seed', type=int, default=-1)

parser.add_argument('--chunk-cache', action='store_true')

parser.add_argument('--n-bases-skip', type=int, default=0)  # number of bases to skip on the left side of the sequence
parser.add_argument('--total-length-from-right', type=int, default=-1)  # number of bases to skip on the left side of the sequence
parser.add_argument('--n-procs', type=int, default=1)  # number of processes over which to parallelize
parser.add_argument('--naivety', default='M', choices=['N', 'M'])
parser.add_argument('--queries')  # restrict to certain query seqs
parser.add_argument('--reco-ids')  # or recombination events
parser.add_argument('--n-max-queries', type=int, default=-1)  # stop after this many queries
parser.add_argument('--only-genes')  # skip all gene matches except for these when parsing the sw output  #'IGHV3-64*04:IGHV1-18*01:IGHV3-23*04:IGHV3-72*01:IGHV5-51*01:IGHD4-23*01:IGHD3-10*01:IGHD4-17*01:IGHD6-19*01:IGHD3-22*01:IGHJ4*02_F:IGHJ5*02_F:IGHJ6*02_F:IGHJ3*02_F:IGHJ2*01_F',

parser.add_argument('--min_observations_to_write', type=int, default=20)  # if we see a gene version fewer times than this, we sum over other alleles, or other versions, etc. (see hmmwriter)

parser.add_argument('--j-subset', default=None)  #'imgt')  # which germline j file to use? NOTE has to correspond to a file <vdjalign-install-dir>/imgt/data/ighj-<j_subset>.fasta
parser.add_argument('--n-max-per-region', default='3:5:2')  # number of best smith-waterman matches (per region, in the order v:d:j) to keep and pass on to the hmm
parser.add_argument('--n-best-events', type=int, default=3)
parser.add_argument('--default-v-fuzz', type=int, default=2)  # TODO play around with these default fuzzes
parser.add_argument('--default-d-fuzz', type=int, default=2)
parser.add_argument('--ighutil-dir', default=os.getenv('HOME') + '/.local')  # this is where '% pip install --user' puts things by default
parser.add_argument('--workdir', default='/tmp/' + os.path.basename(os.getenv('HOME')) + '/hmms/' + str(os.getpid()))

# temporary arguments (i.e. will be removed as soon as they're not needed)
parser.add_argument('--hackey-extra-data-dir', default='data/recombinator')  # dir for tree parameters that I'm not yet inferring. TODO fix that, obviously
parser.add_argument('--tree-parameter-file', default='/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz')
# NOTE command to generate gtr parameter file: [stoat] partis/ > zcat /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/A/04-A-M_gtr_tr-qi-gi.json.gz | jq .independentParameters | grep -v '[{}]' | sed 's/["\:,]//g' | sed 's/^[ ][ ]*//' | sed 's/ /,/' | sort >data/gtr.txt

args = parser.parse_args()
args.only_genes = utils.get_arg_list(args.only_genes)

# ----------------------------------------------------------------------------------------
if args.simulate:
    if args.generate_trees:
        from treegenerator import TreeGenerator, Hist
        treegen = TreeGenerator(args)
        sys.exit(0)
    from recombinator import Recombinator
    assert args.parameter_dir != None and args.outfname != None
    reco = Recombinator(args, total_length_from_right=args.total_length_from_right)
    for ievt in range(args.n_max_queries):
        print ievt,
        sys.stdout.flush()
        reco.combine()
    if not args.no_clean:
        os.rmdir(reco.workdir)
else:
    # assert args.cache_parameters or args.point_estimate or args.partition
    from partitiondriver import PartitionDriver

    args.queries = utils.get_arg_list(args.queries)
    args.reco_ids = utils.get_arg_list(args.reco_ids)
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
