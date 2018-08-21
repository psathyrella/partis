#!/usr/bin/env python
import random
import sys
import time
import argparse
import dendropy
import os

# start = time.time()
# print '      time: %.1f' % (time.time() - start)

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')
sys.path.insert(1, partis_dir + '/packages/baltic')
import utils
import treeutils

# ----------------------------------------------------------------------------------------
def run_lbi(args):
    if args.overwrite:
        raise Exception('not implemented')
    if args.reroot_at_naive:
        assert args.naive_seq_name is not None
        print treeutils.get_ascii_tree(treeutils.get_treestr(args.treefile))
        dendro_tree = treeutils.get_dendro_tree(treefname=args.treefile)
        dendro_tree.reroot_at_node(dendro_tree.find_node_with_taxon_label(args.naive_seq_name), update_bipartitions=True)
        # print dendro_tree.as_ascii_plot(width=100)  # why tf does this show them as all the same depth?
        treestr = dendro_tree.as_string(schema='newick')  #, suppress_rooting=True)
        print treeutils.get_ascii_tree(treestr)
        bio_tree = treeutils.get_bio_tree(treestr=treestr)
    else:
        bio_tree = treeutils.get_bio_tree(treefname=args.treefile)

    treeutils.calculate_LBI(bio_tree, debug=True)

# ----------------------------------------------------------------------------------------
def run_lonr(args):
    if args.lonr_outdir is None:
        raise Exception('have to specify --lonr-outdir')
    if os.path.exists(args.lonr_outdir):
        if args.overwrite:
            utils.prep_dir(args.lonr_outdir, wildlings=['.txt', '.fasta', '.tab', '.phy', '.csv'])
        else:
            print 'output dir exists, not doing anything (override this with --overwrite)'
            return
    else:
        os.makedirs(args.lonr_outdir)

    workdir = '/tmp/%s/%d' % (os.getenv('USER'), random.randint(0,999999))
    os.makedirs(workdir)
    # # installation stuff
    # rcmds = [
    #     'source("https://bioconductor.org/biocLite.R")',
    #     'biocLite("Biostrings")',
    #     'install.packages("seqinr", repos="http://cran.rstudio.com/")',
    # ]
    # utils.run_r(rcmds, workdir)

    r_work_dir = workdir + '/work'
    os.makedirs(r_work_dir)
    rcmds = [
        'source("%s/lonr.R")' % args.lonr_dir,
        'set.seed(1)',  # have only used this for testing a.t.m., but maybe should set the seed to something generally?
        'compute.LONR(method="%s", infile="%s", baseoutdir="%s/", workdir="%s/", outgroup=%s)' % (args.lonr_tree_method, args.seqfile, args.lonr_outdir, r_work_dir,
                                                                                                  ('"%s"' % args.naive_seq_name) if args.reroot_at_naive else 'NULL',
        ),
    ]
    utils.run_r(rcmds, workdir, debug=True)
    os.rmdir(r_work_dir)
    os.rmdir(workdir)

parser = argparse.ArgumentParser()
parser.add_argument('--treefile', help='input tree file name in newick format')
parser.add_argument('--seqfile', help='input fasta file with aligned sequences corresponding to <treefile>')
parser.add_argument('--outfile', help='output file name in yaml format')
parser.add_argument('--naive-seq-name')
parser.add_argument('--reroot-at-naive', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--lonr-tree-method', default='dnapars', choices=['dnapars', 'neighbor'], help='which phylip method should lonr use to infer the tree (maximum parsimony or neighbor-joining)? (their original defaults were dnapars for less than 100 sequences, neighbor for more)')
parser.add_argument('--lonr-outdir', help='directory for the various lonr output files')
parser.add_argument('--lonr-dir', default=partis_dir + '/bin')
args = parser.parse_args()

# run_lbi(args)
run_lonr(args)
