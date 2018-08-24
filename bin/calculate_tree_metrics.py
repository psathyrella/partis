#!/usr/bin/env python
import copy
import csv
import glob
import random
import sys
import time
import argparse
import dendropy
import os
import colored_traceback.always

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
        print '%s --overwrite not implemented for lbi' % utils.color('red', 'warning')
    if args.treefile is None:
        raise Exception('need to set --treefile to run lbi (could instead use tree from lonr, maybe I should implement that?)')

    if args.reroot_at_naive:
        print treeutils.get_ascii_tree(treeutils.get_treestr(args.treefile))
        dendro_tree = treeutils.get_dendro_tree(treefname=args.treefile)
        dendro_tree.reroot_at_node(dendro_tree.find_node_with_taxon_label(args.naive_seq_name), update_bipartitions=True)
        # print dendro_tree.as_ascii_plot(width=100)  # why tf does this show them as all the same depth?
        treestr = dendro_tree.as_string(schema='newick')  #, suppress_rooting=True)
        print treeutils.get_ascii_tree(treestr)
        bio_tree = treeutils.get_bio_tree(treestr=treestr)
    else:
        bio_tree = treeutils.get_bio_tree(treefname=args.treefile)

    treeutils.calculate_LBI(bio_tree)

# ----------------------------------------------------------------------------------------
def parse_lonr(args):
    edgefos = []  # from    to      weight  distance
    with open(args.lonr_outdir + '/' + args.lonr_files['edgefname']) as edgefile:
        reader = csv.DictReader(edgefile, delimiter='\t')
        for line in reader:
            edgefos.append(line)

    # NOTE have to build the tree from the edge file, since the lonr code seems to add nodes that aren't in the newick file (which is just from phylip).
    all_nodes = set([e['from'] for e in edgefos] + [e['to'] for e in edgefos])
    all_nodes.remove('1')
    tns = dendropy.TaxonNamespace(all_nodes)
    dtree = dendropy.Tree(taxon_namespace=tns)
    remaining_nodes = copy.deepcopy(all_nodes)

    root_edgefos = [efo for efo in edgefos if efo['from'] == '1']
    for efo in root_edgefos:
        dtree.seed_node.new_child(label=efo['to'], taxon=tns.get_taxon(efo['to']), edge_length=efo['distance'])
        remaining_nodes.remove(efo['to'])

    while len(remaining_nodes) > 0:
        n_removed = 0  # I think I don't need this any more (it only happened before I remembered to remove the root node), but it doesn't seem like it'll hurt)
        for lnode in dtree.leaf_node_iter():
            children = [efo for efo in edgefos if efo['from'] == lnode.taxon.label]
            if len(children) > 0:
                print '    adding children to %s:' % lnode.taxon.label
            for chfo in children:
                lnode.new_child(label=chfo['to'], taxon=tns.get_taxon(chfo['to']), edge_length=chfo['distance'])
                remaining_nodes.remove(chfo['to'])
                n_removed += 1
                print '              %s' % chfo['to']
        print '  remaining: %d' % len(remaining_nodes)
        if len(remaining_nodes) > 0 and n_removed == 0:  # if there's zero remaining, we're just about to break anyway
            print '  didn\'t remove any, so breaking: %s' % remaining_nodes
            break

    print dtree.as_ascii_plot(width=100)
    sys.exit()

    # lonrfos = []  # mutation,LONR,mutation.type,position,father,son,flag
    # with open(args.lonr_outdir + '/' + args.lonr_files['lonrfname']) as lonrfile:
    #     reader = csv.DictReader(lonrfile)
    #     for line in reader:
    #         lonrfos.append(line)

    # newick_tree = treeutils.get_dendro_tree(treefname=args.lonr_outdir + '/' + args.lonr_files['phy.treefname'])
    # # print newick_tree.as_python_source()
    # for node in newick_tree.postorder_node_iter():
    #     if node.taxon is None:  # internal node, so we have to figure out what phylip called it
    #         continue
    #     print '  %s' % node.taxon.label
    #     parents = [efo['from'] for efo in edgefos if efo['to'] == node.taxon.label]
    #     if len(set(parents)) != 1:
    #         raise Exception('number of parents not equal to 1: %s' % set(parents))
    #     parent_label = parents[0]
    #     if node.parent_node.taxon is None:
    #         node.parent_node.taxon = dendropy.datamodel.taxonmodel.Taxon(label=parent_label)
    #         print '    set parent to %s' % parent_label
    #     elif node.parent_node.taxon.label != parent_label:
    #         raise Exception('inconsistent names for parent node: %s %s' % (node.parent_node.taxon.label, parent_label))

# ----------------------------------------------------------------------------------------
def run_lonr(args):
    if args.treefile is not None:  # TODO
        print 'note: lonr is at the moment still calculating its own trees'

    glob_strs = ['*.txt', '*.fasta', '*.tab', '*.phy', '*.csv', '*.dis', '*.nwk']
    if args.lonr_outdir is None:
        raise Exception('have to specify --lonr-outdir')
    if os.path.exists(args.lonr_outdir):
        if args.overwrite:
            utils.prep_dir(args.lonr_outdir, wildlings=glob_strs)
        else:
            print 'output dir exists, not doing anything (override this with --overwrite)'
            return
    else:
        os.makedirs(args.lonr_outdir)

    workdir = '/tmp/%s/%d' % (os.getenv('USER'), random.randint(0, 999999))
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
        'source("%s")' % args.lonr_code_file,
        'set.seed(1)',  # have only used this for testing a.t.m., but maybe should set the seed to something generally?
        'G.phy.outfname = "%s"'  % args.lonr_files['phy.outfname'],  # this is a pretty shitty way to do this, but the underlying problem is that there's too many files, but I don't want to parse them all into one or two files in R, so I need to pass all of 'em to the calling python script
        'G.phy.treefname = "%s"' % args.lonr_files['phy.treefname'],
        'G.outseqs.fname = "%s"' % args.lonr_files['outseqs.fname'],
        'G.edgefname = "%s"'     % args.lonr_files['edgefname'],
        'G.names.fname = "%s"'   % args.lonr_files['names.fname'],
        'G.lonrfname = "%s"'     % args.lonr_files['lonrfname'],
        'compute.LONR(method="%s", infile="%s", workdir="%s/", outgroup=%s)' % (args.lonr_tree_method, args.seqfile, r_work_dir, ('"%s"' % args.naive_seq_name) if args.reroot_at_naive else 'NULL')
    ]
    utils.run_r(rcmds, workdir, debug=True)
    for fn in args.lonr_files.values():
        os.rename(r_work_dir + '/' + fn, args.lonr_outdir + '/' + fn)
    os.rmdir(r_work_dir)
    os.rmdir(workdir)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
available_metrics = ['lbi', 'lonr']
parser.add_argument('--metrics', default=':'.join(available_metrics), help='colon-separated list of tree metrics to calculate (choose from: %s)' % ' '.join(available_metrics))
parser.add_argument('--reroot-at-naive', action='store_true')
parser.add_argument('--lonr-tree-method', default='dnapars', choices=['dnapars', 'neighbor'], help='which phylip method should lonr use to infer the tree (maximum parsimony or neighbor-joining)? (their original defaults were dnapars for less than 100 sequences, neighbor for more)')
parser.add_argument('--lonr-code-file', default=partis_dir + '/bin/lonr.r')

# input
parser.add_argument('--treefile', help='input tree file name in newick format')
parser.add_argument('--seqfile', help='input fasta file with aligned sequences corresponding to <treefile>')
parser.add_argument('--naive-seq-name', help='uid of inferred naive sequence')

# output
parser.add_argument('--outfile', help='output file name in yaml format')
parser.add_argument('--lonr-outdir', help='directory for the various lonr output files')
parser.add_argument('--overwrite', action='store_true')

args = parser.parse_args()
args.metrics = utils.get_arg_list(args.metrics)
if len(set(args.metrics) - set(available_metrics)) > 0:
    raise Exception('unhandled metric(s): %s (choose from: %s)' % (' '.join(set(args.metrics) - set(available_metrics)), ' '.join(available_metrics)))
if args.reroot_at_naive and args.naive_seq_name is None:
    raise Exception('have to specify --naive-seq-name if --reroot-at-naive is set')
args.lonr_files = {  # this is kind of ugly, but it's the cleanest way I can think of to have both this code and the R code know what they're called
    'phy.outfname' : 'phy_out.txt',
    'phy.treefname' : 'phy_tree.nwk',
    'outseqs.fname' : 'outseqs.fasta',
    'edgefname' : 'edges.tab',
    'names.fname' : 'names.tab',
    'lonrfname' : 'lonr.csv',
}

# ----------------------------------------------------------------------------------------
if 'lbi' in args.metrics:
    run_lbi(args)
if 'lonr' in args.metrics:
    run_lonr(args)
    parse_lonr(args)
