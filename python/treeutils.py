import __builtin__
import operator
import string
import itertools
import copy
import collections
import random
import csv
from cStringIO import StringIO
import subprocess
import tempfile
import os
import numpy
import sys
from distutils.version import StrictVersion
import dendropy
import time
import math
if StrictVersion(dendropy.__version__) < StrictVersion('4.0.0'):  # not sure on the exact version I need, but 3.12.0 is missing lots of vital tree fcns
    raise RuntimeError("dendropy version 4.0.0 or later is required (found version %s)." % dendropy.__version__)

import utils

lb_metrics = collections.OrderedDict(('lb' + let, 'lb ' + lab) for let, lab in (('i', 'index'), ('r', 'ratio')))
affy_keys = {'lbi' : ['affinities'], #, 'relative_affinities'],  # uncomment the relative one to get a separate set of values/plots for affinity relative to the cells alive at the same time as each cell (when sampling ancestors, this lets you isolate the effects of selection at each time point)
             'lbr' : ['affinities']}
default_lb_tau = 0.0025
default_lbr_tau_factor = 20

dummy_str = 'x-dummy-x'

def lb_cons_seq(line):
    return utils.cons_seq(0.01, aligned_seqfos=[{'name' : u, 'seq' : s} for u, s in zip(line['unique_ids'], line['seqs'])], tie_resolver_seq=line['naive_seq'])  # consensus seq fcn for use with lb metrics (defining it here since we need to use it in two different places)
def lb_cons_dist(cons_seq, seq):
    return -utils.hamming_distance(cons_seq, seq)

# ----------------------------------------------------------------------------------------
# NOTE the min lbi is just tau, but I still like doing it this way
lb_bounds = {  # calculated to 17 generations, which is quite close to the asymptote
    400 : {  # seq_len
        0.0030: (0.0030, 0.0331),  # if tau is any bigger than this it doesn't really converge
        0.0025: (0.0025, 0.0176),
        0.0020: (0.0020, 0.0100),
        0.0010: (0.0010, 0.0033),
        0.0005: (0.0005, 0.0015),
    }
}

# ----------------------------------------------------------------------------------------
def normalize_lb_val(metric, lbval, tau, seq_len=400):
    if metric == 'lbr':
        return lbval
    if seq_len not in lb_bounds:
        raise Exception('seq len %d not in cached lb bound values (available: %s)' % (seq_len, lb_bounds.keys()))
    if tau not in lb_bounds[seq_len]:
        raise Exception('tau value %f not in cached lb bound values (available: %s)' % (tau, lb_bounds[seq_len].keys()))
    lbmin, lbmax = lb_bounds[seq_len][tau]
    return (lbval - lbmin) / (lbmax - lbmin)

# ----------------------------------------------------------------------------------------
def get_treestr(treefname):
    with open(treefname) as treefile:
        return '\n'.join(treefile.readlines())

# ----------------------------------------------------------------------------------------
def get_dendro_tree(treestr=None, treefname=None, taxon_namespace=None, schema='newick', ignore_existing_internal_node_labels=False, suppress_internal_node_taxa=False, debug=False):  # specify either <treestr> or <treefname>
    # <ignore_existing_internal_node_labels> is for when you want the internal nodes labeled (which we usually do, since we want to calculate selection metrics for internal nodes), but you also want to ignore the existing internal node labels (e.g. with FastTree output, where they're floats)
    # <suppress_internal_node_taxa> on the other hand is for when you don't want to have taxa for any internal nodes (e.g. when calculating the tree difference metrics, the two trees have to have the same taxon namespace, but since they in general have different internal nodes, the internal nodes can't have taxa)
    assert treestr is None or treefname is None
    if ignore_existing_internal_node_labels and suppress_internal_node_taxa:
        raise Exception('doesn\'t make sense to specify both')
    if treestr is None:
        treestr = get_treestr(treefname)
    if debug:
        print '   getting dendro tree from string:\n     %s' % treestr
        if taxon_namespace is not None:
            print '     and taxon namespace:  %s' % ' '.join([t.label for t in taxon_namespace])
    # dendropy doesn't make taxons for internal nodes by default, so it puts the label for internal nodes in node.label instead of node.taxon.label, but it crashes if it gets duplicate labels, so you can't just always turn off internal node taxon suppression
    dtree = dendropy.Tree.get_from_string(treestr, schema, taxon_namespace=taxon_namespace, suppress_internal_node_taxa=(ignore_existing_internal_node_labels or suppress_internal_node_taxa), preserve_underscores=True, rooting='force-rooted')  # make sure the tree is rooted, to avoid nodes disappearing in remove_dummy_branches() (and proably other places as well)
    label_nodes(dtree, ignore_existing_internal_node_labels=ignore_existing_internal_node_labels, suppress_internal_node_taxa=suppress_internal_node_taxa, debug=debug)  # set internal node labels to any found in <treestr> (unless <ignore_existing_internal_node_labels> is set), otherwise make some up (e.g. aa, ab, ac)

    # # uncomment for more verbosity:
    # check_node_labels(dtree, debug=debug)  # makes sure that for all nodes, node.taxon is not None, and node.label *is* None (i.e. that label_nodes did what it was supposed to, as long as suppress_internal_node_taxa wasn't set)
    # if debug:
    #     print utils.pad_lines(get_ascii_tree(dendro_tree=dtree))

    return dtree

# ----------------------------------------------------------------------------------------
def import_bio_phylo():
    if 'Bio.Phylo' not in sys.modules:
        from Bio import Phylo  # slow af to import
    return sys.modules['Bio.Phylo']

# ----------------------------------------------------------------------------------------
def get_bio_tree(treestr=None, treefname=None, schema='newick'):  # NOTE don't use this in future (all current uses are commented)
    Phylo = import_bio_phylo()
    if treestr is not None:
        return Phylo.read(StringIO(treestr), schema)
    elif treefname is not None:
        with open(treefname) as treefile:
            return Phylo.read(treefile, schema)
    else:
        assert False

# ----------------------------------------------------------------------------------------
def get_leaf_depths(tree, treetype='dendropy'):  # NOTE structure of dictionary may depend on <treetype>, e.g. whether non-named nodes are included (maybe it doesn't any more? unless you return <clade_keyed_depths> at least)
    if treetype == 'dendropy':
        depths = {n.taxon.label : n.distance_from_root() for n in tree.leaf_node_iter()}
    elif treetype == 'Bio':
        clade_keyed_depths = tree.depths()  # keyed by clade, not clade name (so unlabelled nodes are accessible)
        depths = {n.name : clade_keyed_depths[n] for n in tree.find_clades()}
    else:
        assert False

    return depths

# ----------------------------------------------------------------------------------------
def get_n_leaves(tree):
    return len(tree.leaf_nodes())

# ----------------------------------------------------------------------------------------
def get_n_nodes(tree):
    return len(list(tree.preorder_node_iter()))

# ----------------------------------------------------------------------------------------
def collapse_nodes(dtree, keep_name, remove_name, debug=False):  # collapse edge between <keep_name> and <remove_name>, leaving remaining node with name <keep_name>
    # NOTE I wrote this to try to fix the phylip trees from lonr.r, but it ends up they're kind of unfixable... but this fcn may be useful in the future, I guess, and it works
    if debug:
        print '    collapsing %s and %s (the former will be the label for the surviving node)' % (keep_name, remove_name)
        print utils.pad_lines(get_ascii_tree(dendro_tree=dtree))
    keep_name_node, remove_name_node = [dtree.find_node_with_taxon_label(n) for n in (keep_name, remove_name)]  # nodes corresponding to {keep,remove}_name, not necessarily respectively the nodes we keep/remove
    swapped = False
    if keep_name_node in remove_name_node.child_nodes():
        assert remove_name_node not in keep_name_node.child_nodes()
        parent_node = remove_name_node
        parent_node.taxon.label = keep_name  # have to rename it, since we always actually keep the parent
        swapped = True
        child_node = keep_name_node
    elif remove_name_node in keep_name_node.child_nodes():
        assert keep_name_node not in remove_name_node.child_nodes()
        parent_node = keep_name_node
        child_node = remove_name_node
    else:
        print '    node names %s and %s don\'t share an edge:' % (keep_name, remove_name)
        print '        keep node children: %s' % ' '.join([n.taxon.label for n in keep_name_node.child_nodes()])
        print '      remove node children: %s' % ' '.join([n.taxon.label for n in remove_name_node.child_nodes()])
        raise Exception('see above')

    if child_node.is_leaf():
        dtree.prune_taxa([child_node.taxon], suppress_unifurcations=False)
        if debug:
            print '       pruned leaf node %s' % (('%s (renamed parent to %s)' % (remove_name, keep_name)) if swapped else remove_name)
    else:
        found = False
        for edge in parent_node.child_edge_iter():
            if edge.head_node is child_node:
                edge.collapse()  # removes child node (in dendropy language: inserts all children of the head_node (child) of this edge as children of the edge's tail_node (parent)) Doesn't modify edge lengths by default (i.e. collapsed edge should have zero length).
                found = True
                break
        assert found
        if debug:
            print '     collapsed edge between %s and %s' % (keep_name, remove_name)

    if debug:
        print utils.pad_lines(get_ascii_tree(dendro_tree=dtree))
    assert dtree.find_node_with_taxon_label(remove_name) is None

# ----------------------------------------------------------------------------------------
def check_node_labels(dtree, debug=False):
    if debug:
        print 'checking node labels for:'
        print utils.pad_lines(get_ascii_tree(dendro_tree=dtree, width=250))
    for node in dtree.preorder_node_iter():
        if node.taxon is None:
            raise Exception('taxon is None')
        if debug:
            print '    ok: %s' % node.taxon.label
        if node.label is not None:
            raise Exception('node.label not set to None')

# ----------------------------------------------------------------------------------------
# by default, mostly adds labels to internal nodes (also sometimes the root node) that are missing them
def label_nodes(dendro_tree, ignore_existing_internal_node_labels=False, ignore_existing_internal_taxon_labels=False, suppress_internal_node_taxa=False, initial_length=3, debug=False):
    if ignore_existing_internal_node_labels and suppress_internal_node_taxa:
        raise Exception('doesn\'t make sense to specify both')
    if debug:
        print '   labeling nodes'
        # print '    before:'
        # print utils.pad_lines(get_ascii_tree(dendro_tree))
    tns = dendro_tree.taxon_namespace
    initial_names = set([t.label for t in tns])  # should all be leaf nodes, except the naive sequence (at least for now)
    if debug:
        print '           initial taxon labels: %s' % ' '.join(sorted(initial_names))
    potential_names, used_names = None, None
    new_label, potential_names, used_names = utils.choose_new_uid(potential_names, used_names, initial_length=initial_length, shuffle=True)
    skipped_dbg, relabeled_dbg = [], []
    for node in dendro_tree.preorder_node_iter():
        if node.taxon is not None and not (ignore_existing_internal_taxon_labels and not node.is_leaf()):
            skipped_dbg += ['%s' % node.taxon.label]
            assert node.label is None  # if you want to change this, you have to start setting the node labels in build_lonr_tree(). For now, I like having the label in _one_ freaking place
            continue  # already properly labeled

        current_label = node.label
        node.label = None
        if suppress_internal_node_taxa and not node.is_leaf():
            continue

        if current_label is None or ignore_existing_internal_node_labels:
            new_label, potential_names, used_names = utils.choose_new_uid(potential_names, used_names)
        else:
            # turning this off since it's slow, and has been here a while without getting tripped (and I'm pretty sure the tns checks, anyway)
            # if tns.has_taxon_label(current_label):
            #     raise Exception('duplicate node label \'%s\'' % current_label)
            new_label = current_label

        # turning this off since it's slow, and has been here a while without getting tripped (and I'm pretty sure the tns checks, anyway)
        # if tns.has_taxon_label(new_label):
        #     raise Exception('failed labeling internal nodes (chose name \'%s\' that was already in the taxon namespace)' % new_label)

        node.taxon = dendropy.Taxon(new_label)
        tns.add_taxon(node.taxon)
        relabeled_dbg += ['%s' % new_label]

    if debug:
        print '      skipped (already labeled): %s' % ' '.join(sorted(skipped_dbg))
        print '                   (re-)labeled: %s' % ' '.join(sorted(relabeled_dbg))
        # print '   after:'
        # print utils.pad_lines(get_ascii_tree(dendro_tree))

# ----------------------------------------------------------------------------------------
def translate_labels(dendro_tree, translation_pairs, debug=False):
    if debug:
        print get_ascii_tree(dendro_tree=dendro_tree)
    for old_label, new_label in translation_pairs:
        taxon = dendro_tree.taxon_namespace.get_taxon(old_label)
        if taxon is None:
            raise Exception('requested taxon with old name \'%s\' not present in tree' % old_label)
        taxon.label = new_label
        if debug:
            print '%20s --> %s' % (old_label, new_label)
    if debug:
        print get_ascii_tree(dendro_tree=dendro_tree)

# ----------------------------------------------------------------------------------------
def get_mean_leaf_height(tree=None, treestr=None):
    assert tree is None or treestr is None
    if tree is None:
        tree = get_dendro_tree(treestr=treestr, schema='newick')
    heights = get_leaf_depths(tree).values()
    return sum(heights) / len(heights)

# ----------------------------------------------------------------------------------------
def get_ascii_tree(dendro_tree=None, treestr=None, treefname=None, extra_str='', width=200, schema='newick'):
    """
        AsciiTreePlot docs (don't show up in as_ascii_plot()):
            plot_metric : str
                A string which specifies how branches should be scaled, one of:
                'age' (distance from tips), 'depth' (distance from root),
                'level' (number of branches from root) or 'length' (edge
                length/weights).
            show_internal_node_labels : bool
                Whether or not to write out internal node labels.
            leaf_spacing_factor : int
                Positive integer: number of rows between each leaf.
            width : int
                Force a particular display width, in terms of number of columns.
            node_label_compose_fn : function object
                A function that takes a Node object as an argument and returns
                the string to be used to display it.
    """
    if dendro_tree is None:
        assert treestr is None or treefname is None
        if treestr is None:
            treestr = get_treestr(treefname)
        dendro_tree = get_dendro_tree(treestr=treestr, schema=schema)
    if get_mean_leaf_height(dendro_tree) == 0.:  # we really want the max height, but since we only care whether it's zero or not this is the same
        return '%szero height' % extra_str
    # elif: get_n_nodes(dendro_tree) > 1:  # not sure if I really need this if any more (it used to be for one-leaf trees (and then for one-node trees), but the following code (that used to be indented) seems to be working fine on one-leaf, one-node, and lots-of-node trees a.t.m.)

    start_char, end_char = '', ''
    def compose_fcn(x):
        if x.taxon is not None:  # if there's a taxon defined, use its label
            lb = x.taxon.label
        elif x.label is not None:  # use node label
            lb = x.label
        else:
            lb = 'o'
        return '%s%s%s' % (start_char, lb, end_char)
    dendro_str = dendro_tree.as_ascii_plot(width=width, plot_metric='length', show_internal_node_labels=True, node_label_compose_fn=compose_fcn)
    special_chars = [c for c in reversed(string.punctuation) if c not in set(dendro_str)]  # find some special characters that we can use to identify the start and end of each label (could also use non-printable special characters, but it shouldn't be necessary)
    if len(special_chars) >= 2:  # can't color them directly, since dendropy counts the color characters as printable
        start_char, end_char = special_chars[:2]  # NOTE the colors get screwed up when dendropy overlaps labels (or sometimes just straight up strips stuff), which it does when it runs out of space
        dendro_str = dendro_tree.as_ascii_plot(width=width, plot_metric='length', show_internal_node_labels=True, node_label_compose_fn=compose_fcn)  # call again after modiying compose fcn (kind of wasteful to call it twice, but it shouldn't make a difference)
        dendro_str = dendro_str.replace(start_char, utils.Colors['blue']).replace(end_char, utils.Colors['end'] + '  ')
    else:
        print '  %s can\'t color tree, no available special characters in get_ascii_tree()' % utils.color('red', 'note:')
    if get_n_nodes(dendro_tree) == 1:
        extra_str += ' (one node)'
    return_lines = [('%s%s' % (extra_str, line)) for line in dendro_str.split('\n')]
    return '\n'.join(return_lines)

# ----------------------------------------------------------------------------------------
def rescale_tree(new_mean_height, dtree=None, treestr=None, debug=False):
    # NOTE if you pass in <dtree>, it gets modified, but if you pass in <treestr> you get back a new dtree (which is kind of a dumb way to set this up, but I don't want to change it now. Although I guess it returns None if you pass <dtree>, so you shouldn't get in too much trouble)
    # TODO (maybe) switch calls of this to dendro's scale_edges() (but note you'd then have to get the mean depth beforehand, since that just multiplies by factor, whereas this rescales to get a particular new height)
    """ rescale the branch lengths in dtree/treestr by a factor such that the new mean height is <new_mean_height> """
    if dtree is None:
        dtree = get_dendro_tree(treestr=treestr, suppress_internal_node_taxa=True)
    mean_height = get_mean_leaf_height(tree=dtree)
    if debug:
        print '  current mean: %.4f   target height: %.4f' % (mean_height, new_mean_height)
    for edge in dtree.postorder_edge_iter():
        if edge.head_node is dtree.seed_node:  # why tf does the root node have an edge where it's the child?
            continue
        if debug:
            print '     %5s  %7e  -->  %7e' % (edge.head_node.taxon.label if edge.head_node.taxon is not None else 'None', edge.length, edge.length * new_mean_height / mean_height)
        edge.length *= new_mean_height / mean_height  # rescale every branch length in the tree by the ratio of desired to existing height (everybody's heights should be the same... but they never quite were when I was using Bio.Phylo, so, uh. yeah, uh. not sure what to do, but this is fine. It's checked below, anyway)
    dtree.update_bipartitions()  # probably doesn't really need to be done
    if debug:
        print '    final mean: %.4f' % get_mean_leaf_height(tree=dtree)
    if treestr:
        return dtree.as_string(schema='newick').strip()

# ----------------------------------------------------------------------------------------
def get_tree_difference_metrics(region, in_treestr, leafseqs, naive_seq, debug=False):
    taxon_namespace = dendropy.TaxonNamespace()  # in order to compare two trees with the metrics below, the trees have to have the same taxon namespace
    in_dtree = get_dendro_tree(treestr=in_treestr, taxon_namespace=taxon_namespace, suppress_internal_node_taxa=True, debug=debug)
    seqfos = [{'name' : 't%d' % (iseq + 1), 'seq' : seq} for iseq, seq in enumerate(leafseqs)]
    out_dtree = get_fasttree_tree(seqfos, naive_seq=naive_seq, taxon_namespace=taxon_namespace, suppress_internal_node_taxa=True, debug=debug)
    in_height = get_mean_leaf_height(tree=in_dtree)
    out_height = get_mean_leaf_height(tree=out_dtree)
    base_width = 100
    print '  %s: comparing chosen and bppseqgen output trees for' % (utils.color('green', 'full sequence' if region == 'all' else region))
    print '    %s' % utils.color('blue', 'input:')
    print get_ascii_tree(dendro_tree=in_dtree, extra_str='      ', width=base_width)
    print '    %s' % utils.color('blue', 'output:')
    print get_ascii_tree(dendro_tree=out_dtree, extra_str='        ', width=int(base_width*out_height/in_height))
    print '                   heights: %.3f   %.3f' % (in_height, out_height)
    print '      symmetric difference: %d' % dendropy.calculate.treecompare.symmetric_difference(in_dtree, out_dtree)
    print '        euclidean distance: %f' % dendropy.calculate.treecompare.euclidean_distance(in_dtree, out_dtree)
    print '              r-f distance: %f' % dendropy.calculate.treecompare.robinson_foulds_distance(in_dtree, out_dtree)

# ----------------------------------------------------------------------------------------
def get_fasttree_tree(seqfos, naive_seq=None, naive_seq_name='XnaiveX', taxon_namespace=None, suppress_internal_node_taxa=False, debug=False):
    if debug:
        print '    running FastTree on %d sequences plus a naive' % len(seqfos)
    with tempfile.NamedTemporaryFile() as tmpfile:
        if naive_seq is not None:
            tmpfile.write('>%s\n%s\n' % (naive_seq_name, naive_seq))
        for sfo in seqfos:
            tmpfile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))  # NOTE the order of the leaves/names is checked when reading bppseqgen output
        tmpfile.flush()  # BEWARE if you forget this you are fucked
        with open(os.devnull, 'w') as fnull:
            treestr = subprocess.check_output('./bin/FastTree -gtr -nt ' + tmpfile.name, shell=True, stderr=fnull)
    if debug:
        print '      converting FastTree newick string to dendro tree'
    dtree = get_dendro_tree(treestr=treestr, taxon_namespace=taxon_namespace, ignore_existing_internal_node_labels=not suppress_internal_node_taxa, suppress_internal_node_taxa=suppress_internal_node_taxa, debug=debug)
    naive_node = dtree.find_node_with_taxon_label(naive_seq_name)
    if naive_node is not None:
        dtree.reroot_at_node(naive_node, suppress_unifurcations=False, update_bipartitions=True)
    return dtree

# ----------------------------------------------------------------------------------------
# copied from https://github.com/nextstrain/augur/blob/master/base/scores.py
# also see explanation here https://photos.app.goo.gl/gtjQziD8BLATQivR6
def set_lb_values(dtree, tau, only_calc_metric=None, dont_normalize=False, multifo=None, debug=False):
    """
    traverses <dtree> in postorder and preorder to calculate the up and downstream tree length exponentially weighted by distance, then adds them as LBI (and divides as LBR)
    """
    def getmulti(node):  # number of reads with the same sequence
        return multifo.get(node.taxon.label, 1) if multifo is not None else 1  # most all of them should be in there, but for instance I'm not adding the dummy branch nodes

    metrics_to_calc = lb_metrics.keys() if only_calc_metric is None else [only_calc_metric]
    if debug:
        print '    setting %s values with tau %.4f' % (' and '.join(metrics_to_calc), tau)

    initial_labels = set([n.taxon.label for n in dtree.preorder_node_iter()])
    dtree = get_tree_with_dummy_branches(dtree, tau)  # this returns a new dtree, but the old tree is a subtree of the new one (or at least its collection of nodes are), and these nodes get modified by the process (hence the reversal fcn below)

    # calculate clock length (i.e. for each node, the distance to that node's parent)
    for node in dtree.postorder_node_iter():  # postorder vs preorder doesn't matter, but I have to choose one
        if node.parent_node is None:  # root node
            node.clock_length = 0.
        for child in node.child_node_iter():
            child.clock_length = child.distance_from_root() - node.distance_from_root()

    # lbi is the sum of <node.down_polarizer> (downward message from <node>'s parent) and its children's up_polarizers (upward messages)

    # traverse the tree in postorder (children first) to calculate message to parents (i.e. node.up_polarizer)
    for node in dtree.postorder_node_iter():
        node.down_polarizer = 0  # used for <node>'s lbi (this probabably shouldn't be initialized here, since it gets reset in the next loop [at least I think they all do])
        node.up_polarizer = 0  # used for <node>'s parent's lbi (but not <node>'s lbi)
        for child in node.child_node_iter():
            node.up_polarizer += child.up_polarizer
        bl = node.clock_length / tau
        node.up_polarizer *= numpy.exp(-bl)  # sum of child <up_polarizer>s weighted by an exponential decayed by the distance to <node>'s parent
        node.up_polarizer += getmulti(node) * tau * (1 - numpy.exp(-bl))  # add the actual contribution (to <node>'s parent's lbi) of <node>: zero if the two are very close, increasing toward asymptote of <tau> for distances near 1/tau (integral from 0 to l of decaying exponential)

    # traverse the tree in preorder (parents first) to calculate message to children (i.e. child1.down_polarizer)
    for node in dtree.preorder_internal_node_iter():
        for child1 in node.child_node_iter():  # calculate down_polarizer for each of <node>'s children
            child1.down_polarizer = node.down_polarizer  # first sum <node>'s down_polarizer...
            for child2 in node.child_node_iter():  # and the *up* polarizers of any other children of <node>
                if child1 != child2:
                    child1.down_polarizer += child2.up_polarizer  # add the contribution of <child2> to its parent's (<node>'s) lbi (i.e. <child2>'s contribution to the lbi of its *siblings*)
            bl = child1.clock_length / tau
            child1.down_polarizer *= numpy.exp(-bl)  # and decay the previous sum by distance between <child1> and its parent (<node>)
            child1.down_polarizer += getmulti(child1) * tau * (1 - numpy.exp(-bl))  # add contribution of <child1> to its own lbi: zero if it's very close to <node>, increasing to max of <tau> (integral from 0 to l of decaying exponential)

    returnfo = {m : {} for m in metrics_to_calc}
    # go over all nodes and calculate lb metrics (can be done in any order)
    for node in dtree.postorder_node_iter():
        vals = {'lbi' : node.down_polarizer, 'lbr' : 0.}
        for child in node.child_node_iter():
            vals['lbi'] += child.up_polarizer
            vals['lbr'] += child.up_polarizer
        if node.down_polarizer > 0.:
            vals['lbr'] /= node.down_polarizer  # it might make more sense to not include the branch between <node> and its parent in either the numerator or denominator (here it's included in the denominator), but this way I don't have to change any of the calculations above

        if dummy_str in node.taxon.label:
            continue
        if node is dtree.seed_node or node.parent_node is dtree.seed_node:  # second clause is only because of dummy root addition (well, and if we are adding dummy root the first clause doesn't do anything)
            vals['lbr'] = 0.
        for metric in metrics_to_calc:
            returnfo[metric][node.taxon.label] = float(vals[metric]) if dont_normalize else normalize_lb_val(metric, float(vals[metric]), tau)

    if debug:
        max_width = str(max([len(n.taxon.label) for n in dtree.postorder_node_iter()]))
        print ('   %'+max_width+'s %s%s      multi') % ('node', ''.join('     %s' % m for m in metrics_to_calc), 16*' ' if 'lbr' in metrics_to_calc else '')
        for node in dtree.preorder_node_iter():
            if dummy_str in node.taxon.label:
                continue
            multi_str = str(getmulti(node)) if multifo is not None else ''
            lbstrs = ['%8.3f' % returnfo[m][node.taxon.label] for m in metrics_to_calc]
            if 'lbr' in metrics_to_calc:
                lbstrs += [' = %-5.3f / %-5.3f' % (returnfo['lbr'][node.taxon.label] * node.down_polarizer, node.down_polarizer)]
            print ('    %' + max_width + 's  %s    %3s') % (node.taxon.label, ''.join(lbstrs), multi_str)

    # this is maybe time consuming, but I want to leave the tree that was passed in as unmodified as I can (especially since a.t.m. I'm running this fcn twice for lbi/lbr)
    for node in dtree.postorder_node_iter():
        delattr(node, 'clock_length')
        delattr(node, 'up_polarizer')
        delattr(node, 'down_polarizer')

    remove_dummy_branches(dtree, initial_labels)

    return returnfo

# ----------------------------------------------------------------------------------------
def set_multiplicities(dtree, annotation, input_metafo, debug=False):
    def get_multi(uid):
        if input_metafo is None:  # NOTE the input meta file key 'multiplicities' *could* be in the annotation but we *don't* want to use it (at least at the moment, since we haven't yet established rules for precedence with 'duplicates')
            if uid not in annotation['unique_ids']:  # could be from wonky names from lonr.r, also could be from FastTree tree where we don't get inferred intermediate sequences
                return 1
            if 'duplicates' not in annotation: # if 'duplicates' isn't in the annotation, it's probably simulation, but even if not, if there's no duplicate info then assuming multiplicities of 1 should be fine (maybe should add duplicate info to simulation? it wouldn't really make sense though, since we don't collapse duplicates in simulation info)
                return 1
            return len(utils.per_seq_val(annotation, 'duplicates', uid)) + 1
        elif annotation is None:
            if uid not in input_metafo:
                return 1
            return input_metafo[uid]['multiplicity']
        else:
            assert False  # doesn't make sense to set both of 'em

    if annotation is None and input_metafo is None:
        raise Exception('have to get the multiplicity info from somewhere')


    multifo = {}
    for node in dtree.postorder_node_iter():
        multifo[node.taxon.label] = get_multi(node.taxon.label)
    return multifo

# ----------------------------------------------------------------------------------------
def get_tree_with_dummy_branches(old_dtree, tau, n_tau_lengths=10, add_dummy_leaves=False, debug=False): # add long branches above root and/or below each leaf, since otherwise we're assuming that (e.g.) leaf node fitness is zero
    # commenting this since I'm pretty sure I've fixed it, but not removing it since if a similar problem surfaces with dummy branch addition, deep copying is an easy way out
    # zero_length_edges = [e for e in old_dtree.preorder_edge_iter() if e.length == 0 and not e.head_node.is_leaf()]
    # if len(zero_length_edges) > 0:  # rerooting to remove dummy branches screws up the tree in some cases with zero length branches (see comment in that fcn)
    #     old_dtree = copy.deepcopy(old_dtree)  # could maybe do this by default, but it'll probably be really slow on large trees (at least iterating through the trees is; although I suppose maybe deepcopy is smater than that)
    #     print '    %s found %d zero length branches in tree, so deep copying before adding dummy branches (this is probably ok ish, but in general it\'s a bad idea to have zero length branches in your trees): %s' % (utils.color('yellow', 'warning'), len(zero_length_edges), ' '.join([e.head_node.taxon.label for e in zero_length_edges]))

    dummy_edge_length = n_tau_lengths * tau

    new_root_taxon = dendropy.Taxon(dummy_str + '-root')
    old_dtree.taxon_namespace.add_taxon(new_root_taxon)
    new_root_node = dendropy.Node(taxon=new_root_taxon)
    new_dtree = dendropy.Tree(seed_node=new_root_node, taxon_namespace=old_dtree.taxon_namespace, is_rooted=True)

    # then add the entire old tree under this new tree
    new_root_node.add_child(old_dtree.seed_node)
    for edge in new_root_node.child_edge_iter():
        edge.length = dummy_edge_length

    if add_dummy_leaves:  # add dummy child branches to each leaf
        for lnode in new_dtree.leaf_node_iter():
            new_label = '%s-%s' % (dummy_str, lnode.taxon.label)
            tns.add_taxon(dendropy.Taxon(new_label))
            new_child_node = lnode.new_child(taxon=tns.get_taxon(new_label), edge_length=dummy_edge_length)

    zero_len_edge_nodes = [e.head_node for n in new_dtree.preorder_node_iter() for e in n.child_edge_iter() if e.length == 0 and not e.head_node.is_leaf()]  # zero len edges above leaves are fine, since leaves don't count for lbr
    if len(zero_len_edge_nodes) > 0:
        print '    %s found %d zero length internal edges in tree, which means lb ratio may mis-categorize branches: %s' % (utils.color('red', 'warning'), len(zero_len_edge_nodes), ' '.join([n.taxon.label for n in zero_len_edge_nodes]))
        # for node in zero_len_edge_nodes:  # we don't really want to modify the tree this drastically here (and a.t.m. this causes a crash later on), but I'm leaving it as a placeholder for how to remove zero length edges
        #     collapse_nodes(new_dtree, node.taxon.label, node.parent_node.taxon.label)  # keep the child, since it can be a leaf
        # print utils.pad_lines(get_ascii_tree(dendro_tree=new_dtree))

    new_dtree.update_bipartitions(suppress_unifurcations=False)  # not sure if I need this? (suppress_unifurcations is because otherwise it removes the branch between the old and new root nodes)

    if debug:
        print '    added dummy branches to tree:'
        print get_ascii_tree(dendro_tree=new_dtree, extra_str='      ', width=350)

    return new_dtree

# ----------------------------------------------------------------------------------------
def remove_dummy_branches(dtree, initial_labels, add_dummy_leaves=False, debug=False):
    if add_dummy_leaves:
        raise Exception('not implemented (shouldn\'t be too hard, but a.t.m. I don\'t think I\'ll need it)')

    if len(dtree.seed_node.child_nodes()) != 1:
        print '  %s root node has more than one child when removing dummy branches: %s' % (utils.color('yellow', 'warning'), ' '.join([n.taxon.label for n in dtree.seed_node.child_nodes()]))
    new_root_node = dtree.seed_node.child_nodes()[0]
    if debug:
        print '  rerooting at %s' % new_root_node.taxon.label
        print '            current children: %s' % ' '.join([n.taxon.label for n in new_root_node.child_node_iter()])
    # NOTE if the new root has a child separated by a zero-length edge, this reroot call for some reason deletes that child from the tree (both with and without suppress_unifurcations set). After messing around a bunch to try to fix it, the message I'm taking is just that zero length branches (and unifurcations) are a bad idea and I should just forbid them
    # UPDATE I think I was just missing the suppress_unifurcations=False in update_bipartitions(), but leaving these comments here in case there was another problem
    # UPDATE actually the reroot still seems to eat a node sometimes if the tree is unrooted (so adding the extra reroot above)
    # UPDATE this is more or less expectd, from dendropy's perspective; see https://github.com/jeetsukumaran/DendroPy/issues/118
    assert dtree.is_rooted  # make sure it's rooted, to avoid unifurcations getting suppressed (even with the arg set to false)
    dtree.reroot_at_node(new_root_node, suppress_unifurcations=False)  # reroot at old root node
    if debug:
        print '       children after reroot: %s' % ' '.join([n.taxon.label for n in new_root_node.child_node_iter()])
    dtree.prune_taxa_with_labels([dummy_str + '-root'], suppress_unifurcations=False)
    dtree.purge_taxon_namespace()  # I'm sure there's a good reason the previous line doesn't do this
    dtree.update_bipartitions(suppress_unifurcations=False)
    if debug:
        print '        children after purge: %s' % ' '.join([n.taxon.label for n in new_root_node.child_node_iter()])

    final_labels = set([n.taxon.label for n in dtree.preorder_node_iter()])
    if initial_labels != final_labels:  # this was only happening with a zero-length node hanging off root (see above), which probably won't happen any more since I'm now removing zero length (non-leaf) branches in bcr-phylo simulator.py
        print '    %s nodes after dummy branch addition and removal not the same as before:' % utils.color('red', 'error')
        print '       missing: %s' % ' '.join(initial_labels - final_labels)
        print '       extra:   %s' % ' '.join(final_labels - initial_labels)
        print '       tree:'
        print utils.pad_lines(get_ascii_tree(dendro_tree=dtree, width=400))

# ----------------------------------------------------------------------------------------
# check whether 1) node depth and 2) node pairwise distances are super different when calculated with tree vs sequences (not really sure why it's so different sometimes, best guess is fasttree sucks, partly because it doesn't put the root node anywhere near the root of the tree)
def compare_tree_distance_to_shm(dtree, annotation, max_frac_diff=0.5, min_warn_frac=0.25, extra_str=None, debug=False):
    common_nodes = [n for n in dtree.preorder_node_iter() if n.taxon.label in annotation['unique_ids']]
    tdepths, mfreqs, fracs = {}, {}, {}
    for node in common_nodes:
        tdepth = node.distance_from_root()
        mfreq = utils.per_seq_val(annotation, 'mut_freqs', node.taxon.label)
        frac_diff = abs(tdepth - mfreq) / tdepth if tdepth > 0 else 0
        if frac_diff > max_frac_diff:
            key = node.taxon.label
            tdepths[key] = tdepth
            mfreqs[key] = mfreq
            fracs[key] = frac_diff
    if debug or len(fracs) > 0:
        warnstr = utils.color('yellow', 'warning ') if len(fracs) / float(len(common_nodes)) > min_warn_frac else ''
        if debug or warnstr != '':
            print '        %stree depth and mfreq differ by more than %.0f%% for %d/%d nodes%s' % (warnstr, 100*max_frac_diff, len(fracs), len(common_nodes), '' if extra_str is None else ' for %s' % extra_str)
        if debug and len(fracs) > 0:
            print '    tree depth   mfreq    frac diff'
            for key, frac in sorted(fracs.items(), key=operator.itemgetter(1), reverse=True):
                print '      %.4f    %.4f     %.4f     %s' % (tdepths[key], mfreqs[key], frac, key)

    dmatrix = dtree.phylogenetic_distance_matrix()
    dmx_taxa = set(dmatrix.taxon_iter())  # phylogenetic_distance_matrix() seems to only return values for leaves, which maybe I'm supposed to expect?
    tdists, mdists, fracs = {}, {}, {}  # NOTE reusing these names is kind of dangerous
    for n1, n2 in itertools.combinations([n for n in common_nodes if n.taxon in dmx_taxa], 2):
        tdist = dmatrix.distance(n1.taxon, n2.taxon)
        mdist = utils.hamming_fraction(utils.per_seq_val(annotation, 'seqs', n1.taxon.label), utils.per_seq_val(annotation, 'seqs', n2.taxon.label))
        frac_diff = abs(tdist - mdist) / tdist if tdist > 0 else 0
        if frac_diff > max_frac_diff:
            key = (n1.taxon.label, n2.taxon.label)
            tdists[key] = tdist
            mdists[key] = mdist
            fracs[key] = frac_diff
    if debug or len(fracs) > 0:
        warnstr = utils.color('yellow', 'warning ') if len(fracs) / float(len(common_nodes)) > min_warn_frac else ''
        if debug or warnstr != '':
            print '        %spairwise distance from tree and sequence differ by more than %.f%% for %d/%d node pairs%s' % (warnstr, 100*max_frac_diff, len(fracs), 0.5 * len(common_nodes) * (len(common_nodes)-1), '' if extra_str is None else ' for %s' % extra_str)
        if debug and len(fracs) > 0:
            print '          pairwise'
            print '     tree dist  seq dist  frac diff'
            for key, frac_diff in sorted(fracs.items(), key=operator.itemgetter(1), reverse=True):
                print '      %.4f     %.4f    %.4f    %s  %s' % (tdists[key], mdists[key], frac_diff, key[0], key[1])

    if debug:
        print utils.pad_lines(get_ascii_tree(dendro_tree=dtree, width=400))
        utils.print_reco_event(annotation)

# ----------------------------------------------------------------------------------------
def calculate_lb_values(dtree, tau, lbr_tau_factor=None, only_calc_metric=None, dont_normalize=False, annotation=None, input_metafo=None, use_multiplicities=False, extra_str=None, iclust=None, debug=False):
    # if <only_calc_metric> is None, we use <tau> and <lbr_tau_factor> to calculate both lbi and lbr (i.e. with different tau)
    #   - whereas if <only_calc_metric> is set, we use <tau> to calculate only the given metric
    # note that it's a little weird to do all this tree manipulation here, but then do the dummy branch tree manipulation in set_lb_values(), but the dummy branch stuff depends on tau so it's better this way
    # <iclust> is just to give a little more granularity in dbg

    if use_multiplicities:
        print '  %s <use_multiplicities> is turned on in lb metric calculation, which is ok, but you should make sure that you really believe the multiplicity values' % utils.color('red', 'warning')

    # TODO this is too slow
    # if annotation is not None:  # check that the observed shm rate and tree depth are similar (we're still worried that they're different if we don't have the annotation, but we have no way to check it)
    #     compare_tree_distance_to_shm(dtree, annotation, extra_str=extra_str)

    if max(get_leaf_depths(dtree).values()) > 1:  # should only happen on old simulation files
        if annotation is None:
            raise Exception('tree needs rescaling in lb calculation (metrics will be wrong), but no annotation was passed in')
        print '  %s leaf depths greater than 1, so rescaling by sequence length' % utils.color('yellow', 'warning')
        dtree.scale_edges(1. / numpy.mean([len(s) for s in annotation['seqs']]))  # using treeutils.rescale_tree() breaks, it seems because the update_bipartitions() call removes nodes near root on unrooted trees

    if debug:
        print '   calculating %s%s with tree:' % (' and '.join(lb_metrics if only_calc_metric is None else [only_calc_metric]), '' if extra_str is None else ' for %s' % extra_str)
        print utils.pad_lines(get_ascii_tree(dendro_tree=dtree, width=400))

    multifo = None
    if use_multiplicities:
        multifo = set_multiplicities(dtree, annotation, input_metafo, debug=debug)

    treestr = dtree.as_string(schema='newick')  # get this before the dummy branch stuff to make more sure it isn't modified
    normstr = 'unnormalized' if dont_normalize else 'normalized'
    if only_calc_metric is None:
        assert lbr_tau_factor is not None  # has to be set if we're calculating both metrics
        if iclust is None or iclust == 0:
            print '    calculating %s lb metrics%s with tau values %.4f (lbi) and %.4f * %d = %.4f (lbr)' % (normstr, '' if extra_str is None else ' for %s' % extra_str, tau, tau, lbr_tau_factor, tau*lbr_tau_factor)
        lbvals = set_lb_values(dtree, tau, only_calc_metric='lbi', dont_normalize=dont_normalize, multifo=multifo, debug=debug)
        tmpvals = set_lb_values(dtree, tau*lbr_tau_factor, only_calc_metric='lbr', dont_normalize=dont_normalize, multifo=multifo, debug=debug)
        lbvals['lbr'] = tmpvals['lbr']
    else:
        if iclust is None or iclust == 0:
            print '    calculating %s %s with tau %.4f' % (normstr, only_calc_metric, tau)
        lbvals = set_lb_values(dtree, tau, only_calc_metric=only_calc_metric, dont_normalize=dont_normalize, multifo=multifo, debug=debug)
    lbvals['tree'] = treestr

    return lbvals

# ----------------------------------------------------------------------------------------
def set_n_generations(seq_len, tau, n_tau_lengths, n_generations, debug=False):
    if n_generations is None:
        assert n_tau_lengths is not None  # have to specify one or the other
        n_generations = max(1, int(seq_len * tau * n_tau_lengths))
        if debug:
            print '   %d generations = seq_len * tau * n_tau_lengths = %d * %.4f * %d = max(1, int(%.2f))' % (n_generations, seq_len, tau, n_tau_lengths, seq_len * tau * n_tau_lengths)
    else:
        if debug:
            print '   %d generations' % n_generations
    return n_generations

# ----------------------------------------------------------------------------------------
def get_tree_for_lb_bounds(bound, metric, seq_len, tau, n_generations, n_offspring, debug=False):
    dtree = dendropy.Tree(is_rooted=True)  # note that using a taxon namespace while you build the tree is *much* slower than labeling it afterward (and we do need labels when we calculate lb values)
    if bound == 'min':
        leaf_node = dtree.seed_node  # pretty similar to the dummy root stuff
        for igen in range(n_generations):
            leaf_node = leaf_node.new_child(edge_length=1./seq_len)
    elif bound == 'max':
        old_leaf_nodes = [l for l in dtree.leaf_node_iter()]
        assert len(old_leaf_nodes) == 1
        new_leaf_nodes = []
        for igen in range(n_generations):
            for ileaf in range(len(old_leaf_nodes)):
                for ioff in range(n_offspring):
                    new_leaf_nodes += [old_leaf_nodes[ileaf].new_child(edge_length=1./seq_len)]
            old_leaf_nodes = new_leaf_nodes
            new_leaf_nodes = []
    else:
        assert False

    return dtree

# ----------------------------------------------------------------------------------------
def calculate_lb_bounds(seq_len, tau, n_tau_lengths=10, n_generations=None, n_offspring=2, only_metrics=None, btypes=None, debug=False):  # NOTE the min is just tau, but I don't feel like deleting this fcn just to keep clear what the min means
    info = {m : {} for m in lb_metrics}
    n_generations = set_n_generations(seq_len, tau, n_tau_lengths, n_generations, debug=debug)
    for metric in [m for m in lb_metrics if only_metrics is None or m in only_metrics]:
        for bound in [b for b in ['min', 'max'] if btypes is None or b in btypes]:
            if metric == 'lbr' and bound == 'min':  # lbr min is always zero (leaves)
                info[metric][bound] = {metric : 0., 'vals' : None}
                continue
            start = time.time()
            dtree = get_tree_for_lb_bounds(bound, metric, seq_len, tau, n_generations, n_offspring, debug=debug)
            label_nodes(dtree)
            lbvals = calculate_lb_values(dtree, tau, only_calc_metric=metric, dont_normalize=True, debug=debug)
            bfcn = __builtins__[bound]  # min() or max()
            info[metric][bound] = {metric : bfcn(lbvals[metric].values()), 'vals' : lbvals}
            if debug:
                bname, bval = bfcn(lbvals[metric].items(), key=operator.itemgetter(1))
                print '  %s of %d %s values (%.1fs): %s  %.4f' % (bound, len(lbvals[metric]), metric, time.time() - start, bname, bval)

    return info

# ----------------------------------------------------------------------------------------
lonr_files = {  # this is kind of ugly, but it's the cleanest way I can think of to have both this code and the R code know what they're called
    'phy.outfname' : 'phy_out.txt',
    'phy.treefname' : 'phy_tree.nwk',
    'outseqs.fname' : 'outseqs.fasta',
    'edgefname' : 'edges.tab',
    'names.fname' : 'names.tab',
    'lonrfname' : 'lonr.csv',
}

# ----------------------------------------------------------------------------------------
def build_lonr_tree(edgefos, debug=False):
    # NOTE have to build the tree from the edge file, since the lonr code seems to add nodes that aren't in the newick file (which is just from phylip).
    all_nodes = set([e['from'] for e in edgefos] + [e['to'] for e in edgefos])
    effective_root_nodes = set([e['from'] for e in edgefos]) - set([e['to'] for e in edgefos])  # "effective" because it can be in an unrooted tree. Not sure if there's always exactly one node that has no inbound edges though
    if len(effective_root_nodes) != 1:
        raise Exception('too many effective root nodes: %s' % effective_root_nodes)
    root_label = list(effective_root_nodes)[0]  # should be '1' for dnapars
    if debug:
        print '      chose \'%s\' as root node' % root_label
    tns = dendropy.TaxonNamespace(all_nodes)
    root_node = dendropy.Node(taxon=tns.get_taxon(root_label))  # NOTE this sets node.label and node.taxon.label to the same thing, which may or may not be what we want  # label=root_label,    (if you start setting the node labels again, you also have to translate them below)
    dtree = dendropy.Tree(taxon_namespace=tns, seed_node=root_node, is_rooted=True)
    remaining_nodes = copy.deepcopy(all_nodes) - set([root_label])  # a.t.m. I'm not actually using <all_nodes> after this, but I still want to keep them separate in case I start using it

    weight_or_distance_key = 'distance'  # maybe should I be using the 'weight' column? I think they're just proportional though so I guess it shouldn't matter (same thing in the line below) # 
    root_edgefos = [efo for efo in edgefos if efo['from'] == root_label]
    for efo in root_edgefos:
        dtree.seed_node.new_child(taxon=tns.get_taxon(efo['to']), edge_length=efo[weight_or_distance_key])  # label=efo['to'],    (if you start setting the node labels again, you also have to translate them below)
        remaining_nodes.remove(efo['to'])

    while len(remaining_nodes) > 0:
        n_removed = 0  # I think I don't need this any more (it only happened before I remembered to remove the root node), but it doesn't seem like it'll hurt)
        for lnode in dtree.leaf_node_iter():
            children = [efo for efo in edgefos if efo['from'] == lnode.taxon.label]
            if debug > 1 and len(children) > 0:
                print '    adding children to %s:' % lnode.taxon.label
            for chfo in children:
                lnode.new_child(taxon=tns.get_taxon(chfo['to']), edge_length=chfo[weight_or_distance_key])  # label=chfo['to'],   (if you start setting the node labels again, you also have to translate them below)
                remaining_nodes.remove(chfo['to'])
                n_removed += 1
                if debug > 1:
                    print '              %s' % chfo['to']
        if debug > 1:
            print '  remaining: %d' % len(remaining_nodes)
        if len(remaining_nodes) > 0 and n_removed == 0:  # if there's zero remaining, we're just about to break anyway
            if debug > 1:
                print '  didn\'t remove any, so breaking: %s' % remaining_nodes
            break

    return dtree

# ----------------------------------------------------------------------------------------
def parse_lonr(outdir, input_seqfos, naive_seq_name, reco_info=None, debug=False):
    def get_node_type_from_name(name, debug=False):  # internal nodes in simulated trees should be labeled like 'mrca-<stuff>' (has to correspond to what bcr-phylo-benchmark did)
        if 'mrca' in name:
            return 'internal'
        elif 'leaf' in name:
            return 'leaf'
        else:
            if debug:
                print '    not sure of node type for \'%s\'' % name
            return None

    # get lonr names (lonr replaces them with shorter versions, I think because of phylip)
    lonr_names, input_names = {}, {}
    with open(outdir + '/' + lonr_files['names.fname']) as namefile:  # headers: "head	head2"
        reader = csv.DictReader(namefile, delimiter='\t')
        for line in reader:
            if line['head'][0] != 'L' and line['head'] != naive_seq_name:  # internal node
                dummy_int = int(line['head'])  # check that it's just a (string of a) number
                assert line['head2'] == '-'
                continue
            input_names[line['head']] = line['head2']  # head2 is our names
            lonr_names[line['head2']] = line['head']

    def final_name(lonr_name):
        return input_names.get(lonr_name, lonr_name)

    # read edge info (i.e., implicitly, the tree that lonr.r used)
    edgefos = []  # headers: "from    to      weight  distance"
    with open(outdir + '/' + lonr_files['edgefname']) as edgefile:
        reader = csv.DictReader(edgefile, delimiter='\t')
        for line in reader:
            line['distance'] = int(line['distance'])
            line['weight'] = float(line['weight'])
            edgefos.append(line)

    dtree = build_lonr_tree(edgefos, debug=debug)

    # switch leaves to input names
    for node in dtree.leaf_node_iter():
        node.taxon.label = input_names[node.taxon.label]
        assert node.label is None  #   (if you start setting the node labels again, you also have to translate them here)
        # node.label = node.taxon.label  #   (if you start setting the node labels again, you also have to translate them here)

    if debug:
        print utils.pad_lines(get_ascii_tree(dendro_tree=dtree, width=250))

    nodefos = {node.taxon.label : {} for node in dtree.postorder_node_iter()}  # info for each node (internal and leaf), destined for output

    # read the sequences for both leaves and inferred (internal) ancestors
    seqfos = {final_name(sfo['name']) : sfo['seq'] for sfo in utils.read_fastx(outdir + '/' + lonr_files['outseqs.fname'])}
    input_seqfo_dict = {sfo['name'] : sfo['seq'] for sfo in input_seqfos}  # just to make sure lonr didn't modify the input sequences
    for node in dtree.postorder_node_iter():
        label = node.taxon.label
        if label not in seqfos:
            raise Exception('unexpected sequence name %s' % label)
        if node.is_leaf() or label == naive_seq_name:
            if label not in input_seqfo_dict:
                raise Exception('leaf node \'%s\' not found in input seqs' % label)
            if seqfos[label] != input_seqfo_dict[label]:
                print 'input: %s' % input_seqfo_dict[label]
                print ' lonr: %s' % utils.color_mutants(input_seqfo_dict[label], seqfos[label], align=True)
                raise Exception('lonr leaf sequence doesn\'t match input sequence (see above)')
        nodefos[label]['seq'] = seqfos[label]

    # read actual lonr info
    lonrfos = []
    if debug:
        print '     pos  mutation   lonr   syn./a.b.d.    parent   child'
    with open(outdir + '/' + lonr_files['lonrfname']) as lonrfile:  # heads: "mutation,LONR,mutation.type,position,father,son,flag"
        reader = csv.DictReader(lonrfile)
        for line in reader:
            assert len(line['mutation']) == 2
            assert line['mutation.type'] in ('S', 'R')
            assert line['flag'] in ('TRUE', 'FALSE')
            mutation = line['mutation'].upper()  # dnapars has it upper case already, but neighbor has it lower case
            parent_name = final_name(line['father'])
            child_name = final_name(line['son'])
            parent_seq = nodefos[parent_name]['seq']
            pos = int(line['position']) - 1  # switch from one- to zero-indexing
            child_seq = nodefos[child_name]['seq']
            if parent_seq[pos] != mutation[0] or child_seq[pos] != mutation[1]:
                print 'parent: %s' % parent_seq
                print ' child: %s' % utils.color_mutants(parent_seq, child_seq, align=True)
                raise Exception('mutation info (%s at %d) doesn\'t match sequences (see above)' % (mutation, pos))

            lonrfos.append({
                'mutation' : mutation,
                'lonr' : float(line['LONR']),
                'synonymous' : line['mutation.type'] == 'S',
                'position' : pos,
                'parent' : parent_name,
                'child' : child_name,
                'affected_by_descendents' : line['flag'] == 'TRUE',
            })
            if debug:
                lfo = lonrfos[-1]
                print '     %3d     %2s     %5.2f     %s / %s        %4s      %-20s' % (lfo['position'], lfo['mutation'], lfo['lonr'], 'x' if lfo['synonymous'] else ' ', 'x' if lfo['affected_by_descendents'] else ' ', lfo['parent'], lfo['child'])

    # check for duplicate nodes (not sure why lonr.r kicks these, but I should probably collapse them at some point)
    # in simulation, we sample internal nodes, but then lonr.r's tree construction forces these to be leaves, but then frequently they're immediately adjacent to internal nodes in lonr.r's tree... so we try to collapse them
    duplicate_groups = utils.group_seqs_by_value(nodefos.keys(), keyfunc=lambda q: nodefos[q]['seq'])
    duplicate_groups = [g for g in duplicate_groups if len(g) > 1]
    if len(duplicate_groups) > 0:
        n_max = 15
        dbg_str = ',  '.join([' '.join(g) for g in duplicate_groups[:n_max]])  # only print the first 15 of 'em, if there's more
        if len(duplicate_groups) > n_max:
            dbg_str += utils.color('blue', ' [...]')
        print '    collapsing %d groups of nodes with duplicate sequences (probably just internal nodes that were renamed by lonr.r): %s' % (len(duplicate_groups), dbg_str)
    for dgroup in duplicate_groups:
        non_phylip_names = [n for n in dgroup if get_node_type_from_name(n) is not None]
        if len(non_phylip_names) == 0:  # and phylip internal node names are of form str(<integer>), so just choose the first alphabetically, because whatever
            name_to_use = sorted(dgroup)[0]
        elif len(non_phylip_names) == 1:
            name_to_use = non_phylip_names[0]
        else:
            raise Exception('wtf %s (should\'ve been either one or zero non-phylip names)' % non_phylip_names)
        names_to_remove = [n for n in dgroup if n != name_to_use]

        for rname in names_to_remove:  # only info in here a.t.m. is the sequence
            del nodefos[rname]
            # NOTE not collapsing nodes in tree to match <nodefos> (see comment on next line)
            # collapse_nodes(dtree, name_to_use, rname, allow_failure=True, debug=True)  # holy fuckballs this is not worth the effort (it doesn't really work because the tree is too screwed up) [just gave up and added the duplicate info to the return dict]

        for lfo in lonrfos:
            for key in ('parent', 'child'):
                if lfo[key] in names_to_remove:
                    lfo[key] = name_to_use

    return {'tree' : dtree.as_string(schema='newick'), 'nodes' : nodefos, 'values' : lonrfos}

# ----------------------------------------------------------------------------------------
def run_lonr(input_seqfos, naive_seq_name, workdir, tree_method, lonr_code_file=None, phylip_treefile=None, phylip_seqfile=None, seed=1, debug=False):
    if lonr_code_file is None:
        lonr_code_file = os.path.dirname(os.path.realpath(__file__)).replace('/python', '/bin/lonr.r')
    if not os.path.exists(lonr_code_file):
        raise Exception('lonr code file %s d.n.e.' % lonr_code_file)
    if tree_method not in ('dnapars', 'neighbor'):
        raise Exception('unexpected lonr tree method %s' % tree_method)

    # # installation stuff
    # rcmds = [
    #     'source("https://bioconductor.org/biocLite.R")',
    #     'biocLite("Biostrings")',
    #     'install.packages("seqinr", repos="http://cran.rstudio.com/")',
    # ]
    # utils.run_r(rcmds, workdir)

    input_seqfile = workdir + '/input-seqs.fa'
    with open(input_seqfile, 'w') as iseqfile:
        for sfo in input_seqfos:
            iseqfile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))

    existing_phylip_output_str = ''
    if phylip_treefile is not None:  # using existing phylip output, e.g. from cft
        tree = get_dendro_tree(treefname=phylip_treefile)
        edgefos = []
        for node in tree.preorder_node_iter():
            for edge in node.child_edge_iter():
                edgefos.append({'from' : node.taxon.label, 'to' : edge.head_node.taxon.label, 'weight' : edge.length})
        existing_edgefname = workdir + '/edges.csv'
        existing_node_seqfname = workdir + '/infered-node-seqs.fa'
        with open(existing_edgefname, 'w') as edgefile:
            writer = csv.DictWriter(edgefile, ('from', 'to', 'weight'))
            writer.writeheader()
            for line in edgefos:
                writer.writerow(line)
        with open(existing_node_seqfname, 'w') as node_seqfile:
            writer = csv.DictWriter(node_seqfile, ('head', 'seq'))
            writer.writeheader()
            for sfo in utils.read_fastx(phylip_seqfile):
                writer.writerow({'head' : sfo['name'], 'seq' : sfo['seq']})
        existing_phylip_output_str = ', existing.edgefile="%s", existing.node.seqfile="%s"' % (existing_edgefname, existing_node_seqfname)

    rcmds = [
        'source("%s")' % lonr_code_file,
        'set.seed(%d)' % seed,
        'G.phy.outfname = "%s"'  % lonr_files['phy.outfname'],  # this is a pretty shitty way to do this, but the underlying problem is that there's too many files, but I don't want to parse them all into one or two files in R, so I need to pass all of 'em to the calling python script
        'G.phy.treefname = "%s"' % lonr_files['phy.treefname'],
        'G.outseqs.fname = "%s"' % lonr_files['outseqs.fname'],
        'G.edgefname = "%s"'     % lonr_files['edgefname'],
        'G.names.fname = "%s"'   % lonr_files['names.fname'],
        'G.lonrfname = "%s"'     % lonr_files['lonrfname'],
        'compute.LONR(method="%s", infile="%s", workdir="%s/", outgroup="%s"%s)' % (tree_method, input_seqfile, workdir, naive_seq_name, existing_phylip_output_str),
    ]
    outstr, errstr = utils.run_r(rcmds, workdir, extra_str='      ', return_out_err=True, debug=debug)
    if debug:
        print utils.pad_lines(outstr)
        print utils.pad_lines(errstr)

    os.remove(input_seqfile)
    if phylip_treefile is not None:
        os.remove(existing_edgefname)
        os.remove(existing_node_seqfname)

# ----------------------------------------------------------------------------------------
def calculate_liberman_lonr(input_seqfos=None, line=None, reco_info=None, phylip_treefile=None, phylip_seqfile=None, tree_method=None, naive_seq_name='X-naive-X', seed=1, debug=False):
    # NOTE see issues/notes in bin/lonr.r
    if phylip_treefile is not None or phylip_seqfile is not None:
        raise Exception('never got this (passing phylip output files to lonr.r) to work -- lonr.r kept barfing, although if you were running exactly the same phylip commands as lonr.r does, it would probably work.')
    assert input_seqfos is None or line is None
    if input_seqfos is None:
        input_seqfos = [{'name' : line['unique_ids'][iseq], 'seq' : line['seqs'][iseq]} for iseq in range(len(line['unique_ids']))]
        input_seqfos.insert(0, {'name' : naive_seq_name, 'seq' : line['naive_seq']})
    if tree_method is None:
        tree_method = 'dnapars' if len(input_seqfos) < 500 else 'neighbor'

    workdir = utils.choose_random_subdir('/tmp/%s' % os.getenv('USER', default='partis-work'))
    os.makedirs(workdir)

    if debug:
        print '  %s' % utils.color('green', 'lonr:')
    run_lonr(input_seqfos, naive_seq_name, workdir, tree_method, phylip_treefile=phylip_treefile, phylip_seqfile=phylip_seqfile, seed=seed, debug=debug)
    lonr_info = parse_lonr(workdir, input_seqfos, naive_seq_name, reco_info=reco_info, debug=debug)

    for fn in lonr_files.values():
        os.remove(workdir + '/' + fn)
    os.rmdir(workdir)

    return lonr_info

# ----------------------------------------------------------------------------------------
def get_tree_metric_lines(annotations, cpath, reco_info, use_true_clusters, min_overlap_fraction=0.5, debug=False):
    # collect inferred and true events
    lines_to_use, true_lines_to_use = None, None
    if use_true_clusters:  # use clusters from the true partition, rather than inferred one
        assert reco_info is not None
        true_partition = utils.get_true_partition(reco_info)
        print '    using %d true clusters to calculate inferred tree metrics (sizes: %s)' % (len(true_partition), ' '.join(str(l) for l in sorted([len(c) for c in true_partition], reverse=True)))
        if debug:
            print '      choosing    N        N       N         frac       (N chosen)'
            print '       from     true  & chosen = in common  in common   (w/out duplicates)'
        lines_to_use, true_lines_to_use = [], []
        chosen_ustrs = set()  # now that we're using the fraction instead of the raw total, we mostly shouldn't get multiple true clusters corresponding to the same inferred cluster, but maybe it'll still happen occasionally
        for cluster in true_partition:
            true_lines_to_use.append(utils.synthesize_multi_seq_line_from_reco_info(cluster, reco_info))  # note: duplicates (a tiny bit of) code in utils.print_true_events()
            n_max_in_common, max_frac_in_common, ustr_to_use = None, None, None  # look for the inferred cluster that has the most uids in common with this true cluster
            for ustr in set(annotations) - chosen_ustrs:  # order will be different in reco info and inferred clusters
                n_in_common = len(set(utils.uids_and_dups(annotations[ustr])) & set(cluster))  # can't just look for the actual cluster since we collapse duplicates, but bcr-phylo doesn't (but maybe I should throw them out when parsing bcr-phylo output)
                frac_in_common = n_in_common**2 / float(len(utils.uids_and_dups(annotations[ustr])) * len(cluster))  # and have to use frac instead of total to guard against inferred clusters that include several true clusters (reminder: these inferred clusters may have been run with --n-final-clusters 1 or something similar)
                if max_frac_in_common is None or frac_in_common > max_frac_in_common:
                    ustr_to_use = ustr
                    n_max_in_common = n_in_common
                    max_frac_in_common = frac_in_common
            if max_frac_in_common is None:
                raise Exception('cluster \'%s\' not found in inferred annotations (probably because use_true_clusters was set)' % ':'.join(cluster))
            if max_frac_in_common < min_overlap_fraction:
                raise Exception('overlap fraction %.3f too small: for true cluster (size %d), highest was for inferred cluster with size %d (%d including duplicates). Maybe need to set --simultaneous-true-clonal-seqs.' % (max_frac_in_common, len(cluster), len(annotations[ustr_to_use]['unique_ids']), len(utils.uids_and_dups(annotations[ustr_to_use]))))
            if debug:
                print '      %4d     %4d     %4d     %4d        %4.2f        (%d)' % (len(set(annotations) - chosen_ustrs), len(cluster), len(utils.uids_and_dups(annotations[ustr_to_use])), n_max_in_common, max_frac_in_common, len(annotations[ustr_to_use]['unique_ids']))
            if max_frac_in_common < 1:
                print '            note: couldn\'t find an inferred cluster that corresponded exactly to the true cluster (best was %d & %d = %d (frac %.2f), where the inferred includes %d duplicates)' % (len(utils.uids_and_dups(annotations[ustr_to_use])), len(cluster), n_max_in_common, max_frac_in_common, utils.n_dups(annotations[ustr_to_use]))
            if ustr_to_use in chosen_ustrs:
                raise Exception('chose the same inferred cluster to correspond to two different true clusters')
            chosen_ustrs.add(ustr_to_use)
            lines_to_use.append(annotations[ustr_to_use])
    else:  # use clusters from the inferred partition (whether from <cpath> or <annotations>), and synthesize clusters exactly matching these using single true annotations from <reco_info> (to repeat: these are *not* true clusters)
        lines_to_use = annotations.values()  # we used to restrict it to clusters in the best partition, but I'm switching since I think whenever there are extra ones in <annotations> we always actually want their tree metrics (at the moment there will only be extra ones if either --calculate-alternative-annotations or --write-additional-cluster-annotations are set, but in the future it could also be the default)
        if reco_info is not None:
            for line in lines_to_use:
                true_line = utils.synthesize_multi_seq_line_from_reco_info(line['unique_ids'], reco_info)
                true_lines_to_use.append(true_line)

    return lines_to_use, true_lines_to_use

# ----------------------------------------------------------------------------------------
def plot_tree_metrics(base_plotdir, lines_to_use, true_lines_to_use, ete_path=None, workdir=None, only_csv=False, debug=False):
    import plotting
    import lbplotting
    start = time.time()

    # first make inferred plots
    inf_plotdir = base_plotdir + '/inferred-tree-metrics'
    utils.prep_dir(inf_plotdir, wildlings=['*.svg', '*.html'], allow_other_files=True, subdirs=lb_metrics.keys())
    fnames = []
    lbplotting.plot_lb_vs_affinity(inf_plotdir + '/lbi', lines_to_use, 'lbi', lb_metrics['lbi'], only_csv=only_csv, fnames=fnames, is_true_line=False, debug=debug)
    if not only_csv:  # all the various scatter plots are really slow
        for lb_metric in lb_metrics:
            lbplotting.plot_lb_scatter_plots('shm', inf_plotdir, lb_metric, lines_to_use, fnames=fnames, is_true_line=False, colorvar='is_leaf')
        lbplotting.plot_lb_distributions(inf_plotdir, lines_to_use, fnames=fnames, only_overall=True)
        if ete_path is not None:
            lbplotting.plot_lb_trees(lb_metrics.keys(), inf_plotdir, lines_to_use, ete_path, workdir, is_true_line=False)
        subdirs = [d for d in os.listdir(inf_plotdir) if os.path.isdir(inf_plotdir + '/' + d)]
        plotting.make_html(inf_plotdir, fnames=fnames, new_table_each_row=True, htmlfname=inf_plotdir + '/overview.html', extra_links=[(subd, '%s/%s/' % (inf_plotdir, subd)) for subd in subdirs])

    # then make true plots
    if true_lines_to_use is not None:
        if 'affinities' not in true_lines_to_use[0] or all(affy is None for affy in true_lines_to_use[0]['affinities']):  # if it's bcr-phylo simulation we should have affinities for everybody, otherwise for nobody
            # print '  %s no affinity information in this simulation, so can\'t plot lb/affinity stuff' % utils.color('yellow', 'note')
            print '    tree metric plotting time (no true plots)): %.1f sec' % (time.time() - start)
            return
        true_plotdir = base_plotdir + '/true-tree-metrics'
        utils.prep_dir(true_plotdir, wildlings=['*.svg', '*.html'], allow_other_files=True, subdirs=lb_metrics.keys())
        fnames = []
        if not only_csv:
            for lb_metric in ['lbi']: #lb_metrics:
                lbplotting.plot_lb_scatter_plots('shm', true_plotdir, lb_metric, true_lines_to_use, fnames=fnames, is_true_line=True, colorvar='is_leaf', only_overall=True)
                lbplotting.plot_lb_scatter_plots('consensus', true_plotdir, lb_metric, true_lines_to_use, fnames=fnames, is_true_line=True, colorvar='affinity')
            for lb_metric in lb_metrics:
                lbplotting.make_lb_affinity_joyplots(true_plotdir + '/joyplots', true_lines_to_use, lb_metric, fnames=fnames)
        lbplotting.plot_lb_vs_affinity(true_plotdir + '/lbi', true_lines_to_use, 'lbi', lb_metrics['lbi'], is_true_line=True, only_csv=only_csv, fnames=fnames, debug=debug)
        lbplotting.plot_lb_vs_ancestral_delta_affinity(true_plotdir + '/lbr', true_lines_to_use, 'lbr', lb_metrics['lbr'], is_true_line=True, only_csv=only_csv, fnames=fnames, debug=debug)
        if not only_csv:
            lbplotting.plot_lb_distributions(true_plotdir, true_lines_to_use, fnames=fnames, is_true_line=True, only_overall=True)
            if ete_path is not None:
                lbplotting.plot_lb_trees(lb_metrics.keys(), true_plotdir, true_lines_to_use, ete_path, workdir, is_true_line=True)
            # for lb_metric, lb_label in lb_metrics.items():
            #     XXX fnames[-1] += lbplotting.plot_true_vs_inferred_lb(true_plotdir + '/' + lb_metric, true_lines_to_use, lines_to_use, lb_metric, lb_label)
            subdirs = [d for d in os.listdir(true_plotdir) if os.path.isdir(true_plotdir + '/' + d)]
            plotting.make_html(true_plotdir, fnames=fnames, extra_links=[(subd, '%s/%s/' % (true_plotdir, subd)) for subd in subdirs])

    print '    tree metric plotting time: %.1f sec' % (time.time() - start)

# ----------------------------------------------------------------------------------------
def get_tree_for_line(line, treefname=None, cpath=None, annotations=None, use_true_clusters=False, debug=False):
    # figure out how we want to get the inferred tree
    if treefname is not None:
        dtree = get_dendro_tree(treefname=treefname, debug=debug)
        origin = 'treefname'
    elif False:  # use_liberman_lonr_tree:  # NOTE see issues/notes in bin/lonr.r
        lonr_info = calculate_liberman_lonr(line=line, reco_info=reco_info, debug=debug)
        dtree = get_dendro_tree(treestr=lonr_info['tree'])
        # line['tree-info']['lonr'] = lonr_info
        origin = 'lonr'
    elif cpath is not None and not use_true_clusters and line['unique_ids'] in cpath.partitions[cpath.i_best]:  # if <use_true_clusters> is set, then the clusters in <lines_to_use> won't correspond to the history in <cpath>, so this won't work NOTE now that I've added the direct check if the unique ids are in the best partition, i can probably remove the use_true_clusters check, but I don't want to mess with it a.t.m.
        assert annotations is not None
        i_only_cluster = cpath.partitions[cpath.i_best].index(line['unique_ids'])
        cpath.make_trees(annotations=annotations, i_only_cluster=i_only_cluster, get_fasttrees=True, debug=False)
        dtree = cpath.trees[i_only_cluster]  # as we go through the loop, the <cpath> is presumably filling all of these in
        origin = 'cpath'
    else:
        seqfos = [{'name' : uid, 'seq' : seq} for uid, seq in zip(line['unique_ids'], line['seqs'])]
        dtree = get_fasttree_tree(seqfos, naive_seq=line['naive_seq'], debug=debug)
        origin = 'fasttree'

    return {'tree' : dtree, 'origin' : origin}

# ----------------------------------------------------------------------------------------
def check_lb_values(line, lbvals):
    for metric in [m for m in lbvals if m in lb_metrics]:
        if len(set(line['unique_ids']) - set(lbvals[metric])) > 0:  # we expect to get extra ones, for inferred ancestral nodes for which we don't have sequences
            raise Exception('uids in annotation not the same as lb info keys\n    missing: %s\n    extra: %s' % (' '.join(set(line['unique_ids']) - set(lbvals[metric])), ' '.join(set(lbvals[metric]) - set(line['unique_ids']))))

# ----------------------------------------------------------------------------------------
def calculate_tree_metrics(annotations, min_tree_metric_cluster_size, lb_tau, lbr_tau_factor=None, cpath=None, treefname=None, reco_info=None, use_true_clusters=False, base_plotdir=None,
                           ete_path=None, workdir=None, dont_normalize_lbi=False, only_csv=False, debug=False):
    print 'getting tree metrics'
    if reco_info is not None:
        if not use_true_clusters:
            print '    note: getting tree metrics on simulation without setting <use_true_clusters> (i.e. probably without setting --simultaneous-true-clonal-seqs)'
        for tmpline in reco_info.values():
            assert len(tmpline['unique_ids']) == 1  # at least for the moment, we're splitting apart true multi-seq lines when reading in seqfileopener.py

    lines_to_use, true_lines_to_use = get_tree_metric_lines(annotations, cpath, reco_info, use_true_clusters)

    # get tree and calculate metrics for inferred lines
    n_before = len(lines_to_use)
    lines_to_use = sorted([l for l in lines_to_use if len(l['unique_ids']) >= min_tree_metric_cluster_size], key=lambda l: len(l['unique_ids']), reverse=True)
    n_after = len(lines_to_use)  # after removing the small ones
    tree_origin_counts = {n : {'count' : 0, 'label' : l} for n, l in (('treefname', 'read from %s' % treefname), ('cpath', 'made from cpath'), ('fasttree', 'ran fasttree'), ('lonr', 'ran liberman lonr'))}
    print '    calculating tree metrics for %d cluster%s with size%s: %s' % (n_after, utils.plural(n_after), utils.plural(n_after), ' '.join(str(len(l['unique_ids'])) for l in lines_to_use))
    print '      skipping %d smaller than %d' % (n_before - n_after, min_tree_metric_cluster_size)
    n_already_there = 0
    for iclust, line in enumerate(lines_to_use):
        if debug:
            print '  %s sequence cluster' % utils.color('green', str(len(line['unique_ids'])))
        if 'tree-info' in line:  # NOTE we used to continue here, but now I've decided we really want to overwrite what's there (although I'm a little worried that there was a reason I'm forgetting not to overwrite them)
            if debug:
                print '       %s overwriting tree metric info that was already in <line>' % utils.color('yellow', 'warning')
            n_already_there += 1
        treefo = get_tree_for_line(line, treefname=treefname, cpath=cpath, annotations=annotations, use_true_clusters=use_true_clusters, debug=debug)
        tree_origin_counts[treefo['origin']]['count'] += 1
        line['tree-info'] = {}  # NOTE <treefo> has a dendro tree, but what we put in the <line> (at least for now) is a newick string
        line['tree-info']['lb'] = calculate_lb_values(treefo['tree'], lb_tau, lbr_tau_factor=lbr_tau_factor, annotation=line, dont_normalize=dont_normalize_lbi, extra_str='inf tree', iclust=iclust, debug=debug)
        check_lb_values(line, line['tree-info']['lb'])  # would be nice to remove this eventually, but I keep runnining into instances where dendropy is silently removing nodes

    print '      tree origins: %s' % ',  '.join(('%d %s' % (nfo['count'], nfo['label'])) for n, nfo in tree_origin_counts.items() if nfo['count'] > 0)
    if n_already_there > 0:
        print '    %s overwriting %d / %d that already had tree info' % (utils.color('yellow', 'warning'), n_already_there, n_after)

    # calculate lb values for true lines/trees
    if reco_info is not None:  # note that if <base_plotdir> *isn't* set, we don't actually do anything with the true lb values
        n_true_before = len(true_lines_to_use)
        true_lines_to_use = sorted([l for l in true_lines_to_use if len(l['unique_ids']) >= min_tree_metric_cluster_size], key=lambda l: len(l['unique_ids']), reverse=True)
        n_true_after = len(true_lines_to_use)
        print '    also doing %d true cluster%s with size%s: %s' % (n_true_after, utils.plural(n_true_after), utils.plural(n_true_after), ' '.join(str(len(l['unique_ids'])) for l in true_lines_to_use))
        print '      skipping %d smaller than %d' % (n_true_before - n_true_after, min_tree_metric_cluster_size)
        for iclust, true_line in enumerate(true_lines_to_use):
            true_dtree = get_dendro_tree(treestr=true_line['tree'])
            true_lb_info = calculate_lb_values(true_dtree, lb_tau, lbr_tau_factor=lbr_tau_factor, annotation=true_line, dont_normalize=dont_normalize_lbi, extra_str='true tree', iclust=iclust, debug=debug)
            true_line['tree-info'] = {'lb' : true_lb_info}
            check_lb_values(true_line, true_line['tree-info']['lb'])  # would be nice to remove this eventually, but I keep runnining into instances where dendropy is silently removing nodes

    if base_plotdir is not None:
        assert ete_path is None or workdir is not None  # need the workdir to make the ete trees
        plot_tree_metrics(base_plotdir, lines_to_use, true_lines_to_use, ete_path=ete_path, workdir=workdir, only_csv=only_csv, debug=debug)

# ----------------------------------------------------------------------------------------
def calculate_non_lb_tree_metrics(metric_method, true_lines, min_tree_metric_cluster_size, base_plotdir=None, ete_path=None, workdir=None, lb_tau=None, only_csv=False, debug=False):  # well, not necessarily really using a tree, but they're analagous to the lb metrics
    # NOTE doesn't allow use of relative affinity atm
    # NOTE these true clusters should be identical to the ones in <true_lines_to_use> in the lb metric fcn, but I guess that depends on reco_info and synthesize_multi_seq_line_from_reco_info() and whatnot behaving properly
    n_before = len(true_lines)
    true_lines = sorted([l for l in true_lines if len(l['unique_ids']) >= min_tree_metric_cluster_size], key=lambda l: len(l['unique_ids']), reverse=True)
    n_after = len(true_lines)
    print '      getting non-lb metric %s for %d true cluster%s with size%s: %s' % (metric_method, n_after, utils.plural(n_after), utils.plural(n_after), ' '.join(str(len(l['unique_ids'])) for l in true_lines))
    print '        skipping %d smaller than %d' % (n_before - n_after, min_tree_metric_cluster_size)
    for iclust, true_line in enumerate(true_lines):
        assert 'tree-info' not in true_line  # could handle it, but don't feel like thinking about it a.t.m.
        if metric_method == 'shm':
            metric_info = {u : -utils.per_seq_val(true_line, 'n_mutations', u) for u in true_line['unique_ids']}
            true_line['tree-info'] = {'lb' : {metric_method : metric_info}}
        elif metric_method == 'fay-wu-h':  # NOTE this isn't actually tree info, but I"m comparing it to things calculated with a tree, so putting it in the same place at least for now
            pass
        elif metric_method == 'consensus':
            cseq = lb_cons_seq(true_line)
            true_line['tree-info'] = {'lb' : {metric_method : {u : lb_cons_dist(cseq, s) for u, s in zip(true_line['unique_ids'], true_line['seqs'])}}}
        elif metric_method == 'delta-lbi':  # it would be nice to not calculate lbi here, but atm we're not rewriting the simulation file with true lb info, so it's only in memory even if we've already run get-tree-metrics with the regular lb metrics (anyway, since we already have the tree, it should be really fast)
            dtree = get_dendro_tree(treestr=true_line['tree'])
            lbi_info = calculate_lb_values(dtree, lb_tau, only_calc_metric='lbi', annotation=true_line, extra_str='true tree', iclust=iclust)['lbi']
            delta_lbfo = {}
            for uid in true_line['unique_ids']:
                node = dtree.find_node_with_taxon_label(uid)
                if node is dtree.seed_node:
                    continue  # maybe I should add it as something? not sure
                delta_lbfo[uid] = lbi_info[uid] - lbi_info[node.parent_node.taxon.label]  # I think the parent should always be in here, since I think we should calculate lbi for every node in the tree
            true_line['tree-info'] = {'lb' : {metric_method : delta_lbfo}}
        elif metric_method == 'lbi-cons':  # it would also be nice to not calculate lbi here
            cseq = lb_cons_seq(true_line)
            lbfo = {'consensus' : {u : lb_cons_dist(cseq, s) for u, s in zip(true_line['unique_ids'], true_line['seqs'])}}
            dtree = get_dendro_tree(treestr=true_line['tree'])
            lbi_info = calculate_lb_values(dtree, lb_tau, only_calc_metric='lbi', annotation=true_line, extra_str='true tree', iclust=iclust)['lbi']
            lbfo['lbi'] = {u : lbi_info[u] for u in true_line['unique_ids']}  # remove the ones that aren't in <true_line> (since we don't have sequences for them, so also no consensus distance)
            for lbm in lbfo:  # normalize to z score
                mean, std = numpy.mean(lbfo[lbm].values()), numpy.std(lbfo[lbm].values(), ddof=1)
                lbfo[lbm] = {u : (v - mean) / std for u, v in lbfo[lbm].items()}
            true_line['tree-info'] = {'lb' : {metric_method : {u : (lbfo['lbi'][u] + lbfo['consensus'][u]) / math.sqrt(2) for u in true_line['unique_ids']}}}
        else:
            assert False

    if base_plotdir is not None:
        assert ete_path is None or workdir is not None  # need the workdir to make the ete trees
        import plotting
        import lbplotting
        if 'affinities' not in true_lines[0] or all(affy is None for affy in true_lines[0]['affinities']):  # if it's bcr-phylo simulation we should have affinities for everybody, otherwise for nobody
            return
        true_plotdir = base_plotdir + '/true-tree-metrics'
        utils.prep_dir(true_plotdir, wildlings=['*.svg', '*.html'], allow_other_files=True, subdirs=[metric_method])
        fnames = []
        if metric_method == 'delta-lbi':
            lbplotting.plot_lb_vs_ancestral_delta_affinity(true_plotdir+'/'+metric_method, true_lines, metric_method, lbplotting.mtitlestr('per-seq', metric_method), is_true_line=True, only_csv=only_csv, fnames=fnames, debug=debug)
        else:
            lbplotting.plot_lb_vs_affinity(true_plotdir+'/'+metric_method, true_lines, metric_method, metric_method.upper(), is_true_line=True, affy_key='affinities', only_csv=only_csv, fnames=fnames)
        # lbplotting.plot_lb_distributions(true_plotdir, true_lines, fnames=fnames, is_true_line=True, metric_method=metric_method) #, only_overall=True)
        # if ete_path is not None:
        #     lbplotting.plot_lb_trees([metric_method], true_plotdir, true_lines, ete_path, workdir, is_true_line=True)
        if not only_csv:
            plotting.make_html(true_plotdir, fnames=fnames, extra_links=[(subd, '%s/%s/' % (true_plotdir, subd)) for subd in [metric_method]])

# ----------------------------------------------------------------------------------------
def run_laplacian_spectra(treestr, workdir=None, plotdir=None, plotname=None, title=None, debug=False):
    #  - https://www.ncbi.nlm.nih.gov/pubmed/26658901/
    #  - instructions here: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12526
    # I think this is what ended up working (thought probably not in docker):
    #  apt-get install libgmp-dev libmpfr-dev
    #  > install.packages("RPANDA",dependencies=TRUE)
    #  ok but then I needed to modify the code, so downloaded the source from cran, and swapped out for the spectR.R that eric sent, then installed with:
    # R CMD INSTALL -l packages/RPANDA/lib packages/RPANDA/  # NOTE needs to happen whenever you modify the R source
    # condensation of docs from the above paper:
    #  - > res<-spectR(Phyllostomidae)  # compute eigenvalues (and some metrics describing the distribution, e.g. skewness, kurtosis, eigengap)
    #  - > plot_spectR(res)  # make plots for eigenvalue spectrum
    #  - if eigengap (largest gap between sorted eigenvalues) is e.g. between 3 and 4, then the tree can be separated into three regions, and you use the BIC stuff to find those regions
    #    - > res<-BICompare(Phyllostomidae,3)
    #    - > plot_BICompare(Phyllostomidae,res)
    #  - > res<-JSDtree(Phyllostomidae_genera)  # pairwise jensen-shannon distances between the 25 phylogenies
    #  - > JSDtree_cluster(res)  # plots heatmap and hierarchical cluster

    if debug:
        print utils.pad_lines(get_ascii_tree(treestr=treestr))
        print treestr

    if workdir is None:
        workdir = utils.choose_random_subdir('/tmp/%s' % os.getenv('USER', default='partis-work'))
    eigenfname = '%s/eigenvalues.txt' % workdir
    os.makedirs(workdir)

    cmdlines = [
        'library(ape, quiet=TRUE)',
        # 'library(RPANDA, quiet=TRUE)',  # old way, before I had to modify the source code because the CRAN version removes all eigenvalues <1 (for method="standard" -- with method="normal" it's <0, which is probably better, but it also seems to smoosh all the eigenvalues to be almost exactly 1)
        'library("RPANDA", lib.loc="%s/packages/RPANDA/lib", quiet=TRUE)' % os.path.dirname(os.path.realpath(__file__)).replace('/python', ''),
        'tree <- read.tree(text = "%s")' % treestr,
        # 'print(tree)',
        'specvals <- spectR(tree, method=c("standard"))',  # compute eigenvalues (and some metrics describing the distribution, e.g. skewness, kurtosis, eigengap)
        # 'print(specvals)',
        'capture.output(specvals$eigenvalues, file="%s")' % eigenfname,
    ]

    outstr, errstr = utils.run_r(cmdlines, workdir, return_out_err=True)  # if it crashes, call it without return_out_err, so it prints stuff as it goes
    errstr = '\n'.join([l.strip() for l in errstr.split('\n') if 'This is vegan' not in l])
    for oestr in (outstr, errstr):
        if oestr.strip() == '':
            continue
        print utils.pad_lines(outstr)

    eigenvalues = []
    with open(eigenfname) as efile:
        for line in efile:
            for tstr in line.split():
                if '[' in tstr:
                    if int(tstr.strip('[]')) != len(eigenvalues) + 1:
                        raise Exception('couldn\'t process line:\n%s' % line)
                else:
                    eigenvalues.append(float(tstr))

    os.remove(eigenfname)
    os.rmdir(workdir)

    if plotdir is not None:
        import plotting
        plotting.plot_laplacian_spectra(plotdir, plotname, eigenvalues, title)
