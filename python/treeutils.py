import string
import itertools
import copy
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
if StrictVersion(dendropy.__version__) < StrictVersion('4.0.0'):  # not sure on the exact version I need, but 3.12.0 is missing lots of vital tree fcns
    raise RuntimeError("dendropy version 4.0.0 or later is required (found version %s)." % dendropy.__version__)

import baltic
import utils

# ----------------------------------------------------------------------------------------
# two classes to work around the fact baltic doesn't yet support one-leaf trees
class TinyLeaf(object):
    def __init__(self, name, length, height):
        self.name = name
        self.numName = self.name
        self.length = length
        self.height = height
        self.branchType = 'leaf'
class OneLeafTree(object):
    def __init__(self, name, height):
        self.leaves = [TinyLeaf(name, height, height)]
        self.Objects = self.leaves
    def traverse_tree(self):
        self.leaves[0].height = self.leaves[0].length
    def toString(self, numName=None):
        return '%s:%.15f;' % (self.leaves[0].name, self.leaves[0].length)

# ----------------------------------------------------------------------------------------
def get_treestr(treefname):
    with open(treefname) as treefile:
        return '\n'.join(treefile.readlines())

# ----------------------------------------------------------------------------------------
def get_dendro_tree(treestr=None, treefname=None, taxon_namespace=None, schema='newick', ignore_existing_internal_node_labels=False):  # specify either <treestr> or <treefname>
    assert treestr is None or treefname is None
    if treestr is None:
        treestr = get_treestr(treefname)
    # dendropy doesn't make taxons for internal nodes by default, so it puts the label for internal nodes in node.label instead of node.taxon.label, but it crashes if it gets duplicate labels, so you can't just turn off internal node taxon suppression
    dtree = dendropy.Tree.get_from_string(treestr, schema, taxon_namespace=taxon_namespace, suppress_internal_node_taxa=ignore_existing_internal_node_labels)
    # check_node_labels(dtree)
    label_nodes(dtree, ignore_existing_internal_node_labels=ignore_existing_internal_node_labels)  # set internal node labels to any found in <treestr> (unless <ignore_existing_internal_node_labels> is set), otherwise make some up (e.g. aa, ab, ac)
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
def get_baltic_tree(treestr):  # NOTE trying to use dendropy in future, it seems the nicest tree handler
    if treestr.count(':') == 1:  # one-leaf tree
        name, lengthstr = treestr.strip().rstrip(';').split(':')
        tree = OneLeafTree(name, float(lengthstr))
    else:
        tree = baltic.tree()
        baltic.make_tree(treestr, tree, verbose=False)
    tree.traverse_tree()
    return tree

# ----------------------------------------------------------------------------------------
def get_leaf_depths(tree, treetype='dendropy'):  # NOTE structure of dictionary may depend on <treetype>, e.g. whether non-named nodes are included (maybe it doesn't any more? unless you return <clade_keyed_depths> at least)
    if treetype == 'dendropy':
        depths = {n.taxon.label : n.distance_from_root() for n in tree.leaf_node_iter()}
    elif treetype == 'baltic':
        assert False  # l.label doesn't seem to be set. Not sure why tf not
        depths = {l.label : l.height for l in tree.leaves}
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
def check_node_labels(dtree, debug=True):
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
# mostly adds labels to internal nodes, but also sometimes the root node
def label_nodes(dendro_tree, ignore_existing_internal_node_labels=False, debug=False):
    if debug:
        print ' labeling nodes:'
        # print '    before:'
        # print utils.pad_lines(get_ascii_tree(dendro_tree))
    tns = dendro_tree.taxon_namespace
    initial_names = set([t.label for t in tns])  # should all be leaf nodes, except the naive sequence (at least for now)
    potential_names, used_names = None, None
    skipped_dbg, relabeled_dbg = [], []
    for node in dendro_tree.preorder_node_iter():
        if node.taxon is not None:
            skipped_dbg += ['%s' % node.taxon.label]
            assert node.label is None  # if you want to change this, you have to start setting the node labels in build_lonr_tree(). For now, I like having the label in _one_ freaking place
            continue

        if node.label is not None and not ignore_existing_internal_node_labels:  # NOTE results in duplicate labels, if that's what the fasta has
            if tns.has_taxon_label(node.label):
                raise Exception('duplicate node label \'%s\'' % node.label)
            label = node.label
        else:
            label, potential_names, used_names = utils.choose_new_uid(potential_names, used_names)
        node.label = None

        if tns.has_taxon_label(label):  # not sure it really makes sense to check this in two places, but it's nice to see which way is failing I guess
            raise Exception('failed labeling internal nodes (chose a name that was already in the tree)')

        tns.add_taxon(dendropy.Taxon(label))
        node.taxon = tns.get_taxon(label)
        relabeled_dbg += ['%s' % label]

    if debug:
        print '    skipped (already labeled): %s' % '  '.join(skipped_dbg)
        print '                (re-)labeled: %s' % '  '.join(relabeled_dbg)
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
        tree = get_dendro_tree(treestr=treestr, schema='newick')  # if we're calling it with a treestr rather than a tree, it's all from old code that only uses newick
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
    elif get_n_leaves(dendro_tree) > 1:  # if more than one leaf
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
        return_lines = [('%s%s' % (extra_str, line)) for line in dendro_str.split('\n')]
        return '\n'.join(return_lines)

        # Phylo = import_bio_phylo()
        # tmpf = StringIO()
        # if treestr is None:
        #     treestr = dendro_tree.as_string(schema=schema)
        # Phylo.draw_ascii(get_bio_tree(treestr=treestr, schema=schema), file=tmpf, column_width=width)  # this is pretty ugly, but the dendropy ascii printer just puts all the leaves at the same depth
        # return '\n'.join(['%s%s' % (extra_str, line) for line in tmpf.getvalue().split('\n')])
    else:
        return '%sone leaf' % extra_str

# ----------------------------------------------------------------------------------------
def rescale_tree(treestr, new_height, debug=False):  # NOTE assumes newick for now
    """ rescale the branch lengths in <treestr> (newick-formatted) by <factor> """
    baltic_tree = get_baltic_tree(treestr)
    mean_height = get_mean_leaf_height(treestr=treestr)
    if debug:
        print '  current mean: %.4f   target height: %.4f' % (mean_height, new_height)
    for ln in baltic_tree.Objects:
        old_length = ln.length
        ln.length *= new_height / mean_height  # rescale every branch length in the tree by the ratio of desired to existing height (everybody's heights should be the same... but they never quite were when I was using Bio.Phylo, so, uh. yeah, uh. not sure what to do, but this is fine. It's checked below, anyway)
        if debug:
            print '     %5s  %7e  -->  %7e' % (ln.numName if ln.branchType == 'leaf' else ln.branchType, old_length, ln.length)
    baltic_tree.traverse_tree()
    treestr = baltic_tree.toString(numName=True)
    if debug:
        print '    final mean: %.4f' % get_mean_leaf_height(treestr=treestr)
    this_eps = 1e-8
    for leaf in get_baltic_tree(treestr).leaves:  # make sure string conversion (and rescaling) went ok
        if leaf.height < 2 * this_eps or new_height < 2 * this_eps:
            continue
        if not utils.is_normed(leaf.height / new_height, this_eps=this_eps):
            raise Exception('tree not rescaled properly for leaf \'%s\':   %.10f   %.10f    %e' % (leaf.numName, leaf.height, new_height, (leaf.height - new_height) / new_height))
    return treestr

# ----------------------------------------------------------------------------------------
def infer_tree_from_leaves(region, in_treestr, leafseqs, naive_seq, naive_seq_name='XnaiveX', debug=False):  # baltic barfs on (some) dashes
    taxon_namespace = dendropy.TaxonNamespace()  # in order to compare two trees with the metrics below, the trees have to have the same taxon namespace
    with tempfile.NamedTemporaryFile() as tmpfile:
        tmpfile.write('>%s\n%s\n' % (naive_seq_name, naive_seq))
        for iseq in range(len(leafseqs)):
            tmpfile.write('>t%s\n%s\n' % (iseq+1, leafseqs[iseq]))  # NOTE the order of the leaves/names is checked when reading bppseqgen output
        tmpfile.flush()  # BEWARE if you forget this you are fucked
        with open(os.devnull, 'w') as fnull:
            out_treestr = subprocess.check_output('./bin/FastTree -gtr -nt ' + tmpfile.name, shell=True, stderr=fnull)
        out_dtree = get_dendro_tree(treestr=out_treestr, taxon_namespace=taxon_namespace, schema='newick')  # see note above
        out_dtree.reroot_at_node(out_dtree.find_node_with_taxon_label(naive_seq_name), update_bipartitions=True)

    in_dtree = get_dendro_tree(treestr=in_treestr, taxon_namespace=taxon_namespace, schema='newick')  # see note above
    in_height = get_mean_leaf_height(in_dtree)
    out_height = get_mean_leaf_height(out_dtree)
    base_width = 100
    print '  %s trees:' % ('full sequence' if region == 'all' else region)
    print '    %s' % utils.color('blue', 'input:')
    print get_ascii_tree(dendro_tree=in_dtree, extra_str='      ', width=base_width)
    print '    %s' % utils.color('blue', 'output:')
    print get_ascii_tree(dendro_tree=out_dtree, extra_str='        ', width=int(base_width*out_height/in_height))

    if debug:
        print '                   heights: %.3f   %.3f' % (in_height, out_height)
        print '      symmetric difference: %d' % dendropy.calculate.treecompare.symmetric_difference(in_dtree, out_dtree)
        print '        euclidean distance: %f' % dendropy.calculate.treecompare.euclidean_distance(in_dtree, out_dtree)
        print '              r-f distance: %f' % dendropy.calculate.treecompare.robinson_foulds_distance(in_dtree, out_dtree)

# ----------------------------------------------------------------------------------------
def modify_bio_tree_for_lbi(btree, tau, transform, debug=False):
    for node in btree.find_clades(order="postorder"):  # the flu is worrying about which nodes are alive when but a.t.m. we're not
        node.alive = True

    depths = btree.depths()

    # Calculate clock length.
    btree.root.clock_length = 0.0
    for node in btree.find_clades():
        for child in node.clades:
            child.clock_length = depths[child] - depths[node]

    # traverse the tree in postorder (children first) to calculate msg to parents
    for node in btree.find_clades(order="postorder"):
        node.down_polarizer = 0
        node.up_polarizer = 0
        for child in node.clades:
            node.up_polarizer += child.up_polarizer
        bl =  node.clock_length / tau
        node.up_polarizer *= numpy.exp(-bl)
        if node.alive: node.up_polarizer += tau*(1-numpy.exp(-bl))

    # traverse the tree in preorder (parents first) to calculate msg to children
    for node in btree.get_nonterminals():
        for child1 in node.clades:
            child1.down_polarizer = node.down_polarizer
            for child2 in node.clades:
                if child1!=child2:
                    child1.down_polarizer += child2.up_polarizer

            bl =  child1.clock_length / tau
            child1.down_polarizer *= numpy.exp(-bl)
            if child1.alive: child1.down_polarizer += tau*(1-numpy.exp(-bl))

    # go over all nodes and calculate the LBI (can be done in any order)
    max_LBI = 0.0
    for node in btree.find_clades(order="postorder"):
        tmp_LBI = node.down_polarizer
        for child in node.clades:
            tmp_LBI += child.up_polarizer

        node.lbi = transform(tmp_LBI)
        if node.lbi > max_LBI:
            max_LBI = node.lbi

    # Normalize LBI to range [0, 1].
    for node in btree.find_clades():
        node.lbi /= max_LBI

    if debug:
        print '  bio lbi values:'
        for node in btree.find_clades():
            print '    %20s  %8.3f' % (node.name, node.lbi)

# ----------------------------------------------------------------------------------------
def modify_dendro_tree_for_lbi(dtree, tau, transform, debug=False):
    for node in dtree.postorder_node_iter():  # the flu is worrying about which nodes are alive when but a.t.m. we're not
        node.alive = True

    # Calculate clock length.
    for node in dtree.postorder_node_iter():  # postorder shouldn't matter, but I have to choose one or the other when I'm copying from the bio version
        if node.parent_node is None:  # root node
            node.clock_length = 0.
        for child in node.child_node_iter():
            child.clock_length = child.distance_from_root() - node.distance_from_root()

    # traverse the tree in postorder (children first) to calculate msg to parents
    for node in dtree.postorder_node_iter():
        node.down_polarizer = 0
        node.up_polarizer = 0
        for child in node.child_node_iter():
            node.up_polarizer += child.up_polarizer
        bl =  node.clock_length / tau
        node.up_polarizer *= numpy.exp(-bl)
        if node.alive: node.up_polarizer += tau*(1-numpy.exp(-bl))

    # traverse the tree in preorder (parents first) to calculate msg to children
    for node in dtree.preorder_internal_node_iter():
        for child1 in node.child_node_iter():
            child1.down_polarizer = node.down_polarizer
            for child2 in node.child_node_iter():
                if child1!=child2:
                    child1.down_polarizer += child2.up_polarizer

            bl =  child1.clock_length / tau
            child1.down_polarizer *= numpy.exp(-bl)
            if child1.alive: child1.down_polarizer += tau*(1-numpy.exp(-bl))

    # go over all nodes and calculate the LBI (can be done in any order)
    max_LBI = 0.0
    for node in dtree.postorder_node_iter():
        tmp_LBI = node.down_polarizer
        for child in node.child_node_iter():
            tmp_LBI += child.up_polarizer

        node.lbi = transform(tmp_LBI)
        if node.lbi > max_LBI:
            max_LBI = node.lbi

    # Normalize LBI to range [0, 1].
    for node in dtree.postorder_node_iter():  # postorder shouldn't matter, but I have to choose one or the other when I'm copying from the bio version
        node.lbi /= max_LBI

    if debug:
        print '       dendro lbi values:'
        for node in dtree.postorder_node_iter():  # postorder shouldn't matter, but I have to choose one or the other when I'm copying from the bio version
            print '    %20s  %8.3f' % (node.taxon.label, node.lbi)

# ----------------------------------------------------------------------------------------
# copied from https://github.com/nextstrain/augur/blob/master/base/scores.py
def calculate_lbi(treestr=None, treefname=None, naive_seq_name=None, tau=0.4, transform=lambda x:x, extra_str=None, debug=False):  # exactly one of <treestr> or <treefname> should be None
    """
    traverses the tree in postorder and preorder to calculate the up and downstream tree length exponentially weighted
    by distance, then adds them as LBI.
    tree     -- biopython tree for whose node the LBI is being computed
    """

    # reroot at naive sequence, and convert to bio tree
    dtree = get_dendro_tree(treestr=treestr, treefname=treefname)
    if naive_seq_name is not None:  # TODO not sure if I should do this or not
        dtree.reroot_at_node(dtree.find_node_with_taxon_label(naive_seq_name), update_bipartitions=True)

    if debug:
        print '  %s%s' % (utils.color('green', 'lbi'), '' if extra_str is None else ' for %s' % extra_str)
        print '      starting with rerooted tree:'
        print utils.pad_lines(get_ascii_tree(dendro_tree=dtree, width=250))

    # # dendropy makes up new 'id's for each otu when it makes nexml, and puts the existing node/taxon labels as 'label's in the <otu> in the nexml file. And the Bio ignores everything but the 'id', so in order to figure out which stupid node each number corresponds to I'd have to parse the nexml myself to get the translation
    # # so... rewrote the lbi fcn for the dendro tree. Seems to get the exact same values.
    # print dtree.as_string(schema='nexml')
    # btree = get_bio_tree(treestr=dtree.as_string(schema='nexml'))
    # modify_bio_tree_for_lbi(btree, tau, transform, debug=debug)
    modify_dendro_tree_for_lbi(dtree, tau, transform, debug=debug)

    return {'tree' : dtree.as_string(schema='newick'), 'values' : {n.taxon.label : float(n.lbi) for n in dtree.postorder_node_iter()}}

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
    dtree = dendropy.Tree(taxon_namespace=tns, seed_node=root_node)
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
def parse_lonr(outdir, input_seqfos, naive_seq_name, debug=False):
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
        tree = get_dendro_tree(treefname=phylip_treefile, schema='newick')
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
def calculate_lonr(input_seqfos=None, line=None, phylip_treefile=None, phylip_seqfile=None, tree_method=None, naive_seq_name='X-naive-X', seed=1, debug=False):
    assert input_seqfos is None or line is None
    if input_seqfos is None:
        input_seqfos = [{'name' : line['unique_ids'][iseq], 'seq' : line['seqs'][iseq]} for iseq in range(len(line['unique_ids']))]
        input_seqfos.insert(0, {'name' : naive_seq_name, 'seq' : line['naive_seq']})
    if tree_method is None:
        tree_method = 'dnapars' if len(input_seqfos) < 1000 else 'neighbor'

    workdir = utils.choose_random_subdir('/tmp/%s' % os.getenv('USER'))
    os.makedirs(workdir)

    if debug:
        print '  %s' % utils.color('green', 'lonr:')
    run_lonr(input_seqfos, naive_seq_name, workdir, tree_method, phylip_treefile=phylip_treefile, phylip_seqfile=phylip_seqfile, seed=seed, debug=debug)
    lonr_info = parse_lonr(workdir, input_seqfos, naive_seq_name, debug=debug)

    for fn in lonr_files.values():
        os.remove(workdir + '/' + fn)
    os.rmdir(workdir)

    return lonr_info

# ----------------------------------------------------------------------------------------
# interface for calculating tree metrics starting from standard <line> annotations (as opposed to bin/calculate_tree_metrics.py, which is more standalone, e.g. from cft)
def calculate_tree_metrics(annotations, min_tree_metric_cluster_size, reco_info=None, use_true_clusters=False, debug=False):
    if reco_info is not None:
        for tmpline in reco_info.values():
            assert len(tmpline['unique_ids']) == 1  # at least for the moment, we're splitting apart true multi-seq lines when reading in seqfileopener.py

    lines_to_use, true_lines_to_use = None, None
    if use_true_clusters:
        assert reco_info is not None
        true_partition = utils.get_true_partition(reco_info)
        print '    using %d true clusters to calculate inferred tree metrics (sizes: %s)' % (len(true_partition), ' '.join([str(len(c)) for c in true_partition]))
        lines_to_use, true_lines_to_use = [], []
        for cluster in true_partition:
            true_lines_to_use.append(utils.synthesize_multi_seq_line_from_reco_info(cluster, reco_info))  # note: duplicates (a tiny bit of) code in utils.print_true_events()
            max_in_common, ustr_to_use = None, None
            for ustr in annotations:  # order will be different in reco info and inferred clusters
                n_in_common = len(set(ustr.split(':')) & set(cluster))  # can't just look for the actual cluster since we collapse duplicates, but bcr-phylo doesn't (but maybe I should throw them out when parsing bcr-phylo output)
                if max_in_common is None or n_in_common > max_in_common:
                    ustr_to_use = ustr
                    max_in_common = n_in_common
            if max_in_common is None:
                raise Exception('cluster \'%s\' not found in inferred annotations (probably because use_true_clusters was set)' % ':'.join(cluster))
            if max_in_common < len(cluster):
                print '    %s couldn\'t find an inferred cluster that shared all sequences with true cluster (best was %d/%d)' % (utils.color('red', 'note:'), max_in_common, len(cluster))
            lines_to_use.append(annotations[ustr_to_use])
    else:
        lines_to_use = annotations.values()
        if reco_info is not None:
            for line in lines_to_use:
                true_line = utils.synthesize_multi_seq_line_from_reco_info(line['unique_ids'], reco_info)
                true_lines_to_use.append(true_line)


    n_clusters_calculated, n_skipped = 0, 0
    for line in lines_to_use:
        if len(line['unique_ids']) < min_tree_metric_cluster_size:
            n_skipped += 1
            continue
        lonr_info = calculate_lonr(line=line, debug=debug)
        lbi_info = calculate_lbi(treestr=lonr_info['tree'], extra_str='inf tree', debug=debug)
        line['tree-info'] = {'lonr' : lonr_info, 'lbi' : lbi_info}
        n_clusters_calculated += 1

    print '  calculated tree metrics for %d cluster%s (skipped %d smaller than %d)' % (n_clusters_calculated, utils.plural(n_clusters_calculated), n_skipped, min_tree_metric_cluster_size)

    if reco_info is not None:
        import plotting
        plotting.compare_tree_metrics(lines_to_use, reco_info)
        # NOTE not actually doing anything with the true lbi info yet (I think I don't want to mix it in with the inferred info, maybe just make plots with it)
        # for true_line in true_lines_to_use:
        #     true_lbi_info = calculate_lbi(treestr=true_line['tree'], extra_str='true tree', debug=debug)  # NOTE this is the tree we _give_ to bppseqgen, not necessarily the exact one implied by the true mutations (er, ok, I'm not sure if that distinction makes sense)
