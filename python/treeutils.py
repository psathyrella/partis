from cStringIO import StringIO
import subprocess
import tempfile
import os
import numpy
import sys

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
def get_bio_tree(treestr=None, treefname=None):
    if 'Bio.Phylo' not in sys.modules:  # NOTE dendropy seems a lot nicer... use that for new stuff
        from Bio import Phylo
    Phylo = sys.modules['Bio.Phylo']
    if treestr is not None:
        return Phylo.read(StringIO(treestr), 'newick')
    elif treefname is not None:
        with open(treefname) as treefile:
            return Phylo.read(treefile, 'newick')
    else:
        assert False

# ----------------------------------------------------------------------------------------
def get_baltic_tree(treestr):
    if treestr.count(':') == 1:  # one-leaf tree
        name, lengthstr = treestr.strip().rstrip(';').split(':')
        tree = OneLeafTree(name, float(lengthstr))
    else:
        tree = baltic.tree()
        baltic.make_tree(treestr, tree, verbose=False)
    tree.traverse_tree()
    return tree

# ----------------------------------------------------------------------------------------
def get_n_leaves(treestr):
    return len(get_baltic_tree(treestr).leaves)

# ----------------------------------------------------------------------------------------
def get_mean_height(treestr):
    tree = get_baltic_tree(treestr)
    heights = [l.height for l in tree.leaves]
    return sum(heights) / len(heights)

# ----------------------------------------------------------------------------------------
def get_ascii_tree(treestr, extra_str='', width=100):
    if get_mean_height(treestr) == 0.:  # we really want the max height, but since we only care whether it's zero or not this is the same
        return '%szero height' % extra_str
    elif get_n_leaves(treestr) > 1:  # if more than one leaf
        if 'Bio.Phylo' not in sys.modules:  # NOTE dendropy seems a lot nicer... use that for new stuff
            from Bio import Phylo
        Phylo = sys.modules['Bio.Phylo']
        tmpf = StringIO()
        Phylo.draw_ascii(get_bio_tree(treestr=treestr), file=tmpf, column_width=width)
        return '\n'.join(['%s%s' % (extra_str, line) for line in tmpf.getvalue().split('\n')])
    else:
        return '%sone leaf' % extra_str

# ----------------------------------------------------------------------------------------
def rescale_tree(treestr, new_height, debug=False):
    """ rescale the branch lengths in <treestr> (newick-formatted) by <factor> """
    tree = get_baltic_tree(treestr)
    mean_height = get_mean_height(treestr)
    for ln in tree.Objects:
        old_length = ln.length
        ln.length *= new_height / mean_height  # rescale every branch length in the tree by the ratio of desired to existing height (everybody's heights should be the same... but they never quite were when I was using Bio.Phylo, so, uh. yeah, uh. not sure what to do, but this is fine. It's checked below, anyway)
        if debug:
            print '  %5s  %7e  -->  %7e' % (ln.numName if ln.branchType == 'leaf' else ln.branchType, old_length, ln.length)
    tree.traverse_tree()
    treestr = tree.toString(numName=True)
    for leaf in get_baltic_tree(treestr).leaves:  # make sure string conversion (and rescaling) went ok
        if not utils.is_normed(leaf.height / new_height, this_eps=1e-8):
            raise Exception('tree not rescaled properly:   %.10f   %.10f    %e' % (leaf.height, new_height, (leaf.height - new_height) / new_height))
    return treestr

# ----------------------------------------------------------------------------------------
def infer_tree_from_leaves(region, in_tree, leafseqs, naive_seq, naive_seq_name='XnaiveX', debug=False):  # baltic barfs on (some) dashes
    if 'dendropy' not in sys.modules:
        import dendropy
    dendropy = sys.modules['dendropy']

    taxon_namespace = dendropy.TaxonNamespace()
    with tempfile.NamedTemporaryFile() as tmpfile:
        tmpfile.write('>%s\n%s\n' % (naive_seq_name, naive_seq))
        for iseq in range(len(leafseqs)):
            tmpfile.write('>t%s\n%s\n' % (iseq+1, leafseqs[iseq]))  # NOTE the order of the leaves/names is checked when reading bppseqgen output
        tmpfile.flush()  # BEWARE if you forget this you are fucked
        with open(os.devnull, 'w') as fnull:
            out_tree = subprocess.check_output('./bin/FastTree -gtr -nt ' + tmpfile.name, shell=True, stderr=fnull)
        out_dtree = dendropy.Tree.get_from_string(out_tree, 'newick', taxon_namespace=taxon_namespace)
        out_dtree.reroot_at_node(out_dtree.find_node_with_taxon_label(naive_seq_name), update_bipartitions=True)
        out_tree = out_dtree.as_string(schema='newick', suppress_rooting=True)

    in_height = get_mean_height(in_tree)
    out_height = get_mean_height(out_tree)
    base_width = 100
    print '  %s trees:' % ('full sequence' if region == 'all' else region)
    print '    %s' % utils.color('blue', 'input:')
    print get_ascii_tree(in_tree, extra_str='      ', width=base_width)
    print '    %s' % utils.color('blue', 'output:')
    print get_ascii_tree(out_tree, extra_str='        ', width=int(base_width*out_height/in_height))

    in_dtree = dendropy.Tree.get_from_string(in_tree, 'newick', taxon_namespace=taxon_namespace)

    if debug:
        print '                   heights: %.3f   %.3f' % (in_height, out_height)
        print '      symmetric difference: %d' % dendropy.calculate.treecompare.symmetric_difference(in_dtree, out_dtree)
        print '        euclidean distance: %f' % dendropy.calculate.treecompare.euclidean_distance(in_dtree, out_dtree)
        print '              r-f distance: %f' % dendropy.calculate.treecompare.robinson_foulds_distance(in_dtree, out_dtree)

# ----------------------------------------------------------------------------------------
# copied from https://github.com/nextstrain/augur/blob/master/base/scores.py 
def calculate_LBI(tree, attr="lbi", tau=0.4, transform=lambda x:x, **kwargs):
    """
    traverses the tree in postorder and preorder to calculate the
    up and downstream tree length exponentially weighted by distance.
    then adds them as LBI
    tree     -- biopython tree for whose node the LBI is being computed
    attr     -- the attribute name used to store the result
    """
    # Calculate clock length.
    tree.root.clock_length = 0.0
    for node in tree.find_clades():
        for child in node.clades:
            child.clock_length = child.attr['num_date'] - node.attr['num_date']

    # traverse the tree in postorder (children first) to calculate msg to parents
    for node in tree.find_clades(order="postorder"):
        node.down_polarizer = 0
        node.up_polarizer = 0
        for child in node.clades:
            node.up_polarizer += child.up_polarizer
        bl =  node.clock_length / tau
        node.up_polarizer *= numpy.exp(-bl)
        if node.alive: node.up_polarizer += tau*(1-numpy.exp(-bl))

    # traverse the tree in preorder (parents first) to calculate msg to children
    for node in tree.get_nonterminals():
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
    for node in tree.find_clades(order="postorder"):
        tmp_LBI = node.down_polarizer
        for child in node.clades:
            tmp_LBI += child.up_polarizer

        node.attr[attr] = transform(tmp_LBI)
        if node.attr[attr] > max_LBI:
            max_LBI = node.attr[attr]

    # Normalize LBI to range [0, 1].
    for node in tree.find_clades():
        node.attr[attr] /= max_LBI
        setattr(node, attr, node.attr[attr])
