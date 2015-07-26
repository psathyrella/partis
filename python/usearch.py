#!/usr/bin/env python
import sys
import itertools
from Bio import Phylo
from cStringIO import StringIO

tree = Phylo.read('tree.phy', 'newick')

debug = True
n_procs = 3
namehashes = {}
leaf_nodes = list(tree.find_clades(terminal=True))
distances = {}

def print_tree(tree):
    treestr = str(tree)  #tree.format('newick')
    for namehash, realnames in namehashes.items():
        treestr = treestr.replace(namehash, realnames)
    # Phylo.draw_ascii(Phylo.read(StringIO(treestr), 'newick'))
    print treestr

def print_partition():
    for tn in tree.get_terminals():
        print '  ', namehashes.get(tn.name, tn.name)

def adjacent_pairs(tmplf):
    pairs = itertools.tee(tmplf)
    pairs[1].next()
    return [adjpair for adjpair in itertools.izip(pairs[0], pairs[1]) if tree.is_monophyletic(adjpair)]

def closest_leaves():
    min_distance = sorted(distances)[0]
    min_pair = distances[min_distance]
    return min_pair

def join_names(names):
    real_names = []
    for name in names:
        real_names.append(namehashes.get(name, name))
    return ' '.join([n for n in real_names])

def merge_leaves():
    closepair = closest_leaves()
    mrca = tree.common_ancestor(closepair)
    joint_name = join_names(closepair)
    joint_name_hash = str(hash(joint_name))
    namehashes[joint_name_hash] = joint_name
    mrca.name = joint_name_hash
    mrca.branch_length += tree.find_any(closepair[0]).branch_length  # NOTE just using the branch length to the first one... although they're not always exactly the same, they seem to be close
    for leaf in closepair:
        for dist, pair in distances.items():
            if leaf in pair:
                del distances[dist]
        for ln in leaf_nodes:
            if ln.name == leaf:
                leaf_nodes.remove(ln)
        tree.collapse(leaf)
    for ln in leaf_nodes:
        newpair = (ln, mrca)
        if tree.is_monophyletic(newpair):
            newdist = tree.distance(*newpair)
            distances[newdist] = [l.name for l in newpair]
            print '      adding %d   %s   and   %s' % (newdist, ln, namehashes[mrca.name])
    leaf_nodes.append(mrca)

    if debug:
        print '\n\nmerged   %s   and   %s' % tuple([namehashes.get(n, n) for n in closepair])
        print_tree(tree)
        print_partition()

def go():
    if debug:
        print_tree(tree)
    print 'initializing distances'
    for pair in adjacent_pairs(leaf_nodes):
        dist = tree.distance(*pair)
        distances[dist] = [l.name for l in pair]
        print '        %.5f   %s  %s' % (dist, pair[0], pair[1])
    while tree.count_terminals() > n_procs:
        merge_leaves()
    
