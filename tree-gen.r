#
# Generate tree file for later use by Recombinator.py
#
# Tree-generating step is slow, so a bunch of trees and just choose from them in Recombinator.
# Run with:
# > R --slave -f tree-gen.r

library(ape)
set.seed(0)
for (itree in 1:1000) {
    write.tree(rtree(10, br = rexp, rate = 10), 'data/trees.tre', append=TRUE)
}
