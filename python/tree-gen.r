#
# Generate tree file for later use by Recombinator.py
#
# Tree-generating step is slow, so write a bunch of trees and just choose from them in Recombinator.
# Run with:
# > R --slave -f tree-gen.r

#\library(ape)
#\set.seed(0)
#\for (itree in 1:1000) {
#\    write.tree(rtree(10, br = rexp, rate = 10), 'data/trees.tre', append=TRUE)
#\}

library(TreeSim)
n_species <- 10
n_trees <- 1
speciation_rate <- 0.1
extinction_rate <- 0.02
age <- 1
trees <- sim.bd.taxa.age(n_species, n_trees, speciation_rate, extinction_rate, frac = 1, age, mrca = FALSE)
#trees <- sim.bd.age(age=age, n_trees=, speciation_rate, K=0, extinction_rate, frac = 1, mrca = FALSE)

#write.tree(trees[[i]], 'tmp.tre', append=FALSE)
#for (i in 1:(n_trees-1)) {
#    print(i)
#    write.tree(trees[[i]], 'tmp.tre', append=TRUE)
#}
#

png(file="foo.png")
plot.phylo(trees[[1]])
add.scale.bar()
dev.off()
