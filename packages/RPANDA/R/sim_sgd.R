sim_sgd <- function(tau, b, d, nu)
{
	lambert <- Simulate_SGD_Phylogeny(tau, 0, 2, b, d, nu)
	if(length(lambert) == 1)
	{
		newick <- paste("(", Lambert2Newick(lambert), ");", sep="")
	}
	else{
		newick <- paste(Lambert2Newick(lambert), ";", sep="")
	}
	phylo <- read.tree(text = newick)
	phylo$tip.label <- paste("s", 1:Ntip(phylo), sep = "")
	return(phylo)
}
