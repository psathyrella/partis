likelihood_sgd <- function(phylo, tot_time, b, d, nu, f)
{
	lambert <- Phylo2Lambert(phylo)
	lambert[1] <- tot_time
	likel <- LikelihoodSGDFromLambert(lambert, b, d, nu, f)
	return(likel[3]+likel[4])
}
