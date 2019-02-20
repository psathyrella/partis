fit_sgd <- function(phylo, tot_time, par, f=1, meth="Nelder-Mead")
{
	nobs <- Ntip(phylo)
	p<-length(par)
	
	lambert <- Phylo2Lambert(phylo)
	lambert[1] <- tot_time
	
	init<-par

	optim_LH <- function(init)
	{
		likel <- LikelihoodSGDFromLambert(lambert, abs(init[1]), abs(init[2]), abs(init[3]), f)
		LH <- likel[3]+likel[4]
		return(-log(LH))
	}
	
	temp <- optim(init, optim_LH)
	res<-list(model="sgd", LH=-temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1), par=c(birth = abs(temp$par[1]), growth = abs(temp$par[1])-abs(temp$par[2]), mutation = abs(temp$par[3])))
	return(res)
}
