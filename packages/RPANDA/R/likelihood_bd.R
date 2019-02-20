likelihood_bd <- function(phylo,tot_time,f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE,dt=0,cond="crown")
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

  nbtips <- Ntip(phylo)
  log_indLikelihood <- c()
  from_past <- cbind(phylo$edge,node.age(phylo)$ages)
  ages <- rbind(from_past[,2:3],c(nbtips+1,0))
  ages <- ages[order(ages[,1]),]
  age <- max(ages[,2])

  for (j in 1:(nbtips-1))
  {
    node <- (nbtips+j)
    edges <- phylo$edge[phylo$edge[,1]==node,]
    tj <- age-ages[edges[1,1],2]
    Psi_timevar_errap_tj <- .Psi(0,tj,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt)
    log_lik_tj <- log(f.lamb(tj)) + log(Psi_timevar_errap_tj)
    log_indLikelihood <- c(log_indLikelihood,log_lik_tj)
  }
  log_indLikelihood <- c(log_indLikelihood,log(.Psi(0,tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt)))

  log_data_lik <- sum(log_indLikelihood)+nbtips*log(f)

  if (cond==FALSE)
  {
    log_final_lik <- log_data_lik
  }

  else if (cond=="stem")
  {
    Phi <- .Phi(tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt)
    log_final_lik <- log_data_lik - log(1-Phi)
  }

  else if (cond=="crown")
  {
    Phi <- .Phi(tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt)
    log_final_lik <- log_data_lik - log(f.lamb(tot_time)) - 2*log(1-Phi)
  }
  return(log_final_lik)
}
