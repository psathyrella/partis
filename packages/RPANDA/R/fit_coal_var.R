fit_coal_var <-
  function (phylo, lamb0=0.1, alpha=1, mu0=0.01, beta=0,
            meth = "Nelder-Mead", N0=0,
            cst.lamb=FALSE, cst.mu=FALSE,
            fix.eps=FALSE, mu.0=FALSE, pos=TRUE)
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

    Vtimes <- sort(branching.times(phylo))
    ntips<-Ntip(phylo)
    if (N0==0)
    {
      N0<-ntips
    }


  nbobs<-length(Vtimes)-1


  #pure birth constant rates (Model 5)
  if (mu.0==TRUE & cst.lamb==TRUE)
  {
    init<-c(lamb0)
  }

  #birth-death constant rates (Model 3)
  else if (cst.mu==TRUE & cst.lamb==TRUE)
  {
    init <- c(lamb0,mu0)
  }

  #pure birth varying speciation rate (Model 6)
  else if (mu.0==TRUE & cst.lamb==FALSE)
  {
    init<-c(lamb0,alpha)
  }

  #birth-death varying speciation rate (Model 4a)
  else if (cst.mu==TRUE & cst.lamb==FALSE)
  {
    init <- c(lamb0,alpha,mu0)
  }

  #birth-death varying extinction rate (Model 4b)
  else if (cst.mu==FALSE & cst.lamb==TRUE)
  {
    init <- c(lamb0,mu0,beta)
  }

  #birth-death varying speciation rate and constant extinction fraction (Model 4c)
  else if (fix.eps==TRUE)
  {
    init <- c(lamb0,alpha,mu0/lamb0)
  }

  #birth-death varying speciation and extinction rates (Model 4d)
  else
  {
    init = c(lamb0,alpha,mu0,beta)
  }

  nbpar <- length(init)

############################################################

  if (mu.0==TRUE & cst.lamb==TRUE)
  {
    optimLH.coalBD <- function(init)
    {
      lamb0 <- init[1]
      LH <- likelihood_coal_var(Vtimes,ntips,lamb0,alpha=0,mu0=0,beta=0,N0,pos=pos)$res
      return(-LH)
    }
  }

  else if (cst.mu==TRUE & cst.lamb==TRUE)
  {
    optimLH.coalBD <- function(init)
    {
      lamb0 <- init[1]
      mu0 <- init[2]
      LH <- likelihood_coal_var(Vtimes,ntips,lamb0,alpha=0,mu0,beta=0,N0,pos=pos)$res
      return(-LH)
    }
  }

  else if (mu.0==TRUE & cst.lamb==FALSE)
  {
    optimLH.coalBD <- function(init)
    {
      lamb0 <- init[1]
      alpha <- init[2]
      LH <- likelihood_coal_var(Vtimes,ntips,lamb0,alpha,mu0=0,beta=0,N0,pos=pos)$res
      return(-LH)
    }
  }

  else if (cst.mu==TRUE & cst.lamb==FALSE)
  {
    optimLH.coalBD <- function(init)
    {
      lamb0 <- init[1]
      alpha <- init[2]
      mu0 <- init[3]
      LH <- likelihood_coal_var(Vtimes,ntips,lamb0,alpha,mu0,beta=0,N0,pos=pos)$res
      return(-LH)
    }
  }

  else if (cst.mu==FALSE & cst.lamb==TRUE)
  {
    optimLH.coalBD <- function(init)
    {
      lamb0 <- init[1]
      mu0 <- init[2]
      beta<- init[3]
      LH <- likelihood_coal_var(Vtimes,ntips,lamb0,alpha=0,mu0,beta,N0,pos=pos)$res
      return(-LH)
    }
  }

  else if (fix.eps==TRUE)
  {
    optimLH.coalBD <- function(init)
    {
      lamb0 <- init[1]
      alpha <- init[2]
      eps <- init[3]
      LH <- likelihood_coal_var(Vtimes,ntips,lamb0,alpha,mu0=lamb0*eps,beta=alpha,N0,pos=pos)$res
      return(-LH)
    }
  }

  else
  {
    optimLH.coalBD <- function(init)
    {
      lamb0 <- init[1]
      alpha <- init[2]
      mu0 <- init[3]
      beta <- init[4]
      LH <- likelihood_coal_var(Vtimes,ntips,lamb0,alpha,mu0,beta,N0,pos=pos)$res
      return(-LH)
    }
  }

  #######################################################################################

  temp <- suppressWarnings(optim(init, optimLH.coalBD, method = meth,control=list(ndeps=10^(-4))))
  if (mu.0==TRUE & cst.lamb==TRUE)
  {
    if (pos==FALSE)
    {
      res <- list(model = "Pure birth constant speciation (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1])
    }
    else
    {
      res <- list(model = "Pure birth constant speciation (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]))
    }
  }

  else if (cst.mu==TRUE & cst.lamb==TRUE)
  {
    if (pos==FALSE)
    {
      res <- list(model = "Birth-death constant rates (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1],mu0 = temp$par[2])
    }
    else
    {
      res <- list(model = "Birth-death constant rates (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]),mu0 = abs(temp$par[2]))
    }
  }

  else if (mu.0==TRUE & cst.lamb==FALSE)
  {
    if (pos==FALSE)
    {
      res <- list(model = "Pure birth varying speciation (coalescent approx)", LH = -temp$value, aicc = 2 *	temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], alpha = temp$par[2])
    }
    else
    {
      res <- list(model = "Pure birth varying speciation (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), alpha = temp$par[2])
    }
  }

  else if (cst.mu==TRUE & cst.lamb==FALSE)
  {
    if (pos==FALSE)
    {
      res <- list(model = "Birth-death varying speciation constant extinction (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], alpha = temp$par[2], mu0 = temp$par[3])
    }
    else
    {
      res <- list(model = "Birth-death varying speciation constant extinction (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), alpha = temp$par[2], mu0 = abs(temp$par[3]))
    }
  }

  else if (cst.mu==FALSE & cst.lamb==TRUE)
  {
    if (pos==FALSE)
    {
      res <- list(model = "Birth-death constant speciation varying extinction (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], mu0 = temp$par[2], beta = temp$par[3])
    }
    else
    {
      res <- list(model = "Birth-death constant speciation varying extinction (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), mu0 = abs(temp$par[2]), beta = temp$par[3])
    }
  }

  else if (fix.eps==TRUE)
  {
    if (pos==FALSE)
    {
      res <- list(model = "Birth-death constant extinction fraction (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], alpha = temp$par[2], eps = temp$par[3])
    }
    else
    {
      res <- list(model = "Birth-death constant extinction fraction (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), alpha = temp$par[2], eps = abs(temp$par[3]))
    }
  }

  else
  {
    if (pos==FALSE)
    {
      res <- list(model = "Birth-death varying speciation and extinction (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = temp$par[1], alpha = temp$par[2], mu0 = temp$par[3],beta = temp$par[4])
    }
    else
    {
      res <- list(model = "Birth-death varying speciation and extinction (coalescent approx)", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), lamb0 = abs(temp$par[1]), alpha = temp$par[2], mu0 = abs(temp$par[3]),beta = temp$par[4])
    }
  }
  return(res)
}
