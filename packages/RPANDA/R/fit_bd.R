fit_bd <-
 function (phylo, tot_time, f.lamb, f.mu, lamb_par, mu_par, f=1,
           meth = "Nelder-Mead", cst.lamb=FALSE, cst.mu=FALSE,
           expo.lamb=FALSE, expo.mu=FALSE, fix.mu=FALSE,
           dt=0, cond="crown")
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

  nobs <- Ntip(phylo)

  if (fix.mu==FALSE)
  {
    init <- c(lamb_par,mu_par)
    p <- length(init)
    optimLH <- function(init)
    {
      lamb_par <- init[1:length(lamb_par)]
      mu_par <- init[(1+length(lamb_par)):length(init)]
      f.lamb.par <- function(t){abs(f.lamb(t,lamb_par))}
      f.mu.par <- function(t){abs(f.mu(t,mu_par))}
      LH <- likelihood_bd(phylo,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt,cond=cond)
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    mu.par <- temp$par[(1+length(lamb_par)):length(init)]
    f.lamb.par <- function(t){abs(f.lamb(t, lamb.par))}
    f.mu.par <- function(t){abs(f.mu(t, mu.par))}
    res <- list(model = "birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1) , lamb_par=lamb.par, mu_par=mu.par, f.lamb=Vectorize(f.lamb.par), f.mu=Vectorize(f.mu.par))
  }

  else
  {
    init <- c(lamb_par)
    p <- length(init)
    optimLH <- function(init)
    {
      lamb_par <- init[1:length(lamb_par)]
      f.lamb.par <- function(t){abs(f.lamb(t,lamb_par))}
      f.mu.par <- function(t){abs(f.mu(t,mu_par))}
      LH <- likelihood_bd(phylo,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=TRUE,expo.lamb=expo.lamb,dt=dt,cond=cond)
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    f.lamb.par <- function(t){abs(f.lamb(t, lamb.par))}
    f.mu.par <- function(t){abs(f.mu(t, mu_par))}
    res <- list(model = "birth.death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=lamb.par, f.lamb=Vectorize(f.lamb.par))
  }
  class(res) <- "fit.bd"
  return(res)
}
