fit_coal_cst <- function (phylo, tau0=1.e-2, gamma=1, cst.rate=FALSE, meth = "Nelder-Mead", N0=0)
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

    Vtimes <- sort(branching.times(phylo))
    ntips<-Ntip(phylo)
    if (N0==0)
    {
      N0<-ntips
    }

 if (cst.rate==FALSE)
 {
    init<-c(tau0,gamma)
    nbpar<-length(init)
    nbobs<-length(Vtimes)-1

    optimLH <- function(init)
    {
      tau0 <- init[1]
      gamma <- init[2]
      LH <- likelihood_coal_cst(Vtimes,ntips,tau0,gamma,N0)$res
      return(-LH)
    }

    temp<-suppressWarnings(optim(init, optimLH, method = meth))

    res <- list(model = "Equilibrium variable rate", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), tau0 = temp$par[1], gamma=temp$par[2])
    return(res)
 }
 else
 {
    init<-c(tau0)

    nbpar<-length(init)
    nbobs<-length(Vtimes)-1

    optimLH <- function(init)
    {
      tau0 <- init[1]
      LH <- likelihood_coal_cst(Vtimes,ntips,tau0,0.0,N0)$res
      return(-LH)
    }

    temp<-suppressWarnings(optim(init, optimLH, method = meth))

    res <- list(model = "Equilibrium constant rate", LH = -temp$value, aicc = 2 * temp$value + 2*nbpar + 2*nbpar*(nbpar+1)/(nbobs-nbpar-1), tau0 = temp$par[1])

    return(res)
 }
}
