.likelihood_coal_exp_mod <- function(Vtimes,ntips,tau0,gamma,N0)
{
  Ttimes <- diff(Vtimes)
  Vtimes <- Vtimes[2:length(Vtimes)]
  nbint <- length(Ttimes)
  samp <- seq((ntips-2),(ntips-nbint-1),by=-1)
  indLikelihood <- samp*(samp+1)/2*2*tau0/N0*exp(gamma*Vtimes)*exp(-samp*(samp+1)/2*2*tau0/N0*1/gamma*exp(gamma*Vtimes)*(1-exp(-gamma*Ttimes)))
  res <- sum(log(indLikelihood))
  return(list("res"=res,"all"=indLikelihood))
}
