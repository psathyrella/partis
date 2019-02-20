.likelihood_coal_cst_mod <-function(Vtimes,ntips,tau0,N0)
{
  Ttimes <- diff(Vtimes)
  nbint<-length(Ttimes)
  samp<-seq((ntips-2),(ntips-nbint-1),by=-1)
  indLikelihood<-samp*(samp+1)/2*2*tau0/N0*exp(-samp*(samp+1)/2*2*tau0/N0*Ttimes)
  res<-sum(log(indLikelihood))
  return(list("res"=res,"all"=indLikelihood))
}
