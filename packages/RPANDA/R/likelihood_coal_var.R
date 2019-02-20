likelihood_coal_var  <- function(Vtimes,ntips,lamb0,alpha,mu0,beta,N0,pos=TRUE)
{

  Ttimes <- diff(Vtimes)
  Vtimes <- Vtimes[2:length(Vtimes)]
  nbint<-length(Ttimes)
  samp<-seq((ntips-2),(ntips-nbint-1),by=-1)

  times<-c(0,sort(Vtimes))

  if (min(abs(lamb0)*exp(alpha*times)-abs(mu0)*exp(beta*times))<=0)
  {
    indLikelihood<-0*vector(length=length(samp))
    res<-sum(log(indLikelihood))
  }

  # Model 3: \lambda(t) = \lambda_0, \mu(t) = \mu_0 & N(t) = N_0 exp(-(\lambda_0 - \mu_0) * t)
  else if ((alpha==0) & (beta==0))
  {
    if (pos==FALSE)
    {
       r <- lamb0 - mu0
       indLikelihood <- samp*(samp+1)/2*2*lamb0/N0*1/(exp(-r*Vtimes))*exp(-samp*(samp+1)/2*2*lamb0/N0/r*exp(r*Vtimes)*(1-exp(-r*Ttimes)))
    }
    else
    {
       r <- abs(abs(lamb0)-abs(mu0))
       indLikelihood <- samp*(samp+1)/2*2*abs(lamb0)/N0*1/(exp(-r*Vtimes))*exp(-samp*(samp+1)/2*2*abs(lamb0)/N0/r*exp(r*Vtimes)*(1-exp(-r*Ttimes)))
    }
    res <- sum(log(indLikelihood))
  }

  else
  {
    if ((beta==0) & !(alpha==0))
    {
      if (pos==FALSE)
      {
         demfun<-function(t){2*lamb0*exp(alpha*t)/(N0*exp(lamb0/alpha*(1-exp(alpha*t))+mu0*t))}
      }
      else
      {
         demfun<-function(t){2*abs(lamb0)*exp(alpha*t)/(N0*exp(abs(lamb0)/alpha*(1-exp(alpha*t))+abs(mu0)*t))}
      }
  }

  else if ((alpha==0) & !(beta==0))
  {
     if (pos==FALSE)
     {
        demfun<-function(t){2*lamb0/(N0*exp(-lamb0*t-mu0/beta*(1-exp(beta*t))))}
     }
     else
     {
        demfun<-function(t){2*abs(lamb0)/(N0*exp(-abs(lamb0)*t-abs(mu0)/beta*(1-exp(beta*t))))}
     }
  }

  else
  {
    if (pos==FALSE)
    {
       demfun<-function(t){2*lamb0*exp(alpha*t)/(N0*exp(lamb0/alpha*(1-exp(alpha*t))-mu0/beta*(1-exp(beta*t))))}
    }
    else
    {
       demfun<-function(t){2*abs(lamb0)*exp(alpha*t)/(N0*exp(abs(lamb0)/alpha*(1-exp(alpha*t))-abs(mu0)/beta*(1-exp(beta*t))))}
    }
  }

  if (FALSE %in% is.finite(demfun(Vtimes)))
  {
     indLikelihood<-0*vector(length=length(samp))
     res<-sum(log(indLikelihood))
  }
  else
  {
     integrals<-c()
     demfunval<-c()
     for (i in 1:length(Vtimes))
     {
        demfunvali<-demfun(Vtimes[i])
        integrali<-integrate(demfun,(Vtimes[i]-Ttimes[i]),Vtimes[i],stop.on.error=FALSE)$value

        demfunval<-c(demfunval,demfunvali)
        integrals<-c(integrals,integrali)
     }

     indLikelihood<-samp*(samp+1)/2*demfunval*exp(-samp*(samp+1)/2*integrals)
     res<-sum(log(indLikelihood))
  }
 }
 return(list("res"=res,"all"=indLikelihood))
}
