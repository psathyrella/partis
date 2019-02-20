likelihood_coal_cst <- function(Vtimes,ntips,tau0,gamma,N0)
{
  if (gamma == 0)
  {
    return (.likelihood_coal_cst_mod(Vtimes,ntips,tau0,N0))
  }
  else
  {
    return(.likelihood_coal_exp_mod(Vtimes,ntips,tau0,gamma,N0))
  }
}
