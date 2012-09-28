#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>

#include "smctc.hh"
#include "simfunctions.hh"

using namespace std;

///The function corresponding to the log posterior density at specified time and position

///   \param lTime The current time (i.e. the index of the current distribution)
///    \param X     The state to consider  **/
double logDensity(long lTime, const mChain<double> & X)
{
  double lp;

  mElement<double> *x = X.GetElement(0);
  mElement<double> *y = x->pNext;
  //Begin with the density exluding the effect of the potential
  lp = log(gsl_ran_ugaussian_pdf(x->value));

  while(y) {
    lp += log(gsl_ran_ugaussian_pdf(y->value - x->value));
    x = y;
    y = x->pNext;
  }

  //Now include the effect of the multiplicative potential function
  lp -= log(1.0 + exp(-(ALPHA(lTime) * (x->value - THRESHOLD) )));
  return lp;
}

///A function to initialise double type markov chain-valued particles
/// \param pRng A pointer to the random number generator which is to be used
smc::particle<mChain<double> > fInitialise(smc::rng *pRng)
{
  // Create a Markov chain with the appropriate initialisation and then assign that to the particles.
  mChain<double> Mc;

  double x = 0;
  for(int i = 0; i < PATHLENGTH; i++) {
    x += pRng->NormalS();
    Mc.AppendElement(x);
  }

  return smc::particle<mChain<double> >(Mc,0);
}

///A function to select a move randomly
/// \param lTime  The current evolution time of the system
/// \param p      The current position of the particle which is to be moved
/// \param pRng   A pointer to the random number generator which is to be used
long fSelect(long lTime, const smc::particle<mChain<double> > & p, smc::rng *pRng)
{
    return 0;
}

void fMove1(long lTime, smc::particle<mChain<double> > & pFrom, smc::rng *pRng)
{
  // The distance between points in the random grid.
  static double delta = 0.025;
  static double gridweight[2*GRIDSIZE+1], gridws = 0;
  static mChain<double> NewPos[2*GRIDSIZE+1];
  static mChain<double> OldPos[2*GRIDSIZE+1];

  // First select a new position from a grid centred on the old position, weighting the possible choises by the
  // posterior probability of the resulting states.
  gridws = 0;
  for(int i = 0; i < 2*GRIDSIZE+1; i++) {
    NewPos[i] = pFrom.GetValue() + ((double)(i - GRIDSIZE))*delta;
    gridweight[i] = exp(logDensity(lTime,NewPos[i]));
    gridws        = gridws + gridweight[i];
  }

  double dRUnif = pRng->Uniform(0,gridws);
  long j = -1;

  while(dRUnif > 0 && j <= 2*GRIDSIZE) {
    j++;
    dRUnif -= gridweight[j];
  }
  
  pFrom.SetValue(NewPos[j]);

  // Now calculate the weight change which the particle suffers as a result
  double logInc = log(gridweight[j]), Inc = 0; 

  for(int i = 0; i < 2*GRIDSIZE+1; i++) {
    OldPos[i] = pFrom.GetValue() - ((double)(i - GRIDSIZE))*delta;
    gridws = 0;
    for(int k = 0; k < 2*GRIDSIZE+1; k++) {
      NewPos[k] = OldPos[i] + ((double)(k-GRIDSIZE))*delta;
      gridweight[k] = exp(logDensity(lTime, NewPos[k]));
      gridws += gridweight[k];
    }
    Inc += exp(logDensity(lTime-1, OldPos[i])) * exp(logDensity(lTime, pFrom.GetValue())) / gridws;
  }
    logInc -= log(Inc);

  pFrom.SetLogWeight(pFrom.GetLogWeight() + logInc);

  for(int i = 0; i < 2*GRIDSIZE+1; i++)
    {
      NewPos[i].Empty();
      OldPos[i].Empty();
    }

  return;
}
///Another move function
void fMove2(long lTime, smc::particle<mChain<double> > & pFrom, smc::rng *pRng)
{
  pFrom.SetLogWeight(pFrom.GetLogWeight() + logDensity(lTime,pFrom.GetValue()) - logDensity(lTime-1,pFrom.GetValue()));
}

///An MCMC step suitable for introducing sample diversity
int fMCMC(long lTime, smc::particle<mChain<double> > & pFrom, smc::rng *pRng)
{
  static smc::particle<mChain<double> > pTo;
  
  mChain<double> * pMC = new mChain<double>;

  for(int i = 0; i < pFrom.GetValue().GetLength(); i++) 
    pMC->AppendElement(pFrom.GetValue().GetElement(i)->value + pRng->Normal(0, 0.5));
  pTo.SetValue(*pMC);
  pTo.SetLogWeight(pFrom.GetLogWeight());

  delete pMC;

  double alpha = exp(logDensity(lTime,pTo.GetValue()) - logDensity(lTime,pFrom.GetValue()));
  if(alpha < 1)
    if (pRng->UniformS() > alpha) {
      return false;
    }

  pFrom = pTo;
  return true;
}

///A function to be integrated in the path sampling step.
double pIntegrandPS(long lTime, const smc::particle<mChain<double> >& pPos, void* pVoid)
{
  double dPos = pPos.GetValue().GetTerminal()->value;
  return (dPos - THRESHOLD) / (1.0 + exp(ALPHA(lTime) * (dPos - THRESHOLD)));
}

///A function which gives the width distribution for the path sampling step.
double pWidthPS(long lTime, void* pVoid)
{
  if(lTime > 1 && lTime < lIterates)
    return ((0.5)*double(ALPHA(lTime+1.0)-ALPHA(lTime-1.0)));
  else 
    return((0.5)*double(ALPHA(lTime+1.0)-ALPHA(lTime)) +(ALPHA(1)-0.0));
}

//The final state weighting function -- how likely is a random path from this distribution to hit the rare set...
double pIntegrandFS(const mChain<double>& dPos, void* pVoid)
{
  if(dPos.GetTerminal()->value > THRESHOLD) {
    return (1.0 + exp(-FTIME*(dPos.GetTerminal()->value-THRESHOLD)));
  }
  else
    return 0;
}

