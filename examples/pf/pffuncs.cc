#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>

#include "smctc.hh"
#include "pffuncs.hh"

using namespace std;

double var_s0 = 4;
double var_u0 = 1;
double var_s  = 0.02;
double var_u  = 0.001;

double scale_y = 0.1;
double nu_y = 10.0;
double Delta = 0.1;

///The function corresponding to the log likelihood at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the index of the current distribution)
///  \param X     The state to consider 
double logLikelihood(long lTime, const cv_state & X)
{
  return - 0.5 * (nu_y + 1.0) * (log(1 + pow((X.x_pos - y[lTime].x_pos)/scale_y,2) / nu_y) + log(1 + pow((X.y_pos - y[lTime].y_pos)/scale_y,2) / nu_y));
}

///A function to initialise particles

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<cv_state> fInitialise(smc::rng *pRng)
{
  cv_state value;
  
  value.x_pos = pRng->Normal(0,sqrt(var_s0));
  value.y_pos = pRng->Normal(0,sqrt(var_s0));
  value.x_vel = pRng->Normal(0,sqrt(var_u0));
  value.y_vel = pRng->Normal(0,sqrt(var_u0));

  return smc::particle<cv_state>(value,logLikelihood(0,value));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng)
{
  cv_state * cv_to = pFrom.GetValuePointer();

  cv_to->x_pos += cv_to->x_vel * Delta + pRng->Normal(0,sqrt(var_s));
  cv_to->x_vel += pRng->Normal(0,sqrt(var_u));
  cv_to->y_pos += cv_to->y_vel * Delta + pRng->Normal(0,sqrt(var_s));
  cv_to->y_vel += pRng->Normal(0,sqrt(var_u));

  pFrom.AddToLogWeight(logLikelihood(lTime, *cv_to));
}
