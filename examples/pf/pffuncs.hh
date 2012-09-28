#include "smctc.hh"

class cv_state
{
public:
  double x_pos, y_pos;
  double x_vel, y_vel;
};

class cv_obs
{
public:
  double x_pos, y_pos;
};

double logLikelihood(long lTime, const cv_state & X);

smc::particle<cv_state> fInitialise(smc::rng *pRng);
long fSelect(long lTime, const smc::particle<cv_state> & p, 
	     smc::rng *pRng);
void fMove(long lTime, smc::particle<cv_state> & pFrom, 
	   smc::rng *pRng);

extern double nu_x;
extern double nu_y;
extern double Delta;

extern cv_obs * y; 
