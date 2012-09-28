#include "smctc.hh"

#include "markovchains/markovchain.h"

extern long lIterates;
extern long lNumber;
extern long lChainLength;
extern double dSchedule;
extern double dThreshold;

double logDensity(long lTime, const mChain<double> & X);

smc::particle<mChain<double> > fInitialise(smc::rng *pRng);
long fSelect(long lTime, const smc::particle<mChain<double> > & p, smc::rng *pRng);
void fMove1(long lTime, smc::particle<mChain<double> > & pFrom, smc::rng *pRng);
void fMove2(long lTime, smc::particle<mChain<double> > & pFrom, smc::rng *pRng);
int fMCMC(long lTime, smc::particle<mChain<double> > & pFrom, smc::rng *pRng);

double pIntegrandPS(long lTime, const smc::particle<mChain<double> >& pPos, void* pVoid);
double pWidthPS(long lTime, void* pVoid);
double pIntegrandFS(const mChain<double>& dPos, void* pVoid);

///The number of grid elements to either side of the current state for the single state move
#define GRIDSIZE 12
///The value of alpha at the specified time
#define ALPHA(T) (double(T)*double(dSchedule) / double(lIterates))
///The terminal version of alpha
#define FTIME    (ALPHA(lIterates))
///The exceedance threshold which we are interested in.
#define THRESHOLD dThreshold
///The number of steps in the Markov chain
#define PATHLENGTH lChainLength
