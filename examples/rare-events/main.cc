#include <iostream>
#include <cmath>
#include "simfunctions.hh"

///Length of Markov Chain
long lChainLength = 15;
///Number of distributions used
long lIterates;
///Annealing schedule constant
double dSchedule = 30.0;
///Rare event threshold
double dThreshold = 5.0;

int main(int argc, char** argv)
{
  cout << "Number of Particles: ";
  long lNumber;
  cin >> lNumber;
  cout << "Number of Iterations: ";
  cin >> lIterates;
  cout << "Threshold: ";
  cin >> dThreshold;
  cout << "Schedule Constant: ";
  cin >> dSchedule;

  try{
    ///An array of move function pointers
    void (*pfMoves[])(long, smc::particle<mChain<double> > &,smc::rng*) = {fMove1, fMove2};
    smc::moveset<mChain<double> > Moveset(fInitialise, fSelect, sizeof(pfMoves) / sizeof(pfMoves[0]), pfMoves, fMCMC);
    smc::sampler<mChain<double> > Sampler(lNumber, SMC_HISTORY_RAM);

    Sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED,0.5);
    Sampler.SetMoveSet(Moveset);

    Sampler.Initialise();
    Sampler.IterateUntil(lIterates);
      
    ///Estimate the normalising constant of the terminal distribution
    double zEstimate = Sampler.IntegratePathSampling(pIntegrandPS, pWidthPS, NULL) - log(2.0);
    ///Estimate the weighting factor for the terminal distribution
    double wEstimate = Sampler.Integrate(pIntegrandFS, NULL);
      
    cout << zEstimate << " " << log(wEstimate) << " " << zEstimate + log(wEstimate) << endl;
  }
  catch(smc::exception  e)
    {
      cerr << e;
      exit(e.lCode);
    }

  return 0;
}

