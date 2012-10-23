//   SMCTC: sampler.hh
//
//   Copyright Adam Johansen, 2008.
//
//   This file is part of SMCTC.
//
//   SMCTC is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   SMCTC is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with SMCTC.  If not, see <http://www.gnu.org/licenses/>.

//! \file
//! \brief Defines the overall sampler object.
//!
//! This file defines the smc::sampler class which is used to implement entire particle systems.

#ifndef __SMC_SAMPLER_HH

#define __SMC_SAMPLER_HH 1.0

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "rng.hh"
#include "history.hh"
#include "moveset.hh"
#include "particle.hh"
#include "smc-exception.hh"

///Specifiers for various resampling algorithms:
enum ResampleType { SMC_RESAMPLE_MULTINOMIAL = 0,
                    SMC_RESAMPLE_RESIDUAL,
                    SMC_RESAMPLE_STRATIFIED,
                    SMC_RESAMPLE_SYSTEMATIC
                  };

///Storage types for the history of the particle system.
enum HistoryType { SMC_HISTORY_NONE = 0,
                   SMC_HISTORY_RAM
                 };

namespace smc
{

/// A template class for an interacting particle system suitable for SMC sampling
template <class Space>
class sampler
{
private:
    ///A random number generator.
    rng* pRng;

    ///Number of particles in the system.
    long N;
    ///The current evolution time of the system.
    long T;

    ///The resampling mode which is to be employed.
    ResampleType rtResampleMode;
    ///The effective sample size at which resampling should be used.
    double dResampleThreshold;
    ///Structure used internally for resampling.
    double* dRSWeights;
    ///Structure used internally for resampling.
    unsigned int* uRSCount;
    ///Structure used internally for resampling.
    unsigned int* uRSIndices;

    ///The particles within the system.
    particle<Space> *pParticles;
    ///The set of moves available.
    moveset<Space> Moves;

    ///The number of MCMC moves which have been accepted during this iteration
    int nAccepted;
    ///A flag which tracks whether the ensemble was resampled during this iteration
    int nResampled;

    ///A mode flag which indicates whether historical information is stored
    HistoryType htHistoryMode;
    ///The historical process associated with the particle system.
    history<particle<Space> > History;

public:
    ///Create an particle system containing lSize uninitialised particles with the specified mode.
    sampler(long lSize, HistoryType htHistoryMode);
    ///Create an particle system constaining lSize uninitialised particles with the specified mode and random number generator.
    sampler(long lSize, HistoryType htHistoryMode, const gsl_rng_type* rngType, unsigned long nSeed);
    ///Dispose of a sampler.
    ~sampler();
    ///Calculates and Returns the Effective Sample Size.
    double GetESS(void) const;
    ///Returns a pointer to the History of the particle system
    const history<particle<Space> > * GetHistory(void) const { return &History; }
    ///Returns the number of particles within the system.
    long GetNumber(void) const {return N;}
    ///Return the value of particle n
    const Space &  GetParticleValue(int n) { return pParticles[n].GetValue(); }
    ///Return the logarithmic unnormalized weight of particle n
    double GetParticleLogWeight(int n) { return pParticles[n].GetLogWeight(); }
    ///Return the unnormalized weight of particle n
    double GetParticleWeight(int n) { return pParticles[n].GetWeight(); }
    ///Returns the current evolution time of the system.
    long GetTime(void) const {return T;}
    ///Initialise the sampler and its constituent particles.
    void Initialise(void);
    ///Integrate the supplied function with respect to the current particle set.
    double Integrate(double(*pIntegrand)(const Space &, void*), void* pAuxiliary);
    ///Integrate the supplied function over the path path using the supplied width function.
    double IntegratePathSampling(double(*pIntegrand)(long, const particle<Space>&, void*), double(*pWidth)(long, void*), void* pAuxiliary);
    ///Perform one iteration of the simulation algorithm.
    void Iterate(void);
    ///Cancel one iteration of the simulation algorithm.
    void IterateBack(void);
    ///Perform one iteration of the simulation algorithm and return the resulting ess
    double IterateEss(void);
    ///Perform iterations until the specified evolution time is reached
    void IterateUntil(long lTerminate);
    ///Move the particle set by proposing an applying an appropriate move to each particle.
    void MoveParticles(void);
    ///Resample the particle set using the specified resmpling scheme.
    void Resample(ResampleType lMode);
    ///Sets the entire moveset to the one which is supplied
    void SetMoveSet(moveset<Space>& pNewMoveset) { Moves = pNewMoveset; }
    ///Set Resampling Parameters
    void SetResampleParams(ResampleType rtMode, double dThreshold);
    ///Dump a specified particle to the specified output stream in a human readable form
    std::ostream & StreamParticle(std::ostream & os, long n);
    ///Dump the entire particle set to the specified output stream in a human readable form
    std::ostream & StreamParticles(std::ostream & os);
    ///Allow a human readable version of the sampler configuration to be produced using the stream operator.
    /// std::ostream & operator<< (std::ostream& os, sampler<Space> & s);

private:
    ///Duplication of smc::sampler is not currently permitted.
    sampler(const sampler<Space> & sFrom);
    ///Duplication of smc::sampler is not currently permitted.
    sampler<Space> & operator=(const sampler<Space> & sFrom);
};


/// The constructor prepares a sampler for use but does not assign any moves to the moveset, initialise the particles
/// or otherwise perform any sampling related tasks. Its main function is to allocate a region of memory in which to
/// store the particle set and to initialise a random number generator.
///
/// \param lSize The number of particles present in the ensemble (at time 0 if this is a variable quantity)
/// \param htHM The history mode to use: set this to SMC_HISTORY_RAM to store the whole history of the system and SMC_HISTORY_NONE to avoid doing so.
/// \tparam Space The class used to represent a point in the sample space.
template <class Space>
sampler<Space>::sampler(long lSize, HistoryType htHM)
{
    pRng = new rng();
    N = lSize;
    pParticles = new particle<Space>[lSize];

    //Allocate some storage for internal workspaces
    dRSWeights = new double[N];
    ///Structure used internally for resampling.
    uRSCount  = new unsigned[N];
    ///Structure used internally for resampling.
    uRSIndices = new unsigned[N];

    //Some workable defaults.
    htHistoryMode = htHM;
    rtResampleMode = SMC_RESAMPLE_STRATIFIED;
    dResampleThreshold = 0.5 * N;
}

/// The constructor prepares a sampler for use but does not assign any moves to the moveset, initialise the particles
/// or otherwise perform any sampling related tasks. Its main function is to allocate a region of memory in which to
/// store the particle set and to initialise a random number generator.
///
/// \param lSize The number of particles present in the ensemble (at time 0 if this is a variable quantity)
/// \param htHM The history mode to use: set this to SMC_HISTORY_RAM to store the whole history of the system and SMC_HISTORY_NONE to avoid doing so.
/// \param rngType The type of random number generator to use
/// \param rngSeed The seed to use for the random number generator
/// \tparam Space The class used to represent a point in the sample space.
template <class Space>
sampler<Space>::sampler(long lSize, HistoryType htHM, const gsl_rng_type* rngType, unsigned long rngSeed)
{
    pRng = new rng(rngType, rngSeed);
    N = lSize;
    pParticles = new particle<Space>[lSize];

    //Allocate some storage for internal workspaces
    dRSWeights = new double[N];
    ///Structure used internally for resampling.
    uRSCount  = new unsigned[N];
    ///Structure used internally for resampling.
    uRSIndices = new unsigned[N];

    //Some workable defaults.
    htHistoryMode  = SMC_HISTORY_RAM;
    rtResampleMode = SMC_RESAMPLE_STRATIFIED;
    dResampleThreshold = 0.5 * N;
}


template <class Space>
sampler<Space>::~sampler()
{
    delete pRng;

    if(dRSWeights)
        delete [] dRSWeights;
    if(uRSCount)
        delete [] uRSCount;
    if(uRSIndices)
        delete [] uRSIndices;
}

template <class Space>
double sampler<Space>::GetESS(void) const
{
    long double sum = 0;
    long double sumsq = 0;

    for(int i = 0; i < N; i++)
        sum += expl(pParticles[i].GetLogWeight());

    for(int i = 0; i < N; i++)
        sumsq += expl(2.0 * (pParticles[i].GetLogWeight()));

    return expl(-log(sumsq) + 2 * log(sum));
}

/// At present this function resets the system evolution time to 0 and calls the moveset initialisor to assign each
/// particle in the ensemble.
///
/// Note that the initialisation function must be specified before calling this function.
template <class Space>
void sampler<Space>::Initialise(void)
{
    T = 0;

    for(int i = 0; i < N; i++)
        pParticles[i] = Moves.DoInit(pRng);

    if(htHistoryMode != SMC_HISTORY_NONE) {
        while(History.Pop());
        nResampled = 0;
        History.Push(N, pParticles, 0, historyflags(nResampled));
    }

    return;
}

/// This function returns the result of integrating the supplied function under the empirical measure associated with the
/// particle set at the present time. The final argument of the integrand function is a pointer which will be supplied
/// with pAuxiliary to allow for arbitrary additional information to be passed to the function being integrated.
///
/// \param pIntegrand The function to integrate with respect to the particle set
/// \param pAuxiliary A pointer to any auxiliary data which should be passed to the function

template <class Space>
double sampler<Space>::Integrate(double(*pIntegrand)(const Space&, void*), void * pAuxiliary)
{
    long double rValue = 0;
    long double wSum = 0;
    for(int i = 0; i < N; i++) {
        rValue += expl(pParticles[i].GetLogWeight()) * pIntegrand(pParticles[i].GetValue(), pAuxiliary);
        wSum  += expl(pParticles[i].GetLogWeight());
    }

    rValue /= wSum;
    return (double)rValue;
}

/// This function is intended to be used to estimate integrals of the sort which must be evaluated to determine the
/// normalising constant of a distribution obtain using a sequence of potential functions proportional to densities with respect
/// to the initial distribution to define a sequence of distributions leading up to the terminal, interesting distribution.
///
/// In this context, the particle set at each time is used to make an estimate of the path sampling integrand, and a
/// trapezoidal integration is then performed to obtain an estimate of the path sampling integral which is the natural logarithm
/// of the ratio of normalising densities.
///
/// \param pIntegrand  The quantity which we wish to integrate at each time
/// \param pWidth      A pointer to a function which specifies the width of each

template <class Space>
double sampler<Space>::IntegratePathSampling(double(*pIntegrand)(long, const particle<Space> &, void*), double(*pWidth)(long, void*), void* pAuxiliary)
{
    if(htHistoryMode == SMC_HISTORY_NONE)
        throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "The path sampling integral cannot be computed as the history of the system was not stored.");

    History.Push(N, pParticles, nAccepted, historyflags(nResampled));
    double dRes = History.IntegratePathSampling(pIntegrand, pWidth, pAuxiliary);
    History.Pop();
    return dRes;
}

/// The iterate function:
///         -# appends the current particle set to the history if desired
///          -# moves the current particle set
///         -# checks the effective sample size and resamples if necessary
///         -# performs a mcmc step if required
///         -# increments the current evolution time
template <class Space>
void sampler<Space>::Iterate(void)
{
    IterateEss();
    return;
}

template <class Space>
void sampler<Space>::IterateBack(void)
{
    if(htHistoryMode == SMC_HISTORY_NONE)
        throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "An attempt to undo an iteration was made; unforunately, the system history has not been stored.");

    History.Pop(&N, &pParticles, &nAccepted, NULL);
    T--;
    return;
}

template <class Space>
double sampler<Space>::IterateEss(void)
{
    //Initially, the current particle set should be appended to the historical process.
    if(htHistoryMode != SMC_HISTORY_NONE)
        History.Push(N, pParticles, nAccepted, historyflags(nResampled));

    nAccepted = 0;

    //Move the particle set.
    MoveParticles();

    //Normalise the weights to sensible values....
    double dMaxWeight = -std::numeric_limits<double>::infinity();
    for(int i = 0; i < N; i++)
        dMaxWeight = std::max(dMaxWeight, pParticles[i].GetLogWeight());
    for(int i = 0; i < N; i++)
        pParticles[i].SetLogWeight(pParticles[i].GetLogWeight() - (dMaxWeight));


    //Check if the ESS is below some reasonable threshold and resample if necessary.
    //A mechanism for setting this threshold is required.
    double ESS = GetESS();
    if(ESS < dResampleThreshold) {
        nResampled = 1;
        Resample(rtResampleMode);
    } else
        nResampled = 0;
    //A possible MCMC step should be included here.
    for(int i = 0; i < N; i++) {
        if(Moves.DoMCMC(T + 1, pParticles[i], pRng))
            nAccepted++;
    }
    // Increment the evolution time.
    T++;

    return ESS;
}

template <class Space>
void sampler<Space>::IterateUntil(long lTerminate)
{
    while(T < lTerminate)
        Iterate();
}

template <class Space>
void sampler<Space>::MoveParticles(void)
{
    for(int i = 0; i < N; i++) {
        Moves.DoMove(T + 1, pParticles[i], pRng);
        //  pParticles[i].Set(pNew.value, pNew.logweight);
    }
}

///Perform resampling.
///Note: this procedure sets all particle weights to zero after resampling.
///\param lMode The sampling mode for the sampler.
template <class Space>
void sampler<Space>::Resample(ResampleType lMode)
{
    //Resampling is done in place.
    double dWeightSum = 0;
    unsigned uMultinomialCount;

    //First set up dRSWeights to have the weights, which have been normalized so that the maximum element has weight 1.
    //This keeps things from underflowing when likelihoods are small.
    double dMaxLogWeight = std::numeric_limits<double>::min();
    for(int i = 0; i < N; ++i) {
        dMaxLogWeight = std::max(dMaxLogWeight, pParticles[i].GetLogWeight());
    }
    for(int i = 0; i < N; ++i) {
        dRSWeights[i] = exp(pParticles[i].GetLogWeight() - dMaxLogWeight);
    }

    //Obtain a count of the number of children each particle has via the chosen strategy.
    //This will be stored in uRSCount.
    switch(lMode) {
    case SMC_RESAMPLE_MULTINOMIAL:
        //Generate a multinomial random vector with parameters (N,dRSWeights[1:N]) and store it in uRSCount
        pRng->Multinomial(N, N, dRSWeights, uRSCount);
        break;

    case SMC_RESAMPLE_RESIDUAL:
        //Sample from a suitable multinomial vector and add the integer replicate
        //counts afterwards.
        dWeightSum = 0;
        for(int i = 0; i < N; ++i) {
            dWeightSum += dRSWeights[i];
        }

        uMultinomialCount = N;
        for(int i = 0; i < N; ++i) {
            dRSWeights[i] = N * dRSWeights[i] / dWeightSum;
            uRSIndices[i] = unsigned(floor(dRSWeights[i])); //Reuse temporary storage.
            dRSWeights[i] = (dRSWeights[i] - uRSIndices[i]);
            uMultinomialCount -= uRSIndices[i];
        }
        //Generate a multinomial random vector with parameters (uMultinomialCount,dRSWeights[1:N]) and store it in uRSCount
        pRng->Multinomial(uMultinomialCount, N, dRSWeights, uRSCount);
        for(int i = 0; i < N; ++i)
            uRSCount[i] += uRSIndices[i];
        break;


    case SMC_RESAMPLE_STRATIFIED:
    default: {
        // Procedure for stratified sampling
        // See Appendix of Kitagawa 1996, http://www.jstor.org/stable/1390750,
        // or p.290 of the Doucet et al book, an image of which is at:
        // http://cl.ly/image/200T0y473k1d/stratified_resampling.jpg
        dWeightSum = 0;
        double dWeightCumulative = 0;
        // Calculate the normalising constant of the weight vector
        for(int i = 0; i < N; i++)
            dWeightSum += dRSWeights[i];
        //Generate a random number between 0 and 1/N.
        double dRand = pRng->Uniform(0, 1.0 / ((double)N));
        // Clear out uRSCount.
        for(int i = 0; i < N; ++i)
            uRSCount[i] = 0;
        // 0 <= j < N will index the current sampling step, whereas k will index the previous step.
        int j = 0, k = 0;
        // dWeightCumulative is \tilde{\pi}^r from the Doucet book.
        dWeightCumulative = dRSWeights[0] / dWeightSum;
        // Advance j along until dWeightCumulative > j/N + dRand
        while(j < N) {
            while((dWeightCumulative - dRand) > ((double)j) / ((double)N) && j < N) {
                uRSCount[k]++; // Accept the particle k.
                j++;
                dRand = pRng->Uniform(0, 1.0 / ((double)N));
            }
            k++;
            dWeightCumulative += dRSWeights[k] / dWeightSum;
        }
        break;
    }

    case SMC_RESAMPLE_SYSTEMATIC: {
        // Procedure for stratified sampling but with a common RV for each stratum
        dWeightSum = 0;
        double dWeightCumulative = 0;
        // Calculate the normalising constant of the weight vector
        for(int i = 0; i < N; i++)
            dWeightSum += dRSWeights[i];
        //Generate a random number between 0 and 1/N times the sum of the weights
        double dRand = pRng->Uniform(0, 1.0 / ((double)N));

        int j = 0, k = 0;
        for(int i = 0; i < N; ++i)
            uRSCount[i] = 0;

        dWeightCumulative = dRSWeights[0] / dWeightSum;
        while(j < N) {
            while((dWeightCumulative - dRand) > ((double)j) / ((double)N) && j < N) {
                uRSCount[k]++;
                j++;

            }
            k++;
            dWeightCumulative += dRSWeights[k] / dWeightSum;
        }
        break;

    }
    }

    //Map count to indices to allow in-place resampling.
    //Here j represents the next particle from the previous generation that is going to be dropped in the current
    //sample, and thus can get filled by the resampling scheme.
    for(unsigned int i = 0, j = 0; i < N; ++i) {
        if(uRSCount[i] > 0) {
            uRSIndices[i] = i;
            while(uRSCount[i] > 1) {
                while(uRSCount[j] > 0) ++j; // find next free spot
                uRSIndices[j++] = i; // assign index
                --uRSCount[i]; // decrement number of remaining offsprings
            }
        }
    }

    //Perform the replication of the chosen.
    for(int i = 0; i < N ; ++i) {
        if(uRSIndices[i] != i)
            pParticles[i].SetValue(pParticles[uRSIndices[i]].GetValue());
        //Reset the log weight of the particles to be zero.
        pParticles[i].SetLogWeight(0);
    }
}

/// This function configures the resampling parameters, allowing the specification of both the resampling
/// mode and the threshold at which resampling is used.
///
/// \param rtMode The resampling mode to be used.
/// \param dThreshold The threshold at which resampling is deemed necesary.
///
/// The rtMode parameter should be set to one of the following:
/// -# SMC_RESAMPLE_MULTINOMIAL to use multinomial resampling
/// -# SMC_RESAMPLE_RESIDUAL to use residual resampling
/// -# SMC_RESAMPLE_STRATIFIED to use stratified resampling
/// -# SMC_RESAMPLE_SYSTEMATIC to use systematic resampling
///
/// The dThreshold parameter can be set to a value in the range [0,1) corresponding to a fraction of the size of
/// the particle set or it may be set to an integer corresponding to an actual effective sample size.

template <class Space>
void sampler<Space>::SetResampleParams(ResampleType rtMode, double dThreshold)
{
    rtResampleMode = rtMode;
    if(dThreshold < 1)
        dResampleThreshold = dThreshold * N;
    else
        dResampleThreshold = dThreshold;
}

template <class Space>
std::ostream & sampler<Space>::StreamParticle(std::ostream & os, long n)
{
    os << pParticles[n] << std::endl;
    return os;
}

template <class Space>
std::ostream & sampler<Space>::StreamParticles(std::ostream & os)
{
    for(int i = 0; i < N - 1; i++)
        os << pParticles[i] << std::endl;
    os << pParticles[N - 1] << std::endl;

    return os;
}

}

namespace std
{
/// Produce a human-readable display of the state of an smc::sampler class using the stream operator.

/// \param os The output stream to which the display should be made.
/// \param s  The sampler which is to be displayed.
template <class Space>
std::ostream & operator<< (std::ostream & os, smc::sampler<Space> & s)
{
    os << "Sampler Configuration:" << std::endl;
    os << "======================" << std::endl;
    os << "Evolution Time:   " << s.GetTime() << std::endl;
    os << "Particle Set Size:" << s.GetNumber() << std::endl;
    os << std::endl;
    os << "Particle Set:" << std::endl;
    s.StreamParticles(os);
    os << std::endl;
    return os;
}
}
#endif
