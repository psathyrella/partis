//  SMCTC: rng.hh  
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
//

//! \file 
//! \brief Random number generation.
//!
//! This file contains the definitions for the smc::rng and smc::rnginfo class.
//! It wraps the random number generation facilities provided by the GSL and provides a convenient interfaces to access several of its more commonly-used features.

#ifndef __SMC_RNG_HH
#define __SMC_RNG_HH 1.0

extern "C" {
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
}

namespace smc {
  ///A gsl-rng information handling class (not templated)
  class gslrnginfo {
  private:
    ///This is a null terminated array of the available random number generators.
    const gsl_rng_type** typePtArray;
    ///The number of available random number generators.
    int nNumber;

  protected:
    gslrnginfo();

  public:
    ///Returns a reference to the sole static instance of this class.
    static gslrnginfo* GetInstance();

    ///Returns the number of available random number generators.
    int GetNumber();
    ///Returns the name of random number generator number nIndex.
    const char* GetNameByIndex(int nIndex);
    ///Returns a pointer to random number generator nIndex.
    const gsl_rng_type* GetPointerByIndex(int nIndex);
    ///Returns a pointer to the random number generator with name szName (or null if it doesn't exist).
    const gsl_rng_type* GetPointerByName(const char* szName);
  };

  ///The global application instance of the gslrnginfo class:
  extern gslrnginfo rngset;

  ///A random number generator class.

  ///    At present this serves as a wrapper for the gsl random number generation code.
  class rng {
  private:
    ///This is the type of random number generator underlying the class.
    const gsl_rng_type* type;
    ///This is a pointer to the internal workspace of the rng including its current state.
    gsl_rng* pWorkspace;

  public:
    ///Initialise the random number generator using default settings
    rng();
    ///Initialise the random number generator using the default seed for the type
    rng(const gsl_rng_type* Type);
    ///Initialise the random number generator using specified type and seed
    rng(const gsl_rng_type* Type,unsigned long int lSeed);

    ///Free the workspace allocated for random number generation
    ~rng();


    ///Provide access to the raw random number generator
    gsl_rng* GetRaw(void);

    ///Generate a multinomial random vector with parameters (n,w[1:k]) and store it in X
    void Multinomial(unsigned n, unsigned k, const double* w, unsigned* X);
    ///Returns a random integer generated uniformly between the minimum and maximum values specified
    long UniformDiscrete(long lMin, long lMax);

    ///Returns a random number generated from a Beta distribution with the specified parameters.
    double Beta(double da, double db);
    ///Returns a random number generated from a Cauchy distribution with the specified scale parameter.
    double Cauchy(double dScale);
    ///Returns a random number generated from an exponential distribution with the specified mean.
    double Exponential(double dMean);
    ///Return a random number generated from a gamma distribution with shape alpha and scale beta.
    double Gamma(double dAlpha, double dBeta);
    ///Returns a random number generated from a Laplacian distribution with the specified scale.
    double Laplacian(double dScale);
    ///Returns a random number generated from a Lognormal distribution of location mu and scale sigma
    double Lognormal(double dMu, double dSigma);
    ///Return a random number generated from a normal distribution with a specified mean and standard deviation
    double Normal(double dMean, double dStd);
    ///Return a random number generated from a standard normal distribution
    double NormalS(void);
    ///Returns a random number from a normal distribution, conditional upon it exceeding the specified threshold.
    double NormalTruncated(double dMean, double dStd, double dThreshold);
    ///Return a student-t random number generated with a specified number of degrees of freedom
    double StudentT(double dDF);
    ///Return a random number generated uniformly between dMin and dMax
    double Uniform(double dMin, double dMax);    
    ///Returns a random number generated from the standard uniform[0,1) distribution
    double UniformS(void);
  };
}
 
 
#endif
