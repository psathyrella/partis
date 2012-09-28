#include <iostream>
#include <cstring>

//! \file
//! \brief This file contains the untemplated functions used for dealing with random number generation.

#include "rng.hh"
#include "smctc.hh"

namespace smc {
    ///The GSL provides a mechanism for obtaining a list of available random number generators.
    ///
    ///This class provides a wrapper for this mechanism and makes it simple to implement software which allows
    ///the nature of the random number generator to be specified at runtime (and to allow it to be a user-specifiable
    ///parameter).
    ///
    ///For example, gslrnginfo::GetNumber can be used to determine how many RNGs are available and gslrnginfo::GetNameByIndex
    ///can then be used to populate a list box with their names. Once the user has selected on the gslrnginfo::GetPointerByName
    ///function can be used to obtain a pointer to the appropriate type and this can be used to produce a random number generator
    ///of the desired type.
    ///
    ///There should be exactly one instance of this class in any program, and that instance is created by
    ///the library. A singleton DP implementation ensures that any additional attempt to instatiate the class simply returns 
    ///a reference to the existing instance.

    gslrnginfo::gslrnginfo()
    {
	typePtArray = gsl_rng_types_setup();

	const gsl_rng_type** ptIndex = typePtArray;
	nNumber = 0;
	while(ptIndex[0]) {
	    ptIndex++;
	    nNumber++;
	}
	return;
    }

    ///Returns a pointer to a single instance of this class.
    gslrnginfo* gslrnginfo::GetInstance()
    {
	static gslrnginfo ginfo;

	return &ginfo;
    }

    ///This function returns the number of available generators
    int gslrnginfo::GetNumber(void)
    {
	return nNumber;
    }

    ///This function returns the name of the specified generator
    const char* gslrnginfo::GetNameByIndex(int nIndex)
    {
	if(0 <= nIndex && nIndex < nNumber)
	    return typePtArray[nIndex]->name;
	return NULL;
    }

    ///This function returns a pointer to the specified generator type
    const gsl_rng_type* gslrnginfo::GetPointerByIndex(int nIndex)
    {
	if(0 <= nIndex && nIndex < nNumber)
	    return typePtArray[nIndex];
	return NULL;
    }

    ///This function returns a pointer to the specified generator type
    const gsl_rng_type* gslrnginfo::GetPointerByName(const char* szName)
    {
	for(int n = 0; n < nNumber; n++)
	    if(!strcmp(typePtArray[n]->name, szName))
		return typePtArray[n];

	return NULL;
    }


    ///When called without any arguments, the constructor for the smc::rng class simply allocates a buffer for a
    ///random number generator of type gsl_rng_default (something which can be set at run-time via an environment
    ///variable) using its default seed (which again can be over-ridden using an environment variable).
    rng::rng(void)
    {
	gsl_rng_env_setup();
	type = gsl_rng_default;
	pWorkspace = gsl_rng_alloc(gsl_rng_default);
    }

    ///When called with a single argument, the constructor for the smc::rng class allocates a buffer for a
    ///random number generator of the specified type and initialises it with the default seed (which can be set using
    ///and environment variable if one wishes to vary it at run-time).
    ///
    ///\param Type The type of a GSL random number generator
    rng::rng(const gsl_rng_type* Type)
    {
	gsl_rng_env_setup();
	type = Type;
	pWorkspace = gsl_rng_alloc(Type);
    }

    ///When called with a pair of arguments, the constructor for the smc::rng class allocates a buffer for the specified
    ///random number generator type and initialises it with the specified seed (note that zero has special significance and
    ///is used to specify the seed with which the generator was originally used).
    ///
    ///\param Type The type of a GSL random number generator
    ///\param lSeed The value with which the generator is to be seeded
    rng::rng(const gsl_rng_type* Type, unsigned long int lSeed)
    {
	gsl_rng_env_setup();
	type = Type;
	pWorkspace = gsl_rng_alloc(Type);
	gsl_rng_set(pWorkspace, lSeed);
    }

    ///The destructor presently does no more than call the gsl_rng_free function to deallocate the memory which was
    ///previously allocate to the random number generator.
    rng::~rng()
    {
	gsl_rng_free(pWorkspace);    
    }

    ///This function returns a pointer to the underlying GSL random number generator which may be used to provide random
    ///number facilities which are not explicitly provided by the intermediate layer of smc::rng.
    gsl_rng* rng::GetRaw(void)
    {
	return pWorkspace;
    }

    ///This function simply passes the relevant arguments on to gsl_ran_multinomial.
    ///     \param n Number of entities to assign.
    ///     \param k Number of categories.
    ///     \param w Weights of category elements
    ///     \param X Array in which to return the sample values.
    void rng::Multinomial(unsigned n, unsigned k, const double* w, unsigned* X)
    {
	gsl_ran_multinomial(pWorkspace, k, n, w, X);
    }


    ///This function simply calls gsl_rng_uniform_int and shifts the
    /// result as appropriate such that the result is an integer generated uniformly from those between
    /// the two arguments (inclusive of those points).
    ///
    ///     \param lMin The smallest value which can be returned
    ///      \param lMax the largest value which can be returned
    long rng::UniformDiscrete(long lMin, long lMax)
    {
	return gsl_rng_uniform_int(pWorkspace, lMax - lMin + 1) + lMin;
    }

    ///This function simply calls gsl_ran_beta with the specified parameters.
    ///     \param da The parameter associated with "x".
    ///     \paran db The parameter associated with "1-x".
    double rng::Beta(double da, double db)
    {
	return gsl_ran_beta(pWorkspace, da, db);
    }

    ///This function simply calls gsl_ran_cauchy with the specified parameter.
    ///     \param dScale The scale parameter of the distribution.
    double rng::Cauchy(double dScale)
    {
	return gsl_ran_cauchy(pWorkspace, dScale);
    }

    ///This function simply calls gsl_ran_exponential with the specified parameters.
    ///     \param dMean The scale (not rate) (and mean) of the distribution.
    double rng::Exponential(double dMean)
    {
	return gsl_ran_exponential(pWorkspace, dMean);
    }

    ///This function simply calls gsl_ran_gamma with the specified parameters.
    ///     \param dAlpha The shape of the distribution (integers lead to Erlang distributions)
    ///     \param dBeta The scale (not rate) of the distribution.
    double rng::Gamma(double dAlpha, double dBeta)
    {
	return gsl_ran_gamma(pWorkspace, dAlpha, dBeta);
    }

    ///This function simply calls gsl_ran_laplace with the specified parameters.
    ///     \param dScale The scale (not rate) of the distribution.
    ///
    double rng::Laplacian(double dScale)
    {
	return gsl_ran_laplace(pWorkspace, dScale);
    }

    ///This function simply calls gsl_ran_lognormal with the specified parameters.
    ///     \param dMu The location parameter of the distribution.
    ///     \param dSigma The scale parameter of the distribution.
    double rng::Lognormal(double dMu, double dSigma)
    {
	return gsl_ran_lognormal(pWorkspace, dMu, dSigma);
    }

    ///This function simply calls gsl_ran_gaussian with the specified standard deviation and shifts the result.
    ///     \param dMean The mean of the distribution.
    ///     \param dStd  The standard deviation of the distribution
    double rng::Normal(double dMean, double dStd)
    {
	return dMean + gsl_ran_gaussian(pWorkspace, dStd);
    }
  
    ///This function simply calls gsl_ran_ugaussian returns the result. 
    double rng::NormalS(void)
    {
	return gsl_ran_ugaussian(pWorkspace);
    }

    ///This function simply calls gsl_ran_gaussian_tail with the specified parameters and performs appropriate shifting.
    ///     \param dMean The mean of the distribution.
    ///     \param dStd  The standard deviation of the distribution
    ///     \param dThreshold The lower truncation threshold.
    double rng::NormalTruncated(double dMean, double dStd, double dThreshold)
    {
	return dMean + gsl_ran_gaussian_tail(pWorkspace, dThreshold - dMean, dStd);
    }

    ///This function simply calls gsl_ran_tdist with the specified number of degrees of freedom.
    ///     \param dDF The number of degrees of freedom.
    double rng::StudentT(double dDF)
    {
	return gsl_ran_tdist(pWorkspace, dDF);
    }

    ///This function simply calls gsl_rng_uniform and scales and shifts appropriately. 
    ///     \param dMin The lowest value with positive density.
    ///     \param dMax The largest value with positive density.
    double rng::Uniform(double dMin, double dMax)
    {
	double rValue;
 
	rValue = gsl_rng_uniform(pWorkspace);
	rValue *= (dMax - dMin);
	rValue += dMin;

	return rValue;
    }

    ///This function simply calls gsl_rng_uniform. 
    double rng::UniformS(void) {
	return gsl_rng_uniform(pWorkspace); 
    }
}
