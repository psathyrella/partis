//
// File: DiscreteRateAcrossSitesTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: 2005
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _DISCRETERATESACROSSSITESTREELIKELIHOOD_H_
#define _DISCRETERATESACROSSSITESTREELIKELIHOOD_H_

#include "TreeLikelihood.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/ParameterList.h>

namespace bpp
{

/**
 * @brief Interface for rate across sites (RAS) implementation.
 *
 * This interface provides methods for dealing with RAS models.
 */
class DiscreteRatesAcrossSitesTreeLikelihood:
  public virtual TreeLikelihood
{
	public:
		DiscreteRatesAcrossSitesTreeLikelihood() {}
		virtual ~DiscreteRatesAcrossSitesTreeLikelihood() {}

	public:

		/**
		 * @brief Get the rate distribution used for the computation.
		 *
		 * @return A const pointer toward the rate distribution of this instance.
		 */
		virtual const DiscreteDistribution* getRateDistribution() const = 0;

		/**
		 * @brief Get the rate distribution used for the computation.
		 *
		 * @return A pointer toward the rate distribution of this instance.
		 */
		virtual DiscreteDistribution* getRateDistribution() = 0;

		/**
		 * @brief Get the likelihood for a site knowing its rate class.
		 *
		 * @param site      The site index.
		 * @param rateClass The rate class index.
		 * @return The likelihood for the specified site and rate class.
		 */
		virtual double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const = 0;
		
		/**
		 * @brief Get the logarithm of the likelihood for a site knowing its rate class.
		 *
		 * @param site      The site index.
		 * @param rateClass The rate class index.
		 * @return The logarithm of the likelihood for the specified site and rate class.
		 */
		virtual double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const = 0;
	
		/**
		 * @brief Get the likelihood for a site knowing its rate class and its ancestral state.
		 *
		 * @param site      The site index.
		 * @param rateClass The rate class index.
		 * @param state     The ancestral state.
		 * @return The likelihood for the specified site and rate class and ancestral state.
		 */
		virtual double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const = 0;
		
		/**
		 * @brief Get the logarithm of the likelihood for a site knowing its rate class and its ancestral state.
		 *
		 * @param site      The site index.
		 * @param rateClass The rate class index.
		 * @param state     The ancestral state.
		 * @return The logarithm of the likelihood for the specified site and rate class and ancestral state..
		 */
		virtual double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const = 0;

		/**
		 * @brief Get the likelihood for each site and each rate class.
		 *
		 * @return A two-dimension vector with all likelihoods.
		 */
		virtual VVdouble getLikelihoodForEachSiteForEachRateClass() const = 0;
		
		/**
		 * @brief Get the logarithm of the likelihood for each site and each rate class.
		 *
		 * @return A two-dimension vector with all log likelihoods:
		 * <code>V[i][j] =</code> likelihood of site i and rate class j.
		 */
		virtual VVdouble getLogLikelihoodForEachSiteForEachRateClass() const = 0;
		
		/**
		 * @brief Get the likelihood for each site and each rate class and each state.
		 *
		 * @return A three-dimension vector with all likelihoods.
		 */
		virtual VVVdouble getLikelihoodForEachSiteForEachRateClassForEachState() const = 0;
		
		/**
		 * @brief Get the logarithm of the likelihood for each site and each rate class and each state.
		 *
		 * @return A three-dimension vector with all log likelihoods:
		 * <code>V[i][j][k} =</code> likelihood of site i and rate class j and state k.
		 */
		virtual VVVdouble getLogLikelihoodForEachSiteForEachRateClassForEachState() const = 0;

		/**
		 * @brief Get the posterior probability for each site of belonging to a
		 * particular rate class.
		 *
		 * @return A two-dimension vector with all posterior probabilities:
		 * <code>V[i][j] =</code> probablity for site i of belonging to rate class j.
		 */
		virtual VVdouble getPosteriorProbabilitiesOfEachRate() const = 0;
		
		/**
		 * @brief Get the posterior rate class (the one with maximum posterior
		 * probability) for each site.
		 *
		 * @return A vector with all rate classes indexes.
		 */
		virtual std::vector<size_t> getRateClassWithMaxPostProbOfEachSite() const = 0;

		/**
		 * @brief Get the posterior rate (the one with maximum posterior
		 * probability) for each site.
		 *
		 * @return A vector with all rate classes indexes.
		 */
		virtual Vdouble getRateWithMaxPostProbOfEachSite() const = 0;
	
		/**
		 * @brief Get the posterior rate, i.e. averaged over all classes
		 * and weighted with posterior probabilities, for each site.
		 *
		 * @return A vector with all rates.
		 */
		virtual Vdouble getPosteriorRateOfEachSite() const = 0;

		/**
		 * @brief Get the parameters associated to the rate distirbution.
		 *
		 * @return A ParameterList object with all rate distribution parameters.
		 */
		virtual ParameterList getRateDistributionParameters() const = 0;

		/**
		 * @brief Get the number of classes.
		 *
		 * @return The Number of classes.
		 */
		virtual size_t getNumberOfClasses() const = 0;

    /**
     * @brief Retrieves all Pij(t) for a particular branch, defined by the upper node.
     *
     * These intermediate results may be used by other methods.
     *
     * @param nodeId The node defining the branch of interest.
     * @param siteIndex The position in the alignment.
     * @return An array of dimension 3, where a[c][x][y] is the probability of substituting from x to y while being in rate class c.
     */
    virtual VVVdouble getTransitionProbabilitiesPerRateClass(int nodeId, size_t siteIndex) const = 0;
		
};

} //end of namespace bpp.

#endif //_DISCRETERATESACROSSSITESTREELIKELIHOOD_H_

