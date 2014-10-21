//
// File: AbstractDiscreteRatesAcrossSitesTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Wed Jun 15 09:42 2005
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

#ifndef _ABSTRACTDISCRETERATEACROSSSITESTREELIKELIHOOD_H_
#define _ABSTRACTDISCRETERATEACROSSSITESTREELIKELIHOOD_H_

#include "AbstractTreeLikelihood.h"
#include "DiscreteRatesAcrossSitesTreeLikelihood.h"
#include "../Model/SubstitutionModel.h"

namespace bpp
{

/**
 * @brief Partial implementation of the DiscreteRatesAcrossSitesTreeLikelihood interface.
 *
 * It contains a pointer toward a DiscreteDistribution object.
 * This object may be shared by several instances of the class.
 */
class AbstractDiscreteRatesAcrossSitesTreeLikelihood:
	public AbstractTreeLikelihood,
	public virtual DiscreteRatesAcrossSitesTreeLikelihood
{
	protected:
		DiscreteDistribution* rateDistribution_;
		
	public:
		AbstractDiscreteRatesAcrossSitesTreeLikelihood(
			DiscreteDistribution* rDist,
			bool verbose = true
		)	throw (Exception);
		
    AbstractDiscreteRatesAcrossSitesTreeLikelihood(
        const AbstractDiscreteRatesAcrossSitesTreeLikelihood& tl) :
      AbstractTreeLikelihood(tl),
      rateDistribution_(tl.rateDistribution_)
    {}

    AbstractDiscreteRatesAcrossSitesTreeLikelihood& operator=(
        const AbstractDiscreteRatesAcrossSitesTreeLikelihood& tl)
    {
      AbstractTreeLikelihood::operator=(tl);
      rateDistribution_ = tl.rateDistribution_;
      return *this;
    }

		virtual ~AbstractDiscreteRatesAcrossSitesTreeLikelihood() {}
		
	public:
		
		/**
		 * @name The TreeLikelihood interface.
		 *
		 * Other methods are implemented in the AbstractTreeLikelihood class.
		 *
		 * @{
		 */
		double getLikelihoodForASiteForAState (size_t site, int state) const;
		double getLogLikelihoodForASiteForAState(size_t site, int state) const;
		ParameterList getDerivableParameters() const;
		ParameterList getNonDerivableParameters() const;
    VVdouble getTransitionProbabilities(int nodeId, size_t siteIndex) const;
	 		/** @} */

		/**
		 * @name The DiscreteRatesAcrossSites interface implementation:
		 *
		 * @{
		 */
		const DiscreteDistribution* getRateDistribution() const { return rateDistribution_; }
		      DiscreteDistribution* getRateDistribution()       { return rateDistribution_; }
		size_t getNumberOfClasses() const { return rateDistribution_->getNumberOfCategories(); } 
		ParameterList getRateDistributionParameters() const;
		VVdouble getLikelihoodForEachSiteForEachRateClass() const;
		VVdouble getLogLikelihoodForEachSiteForEachRateClass() const;
		VVVdouble getLikelihoodForEachSiteForEachRateClassForEachState() const;
		VVVdouble getLogLikelihoodForEachSiteForEachRateClassForEachState() const;
		VVdouble getPosteriorProbabilitiesOfEachRate() const;
		Vdouble getRateWithMaxPostProbOfEachSite() const;
    std::vector<size_t> getRateClassWithMaxPostProbOfEachSite() const;
		Vdouble getPosteriorRateOfEachSite() const;
		/** @} */

		/**
     * @name Generic tools to deal with likelihood arrays
     *
     * @{
     */
    
    /**
     * @brief Set all conditional likelihoods to 1.
     *
     * @param likelihoodArray the likelihood array.
     */
    static void resetLikelihoodArray(VVVdouble & likelihoodArray);

    /**
     * @brief Print the likelihood array to terminal (debugging tool).
     * 
     * @param likelihoodArray the likelihood array.
     */
		static void displayLikelihoodArray(const VVVdouble & likelihoodArray);

    /** @} */
    
};

} //end of namespace bpp.

#endif //_ABSTRACTDISCRETERATEACROSSSITESTREELIKELIHOOD_H_

