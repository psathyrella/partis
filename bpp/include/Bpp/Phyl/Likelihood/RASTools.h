//
// File: RASTools.h
// Created by: Julien Dutheil
// Created on: May 24 2004
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

Julien.Dutheil@univ-montp2.fr

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

#ifndef _RASTOOLS_H_
#define _RASTOOLS_H_

#include "DiscreteRatesAcrossSitesTreeLikelihood.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{

/**
 * @brief Tools to deal with Rates Across Sites (RAS) models.
 */
class RASTools
{
	public:

		/**
		 * @brief Get the rate distribution estimated from a dataset.
		 *
		 * This methods takes an objet implementing the DiscreteRatesAcroossSites
		 * interface as input and use the posterior probabilities of rates for each site
		 * to generate the corresponding distribution.
		 *
		 * @param treeLikelihood A Likelihood calculation implmenting the RAS model interface.
		 * @return The posterior distribution of rate classes.
		 */
		static DiscreteDistribution * getPosteriorRateDistribution(
				const DiscreteRatesAcrossSitesTreeLikelihood & treeLikelihood);
};

} //end of namespace bpp.

#endif // _RASTOOLS_H_

