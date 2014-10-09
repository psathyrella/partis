//
// File: DRTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Sat Dec 29 19:04 2007
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

#ifndef _DRTREELIKELIHOOD_H_
#define _DRTREELIKELIHOOD_H_

#include "AbstractNonHomogeneousTreeLikelihood.h"
#include "DRASDRTreeLikelihoodData.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{

/**
 * @brief Interface for double-recursive (DR) implementation of the likelihood computation.
 *
 * In the DR implementation, three conditional likelihoods are stored at each node, corresponding to the three connected substrees
 * (in case of multifurcations, there may have even more).
 * The DR implementation hence uses 3x more memory than the simple recursive (R) implementation.
 * However, the likelihood of the tree can be computed independently at each node of the tree,
 * which is very convenient for topology estimation (the likelihood change of each NNI movement can be computed directly),
 * ancestral state reconstruction (all marginal ancestral states can be computed in one pass over the tree), or for
 * substitution mapping.
 *
 * This interface provides
 * - a method to access the DR likelihood data structure,
 * - a method to compute the likelihood array at each node.
 *
 * For now, this interface inherits from DiscreteRatesAcrossSitesTreeLikelihood and not TreeLikelihood,
 * since the data structure available accounts for rate across site variation.
 * This may change in the future.
 *
 * @see DRTreeLikelihoodTools
 */
class DRTreeLikelihood:
  public virtual DiscreteRatesAcrossSitesTreeLikelihood
{
  public:
    DRTreeLikelihood() {}
    virtual ~DRTreeLikelihood() {}

#ifndef NO_VIRTUAL_COV
    DRTreeLikelihood* clone() const = 0;
#endif

  public:
  
    
  public:

    /**
     * @name Get the likelihood data structure associated to this class.
     *
     * @{
     */
    virtual DRASDRTreeLikelihoodData* getLikelihoodData() = 0;
    virtual const DRASDRTreeLikelihoodData* getLikelihoodData() const = 0;
    /** @} */
  
    /**
     * @brief Compute the likelihood array at a given node.
     *
     * @param nodeId The id of the node to consider.
     * @param likelihoodArray The array where to store the results.
     */
    virtual void computeLikelihoodAtNode(int nodeId, VVVdouble& likelihoodArray) const = 0;

};

} //end of namespace bpp.

#endif  //_DRTREELIKELIHOOD_H_

