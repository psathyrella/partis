//
// File: DRTreeLikelihoodTools.h
// Created by: Julien Dutheil
// Created on: Mon Janv 17 09:56 2005
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _DRTREELIKELIHOODTOOLS_H_
#define _DRTREELIKELIHOODTOOLS_H_

#include "TreeLikelihoodTools.h"
#include "DRTreeLikelihood.h"
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>

namespace bpp
{

/**
 * @brief Utilitary methods dealing with objects implementing the DRTreeLikelihood interface.
 */
class DRTreeLikelihoodTools:
  public TreeLikelihoodTools
{

  public:
    /**
     * @brief Compute the posterior probabilities for each state and each rate of each distinct site.
     *
     * @param drl A DR tree likelihood object.
     * @param nodeId The id of the node at which probabilities must be computed.
     * @return A 3-dimensional array, with probabilities for each site, each rate and each state.
     */
    static VVVdouble getPosteriorProbabilitiesForEachStateForEachRate(
        const DRTreeLikelihood& drl,
        int nodeId);

    /**
     * @brief Compute the posterior probabilities for each state for a given node.
     *
     * This method calls the getPosteriorProbabilitiesForEachStateForEachRate function
     * and average the probabilities over all sites and rate classes, resulting in a
     * one-dimensionnal frequency array, with one frequency per model state.
     *
     * @param drl A DR tree likelihood object.
     * @param nodeId The id of the node at which probabilities must be computed.
     * @return vector of double with state frequencies for the given node.
     */
    static Vdouble getPosteriorStateFrequencies(
        const DRTreeLikelihood& drl,
        int nodeId);

};

} //end of namespace bpp.

#endif //_DRTREELIKELIHOODTOOLS_H_

