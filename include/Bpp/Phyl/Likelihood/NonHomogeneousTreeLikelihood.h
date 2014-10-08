//
// File: NonHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Tue Oct 9 16:03 2007
// From file: HomogeneousTreeLikelihood.h
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

#ifndef _NONHOMOGENEOUSTREELIKELIHOOD_H_
#define _NONHOMOGENEOUSTREELIKELIHOOD_H_

#include "TreeLikelihood.h"
#include "../Model/SubstitutionModelSet.h"

namespace bpp
{

/**
 * @brief Specialization of the TreeLikelihood interface for the branch non-homogeneous and non-stationary models.
 *
 * The main difference is that the likelihood depends on the position of the root.
 * The frequencies at the root of the tree are new parameters, which are not necessarily equal to the equilibrium frequencies of the substitution model.
 *
 * This interface further assumes that the substitution model is the same for all sites, for a given branch.
 *
 * @see SubstitutionModelSet.
 */
class NonHomogeneousTreeLikelihood :
	public virtual TreeLikelihood
{
	public:
#ifndef NO_VIRTUAL_COV
    NonHomogeneousTreeLikelihood* clone() const = 0;
#endif

  public:
    const SubstitutionModel* getSubstitutionModel(int nodeId, size_t siteIndex) const throw (NodeNotFoundException)
    {
      return getSubstitutionModelForNode(nodeId);
    }

    SubstitutionModel* getSubstitutionModel(int nodeId, size_t siteIndex) throw (NodeNotFoundException)
    {
      return getSubstitutionModelForNode(nodeId);
    }

    /**
     * @brief Get the substitution model associated to a given node.
     *
     * @param nodeId The id of the request node.
     * @return A pointer toward the corresponding model.
     * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
     */
    virtual const SubstitutionModel* getSubstitutionModelForNode(int nodeId) const throw (NodeNotFoundException) = 0;

    /**
     * @brief Get the substitution model associated to a given node.
     *
     * @param nodeId The id of the request node.
     * @return A pointer toward the corresponding model.
     * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
     */
    virtual SubstitutionModel* getSubstitutionModelForNode(int nodeId) throw (NodeNotFoundException) = 0;

    /**
     * @return The set of substitution models associated to this instance.
     */
    virtual const SubstitutionModelSet* getSubstitutionModelSet() const = 0;

    /**
     * @return The set of substitution models associated to this instance.
     */
    virtual SubstitutionModelSet* getSubstitutionModelSet() = 0;
    
    /**
     * @return Set the substitution models for this instance.
     * @throw Exception If the model could not be set (for instance, because of a wrong alphabet type).
     */
    virtual void setSubstitutionModelSet(SubstitutionModelSet* model) throw (Exception) = 0;

    /**
     * @return The parameters on which the root frequencies depend.
     */
    virtual ParameterList getRootFrequenciesParameters() const = 0;

};

} //end of namespace bpp.

#endif	//_NONHOMOGENEOUSTREELIKELIHOOD_H_

