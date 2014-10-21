//
// File: TreeLikelihoodTools.h
// Created by: Julien Dutheil
// Created on: Tue Jun 30 12:25 2009
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

#ifndef _TREELIKELIHOODTOOLS_H_
#define _TREELIKELIHOODTOOLS_H_

#include "TreeLikelihood.h"

//From the STL:
#include <vector>
#include <map>

namespace bpp
{

/**
 * @brief Utilitary methods that work with TreeLikelihood objects.
 */
class TreeLikelihoodTools
{
  public:
    TreeLikelihoodTools() {}
    virtual ~TreeLikelihoodTools() {}

  public:
    /**
     * @brief Compute the expected ancestral frequencies of all states at all (inner) nodes
     * according to a Markov process defined by a given substitution model.
     *
     * The computation is performed for a given site. If the likelihood object has no
     * site partition, then the method will return the same result for all positions.
     *
     * @param tl          [in] A tree likelihood object.
     * @param site        [in] The site for which the frequencies should be computed.
     * @param frequencies [out] A map where to store the results, as a vector of double (the 
     * size of which being equal to the number of states in the model), and with nodes id as keys.
     * @param alsoForLeaves [opt] Tell if frequencies should also be estimated for terminal nodes.
     * @throw Exception In case something bad happens, like an unvalid model set.
     */
    static void getAncestralFrequencies(
        const TreeLikelihood& tl,
        size_t site,
        std::map<int, std::vector<double> >& frequencies,
        bool alsoForLeaves = false) throw (Exception);

    /**
     * @brief Compute the expected ancestral frequencies of all states at all (inner) nodes
     * according to a Markov process defined by a given substitution model.
     *
     * The computation is averaged over all sites. If the likelihood object has no
     * site partition, then the method will return the same result as all single site numbers.
     *
     * @param tl          [in] A tree likelihood object.
     * @param frequencies [out] A map where to store the results, as a vector of double (the 
     * size of which being equal to the number of states in the model), and with nodes id as keys.
     * @param alsoForLeaves [opt] Tell if frequencies should also be estimated for terminal nodes.
     * @throw Exception In case something bad happens, like an unvalid model set.
     */
    static void getAncestralFrequencies(
        const TreeLikelihood& tl,
        std::map<int, std::vector<double> >& frequencies,
        bool alsoForLeaves = false) throw (Exception);

  private:
    /**
     * @brief Recursive method, for internal use only.
     *
     * @see getAncestralFrequencies()
     */
    static void getAncestralFrequencies_(
        const TreeLikelihood& tl,
        size_t siteIndex,
        int parentId,
        const std::vector<double>& ancestralFrequencies,
        std::map<int, std::vector<double> >& frequencies,
        bool alsoForLeaves) throw (Exception);
 
};

} //end of namespace bpp.

#endif //_TREELIKELIHOODTOOLS_H_

