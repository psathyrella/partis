//
// File: PatternTools.h
// Created by: Julien Dutheil
// Created on: Thu Mar 20 13:36:53 2003
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
 
#ifndef _PATTERNTOOLS_H_
#define _PATTERNTOOLS_H_

#include "Tree.h"

#include <Bpp/Numeric/VectorTools.h>

//From SeqLib:
#include <Bpp/Seq/Site.h>
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <map>

namespace bpp
{

/**
 * @brief Utilitary methods to compute site patterns.
 *
 * Theses methods are mainly designed to save computation in likelihood
 * and parsimony methods.
 */
class PatternTools
{
	public:
    /**
     * @brief Extract the sequences corresponding to a given subtree.
     *
     * @param sequenceSet The container to look in.
     * @param node        The root node of the subtree to check.
     * @return A new site container with corresponding sequences.
     * @throw Exception if an error occured.
     */
		static SiteContainer* getSequenceSubset(const SiteContainer& sequenceSet, const Node& node) throw (Exception);
    /**
     * @brief Extract the sequences corresponding to a given set of names.
     *
     * @param sequenceSet The container to look in.
     * @param names       The names of the sequences to look for.
     * @return A new site container with corresponding sequences.
     * @throw Exception if an error occured.
     */
		static SiteContainer* getSequenceSubset(const SiteContainer& sequenceSet, const std::vector<std::string>& names) throw (Exception);
		/**
     * @brief Compress a site container by removing duplicated sites.
     *
     * @param sequenceSet The container to look in.
     * @return A new site container with unique sites.
     * @throw Exception if an error occured.
     */
    static SiteContainer* shrinkSiteSet(const SiteContainer& sequenceSet) throw (Exception);

		/**
     * @brief Look for the occurence of each site in sequences1 in sequences2 and send the
     * position of the first occurence, or -1 if not found.
     *
     * @param sequences1 First container.
     * @param sequences2 Second container.
     * @return A vecotr of positions.
     */
		static Vint getIndexes(const SiteContainer& sequences1, const SiteContainer& sequences2);
};


} //end of namespace bpp.

#endif	//_PATTERNTOOLS_H_

