//
// File: MaseTools.h
// Created by: Julien Dutheil
// Created on: Tue Apr  1 09:16:59 2003
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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
 
#ifndef _MASETOOLS_H_
#define _MASETOOLS_H_

#include "../Container/SequenceContainer.h"
#include "../Container/OrderedSequenceContainer.h"
#include "../Container/SequenceContainerTools.h"
#include "../Container/SiteContainer.h"
#include "../Container/SiteContainerTools.h"
#include <Bpp/Exceptions.h>

namespace bpp
{

/**
 * @brief Utilitary methods that deal with the Mase format.
 *
 * This class particularily covers the Mase+ format, which allows
 * site and sequence selection.
 * Mase+ tags are in the header of the mase file, which is stored
 * in the 'general comment' section of sequence containers.
 * Most of the methods here hence work on the general comments associated
 * to a container.
 */
class MaseTools
{
	public:
    
    /**
     * @brief Get a site selection from a Mase+ header file.
     *
     * @param maseFileHeader The header of the mase+ file as comments lines.
     * @param setName        The name of the set to retrieve.
     * @throw IOException If the specified set is not found.
     */
		static SiteSelection getSiteSet(const Comments& maseFileHeader, const std::string& setName) throw (IOException);

    /**
     * @brief Get a sequence selection from a Mase+ header file.
     *
     * @param maseFileHeader The header of the mase+ file as comments lines.
     * @param setName        The name of the set to retrieve.
     * @throw IOException If the specified set is not found.
     */
    static SequenceSelection getSequenceSet(const Comments& maseFileHeader, const std::string& setName) throw (IOException);

    /**
     * @brief Create a new container corresponding to a site set given in the mase+ format.
     *
     * A new VectorSiteContainer is created, whose destruction is up to the user.
     * The container passed as argument must have 'general comments' in the mase+ format.
     * This function calls the getSiteSet() function on the comments and then calls for
     * SiteContainerTools::getSelectedSites() on the selection.
     *
     * @param sequences The container to get the sites from.
     * @param setName   The name of the set to retrieve.
     * @throw IOException If the specified set is not found.
     */
    static SiteContainer* getSelectedSites(const SiteContainer& sequences, const std::string& setName) throw (IOException);

    /**
     * @brief Create a new container corresponding to a site set given in the mase+ format.
     *
     * A new VectorSequenceContainer is created, whose destruction is up to the user.
     * The container passed as argument must have 'general comments' in the mase+ format.
     * This function calls the getSequenceSet() function on the comments and then calls for
     * SiteContainerTools::getSelectedSequences() on the selection.
     *
     * @param sequences The container to get the sequence from.
     * @param setName   The name of the set to retrieve.
     * @throw IOException If the specified set is not found.
     */
    static SequenceContainer* getSelectedSequences(const OrderedSequenceContainer& sequences, const std::string & setName) throw (IOException);

    /**
     * @brief Get a list of all available site selections.
     *
     * @param maseHeader Comments as described in the Mase+ format specification.
     * @return A vector of selection names.
     */
    static std::map<std::string, size_t> getAvailableSiteSelections(const Comments & maseHeader);

    /**
     * @brief Get a list of all available sequences selections.
     *
     * @param maseHeader Comments as described in the Mase+ format specification.
     * @return A vector of selection names.
     */
		static std::map<std::string, size_t> getAvailableSequenceSelections(const Comments & maseHeader);

		/**
		 * @brief Get the phase of a given coding region from a mase+ header.
		 *
		 * Look for a /codon_start tag with a phase indice and a site selection with name setName.
		 *
		 * @param maseFileHeader Comments in Mase+ format.
		 * @param setName a cds site selection name.
		 * @return 1,2 or 3.
		 * @throw Exception If no corresponding tag found in file.
		 */
		static size_t getPhase(const Comments & maseFileHeader, const std::string &setName) throw (Exception);

};

} //end of namespace bpp.

#endif	//_MASETOOLS_H_

