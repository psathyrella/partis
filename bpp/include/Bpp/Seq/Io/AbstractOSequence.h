//
// File: AbstractOSequence.h
// Created by: Julien Dutheil
// Created on: ?
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

#ifndef _ABSTRACTOSEQUENCE_H_
#define _ABSTRACTOSEQUENCE_H_

#include "OSequence.h"
#include "../Alphabet/Alphabet.h"
#include "../Container/VectorSequenceContainer.h"

// From the STL:
#include <string>
#include <fstream>

namespace bpp
{

/**
 * @brief Partial implementation of the OSequence and OAlignment interfaces.
 */
class AbstractOSequence:
  public virtual OSequence,
  public virtual OAlignment
{

	public: 
		AbstractOSequence() {}
		virtual ~AbstractOSequence() {}

	public:

		/**
		 * @name OSequence methods:
		 *
		 * @{
		 */ 
		void writeSequences(std::ostream& output, const SequenceContainer& sc) const throw (Exception) = 0;
		void writeSequences(const std::string& path, const SequenceContainer& sc, bool overwrite=true) const throw (Exception)
		{
			// Open file in specified mode
      std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
			writeSequences(output, sc);
			output.close();
		}
		/** @} */
    
    /**
		 * @name OAlignment methods:
		 *
     * As a SiteContainer is a specialization of SequenceContainer, it is assumed that a OSequence
     * object can write aligned sequence just like a OAlignment object.
     * Therefore it implements the OAlignment interface by down-casting the SiteContainer
     * to a SequenceContainer. 
		 * @{
		 */ 
		void writeAlignment(std::ostream& output, const SiteContainer& sc) const throw (Exception)
    {
      writeSequences(output, dynamic_cast<const SequenceContainer&>(sc));
    }
		void writeAlignment(const std::string& path, const SiteContainer& sc, bool overwrite=true) const throw (Exception)
		{
      writeSequences(path, dynamic_cast<const SequenceContainer&>(sc), overwrite);
		}
		/** @} */

};

} //end of namespace bpp.

#endif //_ABSTRACTOSEQUENCE_H_

