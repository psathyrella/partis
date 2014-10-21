//
// File: NexusIOSequence.h
// Created by: Julien Dutheil
// Created on: Wed May 27 16:15 2009
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

#ifndef _NEXUSIOSEQUENCE_H_
#define _NEXUSIOSEQUENCE_H_

#include "AbstractIAlignment.h"
#include "../Sequence.h"
#include "../Container/SequenceContainer.h"
#include "../Container/VectorSequenceContainer.h"
#include "../Container/AlignedSequenceContainer.h"

// From the STL:
#include <iostream>

namespace bpp
{

/**
 * @brief The Nexus format reader for sequences.
 *
 * An AlignedSequenceContainer is used instead of a VectorSequenceContainer.
 *
 * This reader is not supposed to be a full parser of the Nexus files,
 * but only extract the sequence data. Only a basic subset of the options
 * are and will be supported.
 *
 * This format is described in the following paper:
 * Maddison D, Swofford D, and Maddison W (1997), _Syst Biol_ 46(4):590-621
 *
 * @author Julien Dutheil
 */
class NexusIOSequence:
  public AbstractIAlignment,
  public virtual ISequence
{
  protected:

    /**
     * @brief The maximum number of chars to be written on a line.
     */
    unsigned int charsByLine_;

    bool checkNames_;

  public:
    /**
     * @brief Build a new Phylip file reader.
     *
     * @param charsByLine The number of base to display in a row (ignored for now, no writing support).
     * @param checkSequenceNames Tell if the names in the file should be checked for unicity (slower, in o(n*n) where n is the number of sequences).
     */
    NexusIOSequence(unsigned int charsByLine = 100, bool checkSequenceNames = true):
      charsByLine_(charsByLine), checkNames_(checkSequenceNames) {}

    virtual ~NexusIOSequence() {}

  public:

    /**
     * @name The AbstractIAlignment interface.
     *
     * @{
     */
    void appendAlignmentFromStream(std::istream& input, SiteContainer& sc) const throw (Exception);
    /** @} */

    /**
     * @name The ISequence interface.
     *
     * As a SiteContainer is a subclass of SequenceContainer, we hereby implement the ISequence
     * interface by downcasting the interface.
     *
     * @{
     */
    virtual SequenceContainer* readSequences(std::istream& input, const Alphabet* alpha) const throw (Exception) {
      return readAlignment(input, alpha);
    }
    virtual SequenceContainer* readSequences(const std::string& path, const Alphabet* alpha) const throw (Exception) {
      return readAlignment(path, alpha);
    }
    /** @} */

    /**
     * @name The IOSequence interface.
     *
     * @{
     */
    const std::string getFormatName() const;
    const std::string getFormatDescription() const;
    /** @} */

    /**
     * @return true if the names are to be checked when reading sequences from files.
     */
    bool checkNames() const { return checkNames_; }

    /**
     * @brief Tell whether the sequence names should be checked when reading from files.
     *
     * @param yn whether the sequence names should be checked when reading from files.
     */
    void checkNames(bool yn) { checkNames_ = yn; }

    
  private:
    //Reading tools:
    const std::vector<std::string> splitNameAndSequence_(const std::string & s) const throw (Exception); 
};

} //end of namespace bpp.

#endif  //_NEXUSIOSEQUENCE_H_

