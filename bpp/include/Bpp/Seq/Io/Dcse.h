//
// File: DCSE.h
// Created by: Julien Dutheil
// Created on: Wed Mar 3 2004
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

#ifndef _DCSE_H_
#define _DCSE_H_

#include "AbstractIAlignment.h"
#include "../Sequence.h"
#include "../Container/SequenceContainer.h"
#include "../Container/AlignedSequenceContainer.h"

namespace bpp
{

/**
 * @brief Support for the Dedicated Comparative Sequence Editor format.
 *
 * Only the sequence information is retrieved.
 * All structural information is dropped for now.
 * 
 * A description of this format may be found here:
 * http://www.psb.ugent.be/rRNA/help/formats/aliformat.html
 */
class DCSE :
  public AbstractIAlignment,
  public virtual ISequence
{
    
  public: 
    DCSE() {};
    virtual ~DCSE() {};

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
};

} //end of namespace bpp.

#endif // _DCSE_H_

