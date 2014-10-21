//
// File: Stockholm.h
// Authors: Julien Dutheil
// Created: Thu Apr 15 2010
//

/*
Copyright or © or Copr. Bio++ Development Team (2010)

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

#ifndef _STOCKHOLM_H_
#define _STOCKHOLM_H_

#include "AbstractOAlignment.h"
#include "../Sequence.h"
#include "../Container/SequenceContainer.h"
#include "../Container/AlignedSequenceContainer.h"

namespace bpp
{

/**
 * @brief The Stockholm alignment file format.
 *
 * Write to Stockholm files.
 * Only sequence data is read/written, annotation and secondary structures are ignored.
 */
class Stockholm:
  public AbstractOAlignment
{
  private:

    bool checkNames_;

  public:
  
    /**
     * @brief Build a new Stockholm object.
     *
     * @param checkSequenceNames Tell if the names in the file should be checked for unicity (slower, in o(n*n) where n is the number of sequences).
     */
    Stockholm(bool checkSequenceNames = true) : checkNames_(checkSequenceNames) {}

    // Class destructor
    virtual ~Stockholm() {}

  public:

    /**
     * @name The OAlignment interface.
     *
     * @{
     */
    void writeAlignment(std::ostream& output, const SiteContainer& sc) const throw (Exception);
    void writeAlignment(const std::string& path, const SiteContainer& sc, bool overwrite = true) const throw (Exception)
    {
      AbstractOAlignment::writeAlignment(path, sc, overwrite);
    }
    /** @} */
  
    /**
     * @name The IOSequence interface.
     *
     * @{
     */
    const std::string getFormatName() const { return "Stockholm file"; };
    const std::string getFormatDescription() const
    {
      return "See http://en.wikipedia.org/wiki/Stockholm_format";
    }
    /** @} */

    /**
     * @warning This is not used for now, will be when reading is implemented.
     * @return true if the names are to be checked when reading sequences from files.
     */
    bool checkNames() const { return checkNames_; }

    /**
     * @brief Tell whether the sequence names should be checked when reading from files.
     *
     * @warning This is not used for now, will be when reading is implemented.
     * @param yn whether the sequence names should be checked when reading from files.
     */
    void checkNames(bool yn) { checkNames_ = yn; }
};

} //end of namespace bpp.

#endif // _FASTA_H_

