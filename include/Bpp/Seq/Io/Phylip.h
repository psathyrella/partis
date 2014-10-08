//
// File: Phylip.h
// Created by: Julien Dutheil
// Created on: Mon Oct 27 12:22:56 2003
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

#ifndef _PHYLIP_H_
#define _PHYLIP_H_

#include "AbstractIAlignment.h"
#include "AbstractOAlignment.h"
#include "../Sequence.h"
#include "../Container/SequenceContainer.h"
#include "../Container/VectorSequenceContainer.h"
#include "../Container/AlignedSequenceContainer.h"

// From the STL:
#include <iostream>

namespace bpp
{

/**
 * @brief The Phylip & co format.
 *
 * An AlignedSequenceContainer is used instead of a VectorSequenceContainer.
 *
 * This format is described on the Phylip package documentation website:
 * http://evolution.genetics.washington.edu/phylip/doc/sequence.html
 */
class Phylip :
  public AbstractIAlignment,
  public AbstractOAlignment,
  public virtual ISequence
{
  private:

    /* this class allows two kinds of Phylip format:
     * traditional, with names limited to 10 chars,
     * and 'extended', defined by PAML, with names separated from sequences by at least 6 white spaces.
     */
    bool extended_;
    /* tells if sequences are in the seuqential or the interleave format/
     */
    bool sequential_;

    /**
     * @brief The maximum number of chars to be written on a line.
     */
    unsigned int charsByLine_;

    bool checkNames_;

    std::string namesSplit_;
  
  public:
    /**
     * @brief Build a new Phylip file reader.
     *
     * @param extended If true, sequences with names longer than 10 characters are allowed.
     * @param sequential If false, sequences are supposed to be interlaved.
     * @param charsByLine The number of base to display in a row.
     * @param checkSequenceNames Tell if the names in the file should be checked for unicity (slower, in o(n*n) where n is the number of sequences).
     * @param split The string to use to split sequence name from content (only for 'extended' format). This will typically be "  " (two spaces) or "\t" (a tabulation).
     */
    Phylip(bool extended = true, bool sequential = true, unsigned int charsByLine = 100, bool checkSequenceNames = true, const std::string& split = "  "):
      extended_(extended), sequential_(sequential), charsByLine_(charsByLine), checkNames_(checkSequenceNames), namesSplit_(split) {}

    virtual ~Phylip() {}

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
     * @return The number of sequences contained in the specified file.
     *
     * This methods parses the firt line of the phylip file.
     * @param path The path of the file to parse.
     */
    unsigned int getNumberOfSequences(const std::string& path) const throw (IOException);

    /**
     * @name The OSequence interface.
     *
     * @{
     */
    void writeAlignment(std::ostream& output, const SiteContainer& sc) const throw (Exception);
    void writeAlignment(const std::string& path, const SiteContainer& sc, bool overwrite) const throw (Exception)
    {
      AbstractOAlignment::writeAlignment(path, sc, overwrite);
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

    /**
     * @return The string used to split sequence name from content.
     */
    const std::string& getSplit() const { return namesSplit_; }

    /**
     * @param split The string to be used to split sequence name from content.
     */
    void setSplit(const std::string& split) { namesSplit_ = split; }
     
  protected:
    //Reading tools:
    const std::vector<std::string> splitNameAndSequence(const std::string& s) const throw (Exception); 
    void readSequential (std::istream& in, SiteContainer& asc) const throw (Exception);
    void readInterleaved(std::istream& in, SiteContainer& asc) const throw (Exception);
    //Writing tools:
    std::vector<std::string> getSizedNames(const std::vector<std::string>& names) const;
    void writeSequential (std::ostream& out, const SequenceContainer& sc, int charsByLine) const;
    void writeInterleaved(std::ostream& out, const SequenceContainer& sc, int charsByLine) const;
};

} //end of namespace bpp.

#endif  //_PHYLIP_H_

