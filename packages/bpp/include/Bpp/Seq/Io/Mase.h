//
// File: Mase.h
// Created by: Guillaume Deuchst
//             Julien Dutheil
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

#ifndef _MASE_H_
#define _MASE_H_

#include "AbstractISequence.h"
#include "AbstractIAlignment.h"
#include "AbstractOSequence.h"
#include "../Sequence.h"
#include "../Container/SequenceContainer.h"
#include "../Container/VectorSequenceContainer.h"
#include <Bpp/Numeric/Range.h>
#include <Bpp/Utils/MapTools.h>

namespace bpp
{

/**
 * @brief A class to store information from the header of Mase files.
 *
 * @author Julien Dutheil
 */
class MaseHeader
{
  private:
    mutable std::map<std::string, std::string> trees_;
    mutable std::map<std::string, MultiRange<size_t> > siteSelections_;
    mutable std::map<std::string, std::vector<size_t> > sequenceSelections_;

  public:
    MaseHeader(): trees_(), siteSelections_(), sequenceSelections_() {}
    virtual ~MaseHeader() {}

  public:
    size_t getNumberOfTrees() const { return trees_.size(); }
    size_t getNumberOfSiteSelections() const { return siteSelections_.size(); }
    size_t getNumberOfSequenceSelections() const { return sequenceSelections_.size(); }

    std::vector<std::string> getTreeNames() const { return MapTools::getKeys(trees_); }
    std::vector<std::string> getSiteSelectionNames() const { return MapTools::getKeys(siteSelections_); }
    std::vector<std::string> getSequenceSelectionNames() const { return MapTools::getKeys(sequenceSelections_); }

    const std::string& getTree(const std::string& name) const throw (Exception) {
      if (trees_.find(name) != trees_.end()) {
        return trees_[name];
      } else {
        throw Exception("MaseHeader::getTree. No tree with name " + name);
      }
    }
    const MultiRange<size_t>& getSiteSelection(const std::string& name) const throw (Exception) {
      if (siteSelections_.find(name) != siteSelections_.end()) {
        return siteSelections_[name];
      } else {
        throw Exception("MaseHeader::getSiteSelection. No site selection with name " + name);
      }
    }
    const std::vector<size_t>& getSequenceSelection(const std::string& name) const throw (Exception) {
      if (sequenceSelections_.find(name) != sequenceSelections_.end()) {
        return sequenceSelections_[name];
      } else {
        throw Exception("MaseHeader::getSequenceSelection. No sequence selection with name " + name);
      }
    }

    void setTree(const std::string& name, const std::string& tree) {
      trees_[name] = tree;
    }
    void setSiteSelection(const std::string& name, const MultiRange<size_t>& ranges) {
      siteSelections_[name] = ranges;
    }
    void setSequenceSelection(const std::string& name, const std::vector<size_t>& set) {
      sequenceSelections_[name] = set;
    }

};

/**
 * @brief The mase sequence file format.
 *
 * In addition to traditional read and write method, this class offers overloaded method
 * with MaseHeader objects, dedicated to header information storage. If used, then the header
 * of the mase file will be parsed accordingly. Otherwise, the header lines will be stored
 * as general comments.
 *
 * @see MaseTools for alternative way of parsing headers.
 */
class Mase:
  public AbstractISequence,
  public AbstractIAlignment,
  public AbstractOSequence
{

  private:

    /**
     * @brief The maximum number of chars to be written on a line.
     */
    unsigned int charsByLine_;
    bool checkNames_;

  public :
    /**
     * @brief Build a new Mase object.
     *
     * @param charsByLine Number of character per line when writing files.
     * @param checkSequenceNames  Tell if the names in the file should be checked for unicity (slower, in o(n*n) where n is the number of sequences).
     */
    Mase(unsigned int charsByLine = 100, bool checkSequenceNames = true): charsByLine_(charsByLine), checkNames_(checkSequenceNames) {}

    // Class destructor
    virtual ~Mase() {}

  public:

    /**
     * @name Reading method including header:
     *
     * @{
     */
    VectorSequenceContainer* readMeta(std::istream& input, const Alphabet* alpha, MaseHeader& header) const throw (Exception)
    {
      readHeader_(input, header);
      return AbstractISequence::readSequences(input, alpha);
    }
    VectorSequenceContainer* readMeta(std::string& path, const Alphabet* alpha, MaseHeader& header) const throw (Exception)
    {
      std::ifstream input(path.c_str(), std::ios::in);
      VectorSequenceContainer* sc = readMeta(input, alpha, header);
      input.close();
      return sc;
    }
    /** @} */
    
    /**
     * @name The AbstractISequence interface.
     *
     * @{
     */
    void appendSequencesFromStream(std::istream& input, SequenceContainer& sc) const throw (Exception);
    /** @} */

    /**
     * @name The AbstractIAlignment interface.
     *
     * @{
     */
    void appendAlignmentFromStream(std::istream& input, SiteContainer& sc) const throw (Exception) {
      appendSequencesFromStream(input, sc); //This might cast an exception if sequences are not aligned! 
    }
    /** @} */


    /**
     * @name The OSequence interface.
     *
     * @{
     */
    void writeSequences(std::ostream& output, const SequenceContainer& sc) const throw (Exception);
    void writeSequences(const std::string& path, const SequenceContainer& sc, bool overwrite = true) const throw (Exception)
    {
      AbstractOSequence::writeSequences(path, sc, overwrite);
    }
    /** @} */

    /**
     * @name Writing methods including header:
     *
     * @{
     */
    void writeMeta(std::ostream& output, const SequenceContainer& sc, const MaseHeader& header) const throw (Exception)
    {
      writeHeader_(output, header);
      writeSequences(output, sc);
    }
    void writeMeta(const std::string& path, const SequenceContainer& sc, const MaseHeader& header, bool overwrite = true) const throw (Exception)
    {
			// Open file in specified mode
      std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
      writeHeader_(output, header);
			writeSequences(output, sc);
			output.close();
    }
    /** @} */

    /**
     * @name The IOSequence interface.
     *
     * @{
     */
    const std::string getFormatName() const { return "MASE file"; }

    const std::string getFormatDescription() const
    {
      return "Optional file comments (preceeded by ;;), sequence comments (preceeded by ;), sequence name, sequence";
    }
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
    void readHeader_(std::istream& input, MaseHeader& header) const throw (Exception);
    void writeHeader_(std::ostream& output, const MaseHeader& header) const;
};

} //end of namespace bpp.

#endif // _MASE_H_

