//
// File: Fasta.h
// Authors: Guillaume Deuchst
//          Julien Dutheil
//          Sylvain Gaillard
// Created: Tue Aug 21 2003
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

#ifndef _BPP_SEQ_IO_FASTA_H_
#define _BPP_SEQ_IO_FASTA_H_

#include "AbstractISequence.h"
#include "AbstractIAlignment.h"
#include "AbstractOSequence.h"
#include "../Sequence.h"
#include "../Container/SequenceContainer.h"
#include "../Container/VectorSequenceContainer.h"
#include "ISequenceStream.h"
#include "OSequenceStream.h"
#include "SequenceFileIndex.h"

namespace bpp
{

/**
 * @brief The fasta sequence file format.
 *
 * Read and write from/to Fasta files.
 */
class Fasta:
  public AbstractISequence,
  public AbstractIAlignment,
  public AbstractOSequence,
  public virtual ISequenceStream,
  public virtual OSequenceStream
{
  private:

    /**
     * @brief The maximum number of chars to be written on a line.
     */
    unsigned int charsByLine_; // Number of char by line (output only)
    bool checkNames_;          // If names must be checked in container
    bool extended_;            // If using HUPO-PSI extensions
    bool strictNames_;         // If name is between '>' and first space

  public:
  
    /**
     * @brief Build a new Fasta object.
     *
     * @param charsByLine Number of character per line when writing files.
     * @param checkSequenceNames  Tells if the names in the file should be checked for unicity (slower, in o(n*n) where n is the number of sequences).
     * @param extended Tells if we should read general comments and sequence comments in HUPO-PSI format.
     * @param strictSequenceNames Tells if the sequence names should be restricted to the characters between '>' and the first blank one.
     */
    Fasta(unsigned int charsByLine = 100, bool checkSequenceNames = true, bool extended = false, bool strictSequenceNames = false): charsByLine_(charsByLine), checkNames_(checkSequenceNames), extended_(extended), strictNames_(strictSequenceNames) {}

    // Class destructor
    virtual ~Fasta() {}

  public:

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
      appendSequencesFromStream(input, sc); //This may raise an exception if sequences are not aligned!
    }
    /** @} */

    /**
     * @name The OSequence interface.
     *
     * @{
     */
    void writeSequences(std::ostream& output, const SequenceContainer& sc) const throw (Exception);
    
    void writeSequences(const std::string& path, const SequenceContainer& sc, bool overwrite=true) const throw (Exception)
    {
      AbstractOSequence::writeSequences(path, sc, overwrite);
    }
    /** @} */
  
    /**
     * @name The IOSequence interface.
     *
     * @{
     */
    const std::string getFormatName() const { return "FASTA file"; };
    const std::string getFormatDescription() const
    {
      return "Sequence name (preceded by >) in one line, sequence content, no comments";
    }
    /** @} */

    /**
     * @name The ISequenceStream interface.
     *
     * @{
     */
    bool nextSequence(std::istream& input, Sequence& seq) const throw (Exception);
    /** @} */

    /**
     * @name The OSequenceStream interface.
     *
     * @{
     */
    void writeSequence(std::ostream& output, const Sequence& seq) const throw (Exception);
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
     * @return true if the sequence name is restricted to be between '>' and the first space character.
     */
    bool strictNames() const { return strictNames_; }

    /**
     * @brief Tell wethed the sequence name should be restrected to the first non blank characters.
     *
     * @param yn whether the sequence names should be restrected.
     */
    void strictNames(bool yn) { strictNames_ = yn; }

    /**
     * @brief The SequenceFileIndex class for Fasta format
     * @author Sylvain Gaillard
     */
    class FileIndex: SequenceFileIndex {
      public:
        FileIndex(): index_(), fileSize_(0) {}
        ~FileIndex() {}
        void build(const std::string& path) throw (Exception);
        std::streampos getSequencePosition(const std::string& id) const throw (Exception);
        size_t getNumberOfSequences() const throw (Exception) {
          return index_.size();
        }
        /**
         * @brief Read the index from a file
         */
        void read(const std::string& path) throw (Exception);
        /**
         * @brief Write the index to a file
         */
        void write(const std::string& path) throw (Exception);
        /**
         * @brief Get a sequence given its ID
         */
        void getSequence(const std::string& seqid, Sequence& seq, const std::string& path) const;
      private:
        std::map<std::string, std::streampos> index_;
        std::streampos fileSize_;
    };
};

} //end of namespace bpp.

#endif // _BPP_SEQ_IO_FASTA_H_

