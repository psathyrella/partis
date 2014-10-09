//
// File: PhredPhd.h
// Created by: Sylvain Gaillard
// Created on: Wed Nov 5 2008
//

/*
Copyright or © or Copr. CNRS, (November 5, 2008)

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

#ifndef _PHREDPHD_H_
#define _PHREDPHD_H_

#include "ISequenceStream.h"
#include "../Sequence.h"
#include "../SequenceWithQuality.h"

namespace bpp {

  /**
   * @brief The phd sequence file format from phred software.
   *
   * This class read DNA SequenceWithQuality from phd files produced by the
   * phred program from the University of Washington.
   *
   * @par Usage
   *
   * @code
   * // Creating a SequenceWithQuality object
   * DNA alpha;
   * SequenceWithQuality seq(&alpha);
   * std::vector<int> pos;
   *
   * // Create a PhredPhd parser
   * PhredPhd pp;
   *
   * // Opening the file
   * std::ifstream in("my_sequence.phd");
   *
   * // Read the sequence
   * pp.nextSequence(in, seq, pos);
   *
   * // Close the file
   * in.close();
   * @endcode
   *
   * @author Sylvain Gaillard
   */
  class PhredPhd: public ISequenceStream {
    public:

      /**
       * @brief Build a new PhredPhd object.
       */
      PhredPhd() {}

      virtual ~PhredPhd() {}

    public:
      /**
       * @name The ISequenceStream interface.
       *
       * @{
       */
      bool nextSequence(
          std::istream& input,
          Sequence& seq
          ) const throw (Exception);
      /** @} */

      /**
       * @brief Read a SequenceWithQuality from stream and store chromatographic positions
       *
       * A more complete parser that read a SequenceWithQuality and store
       * the position of each base call on the chromatogram in a vector of
       * int.
       *
       * @param input The stram to read.
       * @param seq The sequence to fill.
       * @param pos The vector of positions to fill.
       * @throw Exception IOException and Sequence related exceptions.
       */
      bool nextSequence(
          std::istream& input,
          Sequence& seq,
          std::vector<int>& pos
          ) const throw (Exception);

      /**
       * @name The IOFormat interface.
       *
       * @{
       */
      const std::string getDataType() const { return "SequenceWithQuality"; };
      const std::string getFormatName() const { return "phd file"; };
      const std::string getFormatDescription() const {
        return "Sequences following the phd format as describe in the phred documentation.";
      }
      /** @} */

    private:
      /**
       * @brief Global file parser
       *
       * @param input The stream to read
       * @param name The string to store the sequence name
       * @param sequence The string to store the sequence
       * @param qual The vector to store qualities
       * @param pos The vector to store positions
       */
      bool parseFile_(std::istream& input, std::string& name, std::string& sequence, std::vector<int>& qual, std::vector<int>& pos) const;

      /**
       * @brief Parse the DNA part of the file
       *
       * Read the DNA part until `END_DNA' or EOF.
       *
       * @param input The stream to read
       * @param sequence The string to store the sequence
       * @param qual The vector to store qualities
       * @param pos The vector to store positions
       */
      bool parseDNA_(std::istream& input, std::string& sequence, std::vector<int>& qual, std::vector<int>& pos) const;
  };
} //end of namespace bpp

#endif // _PHREDPHD_H_
