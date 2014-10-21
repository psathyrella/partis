// 
// File:    NucleicAlphabetState.h
// Author:  Sylvain Gaillard
// Created: 31/07/2009
// 

/*
Copyright or © or Copr. Bio++ Development Team, (July 29, 2009)

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

#ifndef _NUCLEICALPHABETSTATE_H_
#define _NUCLEICALPHABETSTATE_H_

// From the STL
#include <string>
#include <bitset>

namespace bpp {
  /**
   * @brief This is the base class to describe states in a NucleicAlphabet.
   *
   * This class store a binary code for each state.
   * This binary code is of length 4, one bit per nucleotid.
   * For DNA alphabet, this code looks like this:
   * <table>
   * <tr><th>Letter</th><th>Binary code</th><th>int value</th></tr>
   * <tr><td>-</td><td>0000</td><td>0</td></tr>
   * <tr><td>A</td><td>0001</td><td>1</td></tr>
   * <tr><td>C</td><td>0010</td><td>2</td></tr>
   * <tr><td>G</td><td>0100</td><td>4</td></tr>
   * <tr><td>T</td><td>1000</td><td>8</td></tr>
   * <tr><td>N</td><td>1111</td><td>15</td></tr>
   * <tr><td>...</td><td></td><td></td></tr>
   * <tr><td>M</td><td>0011</td><td>3</td></tr>
   * <tr><td>W</td><td>1001</td><td>9</td></tr>
   * <tr><td>...</td><td></td><td></td></tr>
   * <tr><td>V</td><td>0111</td><td>7</td></tr>
   * </table>
   * This notation allows the use of bitwize operations like:
   *   - build a generic state from to states
   * @code
   * A | G = R  <=>  0001 | 0100 = 0101
   * @endcode
   *   - extract a state from an unresolved one by subtraction of an other state
   * @code
   * S & ~ C = G  <=>  0110 & ~ 0010 = 0100
   * @endcode
   *
   * The binary code is stored as a char because it's the smallest memory word
   * that can be allocated. A char is 8 bits long allowing the use of this
   * class with Alphabet of at least 8 resolved states (enough for known
   * nucleic alphabets!).
   *
   * @author Sylvain Gaillard
   */
  class NucleicAlphabetState: public AlphabetState {
    private:
      int binCode_;

    public:
      NucleicAlphabetState(int num, const std::string& letter, unsigned char code, const std::string& name):
        AlphabetState(num, letter, name), binCode_(code) {}

      // Class destructor
      virtual ~NucleicAlphabetState() {}

    public:
      NucleicAlphabetState* clone() const {
        return new NucleicAlphabetState(* this);
      }
      /**
       * @brief Get the state's binary representation.
       *
       * @return The state's binary representation.
       */
      int getBinaryCode() const { return binCode_; }
      /**
       * @brief Set the state's binary representation.
       *
       * @param code The state's binary representation.
       */
      void setBinaryCode(int code) { binCode_ = code; }
  };
}

#endif // _NUCLEICALPHABETSTATE_H_

