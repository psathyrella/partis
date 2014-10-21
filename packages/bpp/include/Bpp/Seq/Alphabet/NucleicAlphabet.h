//
// File: NucleicAlphabet.h
// Authors: Guillaume Deuchst
//          Julien Dutheil
//          Sylvain Gaillard
// Created on: Tue Jul 22 2003
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

#ifndef _NUCLEICALPHABET_H_
#define _NUCLEICALPHABET_H_

#include "LetterAlphabet.h"
#include "NucleicAlphabetState.h"

#include <map>
#include <iostream>

namespace bpp
{

/**
 * @brief The abstract base class for nucleic alphabets.
 *
 * This class only implements a few methods, it is mainly designed for methods/classes
 * that will require to work with both RNA and DNA.
 */
class NucleicAlphabet :
  public LetterAlphabet
{
  private:
    std::map<int, unsigned int> binCodes_;
    void updateMaps_(int pos, const NucleicAlphabetState& st) {
      if (binCodes_.find(st.getBinaryCode()) == binCodes_.end())
        binCodes_[st.getBinaryCode()] = pos;
    }

	public:
    NucleicAlphabet(): binCodes_() {}

		virtual ~NucleicAlphabet() {}

  protected:
    /**
     * @name Overloaded methods from AbstractAlphabet
     * @{
     */
    void registerState(const NucleicAlphabetState& st) {
      LetterAlphabet::registerState(st);
      updateMaps_(getNumberOfChars(), st);
    }
    void setState(unsigned int pos, const NucleicAlphabetState& st) {
      LetterAlphabet::setState(pos, st);
      updateMaps_(pos, st);
    }
    const NucleicAlphabetState& getStateAt(unsigned int pos) const
      throw (IndexOutOfBoundsException) {
        return dynamic_cast<const NucleicAlphabetState&>(
            AbstractAlphabet::getStateAt(pos)
            );
      }
    NucleicAlphabetState& getStateAt(unsigned int pos)
      throw (IndexOutOfBoundsException) {
        return dynamic_cast<NucleicAlphabetState&>(
            AbstractAlphabet::getStateAt(pos)
            );
      }
    /** @} */

  public:
    /**
     * @name Overloaded methods from AbstractAlphabet
     * @{
     */
    const NucleicAlphabetState& getState(const std::string& letter) const
      throw (BadCharException) {
        return dynamic_cast<const NucleicAlphabetState&>(
            AbstractAlphabet::getState(letter)
            );
      }
    const NucleicAlphabetState& getState(int num) const
      throw (BadIntException) {
        return dynamic_cast<const NucleicAlphabetState&>(
            AbstractAlphabet::getState(num)
            );
      }
    /** @} */

    /**
     * @name Specific methods
     * @{
     */

    /**
     * @brief Get a state by its binary representation.
     * 
     * @param code The binary representation as an unsigned char.
     * @return The NucleicAlphabetState.
     * @throw BadIntException If the code is not a valide state.
     * @author Sylvain Gaillard
     */
    const NucleicAlphabetState& getStateByBinCode(int code) const
      throw (BadIntException) {
        std::map<int, unsigned int>::const_iterator it = binCodes_.find(code);
      if (it == binCodes_.end())
        throw BadIntException(code, "NucleicAlphabet::getState(unsigned char): Binary code not in alphabet", this);
      return getStateAt(it->second);
    }

    /**
     * @brief Subtract states
     *
     * Get the remaining state when subtracting one state to another.
     *
     * @code
     * int a = alpha->charToInt("A");
     * int m = alpha->charToInt("M");
     * int c = alpha->subtract(m, a);
     * 
     * cout << alpha->intToChar(c) << endl;
     *
     * // should print C because M - A = C
     * @endcode
     *
     * @param s1 the first state as an int
     * @param s2 the second state as an int
     * @throw BadIntException if one of the states is not valide.
     * @return The remaining state as an int
     * @author Sylvain Gaillard
     */
    int subtract(int s1, int s2) const throw (BadIntException) {
      return getStateByBinCode(getState(s1).getBinaryCode() & ~ getState(s2).getBinaryCode()).getNum();
    }

    /**
     * @brief Subtract states
     *
     * Get the remaining state when subtracting one state to another.
     *
     * @code
     * string a = "A";
     * string m = "M";
     * 
     * cout << alpha->subtract(m, a) << endl;
     *
     * // should print C because M - A = C
     * @endcode
     *
     * @param s1 the first state as a string
     * @param s2 the second state as a string
     * @throw BadCharException if one of the states is not valide.
     * @return The remaining state as a string
     * @author Sylvain Gaillard
     */
    std::string subtract(const std::string& s1, const std::string& s2) const throw (BadCharException) {
      return intToChar(subtract(charToInt(s1), charToInt(s2)));
    }

    /**
     * @brief Get the overlap between to states
     *
     * Get the overlapping states between two steps.
     *
     * @code
     * int m = alpha->charToInt("M");
     * int r = alpha->charToInt("R");
     * int a = alpha->getOverlap(m, r);
     *
     * cout << alpha->intToChar(a) << endl;
     *
     * // should print A because M = A/C and R = A/G
     * @endcode
     *
     * @param s1 the first state as an int
     * @param s2 the second state as an int
     * @throw BadIntException if one of the states is not valid
     * @return The overlapping state
     * @author Sylvain Gaillard
     */
    int getOverlap(int s1, int s2) const throw (BadIntException) {
      return getStateByBinCode(getState(s1).getBinaryCode() & getState(s2).getBinaryCode()).getNum();
    }

    /**
     * @brief Get the overlap between to states
     *
     * Get the overlapping states between two steps.
     *
     * @code
     * string m = "M";
     * string r = R;
     *
     * cout << alpha->getOverlap(m, r) << endl;
     *
     * // should print A because M = A/C and R = A/G
     * @endcode
     *
     * @param s1 the first state as a string
     * @param s2 the second state as a string
     * @throw BadCharException if one of the states is not valid
     * @return The overlapping state
     * @author Sylvain Gaillard
     */
    std::string getOverlap(const std::string& s1, const std::string& s2) const throw (BadCharException) {
      return intToChar(getOverlap(charToInt(s1), charToInt(s2)));
    }

    /** @} */
	
	public:
		// return 4 : A, C, G, T (or U)
		unsigned int getSize() const { return 4; }

		// return 15 : gap isn't included, generic unresolved bases (N, X, ?, O, 0) count for one
		unsigned int getNumberOfTypes() const { return 15; }
    
    int getUnknownCharacterCode() const { return 14; }

    bool isUnresolved(int state) const { return state > 3; }
    bool isUnresolved(const std::string& state) const { return charToInt(state) > 3; }

};

} //end of namespace bpp.

#endif // _NUCLEICALPHABET_H_

