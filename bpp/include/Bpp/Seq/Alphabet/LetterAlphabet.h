// 
// File:    LetterAlphabet.h
// Author:  Sylvain Gaillard
// Created: 11/09/2009 14:31:05
// 

/*
Copyright or © or Copr. Bio++ Development Team, (September 11, 2009)

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

#ifndef _LETTERALPHABET_
#define _LETTERALPHABET_

// From the STL
#include <string>
#include <vector>
#include <iostream>

// From Seq
#include "AbstractAlphabet.h"

namespace bpp {
  /**
   * @brief Specialized partial implementation of Alphabet using single letters.
   *
   * @author Sylvain Gaillard
   */
  class LetterAlphabet:
    public AbstractAlphabet
  {
    private:
      static const int LETTER_UNDEF_VALUE;
      std::vector<int> letters_;
      bool caseSensitive_;

    public:
      LetterAlphabet(bool caseSensitive = false): letters_(256, LETTER_UNDEF_VALUE), caseSensitive_(caseSensitive) {}
      virtual ~LetterAlphabet() {}

    public:
      bool isCharInAlphabet(char state) const {
        return letters_[static_cast<unsigned int>(state)] != LETTER_UNDEF_VALUE;
      }
      bool isCharInAlphabet(const std::string& state) const {
        return isCharInAlphabet(state[0]);
      }
      int charToInt(const std::string &state) const throw (BadCharException) {
        if (!isCharInAlphabet(state))
          throw BadCharException(state, "LetterAlphabet::charToInt: Unknown state", this);
        return letters_[static_cast<unsigned int>(state[0])];
      }

    protected:
      void registerState(const AlphabetState& st) {
        AbstractAlphabet::registerState(st);
        if (caseSensitive_) {
          letters_[static_cast<unsigned int>(st.getLetter()[0])] = st.getNum();
        } else {
          letters_[static_cast<unsigned int>(tolower(st.getLetter()[0]))] = st.getNum();
          letters_[static_cast<unsigned int>(toupper(st.getLetter()[0]))] = st.getNum();
        }
      }

      void setState(size_t pos, const AlphabetState& st) throw (IndexOutOfBoundsException) {
        AbstractAlphabet::setState(pos, st);
        if (caseSensitive_) {
          letters_[static_cast<unsigned int>(st.getLetter()[0])] = st.getNum();
        } else {
          letters_[static_cast<unsigned int>(tolower(st.getLetter()[0]))] = st.getNum();
          letters_[static_cast<unsigned int>(toupper(st.getLetter()[0]))] = st.getNum();
        }
      }

  };
}

#endif // _LETTERALPHABET_
