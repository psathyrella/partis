//
// File: CaseMaskedAlphabet.h
// Created by: Julien Dutheil
// Created on: Sun Sep 05 2010
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#ifndef _CASEMASKEDALPHABET_H_
#define _CASEMASKEDALPHABET_H_

#include "LetterAlphabet.h"

//From the STL:
#include <string>

namespace bpp
{

/**
 * @brief Case-sensitive letter alphabet.
 *
 * This alphabet is used for parsing comodity.
 * It takes as input another LetterAlphabet and will create duplicate any aphanumerical and upper case state,
 * by creating a lower case version of the state, also named "masked" state.
 * Helper functions are provided to determine whether a given state is masked or unmasked.
 */
class CaseMaskedAlphabet :
  public LetterAlphabet
{
  public:
    const LetterAlphabet* nocaseAlphabet_;

  public:
    CaseMaskedAlphabet(const LetterAlphabet* nocaseAlphabet);
    CaseMaskedAlphabet(const CaseMaskedAlphabet& cma) :
      LetterAlphabet(cma), nocaseAlphabet_(cma.nocaseAlphabet_) {}
    CaseMaskedAlphabet& operator=(const CaseMaskedAlphabet& cma) {
      LetterAlphabet::operator=(cma);
      nocaseAlphabet_ = cma.nocaseAlphabet_;
      return *this;
    }

  public:
		unsigned int getSize() const { return nocaseAlphabet_->getSize(); }
		unsigned int getNumberOfTypes() const { return nocaseAlphabet_->getNumberOfTypes(); }
    std::string getAlphabetType() const { return "Default alphabet"; }
    int getUnknownCharacterCode() const { return nocaseAlphabet_->getUnknownCharacterCode(); }
    bool isUnresolved(int state) const { return nocaseAlphabet_->isUnresolved(state); }
    bool isUnresolved(const std::string& state) const { return nocaseAlphabet_->isUnresolved(state); }

    const Alphabet* getUnmaskedAlphabet() const { return nocaseAlphabet_; }

    bool isMasked(int state) const { return state >= 100; }
    bool isMasked(const std::string& state) const {
      char c = state.c_str()[0];
      return isMasked(c);
    }
    bool isMasked(char state) const {
      return isalpha(state) && !isupper(state);
    }

    /**
     * @brief Get the masked state equivalent to the input one.
     *
     * If the input state is masked, returns it "as is".
     * @param state The input state.
     * @throw BadIntException if the input state is not supported, or if there is no quivallent masked state.
     */
    int getMaskedEquivalentState(int state) const throw (BadIntException);
    /**
     * @brief Get the masked state equivalent to the input one.
     *
     * If the input state is masked, returns it "as is".
     * @param state The input state.
     * @throw BadCharException if the input state is not supported.
     * @throw BadIntException if there is no equivalent masked state.
     */
    const std::string getMaskedEquivalentState(const std::string& state) const throw (BadCharException, BadIntException);

};

} // end of namespace bpp

#endif //_CASEMASKEDALPHABET_H_

