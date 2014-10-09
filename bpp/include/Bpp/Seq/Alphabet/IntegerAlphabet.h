//
// File: BinaryAlphabet.h
// Author: L Gueguen
// Created on: vendredi 20 septembre 2013, à 23h 01
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

#ifndef _INTEGERALPHABET_H_
#define _INTEGERALPHABET_H_

#include "AbstractAlphabet.h"

namespace bpp
{
/**
 * @brief The Integer Alphabet class, letters are from 0 to a given number, MAX.
 * @author Laurent Gueguen
 *
 */
class IntegerAlphabet :
  public AbstractAlphabet
{
private:

  unsigned int MAX_;
  
protected:
  void registerState(const AlphabetState& st)
  {
   AbstractAlphabet::registerState(*(st.clone()));
  }

public:
  // class constructor
  IntegerAlphabet(unsigned int max);

  // class destructor
  virtual ~IntegerAlphabet() {}

public:
  unsigned int getSize() const { return MAX_+1; }
  unsigned int getNumberOfTypes() const { return MAX_+1; }
  std::string getAlphabetType() const { return "Integer alphabet"; }
  int getUnknownCharacterCode() const { return MAX_; }
  bool isUnresolved(int state) const { return state == (int)MAX_; }
  bool isUnresolved(const std::string& state) const { return false; }
};
} // end of namespace bpp.

#endif // _INTEGERALPHABET_H_

