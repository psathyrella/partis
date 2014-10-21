//
// File: BinaryAlphabet.h
// Author: L Gueguen
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

#ifndef _BINARYALPHABET_H_
#define _BINARYALPHABET_H_

#include "AbstractAlphabet.h"

namespace bpp
{
/**
 * @brief The BinaryAlphabet class, letters are 0 and 1.
 * @author Laurent Gueguen
 *
 */
class BinaryAlphabet :
  public AbstractAlphabet
{
protected:
  void registerState(const AlphabetState& st)
  {
   AbstractAlphabet::registerState(*(st.clone()));
  }

public:
  // class constructor
  BinaryAlphabet();

  // class destructor
  virtual ~BinaryAlphabet() {}

public:
  unsigned int getSize() const { return 2; }
  unsigned int getNumberOfTypes() const { return 2; }
  std::string getAlphabetType() const { return "Binary alphabet"; }
  int getUnknownCharacterCode() const { return 3; }
  bool isUnresolved(int state) const { return state == 3; }
  bool isUnresolved(const std::string& state) const { return false; }
};
} // end of namespace bpp.

#endif // _BINARYALPHABET_H_

