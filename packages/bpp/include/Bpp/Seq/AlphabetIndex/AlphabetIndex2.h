//
// File: AlphabetIndex2.h
// Created by: Julien Dutheil
// Created on: Mon Feb 21 17:42 2005
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

#ifndef _ALPHABETINDEX2_H_
#define _ALPHABETINDEX2_H_

#include "../Alphabet/Alphabet.h"
#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

// From the STL:
#include <string>

namespace bpp
{
/**
 * @brief Two dimensionnal alphabet index interface.
 *
 * Derivatives of this interface implement distances between two states.
 */
class AlphabetIndex2 : public Clonable
{
public:
  AlphabetIndex2() {}
  virtual ~AlphabetIndex2() {}

public:
  virtual AlphabetIndex2* clone() const = 0;

  /**
   * @brief Get the index associated to a pair of states.
   *
   * @param state1 First state to consider, as a int value.
   * @param state2 Second state to consider, as a int value.
   * @return The index associated to the pair of states
   */
  virtual double getIndex(int state1, int state2) const = 0;

  /**
   * @brief Get the index associated to a pair of states.
   *
   * @param state1 First state to consider, as a string value.
   * @param state2 Second state to consider, as a string value.
   * @return The index associated to the pair of states
   */
  virtual double getIndex(const std::string& state1, const std::string& state2) const = 0;

  /**
   * @brief Get the alphabet associated to this index.
   *
   * @return Alphabet The alphabet associated to this index.
   */
  virtual const Alphabet* getAlphabet() const = 0;

  /**
   * @return A matrix object with all indices.
   */
  virtual Matrix<double>* getIndexMatrix() const = 0;
};
} // end of namespace bpp.

#endif // _ALPHABETINDEX2_H_

