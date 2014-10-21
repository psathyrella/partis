//
// File: SimpleScore.h
// Created by: Julien Dutheil
// Created on: Fri May 04 09:35 2007
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

#ifndef _SIMPLESCORE_H_
#define _SIMPLESCORE_H_

// from the STL:
#include <string>

#include "AlphabetIndex2.h"
#include "../Alphabet/Alphabet.h"
#include "../Alphabet/AlphabetExceptions.h"
#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{
/**
 * @brief Simple Substitution Matrix, with match and mismatch penalties.
 */
class SimpleScore :
  public virtual AlphabetIndex2
{
private:
  LinearMatrix<double> distanceMatrix_;
  const Alphabet* alphabet_;

public:
  /**
   * @brief Build a new simpleScore object.
   *
   * @param alphabet The alphabet to use.
   * @param match Matching score.
   * @param mismatch Mismatching penalty.
   */
  SimpleScore(const Alphabet* alphabet, double match, double mismatch);

  SimpleScore(const SimpleScore& sc) :
    distanceMatrix_(sc.distanceMatrix_),
    alphabet_(sc.alphabet_)
  {}

  SimpleScore& operator=(const SimpleScore& sc)
  {
    distanceMatrix_ = sc.distanceMatrix_;
    alphabet_ = sc.alphabet_;
    return *this;
  }

  virtual ~SimpleScore() {}

  SimpleScore* clone() const { return new SimpleScore(*this); }

public:
  /**
   * @name Methods from the AlphabetIndex2 interface.
   *
   * @{
   */
  double getIndex(int state1, int state2) const throw (BadIntException);
  double getIndex(const std::string& state1, const std::string& state2) const throw (BadCharException);
  const Alphabet* getAlphabet() const { return alphabet_; }
  LinearMatrix<double>* getIndexMatrix() const;
  /** @} */
};
} // end of namespace bpp.

#endif // _SIMPLESCORE_H_

