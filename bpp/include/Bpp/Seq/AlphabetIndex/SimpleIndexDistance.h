//
// File: SimpleIndexDistance.h
// Created by: Julien Dutheil
// Created on: Tue Apr 21 2005
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

#ifndef _SIMPLEINDEXDISTANCE_H_
#define _SIMPLEINDEXDISTANCE_H_

// from the STL:
#include <string>

#include "AlphabetIndex1.h"
#include "AlphabetIndex2.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetExceptions.h"
#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{
/**
 * @brief Simple dissimilarity distance.
 *
 * Take a one-dimensional index end return the difference between the
 * indexes of two states.
 */
class SimpleIndexDistance :
  public virtual AlphabetIndex2
{
private:
  std::auto_ptr<AlphabetIndex1> index_;
  bool sym_;

public:
  SimpleIndexDistance(AlphabetIndex1* index) :
    index_(index),
    sym_(false)
  {}

  SimpleIndexDistance(const SimpleIndexDistance& sid) :
    index_(dynamic_cast<AlphabetIndex1*>(sid.index_->clone())),
    sym_(sid.sym_)
  {}

  SimpleIndexDistance& operator=(const SimpleIndexDistance& sid)
  {
    index_.reset(dynamic_cast<AlphabetIndex1*>(sid.index_->clone()));
    sym_ = sid.sym_;
    return *this;
  }

  virtual ~SimpleIndexDistance() {}

public:
  double getIndex(int state1, int state2) const throw (BadIntException)
  {
    double d = index_->getIndex(state2) - index_->getIndex(state1);
    return sym_ ? NumTools::abs<double>(d) : d;
  }

  double getIndex(const std::string& state1, const std::string& state2) const throw (BadCharException)
  {
    double d = index_->getIndex(state2) - index_->getIndex(state1);
    return sym_ ? NumTools::abs<double>(d) : d;
  }

  const Alphabet* getAlphabet() const { return index_->getAlphabet(); }

  SimpleIndexDistance* clone() const { return new SimpleIndexDistance(*this); }

  Matrix<double>* getIndexMatrix() const
  {
    size_t n = index_->getAlphabet()->getSize(); //We should change to "supported ints" there...
    RowMatrix<double>* m = new RowMatrix<double>(n, n);
    for (int i = 0; i < static_cast<int>(n); i++)
    {
      for (int j = 0; j < static_cast<int>(n); j++)
      {
        (*m)(i, j) = getIndex(i, j);
      }
    }
    return m;
  }

public:
  void setSymmetric(bool yn) { sym_ = yn; }
  bool isSymmetric() const { return sym_; }
  /**
   * @return The AlphabetIndex1 object associated to this object.
   */
  const AlphabetIndex1& getAlphabetIndex1() const { return *index_; }
  /**
   * @return The AlphabetIndex1 object associated to this object.
   */
  AlphabetIndex1& getAlphabetIndex1() { return *index_; }
};
} // end of namespace bpp.

#endif // _SIMPLEINDEXDISTANCE_H_

