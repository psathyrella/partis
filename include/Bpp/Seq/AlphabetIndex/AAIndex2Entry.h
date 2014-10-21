//
// File: AAIndex2Entry.h
// Created by: Julien Dutheil
// Created on: Fri Jan 19 17:07 2007
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

#ifndef _AAINDEX2ENTRY_H_
#define _AAINDEX2ENTRY_H_

#include "AlphabetIndex2.h"
#include "../Alphabet/ProteicAlphabet.h"
#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{
/**
 * @brief Create a AlphabetIndex2 object from an AAIndex2 entry.
 */
class AAIndex2Entry :
  public virtual AlphabetIndex2
{
private:
  LinearMatrix<double> property_;
  const ProteicAlphabet* alpha_;

public:
  /**
   * @brief Create a new AAIndex2Entry from an input stream.
   *
   * @param input The input stream to use.
   * @param sym Tell if the matrix is symmetric.
   * This option as an effect only if the matrix is specified as a triangle in the entry.
   * If sym==true, the oher triangle will be built by symmetry.
   * If sym==false, the other triangle will be set to (-) the given triangle.
   * @throw IOException if the stream content does not follow the AAIndex2 database entry format.
   */
  AAIndex2Entry(std::istream& input, bool sym = true) throw (IOException);

  AAIndex2Entry(const AAIndex2Entry& index) :
    property_(index.property_),
    alpha_(index.alpha_)
  {}

  AAIndex2Entry& operator=(const AAIndex2Entry& index)
  {
    property_ = index.property_;
    alpha_ = index.alpha_;
    return *this;
  }

  virtual ~AAIndex2Entry() {}

public:
  const Alphabet* getAlphabet() const { return alpha_; }

  AAIndex2Entry* clone() const { return new AAIndex2Entry(*this); }

  double getIndex(int state1, int state2) const throw (BadIntException)
  {
    if (state1 < 0 || state1 > 19) throw BadIntException(state1, "AAIndex2Entry::getIndex(). Invalid state1.", alpha_);
    if (state2 < 0 || state2 > 19) throw BadIntException(state2, "AAIndex2Entry::getIndex(). Invalid state2.", alpha_);
    return property_(state1, state2);
  }

  double getIndex(const std::string& state1, const std::string& state2) const throw (BadCharException)
  {
    return property_(alpha_->charToInt(state1), alpha_->charToInt(state2));
  }

  LinearMatrix<double>* getIndexMatrix() const { return new LinearMatrix<double>(property_); }
};
} // end of namespace bpp.

#endif // _AAINDEX2ENTRY_H_

