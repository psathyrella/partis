//
// File: DefaultNucleotideScore.h
// Created by: Julien Dutheil
// Created on: Fri Jan 19 10:30 2007
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

#ifndef _DEFAULTNUCLEOTIDESCORE_H_
#define _DEFAULTNUCLEOTIDESCORE_H_

// from the STL:
#include <string>

#include "AlphabetIndex2.h"
#include "../Alphabet/NucleicAlphabet.h"
#include "../Alphabet/AlphabetExceptions.h"
#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{
/**
 * @brief Default Substitution Matrix for nucleotide alignments.
 */
class DefaultNucleotideScore :
  public virtual AlphabetIndex2
{
private:
  LinearMatrix<double> distanceMatrix_;
  const NucleicAlphabet* alpha_;

public:
  /**
   * @brief Build a new DefaultNucleotideScore object.
   *
   * @param alphabet The alphabet to use.
   */
  DefaultNucleotideScore(const NucleicAlphabet* alphabet);

  DefaultNucleotideScore(const DefaultNucleotideScore& dns) :
    distanceMatrix_(dns.distanceMatrix_),
    alpha_(dns.alpha_) {}

  DefaultNucleotideScore& operator=(const DefaultNucleotideScore& dns)
  {
    distanceMatrix_ = dns.distanceMatrix_;
    alpha_ = dns.alpha_;
    return *this;
  }

  virtual ~DefaultNucleotideScore() {}

public:
  /**
   * @name Methods from the AlphabetIndex2 interface.
   *
   * @{
   */
  /**
   * @copydoc bpp::AlphabetIndex2::getIndex()
   *
   * If states are unresolved, takes the best score of all possible matches
   * and divides it by the number of different states.
   * @author Sylvain Gaillard
   */
  double getIndex(int state1, int state2) const throw (BadIntException);
  double getIndex(const std::string& state1, const std::string& state2) const throw (BadCharException);
  const Alphabet* getAlphabet() const { return alpha_; }
  DefaultNucleotideScore* clone() const { return new DefaultNucleotideScore(*this); }
  LinearMatrix<double>* getIndexMatrix() const;
  /** @} */
};
} // end of namespace bpp.

#endif // _DEFAULTNUCLEOTIDESCORE_H_

