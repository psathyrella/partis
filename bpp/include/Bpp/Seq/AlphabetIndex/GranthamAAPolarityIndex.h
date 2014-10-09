//
// File: GranthamAAPolarityIndex.h
// Created by: Julien Dutheil
// Created on: Tue Apr 21 2005
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

#ifndef _GRANTHAMAAPOLARITYINDEX_H_
#define _GRANTHAMAAPOLARITYINDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"

namespace bpp
{
/**
 * @brief Polarity index used in Grantham (1974).
 *
 * @code
 * Database: AAindex1
 * Entry: GRAR740102
 *
 * H GRAR740102
 * D Polarity (Grantham, 1974)
 * R LIT:2004143b PMID:4843792
 * A Grantham, R.
 * T Amino acid difference formula to help explain protein evolution
 * J Science 185, 862-864 (1974)
 * I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
 *      8.1    10.5    11.6    13.0     5.5    10.5    12.3     9.0    10.4     5.2
 *      4.9    11.3     5.7     5.2     8.0     9.2     8.6     5.4     6.2     5.9
 * //
 * @endcode
 */
class GranthamAAPolarityIndex :
  public virtual AlphabetIndex1
{
private:
  std::vector<double> polarity_;

public:
  GranthamAAPolarityIndex() :
    polarity_()
  {
    polarity_.resize(20);
    polarity_[ 0] =  8.1; // A
    polarity_[ 1] = 10.5; // R
    polarity_[ 2] = 11.6; // N
    polarity_[ 3] = 13.0; // D
    polarity_[ 4] =  5.5; // C
    polarity_[ 5] = 10.5; // Q
    polarity_[ 6] = 12.3; // E
    polarity_[ 7] =  9.0; // G
    polarity_[ 8] = 10.4; // H
    polarity_[ 9] =  5.2; // I
    polarity_[10] =  4.9; // L
    polarity_[11] = 11.3; // K
    polarity_[12] =  5.7; // M
    polarity_[13] =  5.2; // F
    polarity_[14] =  8.0; // P
    polarity_[15] =  9.2; // S
    polarity_[16] =  8.6; // T
    polarity_[17] =  5.4; // W
    polarity_[18] =  6.2; // Y
    polarity_[19] =  5.9; // V
  }

  virtual ~GranthamAAPolarityIndex() {}

  GranthamAAPolarityIndex* clone() const { return new GranthamAAPolarityIndex(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "GranthamAAPolarityIndex::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return polarity_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return polarity_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(polarity_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _GRANTHAMAAPOLARITYINDEX_H_

