//
// File: AASEA1030Index.h
// Created by: Bastien Boussau
// Created on: Fri Jan 14 10:31 2011
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

#ifndef _AASEA1030INDEX_H_
#define _AASEA1030INDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetTools.h"

namespace bpp
{
/**
 * @brief Percentage of amino acids having a Solvent Exposed Area between 10 and 30 Angström^2 for each type of amino acid, according to http://prowl.rockefeller.edu/aainfo/access.htm.
 *
 *
 */
class AASEA1030Index :
  public AlphabetIndex1
{
private:
  std::vector<double> sea1030_;

public:
  AASEA1030Index() :
    sea1030_()
  {
    sea1030_.resize(20);
    sea1030_[ 0] =  0.17; // A
    sea1030_[ 1] =  0.11; // R
    sea1030_[ 2] =  0.08; // N
    sea1030_[ 3] =  0.10; // D
    sea1030_[ 4] =  0.14; // C
    sea1030_[ 5] =  0.09; // Q
    sea1030_[ 6] =  0.03; // E
    sea1030_[ 7] =  0.13; // G
    sea1030_[ 8] =  0.15; // H
    sea1030_[ 9] =  0.14; // I
    sea1030_[10] =  0.10; // L
    sea1030_[11] =  0.05; // K
    sea1030_[12] =  0.36; // M
    sea1030_[13] =  0.16; // F
    sea1030_[14] =  0.09; // P
    sea1030_[15] =  0.10; // S
    sea1030_[16] =  0.13; // T
    sea1030_[17] =  0.07; // W
    sea1030_[18] =  0.13; // Y
    sea1030_[19] =  0.10; // V
  }

  virtual ~AASEA1030Index() {}

  AASEA1030Index* clone() const { return new AASEA1030Index(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "AASEA1030Index::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return sea1030_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return sea1030_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(sea1030_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _AASEA1030INDEX_H_


