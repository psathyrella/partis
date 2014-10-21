//
// File: AAChouFasmanTurnIndex.h
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

#ifndef _AACHUFASMANTURNINDEX_H_
#define _AACHUFASMANTURNINDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetTools.h"

namespace bpp
{
/**
 * @brief Turn score for the Chou-Fasman algorithm of secondary structure prediction, according to http://prowl.rockefeller.edu/aainfo/chou.htm
 *
 *
 */
class AAChouFasmanTurnIndex :
  public AlphabetIndex1
{
private:
  std::vector<double> turn_;

public:
  AAChouFasmanTurnIndex() :
    turn_()
  {
    turn_.resize(20);
    turn_[ 0] =  66; // A
    turn_[ 1] =  95; // R
    turn_[ 2] =  156; // N
    turn_[ 3] =  146; // D
    turn_[ 4] =  119; // C
    turn_[ 5] =  98; // Q
    turn_[ 6] =  74; // E
    turn_[ 7] =  156; // G
    turn_[ 8] =  95; // H
    turn_[ 9] =  47; // I
    turn_[10] =  59; // L
    turn_[11] =  101; // K
    turn_[12] =  60; // M
    turn_[13] =  60; // F
    turn_[14] =  152; // P
    turn_[15] =  143; // S
    turn_[16] =  96; // T
    turn_[17] =  96; // W
    turn_[18] =  114; // Y
    turn_[19] =  50; // V
  }

  virtual ~AAChouFasmanTurnIndex() {}

  AAChouFasmanTurnIndex* clone() const { return new AAChouFasmanTurnIndex(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "AAChouFasmanTurnIndex::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return turn_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return turn_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(turn_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _AACHUFASMANTURNINDEX_H_


