//
// File: AAChouFasmanBSheetIndex.h
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

#ifndef _AACHUFASMANBSHEETINDEX_H_
#define _AACHUFASMANBSHEETINDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetTools.h"

namespace bpp
{
/**
 * @brief B-sheet score for the Chou-Fasman algorithm of secondary structure prediction, according to http://prowl.rockefeller.edu/aainfo/chou.htm
 *
 *
 */
class AAChouFasmanBSheetIndex :
  public AlphabetIndex1
{
private:
  std::vector<double> bSheet_;

public:
  AAChouFasmanBSheetIndex() :
    bSheet_()
  {
    bSheet_.resize(20);
    bSheet_[ 0] =  83; // A
    bSheet_[ 1] =  93; // R
    bSheet_[ 2] =  89; // N
    bSheet_[ 3] =  54; // D
    bSheet_[ 4] =  119; // C
    bSheet_[ 5] =  110; // Q
    bSheet_[ 6] =  37; // E
    bSheet_[ 7] =  75; // G
    bSheet_[ 8] =  87; // H
    bSheet_[ 9] =  160; // I
    bSheet_[10] =  130; // L
    bSheet_[11] =  74; // K
    bSheet_[12] =  105; // M
    bSheet_[13] =  138; // F
    bSheet_[14] =  55; // P
    bSheet_[15] =  75; // S
    bSheet_[16] =  119; // T
    bSheet_[17] =  137; // W
    bSheet_[18] =  147; // Y
    bSheet_[19] =  170; // V
  }

  virtual ~AAChouFasmanBSheetIndex() {}

  AAChouFasmanBSheetIndex* clone() const { return new AAChouFasmanBSheetIndex(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "AAChouFasmanBSheetIndex::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return bSheet_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return bSheet_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(bSheet_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _AACHUFASMANBSHEETINDEX_H_


