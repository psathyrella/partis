//
// File: AAChouFasmanAHelixIndex.h
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

#ifndef _AACHUFASMANAHELIXINDEX_H_
#define _AACHUFASMANAHELIXINDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetTools.h"

namespace bpp
{
/**
 * @brief A-Helix score for the Chou-Fasman algorithm of secondary structure prediction, according to http://prowl.rockefeller.edu/aainfo/chou.htm
 *
 *
 */
class AAChouFasmanAHelixIndex :
  public AlphabetIndex1
{
private:
  std::vector<double> aHelix_;

public:
  AAChouFasmanAHelixIndex() :
    aHelix_()
  {
    aHelix_.resize(20);
    aHelix_[ 0] =  142; // A
    aHelix_[ 1] =  98; // R
    aHelix_[ 2] =  67; // N
    aHelix_[ 3] =  101; // D
    aHelix_[ 4] =  70; // C
    aHelix_[ 5] =  111; // Q
    aHelix_[ 6] =  151; // E
    aHelix_[ 7] =  57; // G
    aHelix_[ 8] =  100; // H
    aHelix_[ 9] =  108; // I
    aHelix_[10] =  121; // L
    aHelix_[11] =  114; // K
    aHelix_[12] =  145; // M
    aHelix_[13] =  113; // F
    aHelix_[14] =  57; // P
    aHelix_[15] =  77; // S
    aHelix_[16] =  83; // T
    aHelix_[17] =  108; // W
    aHelix_[18] =  69; // Y
    aHelix_[19] =  106; // V
  }

  virtual ~AAChouFasmanAHelixIndex() {}

  AAChouFasmanAHelixIndex* clone() const { return new AAChouFasmanAHelixIndex(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "AAChouFasmanAHelixIndex::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return aHelix_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return aHelix_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(aHelix_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _AACHUFASMANAHELIXINDEX_H_


