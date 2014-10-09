//
// File: AAVolumeIndex.h
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

#ifndef _AAVOLUMEINDEX_H_
#define _AAVOLUMEINDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetTools.h"

namespace bpp
{
/**
 * @brief Volume (Angström^3) of each amino acid, according to http://www.imb-jena.de/IMAGE_AA.html
 *
 *
 */
class AAVolumeIndex :
  public AlphabetIndex1
{
private:
  std::vector<double> volume_;

public:
  AAVolumeIndex() :
    volume_()
  {
    volume_.resize(20);
    volume_[ 0] =  115; // A
    volume_[ 1] =  225; // R
    volume_[ 2] =  160; // N
    volume_[ 3] =  150; // D
    volume_[ 4] =  135; // C
    volume_[ 5] =  180; // Q
    volume_[ 6] =  190; // E
    volume_[ 7] =  75; // G
    volume_[ 8] =  195; // H
    volume_[ 9] =  175; // I
    volume_[10] =  170; // L
    volume_[11] =  200; // K
    volume_[12] =  185; // M
    volume_[13] =  210; // F
    volume_[14] =  145; // P
    volume_[15] =  115; // S
    volume_[16] =  140; // T
    volume_[17] =  255; // W
    volume_[18] =  230; // Y
    volume_[19] =  155; // V
  }

  virtual ~AAVolumeIndex() {}

  AAVolumeIndex* clone() const { return new AAVolumeIndex(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "AAVolumeIndex::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return volume_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return volume_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(volume_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _AAVOLUMEINDEX_H_


