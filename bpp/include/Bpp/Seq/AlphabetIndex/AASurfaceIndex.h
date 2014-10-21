//
// File: AASurfaceIndex.h
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

#ifndef _AASURFACEINDEX_H_
#define _AASURFACEINDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetTools.h"

namespace bpp
{
/**
 * @brief Surface (Angström^2) of each amino acid, according to http://www.imb-jena.de/IMAGE_AA.html
 *
 *
 */
class AASurfaceIndex :
  public AlphabetIndex1
{
private:
  std::vector<double> surface_;

public:
  AASurfaceIndex() :
    surface_()
  {
    surface_.resize(20);
    surface_[ 0] =  115; // A
    surface_[ 1] =  225; // R
    surface_[ 2] =  160; // N
    surface_[ 3] =  150; // D
    surface_[ 4] =  135; // C
    surface_[ 5] =  180; // Q
    surface_[ 6] =  190; // E
    surface_[ 7] =  75; // G
    surface_[ 8] =  195; // H
    surface_[ 9] =  175; // I
    surface_[10] =  170; // L
    surface_[11] =  200; // K
    surface_[12] =  185; // M
    surface_[13] =  210; // F
    surface_[14] =  145; // P
    surface_[15] =  115; // S
    surface_[16] =  140; // T
    surface_[17] =  255; // W
    surface_[18] =  230; // Y
    surface_[19] =  155; // V
  }

  virtual ~AASurfaceIndex() {}

  AASurfaceIndex* clone() const { return new AASurfaceIndex(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "AASurfaceIndex::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return surface_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return surface_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(surface_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _AASURFACEINDEX_H_


