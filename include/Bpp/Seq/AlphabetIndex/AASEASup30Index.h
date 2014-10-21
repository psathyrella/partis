//
// File: AASEASup30Index.h
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

#ifndef _AASEASUP30INDEX_H_
#define _AASEASUP30INDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetTools.h"

namespace bpp
{
/**
 * @brief Percentage of amino acids having a Solvent Exposed Area above 30 Angström^2 for each type of amino acid, according to http://prowl.rockefeller.edu/aainfo/access.htm
 *
 *
 */
class AASEASup30Index :
  public AlphabetIndex1
{
private:
  std::vector<double> seaSup30_;

public:
  AASEASup30Index() :
    seaSup30_()
  {
    seaSup30_.resize(20);
    seaSup30_[ 0] =  0.48; // A
    seaSup30_[ 1] =  0.84; // R
    seaSup30_[ 2] =  0.82; // N
    seaSup30_[ 3] =  0.81; // D
    seaSup30_[ 4] =  0.32; // C
    seaSup30_[ 5] =  0.81; // Q
    seaSup30_[ 6] =  0.93; // E
    seaSup30_[ 7] =  0.51; // G
    seaSup30_[ 8] =  0.66; // H
    seaSup30_[ 9] =  0.39; // I
    seaSup30_[10] =  0.41; // L
    seaSup30_[11] =  0.93; // K
    seaSup30_[12] =  0.44; // M
    seaSup30_[13] =  0.42; // F
    seaSup30_[14] =  0.78; // P
    seaSup30_[15] =  0.70; // S
    seaSup30_[16] =  0.71; // T
    seaSup30_[17] =  0.49; // W
    seaSup30_[18] =  0.67; // Y
    seaSup30_[19] =  0.40; // V
  }

  virtual ~AASEASup30Index() {}

  AASEASup30Index* clone() const { return new AASEASup30Index(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "AASEASup30Index::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return seaSup30_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return seaSup30_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(seaSup30_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _AASEASUP30INDEX_H_


