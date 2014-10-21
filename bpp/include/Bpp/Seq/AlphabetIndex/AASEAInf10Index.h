//
// File: AASEAInf10Index.h
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

#ifndef _AASEAINF10INDEX_H_
#define _AASEAINF10INDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetTools.h"

namespace bpp
{
/**
 * @brief Percentage of amino acids having a Solvent Exposed Area below 10 Angström^2 for each type of amino acid, according to http://prowl.rockefeller.edu/aainfo/access.htm.
 *
 *
 */
class AASEAInf10Index :
  public AlphabetIndex1
{
private:
  std::vector<double> seaInf10_;

public:
  AASEAInf10Index() :
    seaInf10_()
  {
    seaInf10_.resize(20);
    seaInf10_[ 0] =  0.35; // A
    seaInf10_[ 1] =  0.05; // R
    seaInf10_[ 2] =  0.10; // N
    seaInf10_[ 3] =  0.09; // D
    seaInf10_[ 4] =  0.54; // C
    seaInf10_[ 5] =  0.10; // Q
    seaInf10_[ 6] =  0.04; // E
    seaInf10_[ 7] =  0.36; // G
    seaInf10_[ 8] =  0.19; // H
    seaInf10_[ 9] =  0.47; // I
    seaInf10_[10] =  0.49; // L
    seaInf10_[11] =  0.02; // K
    seaInf10_[12] =  0.20; // M
    seaInf10_[13] =  0.42; // F
    seaInf10_[14] =  0.13; // P
    seaInf10_[15] =  0.20; // S
    seaInf10_[16] =  0.16; // T
    seaInf10_[17] =  0.44; // W
    seaInf10_[18] =  0.20; // Y
    seaInf10_[19] =  0.50; // V
  }

  virtual ~AASEAInf10Index() {}

  AASEAInf10Index* clone() const { return new AASEAInf10Index(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "AASEAInf10Index::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return seaInf10_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return seaInf10_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(seaInf10_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _AASEAINF10INDEX_H_


