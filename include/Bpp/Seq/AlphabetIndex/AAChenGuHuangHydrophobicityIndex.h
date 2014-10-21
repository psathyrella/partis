//
// File: AAChenGuHuangHydrophobicityIndex.h
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

#ifndef _AACHENGUHUANGHYDROPHOBICITYINDEX_H_
#define _AACHENGUHUANGHYDROPHOBICITYINDEX_H_

#include "AlphabetIndex1.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetTools.h"

namespace bpp
{
/**
 * @brief Hydrophobicity of each amino acid, according to Table 1 in Chen, Gu and Huang, BMC Bioinformatics 2006.
 *
 * @code
 * Database: AAChenGuHuangHydrophobicity
 * Entry: CHENH06
 *
 * H FAUJ880111
 * D Hydrophobicity (Chen et al., 2006)
 * R PMCID:PMC1780123 PMID:17217506
 * A Hang Chen, Fei Gu, and Zhengge Huang.
 * T Improved Chou-Fasman method for protein secondary structure prediction
 * J BMC Bioinformatics. 2006; 7(Suppl 4): S14.  (2006)
 * //
 * @endcode
 *
 */
class AAChenGuHuangHydrophobicityIndex :
  public AlphabetIndex1
{
private:
  std::vector<double> hydrophobicity_;

public:
  AAChenGuHuangHydrophobicityIndex() :
    hydrophobicity_()
  {
    hydrophobicity_.resize(20);
    hydrophobicity_[ 0] =  0.87; // A
    hydrophobicity_[ 1] =  0.85; // R
    hydrophobicity_[ 2] =  0.09; // N
    hydrophobicity_[ 3] =  0.66; // D
    hydrophobicity_[ 4] =  1.52; // C
    hydrophobicity_[ 5] =  0.00; // Q
    hydrophobicity_[ 6] =  0.67; // E
    hydrophobicity_[ 7] =  0.00; // G
    hydrophobicity_[ 8] =  0.87; // H
    hydrophobicity_[ 9] =  3.15; // I
    hydrophobicity_[10] =  2.17; // L
    hydrophobicity_[11] =  1.64; // K
    hydrophobicity_[12] =  1.67; // M
    hydrophobicity_[13] =  2.87; // F
    hydrophobicity_[14] =  2.77; // P
    hydrophobicity_[15] =  0.07; // S
    hydrophobicity_[16] =  0.07; // T
    hydrophobicity_[17] =  3.77; // W
    hydrophobicity_[18] =  2.76; // Y
    hydrophobicity_[19] =  1.87; // V
  }

  virtual ~AAChenGuHuangHydrophobicityIndex() {}

  AAChenGuHuangHydrophobicityIndex* clone() const { return new AAChenGuHuangHydrophobicityIndex(); }

public:
  double getIndex(int state) const throw (BadIntException)
  {
    if (state < 0 || state > 19) throw BadIntException(state, "AAChenGuHuangHydrophobicityIndex::getIndex(). Invalid state.", &AlphabetTools::PROTEIN_ALPHABET);
    return hydrophobicity_[state];
  }

  double getIndex(const std::string& state) const throw (BadCharException)
  {
    return hydrophobicity_[AlphabetTools::PROTEIN_ALPHABET.charToInt(state)];
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(hydrophobicity_); }

  const Alphabet* getAlphabet() const { return &AlphabetTools::PROTEIN_ALPHABET; }
};
} // end of namespace bpp.

#endif // _AACHENGUHUANGHYDROPHOBICITYINDEX_H_

