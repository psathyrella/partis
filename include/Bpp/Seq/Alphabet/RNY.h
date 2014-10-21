//
// File: RNY.h
// Created by: Laurent Gueguen
// Created on: Tue Jul 31 2007
//

/*
   Copyright or © or Copr. CNRS, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for sequences analysis.

   This software is governed by the CeCILL license under French law and
   abiding by the rules of distribution of free software. You can use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and rights to copy,
   modify and redistribute granted by the license, users are provided
   only with a limited warranty and the software's author, the holder of
   the economic rights, and the successive licensors have only limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading, using, modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean that it is complicated to manipulate, and that also
   therefore means that it is reserved for developers and experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards
   their requirements in conditions enabling the security of their
   systems and/or data to be ensured and, more generally, to use and
   operate it in the same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */


#ifndef _RNY_H_
#define _RNY_H_

#include "AbstractAlphabet.h"
#include "NucleicAlphabet.h"

#include <string>

/**
 * @brief This alphabet is used to deal quotiented triplet RNY +
 * combinations with "-".
 *
 * More explicitly, a quotiented triplet RNY is defined on the alphabet
 * 
 * @f\{C,T,R\} x \{A,C,G,T\} x \{ A,G,Y\}@f.
 *
 * Triplet "---" is the only triplet defined as gap. Other triplets
 * with "-" are defined as unresolved.
 *
 * @author Laurent Guéguen
 */

namespace bpp
{
class RNY : public AbstractAlphabet
{
private:
  const NucleicAlphabet&  nuclalph_;

public:
  RNY(const NucleicAlphabet&);
  ~RNY() {}

public:
  int charToInt(const std::string& state) const throw (BadCharException);
  std::string intToChar(int state) const throw (BadIntException);

  bool containsGap(const std::string& state) const throw (BadCharException);

  unsigned int getSize() const { return 36; }
  unsigned int getNumberOfTypes() const { return 80; }
  int getUnknownCharacterCode() const { return 350; }
  bool isGap(int state) const;
  std::vector<int   > getAlias(      int state) const throw (BadIntException);
  std::vector<std::string> getAlias(const std::string& state) const throw (BadCharException);
  bool isUnresolved(int state) const;
  bool isUnresolved(const std::string& state) const;

  std::string getAlphabetType() const { return "RNY alphabet"; }

  const NucleicAlphabet& getLetterAlphabet() const;

public:
  /**
   * @brief Get the char code for a triplet given the char code of the three underlying positions.
   *
   * The char code of each position must match the nucleic alphabet specified for this alphabet.
   * NB: This performs pos1 + pos2 + pos3 after checking for each position validity.
   * @param pos1 Char description for position 1.
   * @param pos2 Char description for position 2.
   * @param pos3 Char description for position 3.
   * @return The Char code of the triplet.
   */
  std::string getRNY(const std::string&, const std::string&, const std::string&) const throw (BadCharException);

  int getRNY(int, int, int, const Alphabet&) const throw (BadCharException);
};
}

#endif // _RNY_H_

