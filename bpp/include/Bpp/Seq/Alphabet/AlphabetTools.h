//
// File: AlphabetTools.h
// Created by: Julien Dutheil
// Created on: Fri Oct 10 17:27:39 2003
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

#ifndef _ALPHABETTOOLS_H_
#define _ALPHABETTOOLS_H_

#include "DNA.h"
#include "RNA.h"
#include "ProteicAlphabet.h"
#include "DefaultAlphabet.h"
#include "CodonAlphabet.h"
#include "RNY.h"
#include "BinaryAlphabet.h"
#include <Bpp/Numeric/VectorTools.h>

#include <typeinfo>

namespace bpp
{
/**
 * @brief Utilitary functions dealing with alphabets.
 */
class AlphabetTools
{
public:
  static const DNA DNA_ALPHABET;
  static const RNA RNA_ALPHABET;
  static const ProteicAlphabet PROTEIN_ALPHABET;
  static const DefaultAlphabet DEFAULT_ALPHABET;

public:
  AlphabetTools() {}
  virtual ~AlphabetTools() {}

public:
  /**
   * @brief Character identification method for sequence's alphabet identification
   *
   * Return :
   * - -1 gap
   * - 1 DNA specific (no character!)
   * - 2 RNA specific (U)
   * - 3 Protein specific (characters E, F, I, L, P, Q)
   * - 4 Nucleotide specific (no character)
   * - 5 DNA or Protein specific (T)
   * - 6 RNA or Protein specific (no character)
   * - 7 Any alphabet (A, B, C, D, G, H, J, K, M, N, O, R, S, V, W, X, Y, Z, 0)
   * - 0 Unknown character
   *
   * @param state The character to test.
   * @return The type code.
   */
  static int getType(char state);

  /**
   * @brief This checks that all characters in the alphabet are coded by a string of same length.
   *
   * This method is used when states are coded by more than one character, typically: codon alphabets.
   *
   * @param alphabet The alphabet to check.
   * @return True if all text description have the same length (e.g. 3 for codon alphabet).
   */
  static bool checkAlphabetCodingSize(const Alphabet& alphabet) throw (AlphabetException);

  /**
   * @brief This checks that all characters in the alphabet are coded by a string of same length.
   *
   * This function performs the same work as the previous one, but deals with pointers
   * instead of reference. This may be more convenient since we often use pointers on alphabets.
   *
   * @param alphabet a pointer toward the alphabet to check.
   * @return True if all text description have the same length (e.g. 3 for codon alphabet).
   */
  static bool checkAlphabetCodingSize(const Alphabet* alphabet) throw (AlphabetException);

  /**
   * @brief In case that all states in the given alphabet have a string description of same length,
   * send this length; e.g. 3 for codon alphabets.
   *
   * @param alphabet The alphabet to analyse.
   * @return The common size of all text descriptionif there is one. Else throws an AlphabetException.
   */
  static unsigned int getAlphabetCodingSize(const Alphabet& alphabet) throw (AlphabetException);

  /**
   * @brief In case that all states in the given alphabet have a string description of same length,
   * send this length; e.g. 3 for codon alphabets.
   *
   * This function performs the same work as the previous one, but deals with pointers
   * instead of reference. This may be more convenient since we often use pointers on alphabets.
   *
   * @param alphabet a pointer toward the alphabet to analyse.
   * @return The common size of all text descriptionif there is one. Else throws an AlphabetException.
   */
  static unsigned int getAlphabetCodingSize(const Alphabet* alphabet) throw (AlphabetException);

  /**
   * @return True if the alphabet is an instanciation of the NucleicAlphabet class.
   * @param alphabet The alphabet to check.
   */
  static bool isNucleicAlphabet(const Alphabet* alphabet) { return alphabetInheritsFrom<NucleicAlphabet>(alphabet); }

  /**
   * @return True if the alphabet is an instanciation of the DNA class.
   * @param alphabet The alphabet to check.
   */
  static bool isDNAAlphabet(const Alphabet* alphabet) { return alphabetInheritsFrom<DNA>(alphabet); }

  /**
   * @return True if the alphabet is an instanciation of the RNA class.
   * @param alphabet The alphabet to check.
   */
  static bool isRNAAlphabet(const Alphabet* alphabet) { return alphabetInheritsFrom<RNA>(alphabet); }

  /**
   * @return True if the alphabet is an instanciation of the ProteicAlphabet class.
   * @param alphabet The alphabet to check.
   */
  static bool isProteicAlphabet(const Alphabet* alphabet) { return alphabetInheritsFrom<ProteicAlphabet>(alphabet); }

  /**
   * @return True if the alphabet is an instanciation of the Codon class.
   * @param alphabet The alphabet to check.
   */
  static bool isCodonAlphabet(const Alphabet* alphabet) { return alphabetInheritsFrom<CodonAlphabet>(alphabet); }

  /**
   * @return True if the alphabet is an instanciation of the WordAlphabet class.
   * @param alphabet The alphabet to check.
   */
  static bool isWordAlphabet(const Alphabet* alphabet) { return alphabetInheritsFrom<WordAlphabet>(alphabet); }

  /**
   * @return True if the alphabet is an instanciation of the RNY class.
   * @param alphabet The alphabet to check.
   */
  static bool isRNYAlphabet(const Alphabet* alphabet) { return alphabetInheritsFrom<RNY>(alphabet); }

  /**
   * @return True if the alphabet is an instanciation of the BinaryAlphabet class.
   * @param alphabet The alphabet to check.
   */
  static bool isBinaryAlphabet(const Alphabet* alphabet) { return alphabetInheritsFrom<BinaryAlphabet>(alphabet); }

  /**
   * @return True if the alphabet is an instanciation of the DefaultAlphabet class.
   * @param alphabet The alphabet to check.
   */
  static bool isDefaultAlphabet(const Alphabet* alphabet) { return alphabetInheritsFrom<DefaultAlphabet>(alphabet); }

  /**
   * @brief Tell if two characters match according to a given alphabet.
   *
   * Example (DNA):
   * @verbatim
   * A,T: false
   * A,A: true
   * A,N: true
   * A,Y: false
   * N,Y: true
   * N,N: true
   * @endverbatim
   *
   * @return True if the two characters are identical, or are compatible if at least one of them is a generic character.
   * @param alphabet The alphabet to use.
   * @param i First character to check.
   * @param j Secondt character to check.
   */
  static bool match(const Alphabet* alphabet, int i, int j)
  {
    std::vector<int> a = alphabet->getAlias(i);
    std::vector<int> b = alphabet->getAlias(j);
    std::vector<int> u = VectorTools::vectorIntersection(a, b);
    return u.size() > 0;
  }

private:
  template<class Y>
  static bool alphabetInheritsFrom(const Alphabet* alphabet)
  {
    try
    {
      const Y* t = dynamic_cast<const Y*>(alphabet);
      return t != 0; // Solves strange behavior in new gcc?
    }
    catch (std::exception& e)
    {
      return false;
    }
  }
};
} // end of namespace bpp.

#endif  // _ALPHABETTOOLS_H_

