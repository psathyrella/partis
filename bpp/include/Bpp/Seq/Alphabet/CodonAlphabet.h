//
// File: CodonAlphabet.h
// Created by: Julien Dutheil
// Created on: Sun Oct 12 17:41:56 2003
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#ifndef _CODONALPHABET_H_
#define _CODONALPHABET_H_

#include "WordAlphabet.h"
#include "NucleicAlphabet.h"

// From the STL:
#include <string>


namespace bpp
{

/**
 * @brief Codon alphabet class.
 * @author Laurent Guéguen, Julien Dutheil
 * 
 * A codon alphabet object is a particular case of WordAlphabet with three letters. 
 * Since codons are made of 3 nucleic bases (RNA or DNA), this class
 * has a NucleicAlphabet field used to check char description. This
 * nucleic alphabet is passed to the constructor. This class also adds
 * some methods specific to codon manipulation.
 */
class CodonAlphabet:
  public WordAlphabet
{
public: // Constructor and destructor.
		
  /**
   * @brief Builds a new codon alphabet from a nucleic alphabet.
   * 
   * @param alpha The nucleic alphabet to be used.
   */
  CodonAlphabet(const NucleicAlphabet* alpha) :
    WordAlphabet(alpha, 3) {}
  
  virtual ~CodonAlphabet() {}
  
  std::string getAlphabetType() const
  {
    return "Codon alphabet("+ vAbsAlph_[0]->getAlphabetType() + ")";
  }

  
public:
  
  /**
   * @name Codon specific methods
   *
   * @{
   */

  /**
   * @brief Get the int code for a codon given the int code of the three underlying positions.
   *
   * The int code of each position must match the nucleic alphabet specified for this alphabet.
   * @param pos1 Int description for position 1.
   * @param pos2 Int description for position 2.
   * @param pos3 Int description for position 3.
   * @return The int code of the codon.
   */
  virtual int getCodon(int pos1, int pos2, int pos3) const throw (BadIntException);
  
  /**
   * @brief Get the char code for a codon given the char code of the three underlying positions.
   *
   * The char code of each position must match the nucleic alphabet specified for this alphabet.
   * NB: This performs pos1 + pos2 + pos3 after checking for each position validity.
   * @param pos1 Char description for position 1.
   * @param pos2 Char description for position 2.
   * @param pos3 Char description for position 3.
   * @return The Char code of the codon.
   */
  virtual std::string getCodon(const std::string& pos1, const std::string& pos2, const std::string& pos3) const throw (BadCharException);
  
  /**
   * @brief Get the int code of the first position of a codon given its int description.
   * 
   * @param codon The int description of the codon.
   * @return The int description of the first position of the codon.
   */
  virtual int getFirstPosition(int codon) const throw (BadIntException);
  
  /**
   * @brief Get the int code of the second position of a codon given its int description.
   * 
   * @param codon The int description of the codon.
   * @return The int description of the second position of the codon.
   */
  virtual int getSecondPosition(int codon) const throw (BadIntException);
  
  /**
   * @brief Get the int code of the third position of a codon given its int description.
   * 
   * @param codon The int description of the codon.
   * @return The int description of the third position of the codon.
   */
  virtual int getThirdPosition(int codon) const throw (BadIntException);
  
  /**
   * @brief Get the char code of the first position of a codon given its char description.
   * 
   * @param codon The char description of the codon.
   * @return The char description of the first position of the codon.
   */
  virtual std::string getFirstPosition (const std::string& codon) const throw (BadCharException);
  
  /**
   * @brief Get the char code of the second position of a codon given its char description.
   * 
   * @param codon The char description of the codon.
   * @return The char description of the second position of the codon.
   */
  virtual std::string getSecondPosition(const std::string& codon) const throw (BadCharException);
  
  /**
   * @brief Get the char code of the third position of a codon given its char description.
   * 
   * @param codon The char description of the codon.
   * @return The char description of the third position of the codon.
   */
  virtual std::string getThirdPosition(const std::string& codon) const throw (BadCharException);
  
  /**
   * @return The nucleic alphabet associated to this codon alphabet.
   */
  virtual const NucleicAlphabet* const getNucleicAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(vAbsAlph_[0]);
  }
  
  /** @} */
};

} //end of namespace bpp.

#endif	//_CODONALPHABET_H_

