//
// File: WordAlphabet.h
// Created by: Laurent Gueguen
// Created on: Sun Dec 28 2008
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

#ifndef _WORDALPHABET_H_
#define _WORDALPHABET_H_

#include "AbstractAlphabet.h"

// From the STL:
#include <string>
#include <vector>

#include "../Sequence.h"

namespace bpp
{
/**
 * @brief The base class for word alphabets.
 *
 * These alphabets are compounds of several alphabets. The only
 * constraint on these alphabets is that their words have length one
 * (so it is not possible to make WordAlphabets from other
 * WordAlphabets). The construction is made from a vector of pointers
 * to AbstractAlphabets.
 *
 * The strings of the WordAlphabet are concatenations of the strings
 * of the Alphabets. They are made from the resolved letters of the
 * Alphabets.
 */
class WordAlphabet :
  public AbstractAlphabet
{
protected:
  std::vector<const Alphabet* > vAbsAlph_;

public:
  // Constructor and destructor.
  /**
   * @brief Builds a new word alphabet from a vector of Alphabets.
   *
   * The unit alphabets are not owned by the world alphabet, and won't
   * be destroyed when this instance is destroyed.
   *
   * @param vAlpha The vector of Alphabets to be used.
   */
  WordAlphabet(const std::vector<const Alphabet*>& vAlpha);

  /**
   * @brief Builds a new word alphabet from a pointer to number of
   * Alphabets.
   *
   * @param pAlpha The Pointer to the Alphabet to be used.
   * @param num the length of the words.
   */
  WordAlphabet(const Alphabet* pAlpha, unsigned int num);

  virtual ~WordAlphabet() {}

public:
  /**
   * @name Methods redefined from Alphabet
   *
   * @{
   */
  /**
   * @brief Get the complete name of a state given its string description.
   *
   * In case of undefined characters (i.e. N and X for nucleic alphabets),
   * this method will return the name of the undefined word.
   *
   * @param state The string description of the given state.
   * @return The name of the state.
   * @throw BadCharException When state is not a valid char description.
   */
  std::string getName(const std::string& state) const throw (BadCharException);

  int charToInt(const std::string& state) const throw (BadCharException)
  {
    if (state.size() != vAbsAlph_.size())
      throw BadCharException(state, "WordAlphabet::charToInt", this);
    if (containsUnresolved(state))
      return getSize();
    if (containsGap(state))
      return -1;
    else return AbstractAlphabet::charToInt(state);
  }

  unsigned int getSize() const
  {
    return getNumberOfChars() - 2;
  }

  /** @} */

  /**
   * @brief Returns True if the Alphabet of the letters in the word
   * are the same type.
   *
   */
  bool hasUniqueAlphabet() const;

  /**
   * @brief Returns the length of the word
   *
   */
  unsigned int getLength() const
  {
    return static_cast<unsigned int>(vAbsAlph_.size());
  }


  /**
   * @brief Returns the number of resolved states + one for unresolved
   *
   */
  unsigned int getNumberOfTypes() const
  {
    return getNumberOfChars() - 1;
  }

  std::string getAlphabetType() const;
  int getUnknownCharacterCode() const
  {
    return getSize();
  }

  bool isUnresolved(int state) const { return state == getUnknownCharacterCode(); }
  bool isUnresolved(const std::string& state) const { return charToInt(state) == getUnknownCharacterCode(); }

  std::vector<int> getAlias(int state) const throw (BadIntException);
  std::vector<std::string> getAlias(const std::string& state) const throw (BadCharException);
  int getGeneric(const std::vector<int>& states) const throw (BadIntException);
  std::string getGeneric(const std::vector<std::string>& states) const throw (BadCharException);

private:
  /**
   * @name Inner utilitary functions
   *
   * @{
   */
  bool containsUnresolved(const std::string& state) const throw (BadCharException);
  bool containsGap(const std::string& state) const throw (BadCharException);
  void build_();
  /** @} */

public:
  /**
   * @name Word specific methods
   *
   * @{
   */

  /**
   * @brief Get the pointer to the Alphabet  of the n-position.
   *
   * @param n The position in the word (starting at 0).
   * @return The pointer to the Alphabet of the n-position.
   */
  const Alphabet* getNAlphabet(size_t n) const
  {
    if (n >= vAbsAlph_.size())
      throw IndexOutOfBoundsException("WordAlphabet::getNPosition", n, 0, vAbsAlph_.size());

    return vAbsAlph_[n];
  }

  /**
   * @brief Get the int code for a word given the int code of the underlying positions.
   *
   * The int code of each position must match the corresponding alphabet specified at this position.
   * @param vint description for all the positions.
   * @param pos the start position to match in the vector.
   * @return The int code of the word.
   * @throw IndexOutOfBoundsException In case of wrong position.
   */
  virtual int getWord(const std::vector<int>& vint, size_t pos = 0) const throw (IndexOutOfBoundsException);

  /**
   * @brief Get the char code for a word given the char code of the
   * underlying positions.
   *
   * The char code of each position must match the corresponding alphabet specified at this position.
   * @param vpos vector description for all the positions.
   * @param pos the start position to match in the vector.
   * @return The string of the word.
   * @throw IndexOutOfBoundsException In case of wrong position.
   */
  virtual std::string getWord(const std::vector<std::string>& vpos, size_t pos = 0) const throw (IndexOutOfBoundsException, BadCharException);

  /**
   * @brief Get the int code of the n-position of a word given its int description.
   *
   * @param word The int description of the word.
   * @param n The position in the word (starting at 0).
   * @return The int description of the n-position of the word.
   */
  int getNPosition(int word, size_t n) const throw (BadIntException)
  {
    if (n >= vAbsAlph_.size())
      throw IndexOutOfBoundsException("WordAlphabet::getNPosition", n, 0, vAbsAlph_.size());

    std::string s = intToChar(word);
    return vAbsAlph_[n]->charToInt(s.substr(n, 1));
  }

  /**
   * @brief Get the int codes of each position of a word given its int description.
   *
   * @param word The int description of the word.
   * @return The int description of the positions of the codon.
   */

  std::vector<int> getPositions(int word) const throw (BadIntException)
  {
    std::string s = intToChar(word);
    std::vector<int> positions;
    for (size_t i = 0; i < s.size(); i++)
    {
      positions.push_back(vAbsAlph_[i]->charToInt(s.substr(i, 1)));
    }

    return positions;
  }
  /**
   * @brief Get the char code of the n-position of a word given its char description.
   *
   * @param word The char description of the word.
   * @param n The position in the word (starting at 0).
   * @return The char description of the n-position of the word.
   */
  std::string getNPosition (const std::string& word, size_t n) const throw (BadCharException)
  {
    if (n > vAbsAlph_.size())
      throw BadCharException("", "WordAlphabet::getNPosition", this);
    // Test:
    charToInt(word);

    return "" + word.substr(n, 1);
  }


  /**
   * @brief Get the char codes of each position of a word given its char description.
   *
   * @param word The char description of the word.
   * @return The char description of the three positions of the word.
   */

  std::vector<std::string> getPositions(const std::string& word) const throw (BadCharException)
  {
    charToInt(word);
    std::vector<std::string> positions;
    for (size_t i = 0; i < word.size(); i++)
    {
      positions.push_back(word.substr(i, 1));
    }

    return positions;
  }

  /**
   * @brief Translate a whole sequence from letters alphabet to words alphabet.
   *
   * @param sequence A sequence in letters alphabet.
   * @param pos the start postion (default 0)
   * @return The corresponding sequence in words alphabet.
   * @throw AlphabetMismatchException If the sequence alphabet do not match the source alphabet.
   * @throw Exception                 Other kind of error, depending on the implementation.
   */
  Sequence* translate(const Sequence &sequence, size_t = 0) const throw (AlphabetMismatchException, Exception);

  /**
   * @brief Translate a whole sequence from words alphabet to letters alphabet.
   *
   * @param sequence A sequence in words alphabet.
   * @return The corresponding sequence in letters alphabet.
   * @throw AlphabetMismatchException If the sequence alphabet do not match the target alphabet.
   * @throw Exception                 Other kind of error, depending on the implementation.
   */
  Sequence* reverse(const Sequence& sequence) const throw (AlphabetMismatchException, Exception);

  /** @} */

  /**
   * @name Overloaded AbstractAlphabet methods.
   * @{
   */
  unsigned int getStateCodingSize() const { return static_cast<unsigned int>(vAbsAlph_.size()); }
  /** @} */
};
} // end of namespace bpp.

#endif  // _WORDALPHABET_H_

