//
// File: StringSequenceTools.h
// Created by: Julien Dutheil
// Created on: Sun Nov 30 11:29:07 2003
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

#ifndef _STRINGSEQUENCETOOLS_H_
#define _STRINGSEQUENCETOOLS_H_

#include "Alphabet/Alphabet.h"
#include "Alphabet/AlphabetExceptions.h"
#include "SequenceExceptions.h"
#include <Bpp/Exceptions.h>

// From the STL:
#include <string>
#include <vector>

namespace bpp
{

/**
 * @brief Utilitary methods working on raw string objects.
 *
 * Sequences may be stored as strings, but this approach is not as
 * powerful as using true sequence objects.
 * Consider using the Sequence and SequenceTools classes for more methods.
 *
 * Some of the methods implemented here are internally used by the Sequence object.
 *
 * @see Sequence, Site, SequenceTools, SiteTools
 */
class StringSequenceTools
{
	public:
    StringSequenceTools() {};
    ~StringSequenceTools() {};
	
	public:
    
		/**
     * @brief Get a subsequence.
     *
     * @param sequence The input sequence.
     * @param begin    The begining position (included).
     * @param end      The ending position (included).
     * @return A string with the subsequence.
     * @throw Exception If position does not not match the interval [0, length].
     */
		static std::string subseq(const std::string& sequence, int begin, int end) throw (Exception);

		/**
		 * @brief Set up the size of a sequence from the right side.
		 *
		 * All new characters are filled with gaps.
		 * If the specified size is < to the sequence size, the sequence will be truncated.
		 *
     * @param sequence The input sequence.
		 * @param size The new size of the sequence.
		 */
		static std::string setToSizeR(const std::string& sequence, int size);

    /**
		 * @brief Set up the size of a sequence from the left side.
		 *
		 * All new characters are filled with gaps.
		 * If the specified size is < to the sequence size, the sequence will be truncated.
		 *
     * @param sequence The input sequence.
		 * @param size The new size of the sequence.
		 */
		static std::string setToSizeL(const std::string& sequence, int size);

		/**
     * @brief Delete all occurence of a character in the sequence.
     *
     * @param sequence The sequence to parse.
     * @param chars    The character to remove.
     * @return         The sequence with all specified characters removed.
     */
		static std::string deleteChar(const std::string& sequence, char chars);

    /**
     * @brief Delete all occurence of several characters in the sequence.
     *
     * @param sequence The sequence to parse.
     * @param chars    The characters to remove.
     * @return         The sequence with all specified characters removed.
     */
		static std::string deleteChar(const std::string& sequence, std::string chars);

    /**
     * @brief Reverse the sequence.
     *
     * @param sequence The sequence to reverse.
     * @return The reversed sequence.
     */
    static std::string* reverse(const std::string& sequence);

    /**
     * @brief Get the complement of a sequence.
		 * @deprecated Consider working with sequence objects and translators.
     *
     * For this method, sequence is supposed to be of type DNA.
     *
     * @param sequence The sequence to complement.
     * @return The complementary sequence.
		 * @see Sequence
		 * @see NucleicAcidsReplication
     */
		static std::string* complement(const std::string& sequence);

    /**
     * @brief Calculate the local GC content of a sequence.
     *
		 * GC contents are calculated using a window of specified size around the given position.
     * Note : Calculus for last nucleotides (sequence's size - window size) will return
		 * the last possible rate calculated.
     *
     * @param sequence The sequence to analyse.
     * @param pos      The position where to compute the GC content.
     * @param window   The size of the window to use.
     * @return The GC content as a ratio (# of GC / window).
     * @throw BadIntException If the sequence is not of type DNA or RNA.
     * @throw Exception       Error in calculus (if the sequence contain gaps for instance).
     */
		static double getGCcontent(const std::string& sequence, size_t pos, size_t window) throw (BadIntegerException, Exception);

    /**
     * @brief Convert a string sequence to a vector of int.
     *
		 * This method is used in the sequence constructor.
		 * This deals with the most simple cases: the sequence (string) contains
		 * a succession of all characters. This is indeed the case for simple
		 * alphabets, but may be more complicated if the alphabet is coded with
		 * variable code length.
     *
     * @param sequence The sequence to parse.
     * @param alphabet The alphabet to use to code the sequence.
     * @return A vector of int codes.
     * @throw BarCharException If some character does not match the specified alphabet.
     */
		static std::vector<int> codeSequence(const std::string& sequence, const Alphabet* alphabet) throw (BadCharException);

    /**
     * @brief Convert a sequence to its string representation.
     *
     * @param sequence The sequence object to convert.
     * @param alphabet The alphabet to use to decode the sequence.
     * @return A string representation of the sequence.
     * @throw BarIntException If some value does not match the specified alphabet.
     */
    static std::string decodeSequence(const std::vector<int>& sequence, const Alphabet* alphabet) throw (BadIntException);

		/**
		 * @brief Parse a sequence and try to guess the correct alphabet to use.
     *
		 * @param sequence The sequence to parse.
     * @return A pointer toward a new Alphabet object.
     * @throw EmptySequenceException if the sequence is empty.
     * @throw SequenceException if their is an ambiguity between several alphabet.
     * @throw AlphabetException if the sequence does not match any alphabet.
		 */
		static Alphabet* getAlphabetFromSequence(const std::string& sequence) throw (EmptySequenceException, SequenceException, AlphabetException);

};

} //end of namespace bpp.

#endif	//_STRINGSEQUENCETOOLS_H

