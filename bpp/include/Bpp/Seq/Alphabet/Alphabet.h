//
// File: Alphabet.h
// Created by: Guillaume Deuchst
//             Julien Dutheil
// Created on: Tue Jul 22 2003
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

#ifndef _ALPHABET_H_
#define _ALPHABET_H_

#include <string>
#include <vector>

#include "AlphabetExceptions.h"
#include "AlphabetState.h"

/**
 * @mainpage
 *
 * @par
 * This library provides classes to store and analyse biological sequences.
 * Each position in a sequences is coded by a int. An object implementing the bpp::Alphabet
 * interface is used to make the relation between the int code and its more common character
 * representation. Support for DNA, RNA, protein and codon sequences is provided. The
 * bpp::AlphabetTools provides tools to deal with Alphabet objects.
 * The basic bpp::Sequence class contains the code sequence, a name for the sequence
 * and optionally comments. More elaborated classes can be built by inheriting this class.
 * The bpp::SequenceTools static class provides simple analysis tools, like base frequencies measures,
 * concatenation, etc. 
 * 
 * @par
 * SeqLib also provides tools to perform <i>in silico</i> molecular biology, like complementation,
 * transcription, translation, etc. All these methods are particular cases of alphabet translation, and are implemented
 * via the interface bpp::Translator. Of particular interest are the classes bpp::NucleicAcidsReplication,
 * bpp::DNAToRNA and bpp::GeneticCode + derivatives.
 *
 * @par
 * Sequence collections are stored as <em>containers</em>. The simplest container implements the
 * bpp::SequenceContainer interface, providing access to sequences by their name.
 * The bpp::OrderedSequenceContainer adds access by position in the container.
 * The simplest implementation of this interface is the bpp::VectorSequenceContainer, which stores the sequences
 * as a vector of bpp::Sequence objects (or instances inheriting from this class).
 * Input/output from various file formats is provided, including fasta (bpp::Fasta), GenBank (bpp::GenBank) and Mase (bpp::Mase).
 * Tools dealing with containers can be found in the bpp::SequenceContainerTools static class.
 *
 * @par
 * Support for alignments is provided via the bpp::SiteContainer interface, which enables site access.
 * Sites are stored as a distinct class, similar to a "vertical" sequence, called bpp::Site.
 * It shares several methods with the bpp::Sequence object, although it does not contain a name but
 * a position attribute. This attribute can be used to track the position of sites when handling
 * alignments (for instance after removing all gap-containing sites).
 * There are currently two implementations of the bpp::SiteContainer interface:
 * - bpp::AlignedSequenceContainer, inheriting from bpp::VectorSequenceContainer and adding the site access.
 *   Sequence access is hence in @f$o(1)@f$ and site access in @f$o(n)@f$, n being the total number of sites.
 * - bpp::VectorSiteContainer is a symmetric implementation, storing the data as a vector of bpp::Site objects,
 *   providing site access in @f$o(1)@f$ but sequence access in @f$o(m)@f$, m being the total number of sequences.
 * The static classes bpp::SiteTools and bpp::SiteContainerTools provide some tools (for instance: pairwise alignment) to deal respectively with
 * bpp::Site and bpp::SiteContainer objects. I/O is provided for several formats, including Clustal (bpp::Clustal) and Phylip (bpp::Phylip).
 *
 * @par
 * Bio++ SeqLib also contains support for sequence properties, like amino-acids biochemical properties.
 * The interfaces bpp::AlphabetIndex1 and bpp::AlphabetIndex2 provides methods to deal with indices in 1 and 2 dimensions
 * respectively. Several basic properties are provided, together with input from the AAIndex databases.
 */

namespace bpp
{

/**
 * @brief The Alphabet interface.
 * 
 * An alphabet object defines all the states allowed for a particular type of
 * sequence. These states are coded as a string and an integer.
 * The string description is the one found in the text (human comprehensive)
 * description of sequences, typically in sequence files.
 * However, for computionnal needs, this is often more efficient to store the sequences as
 * a vector of integers.
 * The link between the two descriptions is made <i>via</i>
 * the Alphabet classes, and the two methods intToChar() and charToInt().
 * The Alphabet interface also provides other methods, like getting the full name
 * of the states and so on.
 * 
 * The alphabet objects may throw several exceptions derived of the AlphabetException
 * class.
 *
 * @see AlphabetException, BadCharException, BadIntException 
 */
class Alphabet
{
  public:
    Alphabet() {}
    virtual ~Alphabet() {}
  
  public:
    
    /**
     * @brief Get the complete name of a state given its int description.
     *
     * In case of several states with identical number (i.e. N and X for nucleic alphabets),
     * this method returns the name of the first found in the vector.
     *
     * @param state The int description of the given state.
     * @return The name of the state.
     * @throw BadIntException When state is not a valid integer.
     */
    virtual std::string getName(int state) const throw (BadIntException)  = 0;
        
    /**
     * @brief Get the complete name of a state given its string description.
     *
     * In case of several states with identical number (i.e. N and X for nucleic alphabets),
     * this method will return the name of the first found in the vector.
     *
     * @param state The string description of the given state.
     * @return The name of the state.
     * @throw BadCharException When state is not a valid char description.
     */
    virtual std::string getName(const std::string& state) const throw (BadCharException) = 0;

    /**
     * @name = Tests
     *
     * @{
     */

    /**
     * @brief Tell if a state (specified by its int description) is allowed by the
     * the alphabet.
     *
     * @param state The int description.
     * @return 'true' if the state in known.
     */
    virtual bool isIntInAlphabet(int state) const = 0;
    
    /**
     * @brief Tell if a state (specified by its string description) is allowed by the
     * the alphabet.
     *
     * @param state The string description.
     * @return 'true' if the state in known.
     */
    virtual bool isCharInAlphabet(const std::string& state) const = 0;
    /** @} */

    /**
     * @name State access
     *
     * @{
     */

    /**
     * @brief Get a state given its int description.
     *
     * @param state The int description.
     * @return The AlphabetState.
     * @throw BadIntException When state is not a valid integer.
     */
    virtual const AlphabetState& getState(int state) const throw (BadIntException) = 0;

    /**
     * @brief Get a state given its string description.
     *
     * @param state The string description.
     * @return The AlphabetState.
     * @throw BadCharException When state is not a valid string.
     */
    virtual const AlphabetState& getState(const std::string& state) const throw (BadCharException) = 0;

    /** @} */
        
    /**
     * @name Conversion methods
     *
     * @{
     */

    /**
     * @brief Give the string description of a state given its int description.
     *
     * @param state The int description.
     * @return The string description.
     * @throw BadIntException When state is not a valid integer.
     */
    virtual std::string intToChar(int state) const throw (BadIntException) = 0;
        
    /**
     * @brief Give the int description of a state given its string description.
     *
     * @param state The string description.
     * @return The int description.
     * @throw BadCharException When state is not a valid char description.
     */
    virtual int charToInt(const std::string& state) const throw (BadCharException) = 0;
    /** @} */
        
    /**
     * @name Sizes.
     *
     * @{
     */
    
    /**
     * @brief Get the number of supported characters in this alphabet,
     * including generic characters (e.g. return 20 for DNA alphabet).
     *
     * @return The total number of supported character descriptions.
     */
    virtual unsigned int getNumberOfChars() const = 0;
        
    /**
     * @brief Get the number of <strong>distinct</strong> states in alphabet (e.g. return 15 for DNA alphabet).
     * This is the number of integers used for state description.
     *
     * @return The number of distinct states.
     */
    virtual unsigned int getNumberOfTypes() const = 0;
        
    /**
     * @brief Get the number of <strong>resolved</strong> states in the alphabet (e.g. return 4 for DNA alphabet).
     * This is the method you'll need in most cases.
     *
     * @return The number of resolved states.
     */
    virtual unsigned int getSize() const = 0;
    
    /** @} */
        
    /**
     * @name Utilitary methods
     *
     * @{
     */
    
    /**
     * @brief Get all resolved states that match a generic state.
     *
     * If the given state is not a generic code then the output vector will contain this unique code.
     *
     * @param state The alias to resolve.
     * @return A vector of resolved states.
     * @throw BadIntException When state is not a valid integer.
     */
    virtual std::vector<int> getAlias(int state) const throw (BadIntException) = 0;
    
    /**
     * @brief Get all resolved states that match a generic state.
     *
     * If the given state is not a generic code then the output vector will contain this unique code.
     *
     * @param state The alias to resolve.
     * @return A vector of resolved states.
     * @throw BadCharException When state is not a valid char description.
     */
    virtual std::vector<std::string> getAlias(const std::string& state) const throw (BadCharException) = 0;

    /**
     * @brief Get the generic state that match a set of states.
     *
     * If the given states contain generic code, each generic code is first resolved and then the new generic state is returned.
     * If only a single resolved state is given the function return this state.
     *
     * @param states A vector of states to resolve.
     * @return A int code for the computed state.
     * @throw BadIntException When a state is not a valid integer.
     */
    virtual int getGeneric(const std::vector<int>& states) const throw (BadIntException) = 0;

    /**
     * @brief Get the generic state that match a set of states.
     *
     * If the given states contain generic code, each generic code is first resolved and then the new generic state is returned.
     * If only a single resolved state is given the function return this state.
     *
     * @param states A vector of states to resolve.
     * @return A string code for the computed state.
     * @throw BadCharException when a state is not a valid char description.
     * @throw CharStateNotSupportedException when the alphabet does not support Char state for unresolved state.
     */
    virtual std::string getGeneric(const std::vector<std::string>& states) const throw (AlphabetException) = 0;

    /**
     * @return A list of all supported int codes.
     *
     * Note for developers of new alphabets:
     * we return a const reference here since the list is supposed to be
     * stored within the class and should not be modified outside the class.
     */
    virtual const std::vector<int>& getSupportedInts() const = 0;
    
    /**
     * @return A list of all supported character codes.
     *
     * Note for developers of new alphabets:
     * we return a const reference here since the list is supposed to be
     * stored within the class and should not be modified outside the class.
     */
    virtual const std::vector<std::string>& getSupportedChars() const = 0;
    
    /**
     * @return The int code for unknown characters.
     */
    virtual int getUnknownCharacterCode() const = 0;

    /**
     * @return The int code for gap characters.
     */
    virtual int getGapCharacterCode() const = 0;

    /**
     * @param state The state to test.
     * @return 'True' if the state is a gap.
     */
    virtual bool isGap(int state) const = 0;

    /**
     * @param state The state to test.
     * @return 'True' if the state is a gap.
     */
    virtual bool isGap(const std::string& state) const = 0;

    /**
     * @param state The state to test.
     * @return 'True' if the state is unresolved.
     */
    virtual bool isUnresolved(int state) const = 0;

    /**
     * @param state The state to test.
     * @return 'True' if the state is unresolved.
     */
    virtual bool isUnresolved(const std::string& state) const = 0;

    /** @} */

    /**
     * @brief Identification method.
     *
     * Used to tell if two alphabets describe the same type of sequences.
     * For instance, this method is used by sequence containers to compare two alphabets and
     * allow or deny addition of sequences.
     *
     * @return A text describing the alphabet.
     */
    virtual std::string getAlphabetType() const = 0;

    /**
     * @brief Get the size of the string coding a state.
     * @return The size of the tring coding each states in the Alphabet.
     * @author Sylvain Gaillard
     */
    virtual unsigned int getStateCodingSize() const = 0;
};

} //end of namespace bpp.

#endif // _ALPHABET_H_

