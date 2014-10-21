//
// File: Transliterator.h
// Created by: Julien Dutheil
// Created on: Sun Oct 12 14:25:25 2003
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

#ifndef _TRANSLITERATOR_H_
#define _TRANSLITERATOR_H_

#include "Alphabet/Alphabet.h"
#include "Sequence.h"

namespace bpp
{

/**
 * @brief This interface is used when translating a sequence from an alphabet to another: it gives the translation rules, eg: RNA -> DNA.
 */
class Transliterator
{
	public:
		Transliterator() {}
		virtual ~Transliterator() {}
	
	public:

		/**
		 * @brief Get the source alphabet.
		 *
		 * @return The source alphabet.
		 */
		virtual const Alphabet * getSourceAlphabet() const = 0;
		
		/**
		 * @brief Get the target alphabet.
		 *
		 * @return The target alphabet.
		 */
		virtual const Alphabet * getTargetAlphabet() const = 0;
		
		/**
		 * @brief Translate a given state coded as a int from source alphabet to target alphabet.
		 *
		 * @param state A state in source alphabet.
		 * @return The corresponding state in target alphabet.
		 * @throw BadIntException If the state is not a proper state for source alphabet.
		 * @throw Exception       Other kind of error, depending on the implementation.
		 */
		virtual int translate(int state) const throw (BadIntException, Exception) = 0;		
		
		/**
		 * @brief Translate a given state coded as a string from source alphabet to target alphabet.
		 *
		 * @param state A state in source alphabet.
		 * @return The corresponding state in target alphabet.
		 * @throw BadCharException If the state is not a proper state for source alphabet.
		 * @throw Exception        Other kind of error, depending on the implementation.
		 */
		virtual std::string translate(const std::string & state) const throw (BadCharException, Exception) = 0;
		
		/**
		 * @brief Translate a whole sequence from source alphabet to target alphabet.
		 *
		 * @param sequence A sequence in source alphabet.
		 * @return The corresponding sequence in target alphabet.
		 * @throw AlphabetMismatchException If the sequence alphabet do not match the source alphabet.
		 * @throw Exception                 Other kind of error, depending on the implementation.
		 */	
		virtual Sequence * translate(const Sequence & sequence) const throw (AlphabetMismatchException, Exception) = 0; 
};

/**
 * @brief The same as previous, but can perform the reverse translation, eg: RNA -> DNA and DNA -> RNA;
 */
class ReverseTransliterator:
  public virtual Transliterator
{
	public:
		ReverseTransliterator() {}
		virtual ~ReverseTransliterator() {}
	
	public:
		
		/**
		 * @brief Translate a given state coded as a int from target alphabet to source alphabet.
		 *
		 * @param state A state in target alphabet.
		 * @return The corresponding state in source alphabet.
		 * @throw BadIntException If the state is not a proper state for target alphabet.
		 * @throw Exception       Other kind of error, depending on the implementation.
		 */
		virtual int reverse(int state) const throw (BadIntException, Exception) = 0;
		
		/**
		 * @brief Translate a given state coded as a string from target alphabet to source alphabet.
		 *
		 * @param state A state in target alphabet.
		 * @return The corresponding state in source alphabet.
		 * @throw BadCharException If the state is not a proper state for target alphabet.
		 * @throw Exception        Other kind of error, depending on the implementation.
		 */
		virtual std::string reverse(const std::string & state) const throw (BadCharException, Exception) = 0;			

		/**
		 * @brief Translate a whole sequence from target alphabet to source alphabet.
		 *
		 * @param sequence A sequence in target alphabet.
		 * @return The corresponding sequence in source alphabet.
		 * @throw AlphabetMismatchException If the sequence alphabet do not match the target alphabet.
		 * @throw Exception                 Other kind of error, depending on the implementation.
		 */	
		virtual Sequence * reverse(const Sequence & sequence) const throw (AlphabetMismatchException, Exception) = 0; 
};

/**
 * @brief Partial implementation of the Transliterator interface.
 */
class AbstractTransliterator:
  public virtual Transliterator
{
	public:
		AbstractTransliterator() {}
		virtual ~AbstractTransliterator() {}
	
	public:
		virtual int translate(int state) const throw (BadIntException, Exception) = 0;		
		virtual std::string translate(const std::string & state) const throw (BadCharException, Exception) = 0;
		virtual Sequence * translate(const Sequence & sequence) const throw (AlphabetMismatchException, Exception); 
};

/**
 * @brief Partial implementation of the ReverseTransliterator interface.
 */
class AbstractReverseTransliterator:
  public ReverseTransliterator,
  public AbstractTransliterator
{
	public:
		AbstractReverseTransliterator() {}
		virtual ~AbstractReverseTransliterator() {}
	
	public:
		//These two redeclarations must be here because of the multiple inheritance.
		virtual const Alphabet * getSourceAlphabet() const = 0;
		virtual const Alphabet * getTargetAlphabet() const = 0;
		virtual int reverse(int state) const throw (BadIntException, Exception) = 0;
		virtual std::string reverse(const std::string & state) const throw (BadCharException, Exception) = 0;			
		virtual Sequence * reverse(const Sequence & sequence) const throw (AlphabetMismatchException, Exception);
};

} //end of namespace bpp.

#endif	//_TRANSLITERATOR_H_

