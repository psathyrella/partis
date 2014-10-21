//
// File: SequenceContainer.h
// Created by: Guillaume Deuchst
//             Julien Dutheil
// Created on: Fri Jul 25 2003
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

#ifndef _SEQUENCECONTAINER_H_
#define _SEQUENCECONTAINER_H_

#include "../Alphabet/Alphabet.h"
#include "../Sequence.h"
#include "SequenceContainerExceptions.h"
#include <Bpp/Clonable.h>

// From the STL:
#include <string>

namespace bpp
{

/**
 * @brief The SequenceContainer interface.
 *
 * This interface is the most general one in the container hierarchy.
 * No assumption is made on the sequences in the container (no ordering, no alignment).
 * Sequences may be retrieved using their names, which must be unique.
 *
 * The container is the only one responsible for the allocation/destruction of sequences it
 * contains. This means that any sequence passed to it will be <strong>copied</strong> into the container.
 * The container also provides methods that send const pointers towards these sequences
 * (without performing any copy of the underlying objects).
 *
 * Notes :
 * 1. methods for adding sequences to the container are not declared here
 * (so they can't be used throught this interface),
 * because these methods take sequence container's type specific parameters
 * (i.e. a key for map sequence containers);
 * 2. to delete a sequence from a container, one must use the appropriate method
 * (removeSequence() and deleteSequence()).
 * These methods performs a few check, and properly update pointers.
 * You should never delete a sequence from a container by yourself.
 *
 * @see Sequence
 */

class SequenceContainer:
  public virtual Clonable
{
	public:
		SequenceContainer() {}
		virtual ~SequenceContainer() {}

	public:
		/**
		 * @brief Get sequence container's alphabet.
		 *
		 * @return The alphabet associated to this container.
		 */
		virtual const Alphabet* getAlphabet() const = 0;
		
		/**
		 * @brief Get the content of a sequence.
		 *
		 * @param name The name of the sequence.
		 * @return The content of the sequence as a vector of integers.
		 * @throw SequenceNotFoundException If the name does not match any sequence in the container.
		 */
		virtual const std::vector<int>& getContent(const std::string& name) const throw (SequenceNotFoundException) = 0;  
		
		/**
		 * @brief Convert a particular sequence to a string.
		 *
		 * @param name The name of the sequence.
		 * @return A string describing the content of the sequence.
		 * @throw SequenceNotFoundException If the name does not match any sequence in the container.
		 */
		virtual std::string toString(const std::string& name) const throw (SequenceNotFoundException) = 0;  

		/**
		 * @brief Retrieve a sequence object from the container.
		 *
		 * @param name The name of the sequence.
		 * @return A reference toward the Sequence with corresponding name.
		 * @throw SequenceNotFoundException If the name does not match any sequence in the container.
		 */
		virtual const Sequence& getSequence(const std::string& name) const throw (SequenceNotFoundException) = 0;

		/**
		 * @brief Check if a sequence with a given name is present in the container.
		 *
		 * @param name The name of the sequence.
		 * @return True if a sequence with the given name is present in the container.
		 */
		virtual bool hasSequence(const std::string& name) const = 0;

		/**
		 * @brief Add a sequence to the container.
		 *
		 * @param sequence  The sequence to add.
		 * @param checkName Tell if the container must check if the name of the sequence
		 * is already used in the container before adding it.
		 * @throw Exception Any other kind of exception, if the name of the sequence is
		 * already used, are whatever else depending on the implementation.
		 */
		virtual void addSequence(const Sequence& sequence, bool checkName) throw (Exception) = 0;

		/**
		 * @brief Replace a sequence in the container.
		 *
		 * @param name      The name of the sequence.
		 * @param sequence  The sequence to add.
		 * @param checkName Tell if the container must check if the name of the sequence
		 * is already used in the container before adding it.
		 * @throw SequenceNotFoundException If the name does not match any sequence in
		 * the container.
		 * @throw Exception Any other kind of exception, if the name of the sequence is
		 * already used, are whatever else depending on the implementation.
		 */
		virtual void setSequence(const std::string& name, const Sequence& sequence, bool checkName) throw (Exception) = 0;

		/**
		 * @brief Extract (and remove) a sequence from the container.
		 *
		 * @param name The name of the sequence.
		 * @throw SequenceNotFoundException If the name does not match any sequence in
		 * the container.
		 */
		virtual Sequence* removeSequence(const std::string& name) throw (SequenceNotFoundException, Exception) = 0;
		
		/**
		 * @brief Delete a sequence of the container.
		 *
		 * @param name The name of the sequence.
		 * @throw SequenceNotFoundException If the name does not match any sequence in
		 * the container.
		 */
		virtual void deleteSequence(const std::string& name) throw (SequenceNotFoundException, Exception) = 0;

		/**
		 * @brief Get the number of sequences in the container.
		 *
		 * @return The number of sequences in the container.
		 */
		virtual size_t getNumberOfSequences() const = 0;

		/**
		 * @brief Get all the names of the sequences in the container.
		 *
		 * @return A vector of strings with all sequence names.
		 */
		virtual std::vector<std::string> getSequencesNames() const = 0;
		
		/**
		 * @brief Set all sequence names.
		 *
		 * @param names A vector of strings with all sequence names.
		 * Its size must be strictly equal to the the size of the container (the number of
		 * sequences).
		 * @param checkNames Tell if the container must check if the name of the sequence
		 * is already used in the container before adding it.
		 * @throw Exception If there are redundant names in the input vector.
		 */
		virtual void setSequencesNames(const std::vector<std::string>& names, bool checkNames) throw (Exception) = 0;

		/**
		 * @brief Get comments of a particular sequence.
		 *
		 * @param name The name of the sequence.
		 * @return The comments associated to sequence with name 'name'.
		 * @throw SequenceNotFoundException If the name does not match any sequence in
		 * the container.
		 */
		virtual const Comments& getComments(const std::string& name) const throw (SequenceNotFoundException) = 0;
		
		/**
		 * @brief Set the comments of a particular sequence.
		 *
		 * @param name     The name of the sequence.
		 * @param comments The comments to set to sequence with name 'name'.
		 * @throw SequenceNotFoundException If the name does not match any sequence in
		 * the container.
		 */
		virtual void setComments(const std::string& name, const Comments& comments) throw (SequenceNotFoundException) = 0;
		
		/**
		 * @brief Get the comments of this container.
		 *
		 * @return The comments associated to this container.
		 */
		virtual const Comments& getGeneralComments() const = 0;

		/**
		 * @brief Set the comments of this container.
		 *
		 * @param comments The comments to be associated to this container.
		 */
		virtual void setGeneralComments(const Comments& comments) = 0;
		
		/**
		 * @brief Delete the comments associated to this container.
		 */
		virtual void deleteGeneralComments() = 0;

		/**
		 * @brief Delete all sequences in the container.
		 */
		virtual void clear() = 0;

		/**
		 * @brief Return a copy of this container, but with no sequence inside.
		 *
		 * This method creates a new SequenceContainer objet.
		 * The class of this container depends on the derivative class.
		 *
		 * @return A new empty container, with the same alphabet as this one.
		 */
		virtual SequenceContainer* createEmptyContainer() const = 0;

    /**
     * @name Provide direct access to sequences content.
     *
     * @warning These operators allow you to modifiy the content of the sequences.
     * No checking is performed for your modifications, so use with care, or
     * consider using the setContent() methods.
     *
     * @{
     */
    
    /**
     * @brief Element access function.
     *
     * Allows direct access to the data stored in the container.
     * 
     * @param sequenceName The sequence name.
     * @param elementIndex The element position within the sequence.
     * @throw SequenceNotFoundException If no corresponding sequence is found in the container.
     * @throw IndexOutOfBoundsException If the element position is not valid.
     */
    virtual int& valueAt(const std::string& sequenceName, size_t elementIndex) throw (SequenceNotFoundException, IndexOutOfBoundsException) = 0;

    /**
     * @brief Element access function.
     *
     * Allows direct access to the data stored in the container.
     * 
     * @param sequenceName The sequence name.
     * @param elementIndex The element position within the sequence.
     * @throw SequenceNotFoundException If no corresponding sequence is found in the container.
     * @throw IndexOutOfBoundsException If the element position is not valid.
     */
    virtual const int& valueAt(const std::string& sequenceName, size_t elementIndex) const throw (SequenceNotFoundException, IndexOutOfBoundsException) = 0;

    /**
     * @brief Element access operator.
     *
     * Allows direct access to the data stored in the container.
     * This method is faster then the valueAt function, but input
     * parameters are not checked!
     * 
     * @param sequenceName The sequence name.
     * @param elementIndex The element position within the sequence.
     */
    virtual int& operator()(const std::string& sequenceName, size_t elementIndex) = 0;

    /**
     * @brief Element access operator.
     *
     * Allows direct access to the data stored in the container.
     * This method is faster then the valueAt function, but input
     * parameters are not checked!
     * 
     * @param sequenceName The sequence name.
     * @param elementIndex The element position within the sequence.
     */
    virtual const int& operator()(const std::string& sequenceName, size_t elementIndex) const = 0;
    /** @} */
};

} //end of namespace bpp.

#endif	// _SEQUENCECONTAINER_H_

