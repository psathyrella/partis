//
// File: SequenceContainerExceptions.h
// Created by: Julien Dutheil
// Created on: Mon Nov  3 17:00:05 2003
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#ifndef _SEQUENCECONTAINEREXCEPTIONS_H_
#define _SEQUENCECONTAINEREXCEPTIONS_H_

#include <Bpp/Exceptions.h>

namespace bpp
{

class SequenceContainer;

/**
 * @brief Exception thrown when a sequence is not found The sequence not found exception base class.
 */
class SequenceNotFoundException :
  public Exception
{

	protected:

		/**
		 * @brief The id of the sequence that was to be found.
		 */
		const std::string id;
	
	public:

		/**
		 * @brief Build a new SequenceNotFoundException object.
		 *
		 * @param text  A message to be passed to the exception hierarchy.
		 * @param seqId A the id of the sequence that was to be found.
		 */
		SequenceNotFoundException(const char * text, const char * seqId = "") :
    	Exception("SequenceNotFoundException: " + std::string(text) + "(" + seqId + ")"),
	    id(seqId) {};

		/**
		 * @brief Build a new SequenceNotFoundException object.
		 *
		 * @param text  A message to be passed to the exception hierarchy.
		 * @param seqId A the id of the sequence that was to be found.
		 */
		SequenceNotFoundException(const std::string & text, const std::string & seqId = "") :
    	Exception("SequenceNotFoundException: " + text + "(" + seqId + ")"),
	    id(seqId) {};
	
		// Class destructor
		virtual ~SequenceNotFoundException() throw() {}
	
	public:

		/**
		 * @brief Get the id of the sequence that was to be found.
		 *
		 * @return The id of the sequence that was to be found.
		 */
		virtual const std::string getSequenceId() const { return id; }
};

/**
 * @brief Exception thrown when an empty container is found.
 */
class EmptyContainerException :
  public Exception
{

	private:

		/**
		 * @brief The empty container.
     */
		const SequenceContainer *container_;
	
	public:

		/**
		 * @brief Build a new EmptyContainerException object.
		 *
		 * @param text  A message to be passed to the exception hierarchy.
		 * @param container The empty container.
		 */
		EmptyContainerException(const std::string& text, const SequenceContainer* container) :
    	Exception("EmptyContainerException: " + text),
	    container_(container) {};
      
    EmptyContainerException(const EmptyContainerException& ece):
      Exception(ece), container_(ece.container_) {}
	
    EmptyContainerException& operator=(const EmptyContainerException& ece)
    {
      Exception::operator=(ece);
      container_ = ece.container_;
      return *this;
    }
	
		// Class destructor
		virtual ~EmptyContainerException() throw() {}
	
	public:

		/**
		 * @return The empty container.
		 */
		virtual const SequenceContainer* getContainer() const { return container_; }
};

} //end of namespace bpp.

#endif	//_SEQUENCECONTAINEREXCEPTIONS_H_

