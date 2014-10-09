//
// File: SequenceIterator.h
// Created by: Julien Dutheil
// Created on: Tue Feb 26 14:27 2013
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

#ifndef _SEQUENCEITERATOR_H_
#define _SEQUENCEITERATOR_H_

#include "Sequence.h"
#include "SequenceWithQuality.h"

namespace bpp
{

/**
 * @brief Generic sequence iterator interface, allowing to loop over sequences.
 */
class SequenceIterator
{
	public:
		SequenceIterator() {}
		virtual ~SequenceIterator() {}
	
	public:
		virtual Sequence* nextSequence() = 0;
		virtual bool hasMoreSequences() const = 0;
};

/**
 * @brief Generic const sequence iterator interface, allowing to loop over const sequences.
 */
class ConstSequenceIterator
{
	public:
		ConstSequenceIterator() {}
		virtual ~ConstSequenceIterator() {}
	
	public:
		virtual const Sequence* nextSequence() = 0;
		virtual bool hasMoreSequences() const = 0;
};

/**
 * @brief Generic sequence iterator interface, allowing to loop over sequences with quality scores.
 */
class SequenceWithQualityIterator:
  public virtual SequenceIterator
{
	public:
		SequenceWithQualityIterator() {}
		virtual ~SequenceWithQualityIterator() {}
	
	public:
		virtual SequenceWithQuality* nextSequence() = 0;
};

/**
 * @brief Generic const sequence iterator interface, allowing to loop over const sequences with quality scores.
 */
class ConstSequenceWithQualityIterator:
  public virtual ConstSequenceIterator
{
	public:
		ConstSequenceWithQualityIterator() {}
		virtual ~ConstSequenceWithQualityIterator() {}
	
	public:
		virtual const SequenceWithQuality* nextSequence() = 0;
};

} //end of namespace bpp.

#endif	//_SEQUENCEITERATOR_H_

