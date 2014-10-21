// 
// File:    SequencePositionIterators.h
// Author:  Sylvain Gaillard
// Created: 23/06/2009 10:35:28
// 

/*
Copyright or © or Copr. Bio++ Development Team, (June 23, 2009)

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

#ifndef _SEQUENCEPOSITIONITERATORS_H_
#define _SEQUENCEPOSITIONITERATORS_H_

// from STL
#include <string>

#include "Sequence.h"

namespace bpp
{

  /**
   * @brief Loop over a Sequence.
   *
   * This is the SequencePositionIterator interface.
   *
   * @author Sylvain Gaillard
   */
  class SequencePositionIterator
  {
    public:
      SequencePositionIterator() {}
      virtual ~SequencePositionIterator() {}

    public: 
      /**
       * @brief Get the actual position of the iterator in the Sequence.
       */
      virtual unsigned int getPosition() const = 0;
      /**
       * @brief Set the position of the iterator.
       * @param pos The position on the Sequence
       */
      virtual void setPosition(unsigned int pos) = 0;
      /**
       * @brief Get the numerical value of the Sequence at current position.
       */
      virtual int getValue() const = 0;
      /**
       * @brief Get the textual value of the Sequence at current position.
       */
      virtual std::string getChar() const = 0; 

      virtual bool operator==(const SequencePositionIterator & it) const = 0;
      virtual bool operator!=(const SequencePositionIterator & it) const = 0;
      virtual SequencePositionIterator & operator+=(int i) = 0;
      virtual SequencePositionIterator & operator-=(int i) = 0;
      virtual SequencePositionIterator& operator++() = 0;

      /**
       * @brief Tells if there is more positions in the Sequence.
       * @return true if there is more positions in the Sequence
       */
      virtual bool hasMorePositions() const = 0;
      /**
       * @brief Get the Sequence on which the iterator loops.
       * @return A reference toward the Sequence object.
       */
      virtual const Sequence & getSequence() const = 0;
  };

  /**
   * @brief Partial implementation of the SequencePositionIterator interface.
   *
   * @author Sylvain Gaillard
   */
  class AbstractSequencePositionIterator :
    public virtual SequencePositionIterator
  {
    private:
      const Sequence* sequence_;
      unsigned int currentPosition_;


    public:
      AbstractSequencePositionIterator(const Sequence& seq, unsigned int pos = 0) :
        sequence_(&seq), currentPosition_(pos) {}
      
      AbstractSequencePositionIterator(const AbstractSequencePositionIterator& aspi) :
        sequence_(aspi.sequence_), currentPosition_(aspi.currentPosition_) {}
      
      AbstractSequencePositionIterator& operator=(const AbstractSequencePositionIterator& aspi)
      {
        sequence_ = aspi.sequence_;
        currentPosition_ = aspi.currentPosition_;
        return *this;
      }
      
      virtual ~AbstractSequencePositionIterator() {}

    public:
      
      /**
       * @name Comparison operators
       *
       * @{
       */
      bool operator==(const SequencePositionIterator& it) const;
      bool operator!=(const SequencePositionIterator& it) const;
      /** @} */

      unsigned int getPosition() const;
      void setPosition(unsigned int pos);
      int getValue() const;
      std::string getChar() const;
      const Sequence& getSequence() const;
  };

  /**
   * @brief Loop over all positions in a Sequence
   *
   * This is the simplest implementation of SequencePositionIterator.
   * It just loops over all positions of a Sequence.
   *
   * @code
   * Sequence seq = Sequence("seq1", "ATTCGATCCG-G", &AlphabetTools::DNA_ALPHABET);
   * for (SimpleSequencePositionIterator it(seq) ; it.hasMorePositions() ; ++it) {
   *   cout << it.getPosition() << " : " << it.getValue() << " (" << it.getChar() << ")" << endl;
   * }
   * @endcode
   *
   * @author Sylvain Gaillard
   */
  class SimpleSequencePositionIterator:
    public AbstractSequencePositionIterator
  {
    public:
      /**
       * @name Constructors and destructor
       *
       * @{
       */

      /**
       * @brief General constructor
       *
       * @param seq A reference toward the Sequence object we want to loop over
       * @param pos Optional integer where to start on the Sequence object
       *
       */
      SimpleSequencePositionIterator(const Sequence& seq, unsigned int pos = 0):
        AbstractSequencePositionIterator(seq, pos) {}
      /**
       * @brief Copie constructor.
       *
       * @param it A reference toward a SequencePositionIterator
       */
      SimpleSequencePositionIterator(const SequencePositionIterator& it);
      virtual ~SimpleSequencePositionIterator() {}
       /** @} */

    public:
      /**
       * @name Operators
       *
       * @{
       */
      SimpleSequencePositionIterator & operator++();
      virtual SimpleSequencePositionIterator operator++(int i);
      SimpleSequencePositionIterator & operator+=(int i);
      SimpleSequencePositionIterator & operator-=(int i);
      virtual SimpleSequencePositionIterator operator+(int i) const;
      virtual SimpleSequencePositionIterator operator-(int i) const;
       /** @} */

      bool hasMorePositions() const;
  };

}

#endif //_SEQUENCEPOSITIONITERATORS_H_
