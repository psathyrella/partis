//
// File:       SequenceWithAnnotationTools.h
// Authors:    Julien Dutheil
// Created on: 06 Sep 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (Sep 06, 2010)

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

#ifndef _SEQUENCEWITHANNOTATIONTOOLS_H_
#define _SEQUENCEWITHANNOTATIONTOOLS_H_

#include "SequenceTools.h"
#include "SequenceWithAnnotation.h"
#include <Bpp/Numeric/VectorTools.h>

namespace bpp {

  class SequenceMask :
    public virtual SequenceAnnotation
  {
    private:
      bool removable_;
      std::vector<bool> mask_;

    public:
      static const std::string MASK; 

    public:

      /**
       * @name Constructors
       * @{
       */

      /**
       * @brief Build a new SequenceMask object
       *
       * Build a new SequenceMask object and set the mask to false.
       *
       * @param size The size of the sequence. 
       * @param removable Tell if this listener can be removed by the user.
       */
      SequenceMask(size_t size = 0, bool removable = true) :
        removable_(removable),
        mask_(size, false)
      {}


      /**
       * @brief Build a new SequenceMask object
       *
       * Build a new SequenceMask object and set the mask as a vector of bool.
       *
       * @param mask The boolean mask
       * @param removable Tell if this listener can be removed by the user.
       */
      SequenceMask(const std::vector<bool>& mask, bool removable = true) :
        removable_(removable),
        mask_(mask)
      {}

      /** @} */

      /**
       * @name Destructor
       * @{
       */
      virtual ~SequenceMask() {}
      /** @} */

      /**
       * @name The Clonable interface
       * @{
       */
#ifdef NO_VIRTUAL_COV
      Clonable*
#else
      SequenceMask*
#endif
      clone() const { return new SequenceMask(*this); }
      /** @} */

    public:
      void init(const Sequence& seq)
      {
        mask_.resize(seq.size());
        std::fill(mask_.begin(), mask_.end(), false);
      }

      const std::string& getType() const { return MASK; }

      bool isValidWith(const SequenceWithAnnotation& sequence, bool throwException = true) const
      {
        if (throwException && mask_.size() != sequence.size()) throw Exception("SequenceMask. The mask size must match the sequence size.");
        return (mask_.size() == sequence.size());
      }

      bool isRemovable() const { return removable_; }
      bool isShared() const { return false; }
      void beforeSequenceChanged(const SymbolListEditionEvent& event) {}
      void afterSequenceChanged(const SymbolListEditionEvent& event);
      void beforeSequenceInserted(const SymbolListInsertionEvent& event) {}
      void afterSequenceInserted(const SymbolListInsertionEvent& event);
      void beforeSequenceDeleted(const SymbolListDeletionEvent& event) {}
      void afterSequenceDeleted(const SymbolListDeletionEvent& event);
      void beforeSequenceSubstituted(const SymbolListSubstitutionEvent& event) {}
      void afterSequenceSubstituted(const SymbolListSubstitutionEvent& event) {}

      size_t getSize() const { return mask_.size(); }

      const bool operator[](size_t i) const { return mask_[i]; }

      void setMask(const std::vector<bool>& mask) {
        if (mask.size() != mask_.size())
          throw DimensionException("SequenceMask::setMask. Trying to replace mask by a vector with different length.", mask.size(), mask_.size());
        mask_ = mask;
      }

      /**
       * @return The mask as a vector.
       */
      const std::vector<bool>& getMask() const { return mask_; }

      void setMask(size_t pos, bool mask) {
        if (pos >= mask_.size())
          throw Exception("SequenceMask::setMask. Vector overflow. Scores number: " + TextTools::toString(mask_.size()) + ", but trying to insert mask at position " + TextTools::toString(pos) + ".");
        mask_[pos] = mask;
      }
      
      void setMask(size_t pos, const std::vector<bool>& mask) {
        if (pos + mask.size() > mask_.size())
          throw Exception("SequenceMask::setMask. Vector overflow. Scores number: " + TextTools::toString(mask_.size()) + ", but trying to insert " + TextTools::toString(mask.size()) + " scores at position " + TextTools::toString(pos) + ".");
        std::copy(mask.begin(), mask.end(), mask_.begin() + pos); 
      }

      bool merge(const SequenceAnnotation& anno) {
        try {
          const SequenceMask* mask = & dynamic_cast<const SequenceMask&>(anno);
          VectorTools::append(mask_, mask->getMask());
          return true;
        } catch (std::exception& e) {
          return false;
        }
      }
 
      SequenceAnnotation* getPartAnnotation(size_t pos, size_t len) const throw (Exception) {
        return new SequenceMask(std::vector<bool>(mask_.begin() + pos, mask_.begin() + pos + len), removable_);
      }
  };
  
  /**
   * @brief The SequenceWithAnnotationTools static class
   *
   * Implement methods to manipulate SequencesWithAnnotation
   *
   * @author Julien Dutheil
   */

  class SequenceWithAnnotationTools {

    public:
      /**
       * @brief Parse a sequence with a CaseMaskedAlphabet and creates a new SequenceWithAnnotation object with original alphabet and a mask annotation.
       *
       * @param seq The sequence to parse.
       * @return A new SequenceWithAnnotation object.
       * @throw AlphabetException if the input sequence does not have a CaseMaskedAlphabet.
       */
      SequenceWithAnnotation* createMaskAnnotation(const Sequence& seq) throw (AlphabetException);

  };
}

#endif // _SEQUENCEWITHANNOTATIONTOOLS_H_
