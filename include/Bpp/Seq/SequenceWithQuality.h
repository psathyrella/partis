// 
// File:    SequenceWithQuality.h
// Authors: Sylvain Gaillard
//          Vincent Cahais
//          Julien Dutheil
// Created: 19/01/2010 16:01:20
// 

/*
Copyright or © or Copr. Bio++ Development Team, (January 19, 2010)

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

#ifndef _SEQUENCEQUALITY_H_
#define _SEQUENCEQUALITY_H_

#include "SequenceWithAnnotation.h"
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/VectorExceptions.h>

// From the STL

#include <string>
#include <vector>

namespace bpp {
  /**
   * @brief The SequenceQuality class
   *
   * This is a sequence with quality score associated to each element.
   * The score is a signed int value that can represent the phred or the
   * Solexa quality score for nucleic sequence.
   *
   * @author Sylvain Gaillard, Vincent Cahais, Julien Dutheil
   */
  class SequenceQuality :
    public virtual SequenceAnnotation
  {
    private:
      bool removable_;
      std::vector<int> qualScores_;

    public:
      static const std::string QUALITY_SCORE; 
      static const int DEFAULT_QUALITY_VALUE;

    public:

      /**
       * @name Constructors
       * @{
       */

      /**
       * @brief Build a new SequenceQuality object
       *
       * Build a new SequenceQuality object and set the quality scores to
       * the default value DEFAULT_QUALITY_VALUE.
       *
       * @param size The size of the sequence. 
       * @param removable Tell if this listener can be removed by the user.
       */
      SequenceQuality(size_t size = 0, bool removable = true) :
        removable_(removable),
        qualScores_(size, DEFAULT_QUALITY_VALUE) {}


      /**
       * @brief Build a new SequenceQuality object
       *
       * Build a new SequenceQuality and assign quality scores from
       * a vector of int.
       *
       * @param quality The quality scores
       * @param removable Tell if this listener can be removed by the user.
       */
      SequenceQuality(const std::vector<int>& quality, bool removable = true) :
        removable_(removable),
        qualScores_(quality)
      {
      //    if (size() != qualScores_.size())
      //      throw DimensionException("SequenceWithQuality constructor: sequence and quality must have the same length", qualScores_.size(), size());
      }

      /** @} */

      /**
       * @name Destructor
       * @{
       */
      virtual ~SequenceQuality() {}
      /** @} */

      /**
       * @name The Clonable interface
       * @{
       */
#ifdef NO_VIRTUAL_COV
      Clonable*
#else
      SequenceQuality*
#endif
      clone() const { return new SequenceQuality(*this); }
      /** @} */

    public:
      void init(const Sequence& seq)
      {
        qualScores_.resize(seq.size());
        std::fill(qualScores_.begin(), qualScores_.end(), DEFAULT_QUALITY_VALUE);
      }

      const std::string& getType() const { return QUALITY_SCORE; }

      bool isValidWith(const SequenceWithAnnotation& sequence, bool throwException = true) const
      {
        if (throwException && qualScores_.size() != sequence.size()) throw Exception("SequenceQuality. Quality scores must match the sequence size.");
        return (qualScores_.size() == sequence.size());
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

      size_t getSize() const { return qualScores_.size(); }

      const int& operator[](size_t i) const { return qualScores_[i]; }
      int& operator[](size_t i) { return qualScores_[i]; }

      void setScores(const std::vector<int>& scores) {
        if (scores.size() != qualScores_.size())
          throw DimensionException("SequenceQuality::setScores. Trying to replace score by a vector with different length.", scores.size(), qualScores_.size());
        qualScores_ = scores;
      }

      /**
       * @return All scores as a vector.
       */
      const std::vector<int>& getScores() const { return qualScores_; }

      void setScore(size_t pos, int score) {
        if (pos >= qualScores_.size())
          throw Exception("SequenceQuality::setScore. Vector overflow. Scores number: " + TextTools::toString(qualScores_.size()) + ", but trying to insert score at position " + TextTools::toString(pos) + ".");
        qualScores_[pos] = score;
      }
      
      void setScores(size_t pos, const std::vector<int>& scores) {
        if (pos + scores.size() > qualScores_.size())
          throw Exception("SequenceQuality::setScores. Vector overflow. Scores number: " + TextTools::toString(qualScores_.size()) + ", but trying to insert " + TextTools::toString(scores.size()) + " scores at position " + TextTools::toString(pos) + ".");
        std::copy(scores.begin(), scores.end(), qualScores_.begin() + pos); 
      }
    
      bool merge(const SequenceAnnotation& anno) {
        try {
          const SequenceQuality* qual = & dynamic_cast<const SequenceQuality&>(anno);
          VectorTools::append(qualScores_, qual->getScores());
          return true;
        } catch (std::exception& e) {
          return false;
        }
      }

      SequenceQuality* getPartAnnotation(size_t pos, size_t len) const throw (Exception) {
        return new SequenceQuality(std::vector<int>(qualScores_.begin() + pos, qualScores_.begin() + pos + len), removable_);
      }
  };



  /**
   * @brief A SequenceWithAnnotation class with quality scores attached.
   *
   * This classes adds some usefull functions to handle quality scores.
   *
   * @see SequenceQuality
   * @author Sylvain Gaillard, Vincent Cahais, Julien Dutheil
   */
  class SequenceWithQuality :
    public SequenceWithAnnotation
  {
    private:
      SequenceQuality* qualScores_;

    public:

      /**
       * @name Constructors
       * @{
       */

      /**
       * @brief Build a new empty SequenceWithQuality
       *
       * @param alpha    A pointer to an Alphabet
       *
       * @throw BadCharException if a state is not alowed by the Alphabet
       */
      SequenceWithQuality(
          const Alphabet* alpha
          ):
        SequenceWithAnnotation(alpha),
        qualScores_(new SequenceQuality(0, false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality from a std::string
       *
       * Build a new SequenceWithQuality and set the quality scores to
       * the default value DEFAULT_QUALITY_VALUE.
       *
       * @param name     The name of the sequence
       * @param sequence The string representing the sequence
       * @param alpha    A pointer to an Alphabet
       *
       * @throw BadCharException if a state is not alowed by the Alphabet
       */
      SequenceWithQuality(
          const std::string& name,
          const std::string& sequence,
          const Alphabet* alpha
          ) throw (BadCharException):
        SequenceWithAnnotation(name, sequence, alpha),
        qualScores_(new SequenceQuality(sequence.size(), false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality from a std::string
       *
       * Build a new SequenceWithQuality and set the quality scores to
       * the default value DEFAULT_QUALITY_VALUE.
       *
       * @param name     The name of the sequence
       * @param sequence The string representing the sequence
       * @param comments Comments to add to the sequence
       * @param alpha    A pointer to an Alphabet
       *
       * @throw BadCharException if a state is not alowed by the Alphabet
       *
       * @author Vincent Cahais
       */
      SequenceWithQuality(
          const std::string& name,
          const std::string& sequence,
          const Comments& comments,
          const Alphabet* alpha
          ) throw (BadCharException):
        SequenceWithAnnotation(name, sequence, comments, alpha),
        qualScores_(new SequenceQuality(sequence.size(), false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality from a std::string
       *
       * Build a new SequenceWithQuality and assign quality scores from
       * a vector of int.
       *
       * @param name The name of the sequence
       * @param sequence The string representing the sequence
       * @param quality The quality scores
       * @param alpha A pointer to an alphabet
       *
       * @throw BadCharException if a state is not alowed by the Alphabet
       * @throw DimensionException if the number of quality values is not equal
       * to the number of sequence states
       */
      SequenceWithQuality(
          const std::string& name,
          const std::string& sequence,
          const std::vector<int>& quality,
          const Alphabet* alpha)
        throw (BadCharException, DimensionException):
        SequenceWithAnnotation(name, sequence, alpha),
        qualScores_(new SequenceQuality(quality, false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality from a std::string
       *
       * Build a new SequenceWithQuality and assign quality scores from
       * a vector of int.
       *
       * @param name The name of the sequence
       * @param sequence The string representing the sequence
       * @param quality The quality scores
       * @param comments Comments to add to the sequence
       * @param alpha A pointer to an alphabet
       *
       * @throw BadCharException if a state is not alowed by the Alphabet
       * @throw DimensionException if the number of quality values is not equal
       * to the number of sequence states
       *
       * @author Vincent Cahais
       */
      SequenceWithQuality(
          const std::string& name,
          const std::string& sequence,
          const std::vector<int>& quality,
          const Comments& comments,
          const Alphabet* alpha)
        throw (BadCharException, DimensionException):
        SequenceWithAnnotation(name, sequence, comments, alpha),
        qualScores_(new SequenceQuality(quality, false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality from a std::vector<int>
       *
       * Build a new SequenceWithQuality and set the quality scores to
       * the default value DEFAULT_QUALITY_VALUE.
       *
       * @param name The name of the sequence
       * @param sequence The sequence in int
       * @param alpha A pointer to an Alphabet
       *
       * @throw BadIntException if a state is not alowed by the Alphabet
       */
      SequenceWithQuality(
          const std::string& name,
          const std::vector<int>& sequence,
          const Alphabet* alpha)
        throw (BadIntException):
        SequenceWithAnnotation(name, sequence, alpha),
        qualScores_(new SequenceQuality(sequence.size(), false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality from a std::vector<int>
       *
       * Build a new SequenceWithQuality and set the quality scores to
       * the default value DEFAULT_QUALITY_VALUE.
       *
       * @param name The name of the sequence
       * @param sequence The sequence in int
       * @param comments Comments to add to the sequence
       * @param alpha A pointer to an Alphabet
       *
       * @throw BadIntException if a state is not alowed by the Alphabet
       *
       * @author Vincent Cahais
       */
      SequenceWithQuality(
          const std::string& name,
          const std::vector<int>& sequence,
          const Comments& comments,
          const Alphabet* alpha)
        throw (BadIntException):
        SequenceWithAnnotation(name, sequence, comments, alpha),
        qualScores_(new SequenceQuality(sequence.size(), false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality from a std::vector<int>
       *
       * Build a new SequenceWithQuality and assign quality scores from
       * a vector of int.
       *
       * @param name The name of the sequence
       * @param sequence The sequence in int
       * @param quality The quality scores
       * @param alpha A pointer to an Alphabet
       *
       * @throw BadIntException if a state is not alowed by the Alphabet
       * @throw DimensionException if the number of quality values is not equal
       * to the number of sequence states
       */
      SequenceWithQuality(
          const std::string& name,
          const std::vector<int>& sequence,
          const std::vector<int>& quality,
          const Alphabet* alpha)
        throw (BadIntException, DimensionException):
        SequenceWithAnnotation(name, sequence, alpha),
        qualScores_(new SequenceQuality(quality, false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality from a std::vector<int>
       *
       * Build a new SequenceWithQuality and assign quality scores from
       * a vector of int.
       *
       * @param name The name of the sequence
       * @param sequence The sequence in int
       * @param quality The quality scores
       * @param comments Comments to add to the sequence
       * @param alpha A pointer to an Alphabet
       *
       * @throw BadIntException if a state is not alowed by the Alphabet
       * @throw DimensionException if the number of quality values is not equal
       * to the number of sequence states
       *
       * @author Vincent Cahais
       */
      SequenceWithQuality(
          const std::string& name,
          const std::vector<int>& sequence,
          const std::vector<int>& quality,
          const Comments& comments,
          const Alphabet* alpha)
        throw (BadIntException, DimensionException):
        SequenceWithAnnotation(name, sequence, comments, alpha),
        qualScores_(new SequenceQuality(quality, false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality
       *
       * Build a new SequenceWithQuality from a Sequence object and set the
       * quality scores to the default value DEFAULT_QUALITY_VALUE.
       *
       * @param s The Sequence object
       */
      SequenceWithQuality(const Sequence& s) :
        SequenceWithAnnotation(s), qualScores_(new SequenceQuality(s.size(), false))
      {
        addAnnotation(qualScores_);
      }

      /**
       * @brief Build a new SequenceWithQuality
       *
       * Build a new SequenceWithQuality from a Sequence object and set the
       * quality scores from a vector of int.
       *
       * @param s The Sequence object
       * @param sc The quality scores
       *
       * @throw DimensionException if the number of quality values is not equal
       * to the number of sequence states
       */
      SequenceWithQuality(
          const Sequence& s,
          const std::vector<int>& sc)
        throw (DimensionException):
        SequenceWithAnnotation(s),
        qualScores_(new SequenceQuality(sc, false))
      {
        addAnnotation(qualScores_);
      }

      /** @} */

      /**
       * @name Destructor
       * @{
       */
      virtual ~SequenceWithQuality() {}
      /** @} */

      SequenceWithQuality(const SequenceWithQuality& sequence) : 
        SequenceWithAnnotation(sequence), qualScores_(0)
      {
        qualScores_ = dynamic_cast<SequenceQuality*>(&getAnnotation(SequenceQuality::QUALITY_SCORE));                  
      }

      SequenceWithQuality& operator=(const SequenceWithQuality& sequence)
      { 
        SequenceWithAnnotation::operator=(sequence);
        qualScores_ = dynamic_cast<SequenceQuality*>(&getAnnotation(SequenceQuality::QUALITY_SCORE));
        return *this;
      }

      /**
       * @name The Clonable interface
       * @{
       */
#ifdef NO_VIRTUAL_COV
      Clonable*
#else
      SequenceWithQuality*
#endif
      clone() const { return new SequenceWithQuality(*this); }
      /** @} */

      /**
       * @name Dealing with quality
       * @{
       */

      /**
       * @brief Set the quality score
       *
       * @param pos The position where the quality must be set
       * @param quality The quality value
       *
       * @throw IndexOutOfBoundsException if pos is greater than the
       * sequence size
       */
      void setQuality(size_t pos, int quality) throw (IndexOutOfBoundsException) {
        //if (pos >= qualScores_->getSize())
        //  throw IndexOutOfBoundsException("SequenceWithQuality::setQuality: pos out of bounds", pos, 0, qualScores_->getSize() - 1);
        //qualScores_[pos] = quality;
        qualScores_->setScore(pos, quality);
      }
      
      /**
       * @brief Get the quality score
       *
       * @param pos The position where the quality is read
       *
       * @return The quality score
       *
       * @throw IndexOutOfBoundsException if pos is greater than the
       * sequence size
       */
      int getQuality(size_t pos) const throw (IndexOutOfBoundsException) {
        if (pos >= qualScores_->getSize())
          throw IndexOutOfBoundsException("SequenceWithQuality::getQuality: pos out of bounds", pos, 0, qualScores_->getSize() - 1);
        return (*qualScores_)[pos];
      }

      /**
       * @brief Set the whole quality scores
       *
       * @param quality The vector of quality scores
       *
       * @throw DimensionException if the quality vector does not feet the
       * sequence size
       */
      void setQualities(const std::vector<int>& quality) throw (DimensionException) {
        if (quality.size() != qualScores_->getSize())
          throw DimensionException("SequenceWithQuality::setQualities: quality must fit sequence size", quality.size(), qualScores_->getSize());
        qualScores_->setScores(quality);
      }

      /**
       * @brief Get the whole quality scores
       *
       * @return A reference to the quality vector
       */
      const std::vector<int>& getQualities() const {
        return qualScores_->getScores();
      }

      void append(const std::vector<int>& content)
        throw (BadIntException)
      {
        SequenceWithAnnotation::append(content);
      }

      /**
       * @brief Append content with quality
       *
       * @param content A vector of int to append to the sequence
       * @param qualities A vector of int to append to the qualities
       *
       * @throw BadIntException if one of the content int is not in the
       * Alphabet
       * @throw DimensionException if qualities does not have the same size as
       * content
       */
      void append(
          const std::vector<int>& content,
          const std::vector<int>& qualities)
        throw (BadIntException, DimensionException)
      {
        if (content.size() != qualities.size())
          throw DimensionException("SequenceWithQuality::append: qualities must fit content size", qualities.size(), content.size());
        size_t pos = qualScores_->getSize();
        append(content); //This automatically extend scores array with default values through the listener
        //Update scores:
        qualScores_->setScores(pos, qualities);
      }

      void append(const std::vector<std::string>& content)
        throw (BadCharException)
      {
        SequenceWithAnnotation::append(content);
      }

      /**
       * @brief Append content with quality
       *
       * @param content A vector of string to append to the sequence
       * @param qualities A vector of int to append to the qualities
       *
       * @throw BadCharException if one of the content string is not in the
       * Alphabet
       * @throw DimensionException if qualities does not have the same size as
       * content
       */
      void append(
          const std::vector<std::string>& content,
          const std::vector<int>& qualities)
        throw (BadCharException, DimensionException)
      {
        if (content.size() != qualities.size())
          throw DimensionException("SequenceWithQuality::append: qualities must fit content size", qualities.size(), content.size());
        size_t pos = qualScores_->getSize();
        SequenceWithAnnotation::append(content); //This automatically extend scores array with default values through the listener
        //Update scores:
        qualScores_->setScores(pos, qualities);
      }

      void append(const std::string& content)
        throw (BadCharException)
      {
        SequenceWithAnnotation::append(content);
      }

      /**
       * @brief Append content with quality
       *
       * @param content A string to append to the sequence
       * @param qualities A vector of int to append to the qualities
       *
       * @throw BadCharException if one of the character of the string is not in
       * the Alphabet
       * @throw DimensionException if qualities does not have the same size as
       * content
       */
      void append(
          const std::string& content,
          const std::vector<int>& qualities)
        throw (BadCharException, DimensionException)
      {
        if (content.size() / this->getAlphabet()->getStateCodingSize()
            != qualities.size())
          throw DimensionException("SequenceWithQuality::append: qualities must fit content size", qualities.size(), content.size() / this->getAlphabet()->getStateCodingSize());
        size_t pos = qualScores_->getSize();
        SequenceWithAnnotation::append(content); //This automatically extend scores array with default values through the listener
        //Update scores:
        qualScores_->setScores(pos, qualities);
      }

      void addElement(
          const std::string& c)
        throw (BadCharException)
      {
        SequenceWithAnnotation::addElement(c);
      }

      /**
       * @brief Add a character to the end of the list with quality
       *
       * @param c The element to add to the sequence
       * @param q The quality of this element
       *
       * @throw BadCharException if one of the character of the string is not in
       * the Alphabet
       */
      void addElement(
          const std::string& c, int q)
        throw (BadCharException)
      {
        SequenceWithAnnotation::addElement(c);
        qualScores_->setScore(size() - 1, q);
      }

      void addElement(size_t pos, const std::string& c)
        throw (BadCharException, IndexOutOfBoundsException)
      {
        SequenceWithAnnotation::addElement(pos, c);
      }

      /**
       * @brief Add a character to a certain position in the list with quality
       *
       * @param pos The position where the element will be inserted
       * @param c The element to add to the sequence
       * @param q The quality of this element
       *
       * @throw BadCharException if one of the character of the string is not in
       * the Alphabet
       * @throw IndexOutOfBoundsException if pos is greater than the sequence
       * size
       */
      void addElement(
          size_t pos, const std::string& c, int q)
        throw (BadCharException, IndexOutOfBoundsException)
      {
        SequenceWithAnnotation::addElement(pos, c);
        qualScores_->setScore(pos, q);
      }

      void addElement(int v)
        throw (BadIntException)
      {
        SequenceWithAnnotation::addElement(v);
      }

      /**
       * @brief Add a character to the end of the list with quality
       *
       * @param v The element to add to the sequence
       * @param q The quality of this element
       *
       * @throw BadIntException if the value does not match the current Alphabet
       */
      void addElement(int v, int q)
        throw (BadIntException)
      {
        SequenceWithAnnotation::addElement(v);
        qualScores_->setScore(size() - 1, q);
      }

      void addElement(size_t pos, int v)
        throw (BadIntException, IndexOutOfBoundsException)
      {
        SequenceWithAnnotation::addElement(pos, v);
      }

      /**
       * @brief Add a character to a certain position in the list with quality
       *
       * @param pos The position where the element will be inserted
       * @param v The element to add to the sequence
       * @param q The quality of this element
       *
       * @throw BadIntException if the value does not match the current Alphabet
       * @throw IndexOutOfBoundsException if pos is greater than the sequence
       * size
       */
      void addElement(size_t pos, int v, int q)
        throw (BadCharException, IndexOutOfBoundsException)
      {
        SequenceWithAnnotation::addElement(pos, v);
        qualScores_->setScore(pos, q);
      }

      /** @} */

  };

} // end of namespace bpp.

#endif // _SEQUENCEWITHQUALITY_H_

