//
// File: StateMap.h
// Created by: Julien Dutheil
// Created on: Wed Jun 13 15:03 2012
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

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

#ifndef _STATEMAP_H_
#define _STATEMAP_H_

#include <Bpp/Clonable.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>
#include <Bpp/Numeric/VectorTools.h>

//From the STL:
#include <vector>
#include <string>

namespace bpp
{

  /**
   * @brief Map the states of a given alphabet which have a model state.
   */
  class StateMap:
    public virtual Clonable
  {
    public:
      virtual ~StateMap() {}
      virtual StateMap* clone() const = 0;

    public:
      /**
       * @return The associated alphabet.
       */
      virtual const Alphabet* getAlphabet() const = 0;

      /**
       * @return The number of states supported by the model.
       */
      virtual size_t getNumberOfStates() const = 0;

      /**
       * @param index The model state.
       * @return The corresponding alphabet character as int code.
       */
      virtual int stateAsInt(size_t index) const = 0;

      /**
       * @param index The model state.
       * @return The corresponding alphabet character as character code.
       */
      virtual std::string stateAsChar(size_t index) const = 0;

      /**
       * @param code The int code of the character to check.
       * @return The corresponding model state, is any.
       */
      virtual size_t whichState(int code) const = 0;
      
      /**
       * @param code The character code of the character to check.
       * @return The corresponding model state, is any.
       */
      virtual size_t whichState(const std::string& code) const = 0;
  };

  /**
   * @brief A convenience partial implementation of the StateMap interface.
   *
   * Model states are stored as their corresponding int codes, stored in a vector 'states_'.
   * This vector has to be initialized and filled by the derived class.
   */
  class AbstractStateMap:
    public virtual StateMap
  {
    protected:
      const Alphabet* alphabet_;
      std::vector<int> states_;

    public:
      AbstractStateMap(const Alphabet* alphabet):
        alphabet_(alphabet),
        states_()
      {}

      AbstractStateMap(const AbstractStateMap& absm):
        alphabet_(absm.alphabet_),
        states_(absm.states_)
      {}

      AbstractStateMap& operator=(const AbstractStateMap& absm)
      {
        alphabet_ = absm.alphabet_;
        states_ = absm.states_;
        return *this;
      }

    public:
      virtual const Alphabet* getAlphabet() const { return alphabet_; }
      virtual size_t getNumberOfStates() const { return states_.size(); }
      virtual int stateAsInt(size_t index) const { return states_[index]; }
      virtual std::string stateAsChar(size_t index) const { return alphabet_->intToChar(states_[index]); }
      virtual size_t whichState(int code) const throw (Exception) {
        try { return VectorTools::which(states_, code); }
        catch (ElementNotFoundException<int>& ex) { throw Exception("AbstractStateMap::whichState. Unsupported alphabet char: " + code); }
      }
      virtual size_t whichState(const std::string& code) const throw (Exception) {
        try { return VectorTools::which(states_, alphabet_->charToInt(code)); }
        catch (ElementNotFoundException<int>& ex) { throw Exception("AbstractStateMap::whichState. Unsupported alphabet char: " + code); }
      }

  };

  /**
   * @brief This class implements a state map where all resolved states are modeled.
   *
   * For nucleotides, the underlying states are for instance: A (0), C (1), G (2), T/U (3).
   * Optionally, gaps can be modeled.
   */
  class CanonicalStateMap:
    public AbstractStateMap
  {
    public:
      CanonicalStateMap(const Alphabet* alphabet, bool includeGaps);
      virtual CanonicalStateMap* clone() const { return new CanonicalStateMap(*this); }

  };

}// end of namespace bpp

#endif //_STATEMAP_H_

