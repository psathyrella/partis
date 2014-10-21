// 
// File:    AlphabetState.h
// Author:  Sylvain Gaillard
// Created: 29/07/2009 13:56:01
// 

/*
Copyright or © or Copr. Bio++ Development Team, (July 29, 2009)

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

#ifndef _ALPHABETSTATE_H_
#define _ALPHABETSTATE_H_

#include <Bpp/Clonable.h>

// From the STL
#include <string>

namespace bpp {
  /**
   * @brief This is the base class to describe states in an Alphabet.
   *
   * @author Sylvain Gaillard
   */
  class AlphabetState: public virtual Clonable {
    private:
      int num_;
      std::string letter_;
      std::string name_;

    public:
      AlphabetState(int num, const std::string& letter, const std::string& name): num_(num), letter_(letter), name_(name) {}

      // Class destructor
      virtual ~AlphabetState() {}

    public:
      /**
       * @name The Clonable interface.
       * @{
       */
#ifdef NO_VIRTUAL_COV
      Clonable*
#else
      AlphabetState*
#endif
      clone() const { return new AlphabetState(* this); }
      /** @} */
      /**
       * @brief Get the state's number.
       *
       * @return The state's number (i.e. -1 for gap (-)).
       */
      int getNum() const { return num_; }
      /**
       * @brief Set the state's number.
       *
       * @param num The state's number.
       */
      void setNum(int num) { num_ = num; }

      /** 
       * @brief Get the letter(s) corresponding to the state.
       *
       * The letter is a string because it may more than one char
       * (for instance: codon).
       *
       * @return The state's letter.
       */
      const std::string& getLetter() const { return letter_; }
      /**
       * @brief Set the letter(s) of the state.
       *
       * @param letter The state's letter.
       */
      void setLetter(const std::string& letter) { letter_ = letter; }

      /**
       * @brief Get the name of the state.
       *
       * @return The full name of the state (i.e. Adenine).
       */
      const std::string& getName() const { return name_; }
      /**
       * @brief Set the name of the state.
       *
       * @param name The state's name
       */
      void setName(const std::string& name) { name_ = name; }

      /**
       * @brief operator ==
       *
       * Comparison is done on state num
       */
      bool operator == (AlphabetState& l2) {
        return getNum() == l2.getNum();
      }
  };
}

#endif // _ALPHABETSTATE_H_

