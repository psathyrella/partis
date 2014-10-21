// 
// File:    ProteicAlphabetState.h
// Author:  Sylvain Gaillard
// Created: 29/07/2009 13:56:01
// 

/*
Copyright or © or Copr. CNRS, (July 29, 2009)

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

#ifndef _PROTEICALPHABETSTATE_H_
#define _PROTEICALPHABETSTATE_H_

// From the STL
#include <string>

namespace bpp {
  /**
   * @brief This is the base class to describe states in a ProteicAlphabet.
   *
   * @author Sylvain Gaillard
   */
  class ProteicAlphabetState: public AlphabetState {
    private:
      std::string abbr_;

    public:
      ProteicAlphabetState(int num, const std::string & letter, const std::string & abbr, const std::string & name): AlphabetState(num, letter, name), abbr_(abbr) {}

      // Class destructor
      virtual ~ProteicAlphabetState() {}

    public:
      ProteicAlphabetState * clone() const {
        return new ProteicAlphabetState(* this);
      }
      /**
       * @brief Get the state's abbreviation.
       *
       * @return The state's abbreviation.
       */
      const std::string & getAbbreviation() const { return abbr_; }
      /**
       * @brief Set the state's abbreviation.
       *
       * @param abbr The state's abbreviation.
       */
      void setAbbreviation(const std::string & abbr) { abbr_ = abbr; }
  };
}

#endif // _PROTEICALPHABETSTATE_H_

