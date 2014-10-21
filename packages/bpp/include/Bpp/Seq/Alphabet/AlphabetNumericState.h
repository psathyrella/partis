// 
// File:    AlphabetState.h
// Author:  Laurent Guéguen
// Created: 03/2010
// 

/*
  Copyright or © or Copr. Bio++ Development Team, (March, 2010)

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

#ifndef _ALPHABETNUMERICSTATE_H_
#define _ALPHABETNUMERICSTATE_H_

// From the STL
#include <string>

// From bpp-core
#include <Bpp/Clonable.h>

// From bpp-seq

#include "AlphabetState.h"

namespace bpp {
  /**
   * @brief States that do have a double value
   *
   * @author Laurent Guéguen
   */
  class AlphabetNumericState: public AlphabetState {
  private:
    double value_;
    
  public:
    AlphabetNumericState(int num, double value, const std::string& letter, const std::string& name): AlphabetState(num, letter, name), value_(value) {}
    
    // Class destructor
    virtual ~AlphabetNumericState() {}

  public:
    /**
     * @name The Clonable interface.
     * @{
     */
#ifdef NO_VIRTUAL_COV
    Clonable*
#else
    AlphabetNumericState*
#endif
    clone() const { return new AlphabetNumericState(* this); }
    /** @} */
    /**
     * @brief Get the state value.
     *
     * @return The state value
     */
    double getValue() const { return value_; }

    /**
     * @brief Set the state value.
     * @param value The given state value.
     */
    void setValue(double value) { value_ = value; }

  };
}

#endif // _ALPHABETNUMERICSTATE_H_

