//
// File: HmmStateAlphabet.h
// Created by: Julien Dutheil
// Created on: Fri Oct 26 11:07 2007
//

/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _HMMSTATEALPHABET_H_
#define _HMMSTATEALPHABET_H_

#include "HmmExceptions.h"
#include "../Parametrizable.h"
#include "../../Exceptions.h"

//From the STL:
#include <vector>

//using namespace bpp;

namespace bpp {

  class StateListener;
  class StateChangedEvent;

  /**
   * @brief Hidden states alphabet.
   *
   * Implementations of this interface describe the set of hidden states of a Hidden Markov Model.
   */
  class HmmStateAlphabet:
    public virtual Parametrizable
  {
  public:
    HmmStateAlphabet() {}
    virtual ~HmmStateAlphabet() {}

  public:
    /**
     * @param stateIndex The index of a hidden state.
     * @return The corresponding hidden state.
     * @see getNumberOfStates
     */
    virtual const Clonable& getState(unsigned int stateIndex) const throw (HmmBadStateException) = 0;


    virtual unsigned int getNumberOfStates() const = 0;

    /**
     * @brief Tell if this instance can work with the instance of alphabet given as input.
     *
     * In many case, this will return true is the pointer provided as argument refers to this object.
     *
     * @param stateAlphabet The alphabet to check.
     * @return true if the matrix is compatible with the given alphabet.
     */
    virtual bool worksWith(const HmmStateAlphabet* stateAlphabet) const = 0;
 
  };

  class StateListener
  {
  public:
    StateListener() {}
    virtual ~StateListener() {}

  public:
    virtual void stateChanged(StateChangedEvent& event) = 0;
  };

  class StateChangedEvent
  {
  protected:
    std::vector<unsigned int> states_;

  public:
    StateChangedEvent(unsigned int stateIndex): states_(1)
    {
      states_[0] = stateIndex;
    }
    StateChangedEvent(std::vector<unsigned int>& states): states_(states) {}

  public:
    const std::vector<unsigned int>& getStates() const { return states_; }
    std::vector<unsigned int>& getStates() { return states_; }

  };

}
#endif //_HMMSTATEALPHABET_H_

