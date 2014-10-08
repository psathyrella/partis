//
// File: HmmTransitionMatrix.h
// Created by: Julien Dutheil
// Created on: Sat Oct 27 10:24 2007
//

/*
  Copyright or © or Copr. Bio++Development Team, (November 16, 2004)

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

#ifndef _HMMTRANSITIONMATRIX_H_
#define _HMMTRANSITIONMATRIX_H_

#include "HmmStateAlphabet.h"
#include "HmmExceptions.h"
#include "../Matrix/Matrix.h"
#include "../Parametrizable.h"

//Fron the STL:
#include <vector>

namespace bpp
{

/**
 * @brief Describe the transition probabilities between hidden states of a Hidden Markov Model.
 *
 * This class is part of the HMM framework.
 */
  class HmmTransitionMatrix:
    public virtual Parametrizable
  {
  public:

    /**
     * @return The hidden alphabet associated to this model.
     */

    virtual const HmmStateAlphabet* getHmmStateAlphabet() const = 0;

    /**
     * @brief Set the new hidden state alphabet.
     * @param stateAlphabet The new state alphabet.
     * @throw UnvalidStateAlphabetException if the new alphabet is uncorrect (for instance is NULL pointer).
     */
    virtual void setHmmStateAlphabet(const HmmStateAlphabet* stateAlphabet) throw (HmmUnvalidAlphabetException) = 0;
   
    /**
     * @return The number of states in the model.
     */
    virtual unsigned int getNumberOfStates() const = 0;

    /**
     * @brief Get the transition probability between two states.
     *
     * @param i initial state.
     * @param j final state.
     * @return the transition probability between the two states.
     */
    virtual double Pij(unsigned int i, unsigned int j) const = 0;

    /**
     * @brief Get all transition probabilities as a matrix.
     *
     * @return A n*n matrix will all transition probabilities (n being the number of hidden states).
     */
    virtual const Matrix<double>& getPij() const = 0;

    /**
     * @return The vector of equilibrium frequencies of the Markov chain described by the matrix.
     */
    virtual const std::vector<double>& getEquilibriumFrequencies() const = 0;

  };

} //end of namespace bpp

#endif //_HMMTRANSITIONMATRIX_H_

