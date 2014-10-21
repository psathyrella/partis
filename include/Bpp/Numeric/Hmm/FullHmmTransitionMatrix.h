//
// File: FullHmmTransitionMatrix.h
// Created by: Laurent Guéguen
// Created on: samedi 21 septembre 2013, à 00h 41
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

#ifndef _FULLHMMTRANSITIONMATRIX_H_
#define _FULLHMMTRANSITIONMATRIX_H_

#include "../Prob/Simplex.h"

#include "AbstractHmmTransitionMatrix.h"

#include "../AbstractParameterAliasable.h"

namespace bpp
{

/**
 * @brief Describe the transition probabilities between hidden states of a Hidden Markov Model.
 *
 *
 * Transition probabilities are described through simplexes.
 *
 *
 * Parameters are denoted \c "I.Simplex.thetaJ" where \c I is the line number
 * of the transition matrix, and \c "thetaJ" is the matching parameter
 * in the Ith Simplex. 
 *
 * @see Simplex
 */
  
class FullHmmTransitionMatrix:
  public virtual AbstractHmmTransitionMatrix,
  public AbstractParametrizable
{
private:
  std::vector<Simplex> vSimplex_;

public:

  FullHmmTransitionMatrix(const HmmStateAlphabet* alph, const std::string& prefix = "");

  FullHmmTransitionMatrix(const FullHmmTransitionMatrix& hptm);

  FullHmmTransitionMatrix& operator=(const FullHmmTransitionMatrix& hptm);

  FullHmmTransitionMatrix* clone() const { return new FullHmmTransitionMatrix(*this);}

  /**
   * @brief Set the matrix of the transition probabilities.
   *
   */
  
  void setTransitionProbabilities(const Matrix<double>& mat);
  
  /**
   * @brief Get the transition probability between two states.
   *
   * @param i initial state.
   * @param j final state.
   * @return the transition probability between the two states.
   */

  double Pij(unsigned int i, unsigned int j) const
  {
    return vSimplex_[i].prob(j);
  }

  /**
   * @brief Get all transition probabilities as a matrix.
   *
   * @return A n*n matrix will all transition probabilities (n being the number of hidden states).
   */

  const Matrix<double>& getPij() const;

  /**
   * @return The vector of equilibrium frequencies of the Markov chain described by the matrix.
   */

  const std::vector<double>& getEquilibriumFrequencies() const;


  /*
   * @brief From AbstractParametrizable interface
   *
   */

  void fireParameterChanged(const ParameterList& parameters);
  
};

} //end of namespace bpp

#endif //_FULLHMMTRANSITIONMATRIX_H_

