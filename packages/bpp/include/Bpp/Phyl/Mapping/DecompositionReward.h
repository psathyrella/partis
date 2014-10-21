//
// File: DecompositionReward.h
// Created by: Laurent Guéguen
// Created on: mercredi 27 mars 2013, à 12h 29
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#ifndef _DECOMPOSITIONREWARD_H_
#define _DECOMPOSITIONREWARD_H_

#include "Reward.h"

#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{

/**
 * @brief Analytical reward using the eigen decomposition method.
 *
 * The codes is adapted from the original R code by Paula Tataru and
 * Asger Hobolth to the formula in the article of Minin & Suchard.
 *
 * Minin, V.N. and Suchard, M.A., 
 * Fast, accurate and simulation-free stochastic mapping
 * Philosophical Transactions of the Royal Society B 2008 363:3985-95.
 *
 * Only reversible models are supported for now.
 *
 * @author Laurent Guéguen
 */
  
class DecompositionReward:
  public AbstractReward
{
private:
  const ReversibleSubstitutionModel* model_;
  size_t nbStates_;
  mutable RowMatrix<double> jMat_, v_, vInv_;
  mutable std::vector<double> lambda_;
  RowMatrix<double> bMatrice_, insideProduct_;
  mutable RowMatrix<double> rewards_;
  mutable double currentLength_;
	
public:
  DecompositionReward(const SubstitutionModel* model, AlphabetIndex1* alphIndex);
		
  DecompositionReward(const DecompositionReward& dr) :
    AbstractReward(dr),
    model_(dr.model_),
    nbStates_(dr.nbStates_),
    jMat_(dr.jMat_),
    v_(dr.v_),
    vInv_(dr.vInv_),
    lambda_(dr.lambda_),
    bMatrice_(dr.bMatrice_),
    insideProduct_(dr.insideProduct_),
    rewards_(dr.rewards_),
    currentLength_(dr.currentLength_)
  {}
    
  DecompositionReward& operator=(const DecompositionReward& dr)
  {
    AbstractReward::operator=(dr);
    model_          = dr.model_;
    nbStates_       = dr.nbStates_;
    jMat_           = dr.jMat_;
    v_              = dr.v_;
    vInv_           = dr.vInv_;
    lambda_         = dr.lambda_;
    bMatrice_       = dr.bMatrice_;
    insideProduct_  = dr.insideProduct_;
    rewards_        = dr.rewards_;
    currentLength_  = dr.currentLength_;
    return *this;
  }				
		
  virtual ~DecompositionReward() {}
		
  DecompositionReward* clone() const { return new DecompositionReward(*this); }

public:
  double getReward(int initialState, int finalState, double length) const;

  Matrix<double>* getAllRewards(double length) const;
    
  /**
   * @brief Set the substitution model.
   *
   * @param model A pointer toward the substitution model to use. Only
   * reversible models are currently supported. Setting a
   * non-reversible model will throw an exception.
   *
   */
  void setSubstitutionModel(const SubstitutionModel* model);

protected:
  void computeRewards_(double length) const;
  void jFunction_(const std::vector<double>& lambda, double t, RowMatrix<double>& result) const;
  void alphabetIndexHasChanged() throw (Exception);

private:
  void resetStates_();
  void computeBMatrice_();
  void computeEigen_();
  void computeProducts_();
};

} //end of namespace bpp.

#endif // _DECOMPOSITIONREWARD_H_

