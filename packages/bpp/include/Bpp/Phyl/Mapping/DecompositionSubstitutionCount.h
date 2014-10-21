//
// File: DecompositionSubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Thu Mar 17 16:08 2011
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

#ifndef _DECOMPOSITIONSUBSTITUTIONCOUNT_H_
#define _DECOMPOSITIONSUBSTITUTIONCOUNT_H_

#include "WeightedSubstitutionCount.h"

#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{

/**
 * @brief Analytical substitution count using the eigen decomposition method.
 *
 * The codes is adapted from the original R code by Paula Tataru and Asger Hobolth.
 * Only reversible models are supported for now.
 *
 * @author Julien Dutheil
 */
class DecompositionSubstitutionCount:
  public AbstractSubstitutionCount,
  public AbstractWeightedSubstitutionCount
{
	private:
		const ReversibleSubstitutionModel* model_;
    size_t nbStates_;
		mutable RowMatrix<double> jMat_, v_, vInv_;
    mutable std::vector<double> lambda_;
    std::vector< RowMatrix<double> > bMatrices_, insideProducts_;
    mutable std::vector< RowMatrix<double> > counts_;
    mutable double currentLength_;
	
	public:
		DecompositionSubstitutionCount(const ReversibleSubstitutionModel* model, SubstitutionRegister* reg, const AlphabetIndex2* weights = 0);
		
    DecompositionSubstitutionCount(const DecompositionSubstitutionCount& dsc) :
      AbstractSubstitutionCount(dsc),
      AbstractWeightedSubstitutionCount(dsc),
      model_(dsc.model_),
      nbStates_(dsc.nbStates_),
      jMat_(dsc.jMat_),
      v_(dsc.v_),
      vInv_(dsc.vInv_),
      lambda_(dsc.lambda_),
      bMatrices_(dsc.bMatrices_),
      insideProducts_(dsc.insideProducts_),
      counts_(dsc.counts_),
      currentLength_(dsc.currentLength_)
    {}				
    
    DecompositionSubstitutionCount& operator=(const DecompositionSubstitutionCount& dsc)
    {
      AbstractSubstitutionCount::operator=(dsc);
      AbstractWeightedSubstitutionCount::operator=(dsc);
      model_          = dsc.model_;
      nbStates_       = dsc.nbStates_;
      jMat_           = dsc.jMat_;
      v_              = dsc.v_;
      vInv_           = dsc.vInv_;
      lambda_         = dsc.lambda_;
      bMatrices_      = dsc.bMatrices_;
      insideProducts_ = dsc.insideProducts_;
      counts_         = dsc.counts_;
      currentLength_  = dsc.currentLength_;
      return *this;
    }				
		
    virtual ~DecompositionSubstitutionCount() {}
		
    DecompositionSubstitutionCount* clone() const { return new DecompositionSubstitutionCount(*this); }

	public:
		double getNumberOfSubstitutions(int initialState, int finalState, double length, size_t type = 1) const;

    Matrix<double>* getAllNumbersOfSubstitutions(double length, size_t type = 1) const;
    
    std::vector<double> getNumberOfSubstitutionsForEachType(int initialState, int finalState, double length) const;
   
    /**
     * @brief Set the substitution model.
     *
     * @param model A pointer toward the substitution model to use. Only reversible models are currently supported. Setting a non-reversible model will throw an exception.
     */
    void setSubstitutionModel(const SubstitutionModel* model);

  protected:
    void computeCounts_(double length) const;
    void jFunction_(const std::vector<double>& lambda, double t, RowMatrix<double>& result) const;
    void substitutionRegisterHasChanged() throw (Exception);
    void weightsHaveChanged() throw (Exception);

  private:
    void resetStates_();
    void resetBMatrices_();
    void initBMatrices_();
    void fillBMatrices_();
    void computeEigen_();
    void computeProducts_();
};

} //end of namespace bpp.

#endif // _DECOMPOSITIONSUBSTITUTIONCOUNT_H_

