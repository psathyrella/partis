//
// File: OneJumpSubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Tue Nov 24 17:00 2009
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

#ifndef _ONEJUMPSUBSTITUTIONCOUNT_H_
#define _ONEJUMPSUBSTITUTIONCOUNT_H_

#include "SubstitutionCount.h"
#include "../Model/SubstitutionModel.h"

namespace bpp
{

/**
 * @brief Computes the probability that at least one jump occured on a branch, given the initial and final state.
 *
 * This probability is defined as
 * @f[
 * p_{x,y}(l) = \left\{\begin{array}{ll}1 & \mathrm{if} x \neq y \\ 1 - \exp{\left(Q \cdot t\right)}_{x,y} & \mathrm{otherwise.}\end{array}\right.
 * @f]
 *
 * @author Julien Dutheil
 */
class OneJumpSubstitutionCount:
  public AbstractSubstitutionCount
{
	private:
		const SubstitutionModel* model_;
		mutable RowMatrix<double> tmp_;
	
	public:
		OneJumpSubstitutionCount(const SubstitutionModel* model) :
      AbstractSubstitutionCount(new TotalSubstitutionRegister(model->getAlphabet())),
      model_(model), tmp_() {}
		
    OneJumpSubstitutionCount(const OneJumpSubstitutionCount& ojsc) :
      AbstractSubstitutionCount(ojsc),
      model_(ojsc.model_), tmp_(ojsc.tmp_) {}
				
    OneJumpSubstitutionCount& operator=(const OneJumpSubstitutionCount& ojsc)
    {
      AbstractSubstitutionCount::operator=(ojsc),
      model_    = ojsc.model_;
      tmp_      = ojsc.tmp_;
      return *this;
    }
				
		virtual ~OneJumpSubstitutionCount() {}
		
    virtual OneJumpSubstitutionCount* clone() const { return new OneJumpSubstitutionCount(*this); }

	public:
		double getNumberOfSubstitutions(int initialState, int finalState, double length, size_t type = 1) const
    {
      if (finalState != initialState) return 1.;
      else return 1. - model_->Pij_t(initialState, finalState, length);
    }

    Matrix<double>* getAllNumbersOfSubstitutions(double length, size_t type = 1) const;
    
    std::vector<double> getNumberOfSubstitutionsForEachType(int initialState, int finalState, double length) const
    {
      std::vector<double> v(0);
      v[0] = getNumberOfSubstitutions(initialState, finalState, length, 0);
      return v;
    }
    
    void setSubstitutionModel(const SubstitutionModel* model) { model_ = model; }

  /*
   *@param reg pointer to a SubstitutionRegister
   *
   */
        
    void setSubstitutionRegister(SubstitutionRegister* reg) throw (Exception) {
      throw Exception("OneJumpSubstitutionCount::setSubstitutionRegister. This SubstitutionsCount only works with a TotalSubstitutionRegister.");
    }

  private:
    void substitutionRegisterHasChanged() {}

};

} //end of namespace bpp.

#endif //_ONEJUMPSUBSTITUTIONCOUNT_H_

