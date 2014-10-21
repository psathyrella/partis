//
// File: GlobalClockTreeLikelihoodFunctionWrapper.h
// Created by: Julien Dutheil
// Created on: Thu Jul 14 10:53 2011
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

#ifndef _GLOBALCLOCKTREELIKELIHOODFUNCTIONWRAPPER_H_
#define _GLOBALCLOCKTREELIKELIHOODFUNCTIONWRAPPER_H_

#include "TreeLikelihood.h"

namespace bpp
{

class GlobalClockTreeLikelihoodFunctionWrapper:
  public virtual DerivableSecondOrder,
  public AbstractParametrizable
{
  private:
    TreeLikelihood* tl_;

  public:
    GlobalClockTreeLikelihoodFunctionWrapper(TreeLikelihood* tl):
      AbstractParametrizable(""),
      tl_(tl)
    {
      initParameters_();
    }

    GlobalClockTreeLikelihoodFunctionWrapper(const GlobalClockTreeLikelihoodFunctionWrapper& gctlfw):
      AbstractParametrizable(gctlfw), tl_(gctlfw.tl_)
    {}
    
    GlobalClockTreeLikelihoodFunctionWrapper& operator=(const GlobalClockTreeLikelihoodFunctionWrapper& gctlfw) {
      AbstractParametrizable::operator=(gctlfw);
      tl_ = gctlfw.tl_;
      return *this;
    }

    GlobalClockTreeLikelihoodFunctionWrapper* clone() const { return new GlobalClockTreeLikelihoodFunctionWrapper(*this); }

  public:
    void setParameters(const ParameterList& pl) throw (Exception) {
      //For now we go the hard way and recompute everything:
      matchParametersValues(pl);
    }

    double getValue() const throw (Exception) { return tl_->getValue(); }

    void fireParameterChanged(const bpp::ParameterList& pl);

    void enableSecondOrderDerivatives(bool yn) { tl_->enableSecondOrderDerivatives(yn); }
    bool enableSecondOrderDerivatives() const { return tl_->enableSecondOrderDerivatives(); }
    void enableFirstOrderDerivatives(bool yn) { tl_->enableFirstOrderDerivatives(yn); }
    bool enableFirstOrderDerivatives() const { return tl_->enableFirstOrderDerivatives(); }
    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const throw (Exception) { return tl_->getSecondOrderDerivative(variable1, variable2); }
    double getSecondOrderDerivative(const std::string& variable) const throw (Exception) { return tl_->getSecondOrderDerivative(variable); }
    double getFirstOrderDerivative(const std::string& variable) const throw (Exception) { return tl_->getFirstOrderDerivative(variable); } 

    ParameterList getHeightParameters() const;

  private:
    void initParameters_();
    void computeBranchLengthsFromHeights_(const Node* node, double height, ParameterList& brlenPl) throw (Exception);

};

} // end of namespace bpp.

#endif //_GLOBALCLOCKTREELIKELIHOODFUNCTIONWRAPPER_H_

