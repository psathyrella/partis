//
// File: InvariantMixedDiscreteDistribution.h
// Created by: Julien Dutheil
// Created on: Mon Dec 24 12:02 2007
//

/*
  Copyright or © or Copr. CNRS, (November 17, 2004)

  This software is a computer program whose purpose is to provide classes
  for numerical calculus.

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

#ifndef _INVARIANTMIXEDDISCRETEDISTRIBUTION_H_
#define _INVARIANTMIXEDDISCRETEDISTRIBUTION_H_

#include "AbstractDiscreteDistribution.h"

namespace bpp
{

  /**
   * @brief Discrete mixed distribution, with a one-category fixed value (called "invariant") and a user-specified multi-categories distribution.
   *
   * The term "invariant" comes from the use of such distributions in phylogenetics:
   * the fixed category corresponds to a value of 0 and describes invariant positions in an alignment.
   */
  class InvariantMixedDiscreteDistribution:
    public AbstractDiscreteDistribution
  {
  private:
    DiscreteDistribution* dist_;
    double invariant_, p_;
    std::string nestedPrefix_;

  public:
    /**
     * @brief Build a new InvariantMixedDiscreteDistribution object.
     *
     * @param dist            The distribution to use. The mixed distribution will "own" this distribution object.
     * This means that the distribution will be cloned in case of copy of this instance, and will be
     * deleted with this instance.
     * @param p               The probability of being in the invariant category.
     * @param invariant       The value of the invariant category (typically 0, but other values may be specified). 
     */
    InvariantMixedDiscreteDistribution(DiscreteDistribution* dist, double p, double invariant = 0.);
    
    virtual ~InvariantMixedDiscreteDistribution()
    {
      delete dist_;
    }
    
    InvariantMixedDiscreteDistribution(const InvariantMixedDiscreteDistribution& imdd):
      AbstractParameterAliasable(imdd),
      AbstractDiscreteDistribution(imdd),
      dist_(dynamic_cast<DiscreteDistribution *>(imdd.dist_->clone())),
      invariant_(imdd.invariant_), p_(imdd.p_),
      nestedPrefix_(imdd.nestedPrefix_)
    {}
    
    InvariantMixedDiscreteDistribution& operator=(const InvariantMixedDiscreteDistribution& imdd)
    {
      AbstractParameterAliasable::operator=(imdd);
      AbstractDiscreteDistribution::operator=(imdd);
      dist_         = dynamic_cast<DiscreteDistribution *>(imdd.dist_->clone());
      invariant_    = imdd.invariant_;
      p_            = imdd.p_;
      nestedPrefix_ = imdd.nestedPrefix_;
      return *this;
    }

    InvariantMixedDiscreteDistribution* clone() const { return new InvariantMixedDiscreteDistribution(*this); }

  public:

    std::string getName() const {return("Invariant");}

    void fireParameterChanged(const ParameterList & parameters);

    void setNamespace(const std::string& prefix);

    /*
     *@ brief Sets the number of categories of the embedded
     *distribution. The number of categories of this class equals this
     *plus one if the invariant_ is not in the embedded distribution
     *values.
     *
     */
    
    void setNumberOfCategories(unsigned int nbClasses) {
      dist_->setNumberOfCategories(nbClasses);
      updateDistribution();
    }
    
    /**
     * @return The nested, conditional, sub-distribution.
     */
    const DiscreteDistribution * getVariableSubDistribution() const { return dist_; }

    double qProb(double x) const{
      return (x>=p_+(1-p_)*dist_->pProb(invariant_))?dist_->qProb((x-p_)/(1-p_)):dist_->qProb(x/(1-p_));
    }
     
    double pProb(double x) const{
      return (1-p_)*dist_->pProb(x)+(x<invariant_)?0:p_;
    }
                                
    double Expectation(double a) const{
      return (1-p_)*dist_->Expectation(a)+(a<invariant_?0:p_);
    }

    void setMedian(bool median){
      if (median_!=median){
        median_=median;
        dist_->setMedian(median);
        updateDistribution();
      }
    }
  
    void discretize(){
      dist_->discretize();
      updateDistribution();
    }

    void restrictToConstraint(const Constraint& c);

  protected:
    void updateDistribution();


  };

} //end of namespace bpp.

#endif //_INVARIANTMIXEDDISCRETEDISTRIBUTION_H_

