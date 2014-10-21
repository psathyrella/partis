//
// File: AbstractDiscreteDistribution.h
// Created by: Julien Dutheil
// Created on: ?
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 19, 2004)

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

#ifndef _ABSTRACTDISCRETEDISTRIBUTION_H_
#define _ABSTRACTDISCRETEDISTRIBUTION_H_

#include "DiscreteDistribution.h"
#include "../Constraints.h"
#include "../AbstractParameterAliasable.h"

#include <map>

namespace bpp
{

  /**
   * @brief Partial implementation of the DiscreteDistribution interface.
   *
   * This class uses a map to store the cateogry values as keys and probabilities as values.
   * It uses its own comparator class to deal with double precision.
   * By default, category values that differ less than 10E-9 will be considered identical.
   */
  class AbstractDiscreteDistribution :
    public virtual DiscreteDistribution,
    public virtual AbstractParameterAliasable
  {

  public:
    
    /**
     * @brief Comparator class for AbstractDiscreteDistribution.
     */
    class Order
    {
    private:
      double precision_;
      
    public:

      Order(double prec=NumConstants::TINY()):
        precision_(prec)
      {}

      Order(const Order& ord) :
        precision_(ord.precision_)
      {}

      Order& operator=(const Order& ord)
      {
        precision_=ord.precision_;
        return *this;
      }

      double precision() const {
        return precision_;
      }

      void setPrecision(double prec) {
        precision_=prec;
      }
      
      bool operator() (double l1, double l2) const
      {
        return (l1 < l2 - precision_); 
      }

    };

  protected:

    /*
     * The number of categories
     */

    size_t numberOfCategories_;  
    /**
     * These fields must be initialized in the constructor of the derived classes.
     */
    std::map<double, double, Order> distribution_;
    
    std::vector<double> bounds_;

    /**
     * @brief the interval where the distribution is defined/restricted.
     *
     */
    
    IntervalConstraint intMinMax_;    

    /**
     * Tells if the values in the classes is associated to the median or not (default: false)
     *
     */

    bool median_;
    
  public:
    AbstractDiscreteDistribution(size_t nbClasses, const std::string& prefix = ""); 

    /**
     * With additional precision value to discriminate categories (default 1e-12)
     *
     */
    AbstractDiscreteDistribution(size_t nbClasses, double precision, const std::string& prefix = "");
    
    AbstractDiscreteDistribution(const AbstractDiscreteDistribution& adde);

    AbstractDiscreteDistribution& operator=(const AbstractDiscreteDistribution& adde);
    
    virtual ~AbstractDiscreteDistribution() {}

    
  public:

    /**
     * @name The DiscreteDistribution interface.
     *
     * @{
     */
    size_t getNumberOfCategories() const;
    void setNumberOfCategories(size_t nbClasses);
    double getCategory(size_t categoryIndex) const;
    double getProbability(size_t categoryIndex) const;
    double getProbability(double category) const;
    Vdouble getCategories() const;
    Vdouble getProbabilities() const;
    double getValueCategory(double value) const;
    void set(double category, double probability);
    void add(double category, double probability);
    double getInfCumulativeProbability(double category) const;
    double getIInfCumulativeProbability(double category) const;
    double getSupCumulativeProbability(double category) const;
    double getSSupCumulativeProbability(double category) const;
    double rand() const;
    double randC() const throw (Exception) { throw Exception("AbstractDiscreteDistribution::randC. No continuous version available for this distribution."); }

    /*
     *@return value of the internal bound
     *
     */
    
    double getBound(size_t i) const throw (IndexOutOfBoundsException)
    {
      if (i >= numberOfCategories_ - 1)
        throw IndexOutOfBoundsException("AbstractDiscreteDistribution::getBound(i)", i , 0, numberOfCategories_-1);
      return bounds_[i];
    }  


    /*
     *@brief Information about the range of the distribution
     *
     */
     

    double getLowerBound() const
    {
      return intMinMax_.getLowerBound();
    }

    double getUpperBound() const
    {
      return intMinMax_.getUpperBound();
    }

    bool strictLowerBound() const
    {
      return intMinMax_.strictLowerBound();
    }

    bool strictUpperBound() const
    {
      return intMinMax_.strictUpperBound();
    }

    Vdouble getBounds() const;
    
    void print(OutputStream& out) const;

    double precision() const { return distribution_.key_comp().precision();}
      
    void setMedian(bool median) {
      if (median_ != median) {
        median_ = median;
        discretize();
      }
    }
    
    virtual void discretize();

    /** @} */

    /**
     * @brief Restricts the distribution to the domain where the
     * constraint is respected, in addition of other predefined
     * constraints.
     *
     * @param c The Constraint to respect.
     *
     */
    
    virtual void restrictToConstraint(const Constraint& c);
      

  };

} //end of namespace bpp.

#endif  //_ABSTRACTDISCRETEDISTRIBUTION_H_

