//
// File: DiscreteDistribution.h
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

#ifndef _DISCRETEDISTRIBUTION_H_
#define _DISCRETEDISTRIBUTION_H_

#include "../VectorTools.h"
#include "../ParameterAliasable.h"
#include "../NumConstants.h"
#include "../../Exceptions.h"
#include "../../Io/OutputStream.h"

namespace bpp
{

  /**
   * @brief Interface for discrete distribution objects.
   *
   * A discrete distribution usually contains a finite set of
   * categories and a probability associated to each.
   * 
   * Each category (or class) is defined by two bounds, and sometimes by
   * a mean or a median value.
   *
   * A discrete distribution may contain one or several parameters.
   * The probabilities associated to each class usually depends on
   * the parameter values.
   * In some cases, the number and/or bounds of the classes may also 
   * depend on the parameters values, depending on the kind of
   * discretization used.
   */
  class DiscreteDistribution:
    public virtual ParameterAliasable
  {
  public:
    DiscreteDistribution() {}
    virtual ~DiscreteDistribution() {}

  
#ifndef NO_VIRTUAL_COV
    DiscreteDistribution* clone() const = 0;
#endif

  public:
		
    /**
     * @brief Get the name of the distribution.
     *
     * @return The name of this distribution.
     */
    virtual std::string getName() const = 0;

    /**
     * @return The number of categories.
     */
    virtual size_t getNumberOfCategories() const = 0;

    
    /**
     * @brief sets the number of categories and discretizes if there
     * is a change in this number.
     */
    
    virtual void setNumberOfCategories(size_t nbClasses) = 0;
		
    /**
     * @param value 
     * @return The value of the category the value is in. Throws a
     * OutOfRangeException if the value is off the domain of the
     * DiscreteDistribution.
     */
    
    virtual double getValueCategory(double value) const = 0 ;
    
    /**
     * @param categoryIndex Class index.
     * @return The value associated to a given class.
     */

    virtual double getCategory(size_t categoryIndex) const = 0;
		
    /**
     * @param categoryIndex Class index.
     * @return The probability associated to a given class.
     */
    virtual double getProbability(size_t categoryIndex) const = 0;

    /**
     * @param category The value associated to the class.
     * @return The probability associated to a given class.
     */
    virtual double getProbability(double category) const = 0;

    /**
     * @return A vector with all classes values.
     */
    virtual Vdouble getCategories() const = 0;
    /**
     * @return A vector with all probabilities.
     */
    virtual Vdouble getProbabilities() const = 0;

    /**
     * @brief Set the probability associated to a class.
     *
     * If the category does not exist, a new category is created
     * with the corresponding probability.
     * If the category already exist, its probability is set to 'probability'.
     * The sum of all probabilities is not checked.
     * 
     * @param category The class value.
     * @param probability The class probability.
     */
    virtual void set(double category, double probability) = 0;
		
    /**
     * @brief Modify the probability associated to a class.
     *
     * If the category does not exist, a new category is created
     * with the corresponding probability.
     * if the category exists, add 'probability' to the existing probability.
     * The sum of all probabilities is not checked.
     * 
     * @param category The class value.
     * @param probability The class probability.
     */
    virtual void add(double category, double probability) = 0;

    /**
     * @return \f$Pr(x < \mbox{category})\f$.
     * @param category The class value.
     */
    virtual double getInfCumulativeProbability(double category) const = 0;
    /**
     * @return \f$Pr(x \leq \mbox{category})\f$.
     * @param category The class value.
     */
    virtual double getIInfCumulativeProbability(double category) const = 0;
    /**
     * @return \f$Pr(x > \mbox{category})\f$.
     * @param category The class value.
     */
    virtual double getSupCumulativeProbability(double category) const = 0;
    /**
     * @return \f$Pr(x \geq \mbox{category})\f$.
     * @param category The class value.
     */
    virtual double getSSupCumulativeProbability(double category) const = 0;
	
    /**
     * @brief Draw a random number from this distribution.
     *
     * This number will be one of the class values, drawn according
     * to the class probabilities.
     * 
     * @return A random number according to this distribution.
     */
    virtual double rand() const = 0;

    /**
     * @brief Draw a random number from the continuous version of this distribution, if it exists.
     *
     * Uses the continuous version of this distribution to draw a random number.
     * 
     * @return A random number according to this distribution.
     * @throw Exception If there is no continuous version of this distribution.
     */
    virtual double randC() const throw (Exception) = 0;

    /**
     * @brief Return the quantile of the continuous version of the
     * distribution, ie y such that @f$ Prob(X<y)=x @f$
     * 
     **/
     
    virtual double qProb(double x) const  = 0;

    /**
     * @brief Return the cumulative quantile of the continuous version
     * of the distribution, ie @f$ Prob(X<x) @f$.
     *
     **/
     
    virtual double pProb(double x) const  = 0;

    /**
     * @brief Return a primitive function used for the expectation
     * of the continuous version of the distribution, ie
     * @f$ E(X|...<=X<a) = \int_{...}^a t.dProb(t) dt @f$.
     *
     **/
     
    virtual double Expectation(double a) const  = 0;

    /**
     *@brief Sets the median value to true to say that the value in a
     * class is proportional to the median value of the class, the
     * proportionality factor being such that the sum of the values
     * equals the expectation of the distribution. If it is set to
     * false, the value is the mean value in the class.
     *
     * If the median value is modified, the discretization process is
     * launched.
     *
     * @param median tells how the value associated to each class is
     * computed. 
     */
    
    virtual void setMedian(bool median) = 0;
    
    /**
     * @brief Discretizes the distribution in equiprobable classes.
     */

    virtual void discretize() = 0;

    /**
     *@return A vector of all the bounds
     */
    virtual Vdouble getBounds() const = 0;
    
    /**
     * @return the i th internal bound
     */
    virtual double getBound(size_t) const = 0;

    /**
     *@brief methods about the range of the definition
     *
     **/
     
    /**
     * @return The lowest value.
     */
    virtual double getLowerBound() const {
      return (-NumConstants::VERY_BIG());
    }
  
    /**
     * @return The highest value.
     */
    virtual double getUpperBound() const {
      return (NumConstants::VERY_BIG());
    }
  
    /**
     * @return The lowest value.
     */
    virtual bool strictLowerBound() const {
      return (true);
    }
  
    /**
     * @return The highest value.
     */
    virtual bool strictUpperBound() const {
      return (true);
    }
  
    /**
     * @brief Restricts the distribution to the domain where the
     * constraint is respected, in addition of other predefined
     * constraints.
     *
     * If the domain interval is modified, the discretization process
     * is launched.
     *
     * @param c The Constraint to respect.
     *
     */
    
    virtual void restrictToConstraint(const Constraint& c) = 0;
          
    /**
     * @brief Print the distribution (categories and corresponding probabilities) to a stream.
     *
     * @param out The outstream where to print the distribution.
     */
    virtual void print(OutputStream& out) const = 0;

  };

} //end of namespace bpp.

#endif  //_DISCRETEDISTRIBUTION_H_

