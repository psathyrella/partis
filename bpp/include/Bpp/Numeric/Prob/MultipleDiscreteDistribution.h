//
// File: MultipleDiscreteDistribution.h
// Created by: Laurent Guéguen
// Created on: mardi 20 juillet 2010, à 14h 52
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

#ifndef _MULTIPLEDISCRETEDISTRIBUTION_H_
#define _MULTIPLEDISCRETEDISTRIBUTION_H_

#include "../VectorTools.h"
#include "../ParameterAliasable.h"
#include "../../Exceptions.h"
#include "../../Io/OutputStream.h"

namespace bpp
{

  /**
   * @brief Interface for multiple discrete distribution objects.
   *
   * A multiple discrete distribution usually contains a vector of
   * finite set of categories and a probability associated to each.
   * The size of the vector is the dimension of the distribution.
   *
   * Each category (or class) is defined by two bounds, and sometimes
   * by a mean or a median value.
   *
   * A multiple discrete distribution may contain one or several
   * parameters. The probabilities associated to each class usually
   * depends on the parameter values. In some cases, the number and/or
   * bounds of the classes may also depend on the parameters values,
   * depending on the kind of discretization used.
   */

  class MultipleDiscreteDistribution:
    public virtual ParameterAliasable
  {
  public:
    MultipleDiscreteDistribution() {}
    
    virtual ~MultipleDiscreteDistribution() {}

#ifndef NO_VIRTUAL_COV
    MultipleDiscreteDistribution * clone() const = 0;
#endif

  public:

    /**
     * @return The number of categories 
     */
    virtual size_t getNumberOfCategories() const = 0;

    /**
     * @param Vvalue The vector of values to check.
     * @return The vector of categories of the classes the value is
     * in. Throws a ConstraintException if the value is off the domain
     * of the MultipleDiscreteDistribution.
     */
    
    virtual Vdouble getValueCategory(Vdouble& Vvalue) const = 0;
    
    /**
     * @param category The vector of values associated to the class.
     * @return The probability associated to a given class.
     */
    virtual double getProbability(Vdouble& category) const = 0;


  public:

    /**
     * @brief Draw a random vector from this distribution.
     *
     * This vector will be one of the class values, drawn according
     * to the class probabilities.
     * 
     * @return A random number according to this distribution.
     */
    virtual Vdouble rand() const = 0;

    /**
     * @brief Draw a random vector from the continuous version of this distribution, if it exists.
     *
     * Uses the continuous version of this distribution to draw a random vector.
     * 
     * @return A random vector according to this distribution.
     * @throw Exception If there is no continuous version of this distribution.
     */
    virtual Vdouble randC() const = 0;

    /**
     * @brief Checks if the Parameters can respect the given 
     * Constraints (one per dimension) and optionnaly tries to modify
     * their own Constraints.
     *
     * @param vc The vector of Constraint to respect.
     * @param f boolean flag to say if the Constraints must be changed
     * (if possible) (default: true)
     *
     * @return true if the Parameters Constraints are adapted to the
     * given Constraints, false otherwise.
     */
    //virtual bool adaptToConstraint(const std::vector<Constraint&>& vc, bool f=true) = 0;
    
    /**
     * @brief Print the distribution (categories and corresponding probabilities) to a stream.
     *
     * @param out The outstream where to print the distribution.
     */
    //    virtual void print(OutputStream& out) const = 0;

  };

} //end of namespace bpp.

#endif  //_MULTIPLEDISCRETEDISTRIBUTION_H_

