//
// File: TruncatedExponentialDiscreteDistribution.h
// Created by: Julien Dutheil
// Created on: Fri Jan 25 15:24 2008
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#ifndef _TRUNCATEDEXPONENTIALDISCRETEDISTRIBUTION_H_
#define _TRUNCATEDEXPONENTIALDISCRETEDISTRIBUTION_H_

#include "ExponentialDiscreteDistribution.h"
#include "../Constraints.h"
#include "../Random/RandomTools.h"

namespace bpp
{
/**
 * @brief Discretized Truncated (on the right) Exponential
 * distribution, where the probabilities are given the exponential,
 * conditioned by the upper limit.
 *
 * This distribution has two parameters: the traditional exponential law parameter,
 * and the abscissa of the truncation. The distribution will be truncated on the right
 * of this point.
 */
  
class TruncatedExponentialDiscreteDistribution :
  public AbstractDiscreteDistribution
{
protected:
  double lambda_;

  double tp_;

  /*
   * Probability of the condition given x < tp_.
   *
   */
   
  double cond_;
  
public:
  /**
   * @brief Build a new truncated exponential discrete distribution.
   * @param n the number of categories to use.
   * @param lambda The lambda parameter
   * @param truncationPoint The truncation point
   *
   * The Parameters are: lambda @f$ \in [0.000001;\infty[ @f$ and tp .@f$ \in [0;\infty[ @f$
   *
   */

  TruncatedExponentialDiscreteDistribution(size_t n, double lambda = 1., double truncationPoint = 10);

  TruncatedExponentialDiscreteDistribution(const TruncatedExponentialDiscreteDistribution& dist) :
    AbstractParameterAliasable(dist),
    AbstractDiscreteDistribution(dist),
    lambda_(dist.lambda_),
    tp_(dist.tp_),
    cond_(dist.cond_)
  {}

  TruncatedExponentialDiscreteDistribution& operator=(const TruncatedExponentialDiscreteDistribution& dist)
  {
    AbstractParameterAliasable::operator=(dist);
    AbstractDiscreteDistribution::operator=(dist);
    lambda_= dist.lambda_;
    tp_ = dist.tp_;
    cond_ = dist.cond_;
    return *this;
  }

  ~TruncatedExponentialDiscreteDistribution(){};

  TruncatedExponentialDiscreteDistribution* clone() const { return new TruncatedExponentialDiscreteDistribution(*this); }

public:

  std::string getName() const {return("TruncExponential");}

  void fireParameterChanged(const ParameterList& parameters);

  double randC() const throw (Exception)
  {
    double x = RandomTools::randExponential(1. / getParameterValue("lambda"));
    while (!intMinMax_.isCorrect(x))
      x = RandomTools::randExponential(1. / getParameterValue("lambda"));

    return x;
  }

  double pProb(double x) const
  {
    if (x>=tp_)
      return 1.;
    else
      return (1. - exp(-lambda_ * x))/cond_;
  }

  double qProb(double x) const
  {
    if (x==1)
      return tp_;
    else
      return -log(1. - cond_*x) / lambda_;
  }

  double Expectation(double a) const
  {
    if (a<tp_)
      return (1. / lambda_ - exp(-a * lambda_) * (a + 1. / lambda_))/cond_;
    else
      return (1. / lambda_ - exp(-tp_ * lambda_) * (tp_ + 1. / lambda_))/cond_;
      
  }

  void restrictToConstraint(const Constraint& c);
};
} //end of namespace bpp.

#endif  //_TRUNCATEDEXPONENTIALDISCRETEDISTRIBUTION_H_

