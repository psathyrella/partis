//
// File: ExponentialDiscreteDistribution.h
// Created by: Julien Dutheil
// Created on: Tue Nov 13 12:37 2007
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

#ifndef _EXPONENTIALDISCRETEDISTRIBUTION_H_
#define _EXPONENTIALDISCRETEDISTRIBUTION_H_

#include "AbstractDiscreteDistribution.h"
#include "../Constraints.h"
#include "../Random/RandomTools.h"

namespace bpp
{
/**
 * @brief Discretized Exponential distribution.
 *
 * @author Julien Dutheil
 */
class ExponentialDiscreteDistribution :
  public AbstractDiscreteDistribution
{
protected:
  double lambda_;

public:
  /**
   * @brief Build a new discretized exponential distribution.
   * @param n the number of categories to use.
   * @param lambda The lambda parameter.
   *
   * The Parameter is: lambda @f$ \in [0;\infty[ @f$.
   *
   */

  ExponentialDiscreteDistribution(size_t n, double lambda = 1.);

  ExponentialDiscreteDistribution(const ExponentialDiscreteDistribution& dist) :
    AbstractParameterAliasable(dist),
    AbstractDiscreteDistribution(dist),
    lambda_(dist.lambda_)
  {}

  ExponentialDiscreteDistribution& operator=(const ExponentialDiscreteDistribution& dist)
  {
    AbstractParameterAliasable::operator=(dist);    
    AbstractDiscreteDistribution::operator=(dist);
    lambda_ = dist.lambda_;
    return *this;
  }

  ~ExponentialDiscreteDistribution(){};

  ExponentialDiscreteDistribution* clone() const { return new ExponentialDiscreteDistribution(*this); }

public:
  std::string getName() const {return("Exponential");}
  
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
    return 1. - exp(-lambda_ * x);
  }

  double qProb(double x) const
  {
    return -log(1. - x) / lambda_;
  }

  double Expectation(double a) const
  {
    return 1. / lambda_ - exp(-a * lambda_) * (a + 1. / lambda_);
  }
};
} //end of namespace bpp.

#endif  //_EXPONENTIALDISCRETEDISTRIBUTION_H_

