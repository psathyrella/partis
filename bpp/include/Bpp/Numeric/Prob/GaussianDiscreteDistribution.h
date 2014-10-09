//
// File: GaussianDiscreteDistribution.h
// Created by: Laurent Guéguen
// Created on: April 2010
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

#ifndef _GAUSSIANDISCRETEDISTRIBUTION_H_
#define _GAUSSIANDISCRETEDISTRIBUTION_H_

#include "AbstractDiscreteDistribution.h"
#include "../Constraints.h"
#include "../Random/RandomTools.h"

namespace bpp
{

/**
 * @brief Discretized Gaussian distribution.
 *
 * @author Laurent Guéguen
 */
class GaussianDiscreteDistribution:
  public AbstractDiscreteDistribution
{
  private:
  double mu_, sigma_;
  
  public:
    /**
     * @brief Build a new discretized normal distribution.
     * @param n the number of categories to use.
     * @param mu The mean parameter.
     * @param sigma The standard deviation parameter.
     *
     * The Parameters are: mu @f$ \in ]-\infty;\infty[ @f$ and sigma
     * @f$ \in ]0;\infty[ @f$.
     *
     */
  GaussianDiscreteDistribution(size_t n, double mu=0., double sigma = 1.);

  GaussianDiscreteDistribution(const GaussianDiscreteDistribution&);

  GaussianDiscreteDistribution& operator=(const GaussianDiscreteDistribution&);
  
  virtual ~GaussianDiscreteDistribution();

  GaussianDiscreteDistribution* clone() const { return new GaussianDiscreteDistribution(*this); }
  
public:

  std::string getName() const {return("Gaussian");}

  void fireParameterChanged(const ParameterList & parameters);
  
  double randC() const throw (Exception)
  {
    return RandomTools::randGaussian(mu_,sigma_);
  }

  double qProb(double x) const;
     
  double pProb(double x) const;

  double Expectation(double a) const;

};
  
} //end of namespace bpp.

#endif  //_GAUSSIANDISCRETEDISTRIBUTION_H_

