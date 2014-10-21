//
// File: GammaDiscreteDistribution.h
// Created by: Julien Dutheil
// Created on: Sun Oct 26 20:36:12 2003
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

#ifndef _GAMMADISCRETEDISTRIBUTION_H_
#define _GAMMADISCRETEDISTRIBUTION_H_

#include "AbstractDiscreteDistribution.h"
#include "../Constraints.h"
#include "../Random/RandomTools.h"

namespace bpp
{

  /**
   * @brief Discretized Gamma distribution.
   *
   * @author Julien Dutheil, Laurent Guéguen, with original code from Tal Pupko and Ziheng Yang.
   */
  
  class GammaDiscreteDistribution:
    public AbstractDiscreteDistribution
  {
  private:

    double alpha_, beta_;

    // To prevent useless computations
    double ga1_;
    
  public:
    std::string getName() const {return("Gamma");}

    /**
     * @brief Build a new discretized gamma distribution.
     * @param n the number of categories to use.
     * @param alpha The alpha parameter (shape)
     * @param beta The beta parameter (rate)
     * @param minimumAlpha The minimum allowed value for parameter alpha.
     * @param minimumBeta The minimum allowed value for parameter beta.
     *
     * The Parameters are: alpha and beta @f$ \in [minimumBound;\infty[ @f$.
     * Small values of alpha and/or beta can lead to discretization issues.
     *
     * If @f$ alpha > 1 @f$, the minimum value of the distribution is
     * set to 1e-12, otherwise it is 0.
     */
    GammaDiscreteDistribution(size_t n, double alpha = 1., double beta = 1., double minimumAlpha = 0.05, double minimumBeta = 0.05);

    GammaDiscreteDistribution(const GammaDiscreteDistribution&);

    GammaDiscreteDistribution& operator=(const GammaDiscreteDistribution&);

    virtual ~GammaDiscreteDistribution();

    GammaDiscreteDistribution* clone() const { return new GammaDiscreteDistribution(*this); }

    void fireParameterChanged(const ParameterList & parameters);
  
    double randC() const throw (Exception)
    {
      double x= RandomTools::randGamma(getParameterValue("alpha"),
                                       getParameterValue("beta"));
      while (!intMinMax_.isCorrect(x))
        x= RandomTools::randGamma(getParameterValue("alpha"),
                                  getParameterValue("beta"));
      return x;
    }

    double qProb(double x) const;
     
    double pProb(double x) const;

    double Expectation(double a) const;

};

} //end of namespace bpp.

#endif  //_GAMMADISCRETEDISTRIBUTION_H_

