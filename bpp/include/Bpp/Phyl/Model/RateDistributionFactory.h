//
// File: RateDistributionFactory.h
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Fri apr 14 11:11 2006
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#ifndef _RATEDISTRIBUTIONFACTORY_H_
#define _RATEDISTRIBUTIONFACTORY_H_

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

// From the STL:
#include <string>

namespace bpp
{

/**
 * @brief Utilitary class for creating rate distributions.
 *
 * @see SubstitutionModelFactory
 */
class RateDistributionFactory
{
public:
  static const std::string CONSTANT_DISTRIBUTION;
  static const std::string GAMMA_DISTRIBUTION;
  static const std::string GAUSSIAN_DISTRIBUTION;
  static const std::string EXPONENTIAL_DISTRIBUTION;

private:
  unsigned int nbClasses_; // For discrete distributions. 
  
public:
  /**
   * @brief Creates a new factory object.
   *
   * Example:
   * @code
   * DiscreteDistribution* dist = RateDistributionFactory().createDiscreteDistribution(RateDistributionFactory::GAMMA_DISTRIBUTION);
   * // or DiscreteDistribution* dist = RateDistributionFactory(10).createDiscreteDistribution(RateDistributionFactory::GAMMA_DISTRIBUTION);
   * // or DiscreteDistribution* dist = RateDistributionFactory().createDiscreteDistribution(RateDistributionFactory::GAMMA_DISTRIBUTION, 10);
   * // dist can be used in any object dealing with rate distributions.
   * @endcode
   */
  RateDistributionFactory(unsigned int nbClasses = 4): nbClasses_(nbClasses) {}
  virtual ~RateDistributionFactory() {}

public:
  /**
   * @brief Get a new dynamically created DiscreteDistribution object.
   *
   * @param distName The name of the distribution to use.
   * @param nbClasses The number of classes to use. This override the value passed to the constructor and is ignored for a constant distribution.
   * @return A pointer toward a new discrete distribution, with default parameter values.
   * @throw Exception If the dist name do not match any available distribution.
   */
  virtual DiscreteDistribution* createDiscreteDistribution(const std::string& distName, unsigned int nbClasses) throw (Exception);
  /**
   * @brief Get a new dynamically created DiscreteDistribution object.
   *
   * @param distName The name of the distribution to use.
   * @return A pointer toward a new discrete distribution, with default parameter values.
   * @throw Exception If the dist name do not match any available distribution.
   */
  virtual DiscreteDistribution* createDiscreteDistribution(const std::string& distName) throw (Exception)
  {
    return createDiscreteDistribution(distName, nbClasses_);
  }

};

} //end of namespace bpp.

#endif //_RATEDISTRIBUTIONFACTORY_H_

