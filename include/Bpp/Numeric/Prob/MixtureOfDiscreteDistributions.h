//
// File: MixtureOfDiscreteDistributions.h
// Created by: Laurent Guéguen
// Created on: mercredi 9 juin 2010, à 22h 44
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

#ifndef _MIXTUREOFDISCRETEDISTRIBUTIONS_H_
#define _MIXTUREOFDISCRETEDISTRIBUTIONS_H_

#include "AbstractDiscreteDistribution.h"

namespace bpp
{
/**
 * @brief A Discrete distribution object defined by a vector of
 * Discrete Distributions and a set of probabilities for these
 * Discrete Distributions.
 *
 * The non-null values of the MixtureOfDiscreteDistributions are all
 * the non-null values of the Discrete Distributions, with
 * probabilities equal to their probabilities in each Discrete
 * Distribution multiplied by the specific probability of this
 * Distribution.
 *
 * Parameters:
 *
 * For the probabilities: they are called \c "theta1",... and defined
 * as @f$ \theta_{i \in 1..\textrm{size-1}} @f$ such that probability
 * of value @f$i @f$ is @f$ (1-\theta_1).(1-\theta_2)...\theta_{i} @f$
 *
 * For the values: they are the parameters of the Discrete
 * Distributions, prefixed by the index in the vector of the Discrete
 * Distributions.
 *
 */
class MixtureOfDiscreteDistributions :
  public AbstractDiscreteDistribution
{
protected:
  std::vector<DiscreteDistribution*> vdd_;

  std::vector<double> probas_;

  std::vector<std::string> vNestedPrefix_;

public:
  /**
   * @brief Builds a new MixtureOfDiscreteDistributions object from a
   * vector of Discrete Distributions and a vector of probabilities.
   * The Discrete Distributions are cloned in the constructor to
   * become attributes.
   *
   * @param distributions The vector of pointers to Discrete
   * Distributions.
   * @param probas The vector of probabilities.
   *
   */

  MixtureOfDiscreteDistributions(const std::vector<DiscreteDistribution*>& distributions, const std::vector<double>& probas);

  ~MixtureOfDiscreteDistributions();

  MixtureOfDiscreteDistributions(const MixtureOfDiscreteDistributions& mdd);

  MixtureOfDiscreteDistributions& operator=(const MixtureOfDiscreteDistributions& mdd);

#if defined(NO_VIRTUAL_COV)
  Clonable* clone() const { return new MixtureOfDiscreteDistributions(*this); }
#else
  MixtureOfDiscreteDistributions* clone() const { return new MixtureOfDiscreteDistributions(*this); }
#endif

public:
  std::string getName() const {return "Mixture"; }

  /**
   * @brief Returns the number of discrete distributions in the
   * mixture.
   *
   */
  size_t getNumberOfDistributions() const {return vdd_.size(); }

  /**
   * @brief Returns a pointer to the n-th discrete distribution in the mixture.
   *
   * @param n tne number of the distribution in the mixture;
   */
  const DiscreteDistribution* getNDistribution(size_t n) const
  {
    return vdd_[n];
  }

  /**
   * @brief Returns the probability of the n-th discrete distribution in the mixture.
   *
   * @param n the number of the distribution in the mixture;
   */
  double getNProbability(size_t n) const
  {
    return probas_[n];
  }

  /**
   * @brief sets the number of categories of EACH submodel to
   * nbClasses, so the number of categories of the mixture is the sum
   * of all these numbers.
   *
   */

  void setNumberOfCategories(size_t nbClasses);

  void fireParameterChanged(const ParameterList& parameters);

  double qProb(double x) const;

  double pProb(double x) const;

  double Expectation(double a) const;

  void setMedian(bool median);

  void restrictToConstraint(const Constraint& c);

  void discretize();

  void setNamespace(const std::string& prefix);

protected:
  void updateDistribution();
};
} // end of namespace bpp.

#endif  // _MIXTUREOFDISCRETEDISTRIBUTIONS_H__

