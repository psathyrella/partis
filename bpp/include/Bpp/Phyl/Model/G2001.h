//
// File: G2001.h
// Created by: Julien Dutheil
// Created on: Mon Aug 07 18:31 2006
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _G2001_H_
#define _G2001_H_

#include "MarkovModulatedSubstitutionModel.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Parameter.h>

namespace bpp
{
/**
 * @brief Galtier's 2001 covarion model.
 *
 * This model is a subclass of the so-called Markov-modulated substitution models,
 * with a Jukes-Cantor rate matrix, of parameter @f$\nu@f$.
 * the original version uses a discrete @f$\Gamma@f$ distribution for rates, but you can
 * use it with virtually any rate distribution.
 *
 * @see MarkovModulatedSubstitutionModel
 *
 * Galtier N., Maximum-likelihood phylogenetic analysis under a covarion-like model (2001).
 * _Molecular Biology and Evolution_, 18:866-73.
 */
class G2001 :
  public MarkovModulatedSubstitutionModel
{
private:
  DiscreteDistribution* rDist_;

  std::string nestedRatePrefix_;

public:
  /**
   * @brief Build a new G2001 substitution model.
   *
   * @param model The substitution model to use. May be of any alphabet type.
   * @param rDist The discrete distribution for rates. The class will own the DiscreteDistribution object,
   * which will be deleted together with this instance.
   * @param nu    The rate matrix parameter.
   * @param normalizeRateChanges Tell if the rate transition matrix should be normalized.
   */
  G2001(ReversibleSubstitutionModel* model, DiscreteDistribution* rDist, double nu = 1., bool normalizeRateChanges = false) :
    MarkovModulatedSubstitutionModel(model, normalizeRateChanges, "G01."),
    rDist_(rDist),
    nestedRatePrefix_("rdist_" + rDist->getNamespace())
  {
    nbRates_ = rDist_->getNumberOfCategories();
    ratesExchangeability_.resize(nbRates_, nbRates_);
    rates_.resize(nbRates_, nbRates_);
    ratesFreq_ = std::vector<double>(nbRates_, 1. / static_cast<double>(nbRates_));
    rDist_->setNamespace(getNamespace() + nestedRatePrefix_);
    addParameters_(rDist_->getIndependentParameters());
    addParameter_(new Parameter("G01.nu", nu, &Parameter::R_PLUS));
    updateRatesModel();
    updateMatrices();
  }

  G2001(const G2001& model) :
    MarkovModulatedSubstitutionModel(model),
    rDist_(dynamic_cast<DiscreteDistribution*>(model.rDist_->clone())),
    nestedRatePrefix_(model.nestedRatePrefix_)
  {}

  G2001& operator=(const G2001& model)
  {
    MarkovModulatedSubstitutionModel::operator=(model);
    rDist_ = dynamic_cast<DiscreteDistribution*>(model.rDist_->clone());
    nestedRatePrefix_ = model.nestedRatePrefix_;
    return *this;
  }

  virtual ~G2001() { delete rDist_; }

  G2001* clone() const { return new G2001(*this); }

public:
  std::string getName() const { return "G01"; }

  /**
   * @brief Re-definition of the super-class method to update the rate distribution too.
   *
   * @param parameters The parameters that have been modified.
   */
  void fireParameterChanged(const ParameterList& parameters)
  {
    rDist_->matchParametersValues(parameters);
    MarkovModulatedSubstitutionModel::fireParameterChanged(parameters);
  }

  /**
   * @return The rate distribution associated to this instance.
   */
  const DiscreteDistribution* getRateDistribution() const { return rDist_; }

  void setNamespace(const std::string& prefix)
  {
    MarkovModulatedSubstitutionModel::setNamespace(prefix);
    // We also need to update the namespace of the nested distribution:
    rDist_->setNamespace(prefix + nestedRatePrefix_);
  }


  double getRate() const {  return 1.; }

  void setRate(double rate) {}

  void addRateParameter() {}

protected:
  void updateRatesModel()
  {
    double nu = getParameterValue("nu");
    for (size_t i = 0; i < nbRates_; i++)
    {
      rates_(i, i) = rDist_->getCategory(i);
      for (size_t j = 0; j < nbRates_; j++)
      {
        if (i == j)
        {
          ratesExchangeability_(i, j) = -static_cast<double>(nbRates_) * nu;
        }
        else
        {
          ratesExchangeability_(i, j) = static_cast<double>(nbRates_) * nu / static_cast<double>(nbRates_ - 1);
        }
      }
    }
  }
};
} // end of namespace bpp.

#endif // _G2001_H_

