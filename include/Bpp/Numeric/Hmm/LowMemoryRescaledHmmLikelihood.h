//
// File: LowMemoryRescaledHmmLikelihood.h
// Created by: Julien Dutheil
// Created on: Wed Dec 16 10:47 2009
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

#ifndef _LOWMEMORYRESCALEDHMMLIKELIHOOD_H_
#define _LOWMEMORYRESCALEDHMMLIKELIHOOD_H_

#include "HmmLikelihood.h"
#include "../AbstractParametrizable.h"
#include "../Matrix/Matrix.h"

// From the STL:
#include <vector>
#include <memory>

namespace bpp
{
/**
 * @brief A modified implementation of the RescaledHmmLikelihood implementation, with lower memory usage.
 *
 * This implementation is similar to the one used in the RescaledHmmLikelihood class,
 * but does not store the full likelihood array. The benefit of it is a significantly reduced
 * memory usage, allowing to compute likelihood for very large data sets.
 *
 * The drawback is that this class can't compute posterior
 * probabilities, neither derivatives of the likelihoods, and can
 * hence only be used to compute likelihoods.
 *
 */
  
class LowMemoryRescaledHmmLikelihood :
  public AbstractHmmLikelihood,
  public AbstractParametrizable
{
private:
  /**
   * @brief The alphabet describing the hidden states.
   */
  std::auto_ptr<HmmStateAlphabet> hiddenAlphabet_;
  std::auto_ptr<HmmTransitionMatrix> transitionMatrix_;
  std::auto_ptr<HmmEmissionProbabilities> emissionProbabilities_;

  /**
   * @brief The likelihood array.
   *
   * Here we use two arrays for the i and i-1 positions
   */
  std::vector<double> likelihood1_;
  std::vector<double> likelihood2_;
  double logLik_;
  size_t maxSize_;

  std::vector<size_t> breakPoints_;

  size_t nbStates_, nbSites_;

public:
  /**
   * @brief Build a new LowMemoryRescaledHmmLikelihood object.
   *
   * @warning the HmmTransitionMatrix and HmmEmissionProbabilities object passed as argument must be non-null
   * and point toward the same HmmStateAlphabet instance. The three object will be copied if needed, and
   * deleted when the hmm likelihood objet is deleted. You should secure a copy before if you don't want them to
   * be destroyed with this object.
   *
   * @param hiddenAlphabet The hidden states alphabet to use.
   * @param transitionMatrix The transition matrix to use.
   * @param emissionProbabilities The emission probabilities to use.
   * @param prefix A namespace for parameter names.
   * @param maxSize the maximum size of the vector of scales. If this size is exceeded, then a temporary likelihood computation is made and stored, and the vector is reset.
   * the size of the vector specify the memory usage of the class. A two low value can lead to numerical precision errors.
   */
  LowMemoryRescaledHmmLikelihood(
    HmmStateAlphabet* hiddenAlphabet,
    HmmTransitionMatrix* transitionMatrix,
    HmmEmissionProbabilities* emissionProbabilities,
    const std::string& prefix,
    size_t maxSize = 1000000) throw (Exception);

  LowMemoryRescaledHmmLikelihood(const LowMemoryRescaledHmmLikelihood& lik) :
    AbstractHmmLikelihood(lik),
    AbstractParametrizable(lik),
    hiddenAlphabet_(dynamic_cast<HmmStateAlphabet*>(lik.hiddenAlphabet_->clone())),
    transitionMatrix_(dynamic_cast<HmmTransitionMatrix*>(lik.transitionMatrix_->clone())),
    emissionProbabilities_(dynamic_cast<HmmEmissionProbabilities*>(lik.emissionProbabilities_->clone())),
    likelihood1_(lik.likelihood1_),
    likelihood2_(lik.likelihood2_),
    logLik_(lik.logLik_),
    maxSize_(lik.maxSize_),
    breakPoints_(lik.breakPoints_),
    nbStates_(lik.nbStates_),
    nbSites_(lik.nbSites_)
  {
    // Now adjust pointers:
    transitionMatrix_->setHmmStateAlphabet(hiddenAlphabet_.get());
    emissionProbabilities_->setHmmStateAlphabet(hiddenAlphabet_.get());
  }

  LowMemoryRescaledHmmLikelihood& operator=(const LowMemoryRescaledHmmLikelihood& lik)
  {
    AbstractHmmLikelihood::operator=(lik);
    AbstractParametrizable::operator=(lik);
    hiddenAlphabet_        = std::auto_ptr<HmmStateAlphabet>(dynamic_cast<HmmStateAlphabet*>(lik.hiddenAlphabet_->clone()));
    transitionMatrix_      = std::auto_ptr<HmmTransitionMatrix>(dynamic_cast<HmmTransitionMatrix*>(lik.transitionMatrix_->clone()));
    emissionProbabilities_ = std::auto_ptr<HmmEmissionProbabilities>(dynamic_cast<HmmEmissionProbabilities*>(lik.emissionProbabilities_->clone()));
    likelihood1_           = lik.likelihood1_;
    likelihood2_           = lik.likelihood2_;
    logLik_                = lik.logLik_;
    maxSize_               = lik.maxSize_;
    breakPoints_           = lik.breakPoints_;
    nbStates_              = lik.nbStates_;
    nbSites_               = lik.nbSites_;

    // Now adjust pointers:
    transitionMatrix_->setHmmStateAlphabet(hiddenAlphabet_.get());
    emissionProbabilities_->setHmmStateAlphabet(hiddenAlphabet_.get());
    return *this;
  }

  virtual ~LowMemoryRescaledHmmLikelihood() {}

#ifndef NO_VIRTUAL_COV
  LowMemoryRescaledHmmLikelihood*
#else
  Clonable*
#endif
  clone() const { return new LowMemoryRescaledHmmLikelihood(*this); }

public:
  const HmmStateAlphabet& getHmmStateAlphabet() const { return *hiddenAlphabet_; }
  HmmStateAlphabet& getHmmStateAlphabet() { return *hiddenAlphabet_; }

  const HmmTransitionMatrix& getHmmTransitionMatrix() const { return *transitionMatrix_; }
  HmmTransitionMatrix& getHmmTransitionMatrix() { return *transitionMatrix_; }

  const HmmEmissionProbabilities& getHmmEmissionProbabilities() const { return *emissionProbabilities_; }
  HmmEmissionProbabilities& getHmmEmissionProbabilities() { return *emissionProbabilities_; }

  void setBreakPoints(const std::vector<size_t>& breakPoints) {
    breakPoints_ = breakPoints;
    computeForward_();
  }

  const std::vector<size_t>& getBreakPoints() const { return breakPoints_; }

  void setParameters(const ParameterList& pl) throw (Exception)
  {
    setParametersValues(pl);
  }

  double getValue() const throw (Exception) { return -logLik_; }

  double getLogLikelihood() const { return logLik_; }

  void fireParameterChanged(const ParameterList& pl);


  double getLikelihoodForASite(size_t site) const
  {
    throw (NotImplementedException("LowMemoryRescaledHmmLikelihood::getLikelihoodForASite. This class can't compute posterior probabilities, use RescaledHmmLikelihood instead."));    
  }


  Vdouble getLikelihoodForEachSite() const
  {
    throw (NotImplementedException("LowMemoryRescaledHmmLikelihood::getLikelihoodForEachSite. This class can't compute posterior probabilities, use RescaledHmmLikelihood instead."));    
  }

  Vdouble getHiddenStatesPosteriorProbabilitiesForASite(size_t site) const
  {
    throw (NotImplementedException("LowMemoryRescaledHmmLikelihood::getHiddenStatesPosteriorProbabilitiesForASite. This class can't compute posterior probabilities, use RescaledHmmLikelihood instead."));    
  }

  void getHiddenStatesPosteriorProbabilities(std::vector< std::vector<double> >& probs, bool append = false) const throw (NotImplementedException)
  {
    throw (NotImplementedException("LowMemoryRescaledHmmLikelihood::getHiddenStatesPosteriorProbabilities. This class can't compute posterior probabilities, use RescaledHmmLikelihood instead."));    
  }

protected:
  void computeForward_();

  void computeDLikelihood_() const
  {
    //    computeDForward_();
  }

  void computeD2Likelihood_() const
  {
    //    computeD2Forward_();
  }
};

}

#endif // _LOWMEMORYRESCALEDHMMLIKELIHOOD_H_

