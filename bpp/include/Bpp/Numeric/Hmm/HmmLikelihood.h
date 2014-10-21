//
// File: HmmLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 26 11:20 2007
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

#ifndef _HMMLIKELIHOOD_H_
#define _HMMLIKELIHOOD_H_


// From NumCalc:
#include "../Function/Functions.h"
#include "../VectorTools.h"

#include "HmmStateAlphabet.h"
#include "HmmTransitionMatrix.h"
#include "HmmEmissionProbabilities.h"

namespace bpp
{

  /**
   * @brief Basal interface for Hidden Markov Models likelihood computation.
   *
   * HmmLikelihood classes compute the probability of data according to parameters (likelihood),
   * using the so-called forward recursion:
   *
   * - Initialisation: @f[ f_0({\cal A}_0), \ldots, f_0({\cal A}_n) @f] The initial frequencies, set to the equilibrium frequencies of the chain.
   * - Recursion (for i=1 to l, the length of the sequence data): @f[ f_i({\cal A}_y) = e_y(D_i) \sum_{x=1}^n f_{i-1}({\cal A}_x) \cdot p_{x,y} @f]
   * - Termination: @f[ \Pr(D) = \sum_{x=1}^n f_l({\cal A}_x) @f]
   * where @f$ {\cal A}_{1..n} @f$ denotes the hidden states of the alphabet, @f$ e_y(D_i) @f$ the
   * probability of the data at position i conditioned on hidden state y (emission probabilities)
   * and @f$ p_{x,y} @f$ is the probability of havving hidden state y at state i+1 knowing there
   * is hidden state x at position i (transition probabilities). These essential elements are given
   * respectively by the HmmEmissionProbabilities and HmmTransitionMatrix object associated to this
   * class. Both objects have to share the same HmmStateAlphabet instance, which describes all
   * allowed hidden states.
   *
   * The HmmLikelihood interface provides essentially two major methods:
   * - A method to retrieve the likelihood value (parameter estimation)
   * - Two methods to retrieve the posterio probabilities of each state using the forward and backward conditionnal likelihoods (posterior decoding).
   */
  
  class HmmLikelihood:
    public virtual DerivableSecondOrder
  {
  public:

#ifndef NO_VIRTUAL_COV
    virtual HmmLikelihood* clone() const = 0;
#endif

    virtual const HmmStateAlphabet& getHmmStateAlphabet() const = 0;
    virtual HmmStateAlphabet& getHmmStateAlphabet() = 0;
    
    virtual const HmmTransitionMatrix& getHmmTransitionMatrix() const = 0;
    virtual HmmTransitionMatrix& getHmmTransitionMatrix() = 0;
    
    virtual const HmmEmissionProbabilities& getHmmEmissionProbabilities() const = 0;
    virtual HmmEmissionProbabilities& getHmmEmissionProbabilities() = 0;

    virtual void getHiddenStatesPosteriorProbabilities(std::vector< std::vector<double> >& probs, bool append) const throw (Exception) = 0;

    virtual Vdouble getHiddenStatesPosteriorProbabilitiesForASite(size_t site) const = 0;

    virtual double getLogLikelihood() const = 0;

    virtual double getDLogLikelihood() const = 0;

    virtual double getD2LogLikelihood() const = 0;

  /**
     * @brief Get the likelihood for a site.
     *
     * @param site The site index to analyse.
     * @return The likelihood for site <i>site</i>.
     */
  
    virtual double getLikelihoodForASite(size_t site) const = 0;

    /**
     * @brief Get the likelihood for each site.
     *
     * @return A vector with all likelihoods for each site.
     */
    virtual Vdouble getLikelihoodForEachSite() const = 0;

    virtual const std::vector<size_t>& getBreakPoints() const = 0;

    virtual void setBreakPoints(const std::vector<size_t>& breakPoints) = 0;

  protected:

    virtual void computeDLikelihood_() const = 0;

    virtual void computeD2Likelihood_() const = 0;

  };


  /*
   * @brief partial impmementation of Hmm Likelihoods.
   *
   */
  
  class AbstractHmmLikelihood:
    public virtual HmmLikelihood
  {
  protected:
    
    mutable double dLogLik_;
    mutable std::string dVariable_;

    mutable double d2LogLik_;
    mutable std::string d2Variable_;

  public:
    AbstractHmmLikelihood();
    
    AbstractHmmLikelihood(const AbstractHmmLikelihood& adhlik);

    AbstractHmmLikelihood& operator=(const AbstractHmmLikelihood& adhlik);

    /* @{
     *
     * @brief From FirstOrder:
     *
     */

    void enableFirstOrderDerivatives(bool yn) {};
    
    bool enableFirstOrderDerivatives() const { return true;}

    double getFirstOrderDerivative(const std::string& variable) const throw (Exception);

    double getDLogLikelihood() const
    {
      return dLogLik_;
    }

    /*
     * @}
     *
     */

    /* @{
     *
     * @brief From SecondOrder:
     *
     */

    void enableSecondOrderDerivatives(bool yn) {};
    
    bool enableSecondOrderDerivatives() const {return true;}

    double getSecondOrderDerivative(const std::string& variable) const throw (Exception);
  
    double getD2LogLikelihood() const
    {
      return d2LogLik_;
    }

    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const throw (Exception) {
      throw (NotImplementedException("AbstractHmmLikelihood::getSecondOrderDerivative is not defined for 2 variables."));
    }
    
    /*
     * @}
     *
     */

  
  };


} //end of namespace bpp.

#endif //_HMMLIKELIHOOD_H_

