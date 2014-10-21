//
// File: RHomogeneousMixedTreeLikelihood.h
// Created by: David Fournier, Laurent Gueguen
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

#ifndef _RHOMOGENEOUSMIXEDTREELIKELIHOOD_H_
#define _RHOMOGENEOUSMIXEDTREELIKELIHOOD_H_

#include "RHomogeneousTreeLikelihood.h"
#include "../Model/SubstitutionModel.h"
#include "../Model/MixedSubstitutionModel.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{
/**
 *@ brief A class to compute the average of several
 *RHomogeneousTreeLikelihood defined from a Mixed Substitution
 *Model.
 *
 * In all the calculs, the average of the likelihoods, probabilities
 * are computed.
 **/

class RHomogeneousMixedTreeLikelihood :
  public RHomogeneousTreeLikelihood
{
private:
  std::vector<RHomogeneousTreeLikelihood*> treeLikelihoodsContainer_;
  std::vector<double> probas_;
  
public:
  /**
   * @brief Build a new RHomogeneousMixedTreeLikelihood object without
   * data.
   *
   * This constructor only initialize the parameters. To compute a
   * likelihood, you will need to call the setData() and the
   * computeTreeLikelihood() methods.
   *
   * @param tree The tree to use.
   * @param model The mixed substitution model to use.
   * @param rDist The rate across sites distribution to use.
   * @param checkRooted Tell if we have to check for the tree to be unrooted.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occured.
   */
  RHomogeneousMixedTreeLikelihood(
    const Tree& tree,
    SubstitutionModel* model,
    DiscreteDistribution* rDist,
    bool checkRooted = true,
    bool verbose = true,
    bool usePatterns = true)
  throw (Exception);

  /**
   * @brief Build a new RHomogeneousMixedTreeLikelihood object with data.
   *
   * This constructor initializes all parameters, data, and likelihood arrays.
   *
   * @param tree The tree to use.
   * @param data Sequences to use.
   * @param model The mixed substitution model to use.
   * @param rDist The rate across sites distribution to use.
   * @param checkRooted Tell if we have to check for the tree to be unrooted.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occured.
   */
  RHomogeneousMixedTreeLikelihood(
    const Tree& tree,
    const SiteContainer& data,
    SubstitutionModel* model,
    DiscreteDistribution* rDist,
    bool checkRooted = true,
    bool verbose = true,
    bool usePatterns = true)
  throw (Exception);

  RHomogeneousMixedTreeLikelihood(const RHomogeneousMixedTreeLikelihood& lik);

  RHomogeneousMixedTreeLikelihood& operator=(const RHomogeneousMixedTreeLikelihood& lik);

  virtual ~RHomogeneousMixedTreeLikelihood();

  RHomogeneousMixedTreeLikelihood* clone() const { return new RHomogeneousMixedTreeLikelihood(*this); }

public:
  /**
   * @name The TreeLikelihood interface.
   *
   * Other methods are implemented in the RHomogeneousTreeLikelihood class.
   *
   * @{
   */
  void setData(const SiteContainer& sites) throw (Exception);

  /** @} */


  /**
   * @name The DiscreteRatesAcrossSites interface implementation:
   *
   * @{
   */
  double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
  double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
  double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
  double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
  /** @} */

public:
  // Specific methods:
  void initialize() throw (Exception);

  void fireParameterChanged(const ParameterList& params);

  void computeTreeLikelihood();

  virtual double getDLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

  virtual void computeTreeDLikelihood(const std::string& variable);

  virtual double getD2LikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

  virtual void computeTreeD2Likelihood(const std::string& variable);

protected:
  /**
   * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
   *
   * @param node The root of the subtree.
   */
  virtual void computeSubtreeLikelihood(const Node* node); // Recursive method.

  virtual void computeDownSubtreeDLikelihood(const Node*);

  virtual void computeDownSubtreeD2Likelihood(const Node*);

  /**
   * @brief This method is used by fireParameterChanged method.
   *
   */
  void computeAllTransitionProbabilities();
  /**
   * @brief This method is used by fireParameterChanged method.
   *
   */
  void computeTransitionProbabilitiesForNode(const Node* node);

  /**
   * @brief This method is mainly for debugging purpose.
   *
   * @param node The node at which likelihood values must be displayed.
   */
  virtual void displayLikelihood(const Node* node);

  virtual void setMinimumBranchLength(double brlen) throw (Exception) {
    RHomogeneousMixedTreeLikelihood::setMinimumBranchLength(brlen);
    for (size_t i = 0; i < treeLikelihoodsContainer_.size(); ++i)
      treeLikelihoodsContainer_[i]->setMinimumBranchLength(brlen);
  }
  virtual void setMaximumBranchLength(double brlen) throw (Exception) {
    RHomogeneousMixedTreeLikelihood::setMaximumBranchLength(brlen);
    for (size_t i = 0; i < treeLikelihoodsContainer_.size(); ++i)
      treeLikelihoodsContainer_[i]->setMaximumBranchLength(brlen);
  }
};
} // end of namespace bpp.

#endif  // _RHOMOGENEOUSMIXEDTREELIKELIHOOD_H_

