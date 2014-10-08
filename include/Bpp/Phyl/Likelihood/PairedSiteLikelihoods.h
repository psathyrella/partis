//
// File: PairedSiteLikelihoods.h
// Created by: Nicolas Rochette
// Created on: January 6, 2011
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

#ifndef _PAIREDSITELLIKELIHOODS_H_
#define _PAIREDSITELLIKELIHOODS_H_

// From the STL:
#include <vector>
#include <string>

// From Bio++
#include "TreeLikelihood.h"
#include <Bpp/Exceptions.h>

namespace bpp
{
/**
 * @brief A container for paired-site likelihoods (likelihoods over
 * the same sites for different models, especially topologies).
 * An instance of this class is, roughly, a list of models, each of
 * them having a name (stored in the <i>modelNames</i> attribute) and
 * a set of site likelihoods (stored in the <i>logLikelihoods</i> attribute).
 */
class PairedSiteLikelihoods
{
private:
  std::vector<std::vector<double> > logLikelihoods_;
  std::vector<std::string> modelNames_;

public:
  PairedSiteLikelihoods();

  /**
   * @brief Build a new object from a site likelihoods array.
   *
   * @param siteLogLikelihoods An nmodels*nsites array of loglikelihoods.
   * @param modelNames <i>(Optional)</i> The names of the models.
   *
   * @throw Exception If the number of sites differ between the models,
   * or if the number of names and site loglikelihood records differ.
   */
  PairedSiteLikelihoods(
    const std::vector<std::vector<double> >& siteLogLikelihoods,
    const std::vector<std::string>& modelNames = std::vector<std::string>()
    ) throw (Exception);

  ~PairedSiteLikelihoods() {}

  /**
   * @brief Append a model.
   *
   * @param siteLogLikelihoods The loglikelihoods of the sites under this model.
   * @param modelName The name of the model.
   *
   * @throw Exception If the number of sites is not the same as in the container.
   */
  void appendModel(
    const std::vector<double>& siteLogLikelihoods,
    const std::string& modelName = "") throw (Exception);

  /**
   * @brief Append a model.
   *
   * @param treelikelihood A TreeLikelihood record.
   *
   * @throw Exception If the number of sites is not the same as in the container.
   */
  void appendModel(const bpp::TreeLikelihood& treelikelihood) throw (Exception);

  /**
   * @brief Append models by concatenation.
   *
   * @param psl the PairedSiteLikelihoods object to append to the caller.
   *
   * @throw Exception If the number of sites in the two object is not equal.
   */
  void appendModels(const PairedSiteLikelihoods& psl) throw (Exception);

  /**
   * @return The site-likelihoods of all models.
   */
  const std::vector<std::vector<double> >& getLikelihoods() const
  {
    return logLikelihoods_;
  }

  /**
   * @return The model names.
   */
  const std::vector<std::string>& getModelNames() const
  {
    return modelNames_;
  }

  /** @brief Get the number of models in the container. */
  size_t getNumberOfModels() const
  {
    return logLikelihoods_.size();
  }

  /**
   * @return The number of sites for each model.
   * @throw Exception If the container is empty.
   */
  std::size_t getNumberOfSites() const throw (Exception)
  {
    try
    {
      return logLikelihoods_.at(0).size();
    }
    catch (std::out_of_range&)
    {
      throw Exception("PairedSiteLikelihoods::nsites: The container is empty, there isn't a number of sites.");
    }
  }

  /**
   * @brief Set the name of a model.
   *
   * @param pos The position of the target model.
   * @param name The new name.
   */
  void setName(std::size_t pos, std::string& name)
  {
    modelNames_.at(pos) = name;
  }

  /**
   * @brief Compute the Expected Likelihood Weights of the models.
   *
   * The weight \f$W_m\f$ of a model is :
   * \f[
   * W_m
   * = \frac{1}{B} \sum_{b \in B} \frac{L_m^{(b)}}{\sum L_k^{(b)}}
   * = \frac{1}{B} \sum_{b \in B} \frac{exp(Y^{(b)}_m - Ymax^{(b)})}{\sum exp(Y_k^{(b)} - Y_{max}^{(b)})}
   * \f]
   * where \f$Y_k^{(b)}\f$ is the loglikelihood of model k for replicate b.
   * @return A pair of vectors containing the names and weights of the models.
   *
   * @param replicates The number of pseudoreplicates over which the weights are to be averaged.
   */
  std::pair< std::vector<std::string>, std::vector<double> > computeExpectedLikelihoodWeights(int replicates = 10000) const;

  /**
   * @brief Draw a nonparametric pseudoreplicate
   *
   * @param length The length of the data.
   * @param scaling The length of the pseudoreplicate, in fraction of
   *  the length of the data.
   *
   * @return A vector of the same length as the data, containing the count
   * of each element in the pseudoreplicate.
   */
  static std::vector<int> bootstrap(std::size_t length, double scaling = 1);
  
};
} // namespace bpp.

#endif  //_PAIREDSITELLIKELIHOODS_H_
