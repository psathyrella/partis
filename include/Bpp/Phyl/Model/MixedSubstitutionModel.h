//
// File: MixedSubstitutionModel.h
// Created by: Laurent Gueguen
//

/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _MIXEDSUBSTITUTIONMODEL_H_
#define _MIXEDSUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"
// #include <Bpp/Seq/Alphabet.all>

#include <vector>
#include <string>
#include <map>
#include <cstring> // C lib for string copy

namespace bpp
{
  /**
   * @brief Interface for Substitution models, defined as a mixture
   * of "simple" substitution models.
   * @author Laurent Guéguen
   *
   */

  class MixedSubstitutionModel :
    public virtual SubstitutionModel
  {
  public:
    MixedSubstitutionModel() {};

    virtual ~MixedSubstitutionModel(){};

    virtual MixedSubstitutionModel* clone() const = 0;

  public:

    /**
     * @brief Returns a specific model from the mixture
     */
    virtual const SubstitutionModel* getNModel(size_t i) const = 0;

    virtual SubstitutionModel* getNModel(size_t i) = 0;

    /**
     * @brief Returns the  probability of a specific model from the mixture
     */
  
    virtual double getNProbability(size_t i) const = 0;

    virtual const std::vector<double>& getProbabilities() const = 0;

    /**
     * @brief Sets the  probability of a specific model from the mixture
     */
  
    virtual void setNProbability(size_t i, double prob) = 0;

    virtual size_t getNumberOfModels() const = 0;

    /**
     * @brief Returns the rates of the submodels.
     */

    virtual const std::vector<double>& getVRates() const = 0;

    /**
     * @brief Returns the rate of a specific submodel.
     */

    virtual double getNRate(size_t i) const = 0;

    /**
     * @brief Sets the rates of the submodels to be proportional to a
     * given vector, and normalizes them so that the mean rate of the
     * mixture equals rate_.
     
     * @param vd a vector of positive values such that the rates of
     * the respective submodels are in the same proportions (ie this
     * vector does not need to be normalized).
     */

    virtual void setVRates(const Vdouble& vd) = 0;

    /**
     * @brief Normalizes the rates of the submodels so that the mean
     * rate of the mixture equals rate_.
     */

    virtual void normalizeVRates() = 0;

    /*
     *@brief Returns the vector of numbers of the submodels in the
     *mixture that match a description.
     *
     */
  
    virtual Vint getSubmodelNumbers(std::string& desc) const = 0;


  };
} // end of namespace bpp.

#endif  // _MIXEDSUBSTITUTIONMODEL_H_
