//
// File: MixtureOfSubstitutionModels.h
// Created by: Laurent Gueguen
// Date: lundi 13 septembre 2010, à 21h 31
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

#ifndef _MIXTUREOFSUBSTITUTIONMODELS_H_
#define _MIXTUREOFSUBSTITUTIONMODELS_H_

// #include <Bpp/Numeric/Prob.all>
#include <Bpp/Numeric/VectorTools.h>
#include "AbstractMixedSubstitutionModel.h"

#include <vector>
#include <string>
#include <map>
#include <cstring> // C lib for string copy

namespace bpp
{
  /**
   * @brief Substitution models defined as a mixture of several
   * substitution models.
   * @author Laurent Guéguen
   *
   * All the models can be of different types (for example T92 or
   * GY94), and each model has a specific probability and rate. 
   *
   *
   * The probabilities and rates of the models are independent
   * parameters, handled directly, under the constraint that the
   * expectation of the rates on the distribution of the models must
   * equal one. 
   *
   * If there are @f$n@f$ models, @f$p_i@f$ is the probability of
   * model i (@f$\sum_{i=1}^{n} p_i = 1@f$) and the probabilities
   * are defined by relative probabilities parameters @f$rp_i@f$
   * (called "relprobai") with:
   * @f[
   * 1 <= i < n, p_i = (1-rp_1)*(1-rp_2)...(1-rp_{i-1})*rp_{i}
   * @f]
   * @f[
   * p_n = (1-rp_1)*(1-rp_2)...(1-rp_{n-1})
   * @f]
   * and
   * @f[
   * \forall 1 <= i < n, rp_i = \frac{p_i}{1-(p_1+...+p_{i-1})}
   * @f]
   * where @f$p_i@f$ stands for the probability of model @f$i@f$.
   *
   *
   * If there are @f$n@f$ models, @f$\rho_i@f$ is the rate and @f$p_i@f$
   * is the probability of model i (@f$\sum_{i=1}^{n} p_i * \rho_i =
   * 1@f$), the rates are defined by relative rates parameters
   * @f$r_i@f$ (called "relratei") with:
   * @f[
   * 1 <= i < n, \rho_i = (1-r_1)*(1-r_2)...(1-r_{i-1})*\frac{r_{i}}{p_i} 
   * @f]
   * @f[
   * \rho_n = \frac{(1-r_1)*(1-r_2)*...*(1-r_{n-1})}{p_n}
   * @f]
   * and
   * @f[
   * \forall 1 <= i < n, r_i = \frac{\rho_i*p_i}{1-(p_1*\rho_1+...+p_{i-1}*\rho_{i-1})} < 1.
   * @f]
   *
   * For example:
   *
   * Mixture(model1=HKY85(kappa=3), model2=T92(theta=0.1),
   *         model2=L95(gamma=2), relrate1=0.2, relrate2=0.9,
   *         relproba1=0.1,relproba2=0.8)
   *
   * define a model as a mixture of 3 different models: HKY85 has
   * probability 0.1 and rate 2, T92 has probability 0.4 and rate 1.8,
   * and L95 has probability 0.5 and rate 0.16.
   *
   *
   * The parameters are named \c "Mixture.relrate1", \c
   * "Mixture.relrate2", \c "Mixture.relproba1", \c
   * "Mixture.relproba2"... in addition to the parameters of the
   * submodels that are prefixed by "Mixture.i_", where i is the order
   * of the model.

   */

  class MixtureOfSubstitutionModels :
    public AbstractMixedSubstitutionModel
  {
  public:

    /*
     *@brief Constructor of a MixtureOfSubstitutionModels, where all
     *the models have rate 1 and equal probability.
     *
     *@param alpha pointer to the Alphabet
     *@param vpModel vector of pointers to SubstitutionModels. All the
     *   SubstitutionModels are owned by the instance.
     */
    
    MixtureOfSubstitutionModels(const Alphabet* alpha,
                                std::vector<SubstitutionModel*> vpModel);
    
    /*
     *@brief Constructor of a MixtureOfSubstitutionModels.
     *
     *@param alpha pointer to the Alphabet
     *@param vpModel vector of pointers to SubstitutionModels. All the
     *   SubstitutionModels are owned by the instance.
     *@param vproba vector of the probabilities of the models
     *@param vrate vector of the rates of the models
     *
     * See above the constraints on the rates and the probabilities of
     * the vectors.
     */
    
    MixtureOfSubstitutionModels(const Alphabet* alpha,
                                std::vector<SubstitutionModel*> vpModel,
                                Vdouble& vproba, Vdouble& vrate);

    MixtureOfSubstitutionModels(const MixtureOfSubstitutionModels&);
    
    MixtureOfSubstitutionModels& operator=(const MixtureOfSubstitutionModels&);
    
    ~MixtureOfSubstitutionModels();
    
    MixtureOfSubstitutionModels* clone() const { return new MixtureOfSubstitutionModels(*this); }

  public:
    std::string getName() const { return "Mixture"; }

    void updateMatrices();
  
    /**
     * @brief Sets the rates of the submodels to follow the constraint
     * that the mean rate of the mixture equals rate_.
     
     * @param vd a vector of positive values such that the rates of
     * the respective submodels are in the same proportions (ie this
     * vector does not need to be normalized).
     */

    virtual void setVRates(const Vdouble& vd);

    /*
     *@brief Returns the vector of numbers of the submodels in the
     *mixture that match a description of the parameters numbers.
     *
     *@param desc is the description of the class indexes of the mixed
     *parameters. Syntax is like: kappa_1,gamma_3,delta_2
     *
     */
  
    Vint getSubmodelNumbers(std::string& desc) const;

    /**
     * @brief applies setFreq to all the models of the mixture and
     * recovers the parameters values.
     *
     **/
  
    void setFreq(std::map<int,double>&);

  };
} // end of namespace bpp.

#endif  // _MIXTUREOFSUBSTITUTIONMODELS_H_
