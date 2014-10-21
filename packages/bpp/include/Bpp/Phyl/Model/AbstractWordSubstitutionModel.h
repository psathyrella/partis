//
// File: AbstractWordSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: Jan 2009
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

#ifndef _ABSTRACTWORDSUBSTITUTIONMODEL_H_
#define _ABSTRACTWORDSUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief Abstract Basal class for words of substitution models.
 * @author Laurent Guéguen
 *
 * Objects of this class are built from several substitution models.
 * Each model corresponds to a position in the word. No model is
 * directly accessible.
 *
 * Only substitutions with one letter changed are accepted.
 *
 * There is one substitution per word per unit of time
 * on the equilibrium frequency, and each position has its specific rate.
 * For example, if there are @f$n@f$ models and \f$\rho_i\f$ is the rate of
 * model i (@f$\sum_{i=1}^{n} \rho_i = 1@f$):
 * @f{eqnarray*}
 * Q_{abc \rightarrow abd} &=& \rho_2 Q^{(2)}_{c \rightarrow d}\\
 * Q_{abc \rightarrow aed} &=& 0\\
 * Q_{abc \rightarrow abc} &=& \rho_0 Q^{(0)}_{a \rightarrow a} + \rho_1 Q^{(1)}_{b \rightarrow b} + \rho_2 Q^{(2)}_{c \rightarrow c})
 * @f}
 *
 * The parameters of this word model are the same as the ones of the
 * models used. Their names have a new suffix, "phi_" where i stands
 * for the position (i.e. the phase) in the word.
 *
 */
class AbstractWordSubstitutionModel :
    public virtual AbstractSubstitutionModel
{
private:
  /**
   * @ brief boolean flag to check if a specific WordAlphabet has been built
   */
  bool new_alphabet_;

protected:
  std::vector<SubstitutionModel*> VSubMod_;
  std::vector<std::string> VnestedPrefix_;

  std::vector<double> Vrate_;

protected:
  static Alphabet* extractAlph(const std::vector<SubstitutionModel*>& modelVector);

protected:
  void updateMatrices();

  /**
   * @brief Called by updateMatrices to handle specific modifications
   * for inheriting classes
   */
  virtual void completeMatrices() = 0;

public:
  /**
   * @brief Build a new AbstractWordSubstitutionModel object from a
   * vector of pointers to SubstitutionModels.
   *
   * @param modelVector the vector of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param st the Namespace.
   */
  AbstractWordSubstitutionModel(
    const std::vector<SubstitutionModel*>& modelVector,
    const std::string& st);

  /**
   * @brief Build a new AbstractWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of
   * desired models.
   *
   * @param pmodel A pointer to the substitution model to use in all
   * the positions. It will be owned by the instance.
   * @param num The number of models involved.
   * @param st the Namespace.
   */
  AbstractWordSubstitutionModel(
    SubstitutionModel* pmodel,
    unsigned int num,
    const std::string& st);

  AbstractWordSubstitutionModel(const AbstractWordSubstitutionModel&);

  AbstractWordSubstitutionModel& operator=(const AbstractWordSubstitutionModel&);

  virtual ~AbstractWordSubstitutionModel();

  void setNamespace(const std::string& prefix);

protected:
  /**
   *@brief Constructor for the derived classes only
   */
  AbstractWordSubstitutionModel(const Alphabet* alph, const std::string&);

public:
  virtual size_t getNumberOfStates() const;

  /**
   * @brief returns the ith model, or Null if i is not a valid number.
   *
   */
  
  const SubstitutionModel* getNModel(size_t i) const {
    if (i< VSubMod_.size())
        return dynamic_cast<const SubstitutionModel*>(VSubMod_[i]);
    else
      return 0;
  }

  size_t getNumberOfModels() const {
    return VSubMod_.size();
  }
  
  /**
   *@brief Estimation of the parameters of the models so that the
   *equilibrium frequencies match the given ones.
   *
   *@param freqs  map of the frequencies
   *
   * When there is one submodel for all the positions, the submodel
   * parameters are fit on the means of the frequencies on each
   * position. Otherwise, each model is fit on the frequencies on its
   * corresponding position in the word.
   *
   **/
  
  virtual void setFreq(std::map<int, double>& freqs);
};
} // end of namespace bpp.

#endif  // ABSTRACTWORDSUBSTITUTIONMODEL_

