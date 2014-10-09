//
// File: WordSubstitutionModel.h
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

#ifndef _WORDSUBSTITUTIONMODEL_H_
#define _WORDSUBSTITUTIONMODEL_H_

#include "AbstractWordSubstitutionModel.h"
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/BppVector.h>

namespace bpp
{
/**
 * @brief Basal class for words of reversible substitution models.
 * @author Laurent Guéguen
 *
 * Only substitutions with one letter changed are accepted. Hence the
 * equilibrium frequency of each word is the product of the
 * equilibrium frequencies of the letters.</p>
 *
 * If there are @f$n@f$ models, @f$\rho_i@f$ is the rate of
 * model i (@f$\sum_{i=1}^{n} \rho_i = 1@f$) and the rates
 * are defined by relative rates parameters @f$r_i@f$
 * (called "relratei") with:
 * @f[
 * 1 <= i < n, \rho_i = (1-r_1).(1-r_2)...(1-r_{i-1}).r_{i}
 * @f]
 * @f[
 * \rho_n = (1-r_1).(1-r_2)...(1-r_{n-1})
 * @f]
 * and
 * @f[
 * \forall 1 <= i < n, r_i = \frac{\rho_i}{1-(\rho_1+...\rho_{i-1})}
 * @f]
 * where @f$\rho_i@f$ stands for the rate of position @f$i@f$.
 */

class WordSubstitutionModel :
  public AbstractWordSubstitutionModel
{
public:
  /**
   * @brief Build a new WordSubstitutionModel object from a
   * Vector of pointers to SubstitutionModels.
   *
   * @param modelVector the Vector of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param st the Namespace.
   */
  WordSubstitutionModel(const std::vector<SubstitutionModel*>& modelVector, const std::string& st = "");

  /**
   * @brief Build a new WordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of
   * desired models.
   *
   * @param pmodel pointer to the substitution model to use in all the
   *  positions. It is owned by the instance.
   * @param num The number of models involved.
   * @param st the Namespace.
   */
  WordSubstitutionModel(SubstitutionModel* pmodel, unsigned int num, const std::string& st = "");

  virtual ~WordSubstitutionModel() {}

  WordSubstitutionModel* clone() const { return new WordSubstitutionModel(*this); }

protected:
  /**
   *@brief Constructor for the derived classes only
   */

  WordSubstitutionModel(const Alphabet* alph, const std::string& = "");

  virtual void updateMatrices();
  virtual void completeMatrices();

public:
  virtual const RowMatrix<double>& getPij_t(double d) const;

  virtual const RowMatrix<double>& getdPij_dt(double d) const;

  virtual const RowMatrix<double>& getd2Pij_dt2(double d) const;

  virtual std::string getName() const;
};
} // end of namespace bpp.

#endif  // _WORDSUBSTITUTIONMODEL

