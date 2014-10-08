//
// File: AbstractBiblioSubstitutionModel.h
// Created by: Laurent Guéguen
// Created on: vendredi 8 juillet 2011, à 20h 17
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

#ifndef _ABSTRACTBIBLIOSUBSTITUTIONMODEL_H_
#define _ABSTRACTBIBLIOSUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"
#include "AbstractSubstitutionModel.h"

#include <Bpp/Numeric/AbstractParameterAliasable.h>

namespace bpp
{
/**
 * @brief Partial implementation of the SubstitutionModel interface
 *   for models that are set for matching the bibliography, and are
 *   only defined through a link to a "real" model.
 *
 */

class AbstractBiblioSubstitutionModel :
  public virtual SubstitutionModel,
  public AbstractParameterAliasable
{
protected:
  /**
   * @brief Tools to make the link between the Parameters of the
   * object and those of pmixmodel_.
   *
   */

  std::map<std::string, std::string> mapParNamesFromPmodel_;

  ParameterList lParPmodel_;

public:
  AbstractBiblioSubstitutionModel(const std::string& prefix);

  AbstractBiblioSubstitutionModel(const AbstractBiblioSubstitutionModel& model);

  AbstractBiblioSubstitutionModel& operator=(const AbstractBiblioSubstitutionModel& model);

  virtual ~AbstractBiblioSubstitutionModel() {}

#ifndef NO_VIRTUAL_COV
  virtual AbstractBiblioSubstitutionModel* clone() const = 0;
#endif

public:
  virtual const SubstitutionModel& getModel() const = 0;

  /*
     *@ brief Methods to supersede SubstitutionModel methods.
   *
   * @{
   */

  const std::vector<int>& getAlphabetChars() const { return getModel().getAlphabetChars(); }

  int getAlphabetChar(size_t i) const { return getModel().getAlphabetChar(i); }

  std::vector<size_t> getModelStates(int i) const { return getModel().getModelStates(i); }

  virtual double freq(size_t i) const { return getModel().freq(i); }

  virtual double Qij(size_t i, size_t j) const { return getModel().Qij(i, j); }

  virtual double Pij_t    (size_t i, size_t j, double t) const { return getModel().Pij_t(i, j, t); }
  virtual double dPij_dt  (size_t i, size_t j, double t) const { return getModel().dPij_dt (i, j, t); }
  virtual double d2Pij_dt2(size_t i, size_t j, double t) const { return getModel().d2Pij_dt2(i, j, t); }

  virtual const Vdouble& getFrequencies() const { return getModel().getFrequencies(); }

  const Matrix<double>& getGenerator() const { return getModel().getGenerator(); }

  const Matrix<double>& getExchangeabilityMatrix() const { return getModel().getExchangeabilityMatrix(); }

  double Sij(size_t i, size_t j) const { return getModel().Sij(i, j); }

  const Matrix<double>& getPij_t(double t) const { return getModel().getPij_t(t); }

  const Matrix<double>& getdPij_dt(double t) const { return getModel().getdPij_dt(t); }

  const Matrix<double>& getd2Pij_dt2(double t) const { return getModel().getd2Pij_dt2(t); }

  void enableEigenDecomposition(bool yn) { getModel().enableEigenDecomposition(yn); }

  bool enableEigenDecomposition() { return getModel().enableEigenDecomposition(); }

  bool isDiagonalizable() const { return getModel().isDiagonalizable(); }

  bool isNonSingular() const { return getModel().isNonSingular(); }

  const Vdouble& getEigenValues() const { return getModel().getEigenValues(); }

  const Vdouble& getIEigenValues() const { return getModel().getIEigenValues(); }

  const Matrix<double>& getRowLeftEigenVectors() const { return getModel().getRowLeftEigenVectors(); }
  const Matrix<double>& getColumnRightEigenVectors() const { return getModel().getColumnRightEigenVectors(); }

  double getRate() const { return getModel().getRate(); }

  void setRate(double rate) { return getModel().setRate(rate); }

  void addRateParameter();

  void setFreqFromData(const SequenceContainer& data, double pseudoCount = 0);

  void setFreq(std::map<int, double>& frequ);

  const Alphabet* getAlphabet() const { return getModel().getAlphabet(); }

  size_t getNumberOfStates() const { return getModel().getNumberOfStates(); }

  double getInitValue(size_t i, int state) const throw (BadIntException) { return getModel().getInitValue(i, state); }

  const FrequenciesSet* getFrequenciesSet() const {return getModel().getFrequenciesSet(); }

  /*
   * @}
   *
   */

  /*
     *@ brief Methods to supersede AbstractSubstitutionModel methods.
   *
   * @{
   */

  /**
   * @brief Tells the model that a parameter value has changed.
   *
   * This updates the matrices consequently.
   */
  virtual void fireParameterChanged(const ParameterList& parameters)
  {
    AbstractParameterAliasable::fireParameterChanged(parameters);
    if (parameters.size() > 1 || (parameters.size() == 1 && parameters[0].getName() != getNamespace() + "rate")) {
      updateMatrices();
    }
  }

protected:
  virtual void updateMatrices();
  virtual SubstitutionModel& getModel() = 0;

public:
  double getScale() const { return getModel().getScale(); }

  void setScale(double scale) { getModel().setScale(scale); }

  /*
   * @}
   */
};
} // end of namespace bpp.


#endif  // _ABSTRACTBIBLIOSUBSTITUTIONMODEL_H_

