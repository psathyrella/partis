//
// File: MarkovModulatedSubstitutionModel.h
// Created by: Julien Dutheil
// Created on: Sat Aug 05 08:21 2006
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

#ifndef _MARKOVMODULATEDSUBSTITUTIONMODEL_H_
#define _MARKOVMODULATEDSUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

namespace bpp
{

  /**
   * @brief Partial implementation of the Markov-modulated class of substitution
   * models.
   *
   * This class wraps a substitution model and provide a matrix describing rate changes.
   * The rate matrix must be initialized by derived classes of this class.
   * Using these matrices, the diagonalization procedure of Galtier and Jean-Marie is used.
   *
   * Such models can be described using two matrices:
   * a substitution matrix, @f$M@f$, with size @f$m@f$, which is a "standard" substitution model of any alphabet type,
   * and a rate matrix @f$G@f$ of size @f$g@f$.
   * The generator of the markov-modulated model, @f$Q@f$ can be written using Kronecker matrix operands:
   * @f[
   * Q=D_R \otimes M + G \otimes I_m,
   * @f]
   * where @f$D_R@f$ is the diagonal matrix of all rates, and @f$I_m@f$ is the identity matrix of size @f$m@f$.
   *
   * This generator is normalized so that branch lengths are measured in unit of mean number of substitutions per site,
   * where susbstitution here means "change of alphabet state".
   * Rate changes are not counted.
   *
   * Galtier N. and Jean-Marie A., Markov-modulated Markov chains and the covarion process of molecular evolution (2004).
   * _Journal of Computational Biology_, 11:727-33.
   */
  class MarkovModulatedSubstitutionModel:
    public virtual ReversibleSubstitutionModel,
    public AbstractParameterAliasable
  {

  protected:
    ReversibleSubstitutionModel* model_;
    size_t nbStates_; //Number of states in model
    size_t nbRates_; //Number of rate classes

    /**
     * @name Rate generator.
     *
     * These variables must be initialized in the constructor of the derived class.
     * @{
     */
    RowMatrix<double> rates_;                //All rates values
    RowMatrix<double> ratesExchangeability_; //All rates transitions
    Vdouble           ratesFreq_;            //All rates equilibrium frequencies
    /**@}*/
    RowMatrix<double> ratesGenerator_;       //All rates transitions
    
    /**
     * @brief The list of supported chars.
     */
    std::vector<int> chars_;

    /**
     * @brief The generator matrix \f$Q\f$ of the model.
     */
    RowMatrix<double> generator_;

    /**
     * @brief The exchangeability matrix \f$S\f$ of the model.
     */
    RowMatrix<double> exchangeability_;

    /**
     * @brief The \f$U\f$ matrix made of left eigen vectors (by row).
     */
    RowMatrix<double> leftEigenVectors_;

    /**
     * @brief The \f$U^-1\f$ matrix made of right eigen vectors (by column).
     */
    RowMatrix<double> rightEigenVectors_;

    /**
     * @brief The vector of real parts of eigen values.
     */
    Vdouble eigenValues_;

    /**
     * @brief The vector of imaginary parts of the eigen values (zero
     * in case of reversible pmodel).
     */
    Vdouble iEigenValues_;

    /**
     * @brief Tell if the eigen decomposition should be performed.
     */
    bool eigenDecompose_;

    /**
     * @brief These ones are for bookkeeping:
     */
    mutable RowMatrix<double> pijt_;
    mutable RowMatrix<double> dpijt_;
    mutable RowMatrix<double> d2pijt_;

    /**
     * @brief The vector of equilibrium frequencies.
     */
    Vdouble freq_;

    bool normalizeRateChanges_;

    std::string nestedPrefix_;

  public:
    /**
     * @brief Build a new MarkovModulatedSubstitutionModel object.
     *
     * @param model The substitution model to use. Can be of any alphabet type, and will be owned by this instance.
     * @param normalizeRateChanges Tells if the branch lengths must be computed in terms of rate and state
     * @param prefix The parameter namespace to be forwarded to the AbstractParametrizable constructor.
     * changes instead of state change only.
     * NB: In most cases, this parameter should be set to false.
     */
    MarkovModulatedSubstitutionModel(ReversibleSubstitutionModel* model, bool normalizeRateChanges, const std::string& prefix) :
      AbstractParameterAliasable(prefix),
      model_(model), nbStates_(model->getNumberOfStates()), nbRates_(0), rates_(), ratesExchangeability_(),
      ratesFreq_(), ratesGenerator_(), chars_(), generator_(), exchangeability_(),
      leftEigenVectors_(), rightEigenVectors_(), eigenValues_(), iEigenValues_(), eigenDecompose_(true), 
      pijt_(), dpijt_(), d2pijt_(), freq_(),
      normalizeRateChanges_(normalizeRateChanges),
      nestedPrefix_("model_" + model->getNamespace())
    {
      model_->setNamespace(prefix + nestedPrefix_);
      addParameters_(model_->getIndependentParameters());
    }
    
    MarkovModulatedSubstitutionModel(const MarkovModulatedSubstitutionModel& model);
    MarkovModulatedSubstitutionModel& operator=(const MarkovModulatedSubstitutionModel& model);

    virtual ~MarkovModulatedSubstitutionModel() { delete model_; }

#ifndef NO_VIRTUAL_COV
    MarkovModulatedSubstitutionModel*
#else
    Clonable*
#endif
    clone() const = 0;

  public:
	  
    const Alphabet* getAlphabet() const { return model_->getAlphabet(); }

    size_t getNumberOfStates() const { return nbStates_ * nbRates_; }

    const Vdouble& getFrequencies() const { return freq_; }
    
    const Matrix<double>& getExchangeabilityMatrix() const { return exchangeability_; }
    
    const Matrix<double>& getGenerator() const { return generator_; }
    
    const Matrix<double>& getPij_t(double t) const;
    const Matrix<double>& getdPij_dt(double t) const;
    const Matrix<double>& getd2Pij_dt2(double t) const;
    
    const Vdouble& getEigenValues() const { return eigenValues_; }
    const Vdouble& getIEigenValues() const { return iEigenValues_; }

    bool isDiagonalizable() const { return true; }
    bool isNonSingular() const { return true; }

    const Matrix<double>& getRowLeftEigenVectors() const { return leftEigenVectors_; }
    const Matrix<double>& getColumnRightEigenVectors() const { return rightEigenVectors_; }
    
    double freq(size_t i) const { return freq_[i]; }
    double Sij(size_t i, size_t j) const { return exchangeability_(i, j); }
    double Qij(size_t i, size_t j) const { return generator_(i, j); }
    
    double Pij_t    (size_t i, size_t j, double t) const { return getPij_t(t)(i, j); }
    double dPij_dt  (size_t i, size_t j, double t) const { return getdPij_dt(t)(i, j); }
    double d2Pij_dt2(size_t i, size_t j, double t) const { return getd2Pij_dt2(t)(i, j); }
    
    double getInitValue(size_t i, int state) const throw (IndexOutOfBoundsException, BadIntException);
    
    void setFreqFromData(const SequenceContainer& data, double pseudoCount = 0)
    {
      model_->setFreqFromData(data, pseudoCount);
      updateMatrices();
    }

    const std::vector<int>& getAlphabetChars() const
    {
      return chars_;
    }

    int getAlphabetChar(size_t i) const
    {
      return chars_[i]; 
    }
   
    std::vector<size_t> getModelStates(int i) const
    {
      std::vector<size_t> states(nbRates_ * nbStates_);
      std::vector<size_t> nestedStates = model_->getModelStates(i);
      for(size_t j = 0; j < nbRates_; j++)
        for(size_t k = 0; k < nestedStates.size(); k++)
          states.push_back(j * nbRates_ + states[k]);
      return states;
    }

    const ReversibleSubstitutionModel* getNestedModel() const { return model_; }

    /**
     * @brief Get the rate category corresponding to a particular state in the compound model.
     *
     * @param i The state.
     * @return The corresponding rate category.
     * @see getState;
     */
    size_t getRate(size_t i) const
    {
      return i / nbStates_; 
    }

    double getRate() const { return model_->getRate(); }

    void setRate(double rate) { model_->setRate(rate); }

    double getScale() const
    {
      std::vector<double> v;
      MatrixTools::diag(generator_, v);
      return -VectorTools::scalar<double, double>(v, freq_);
    }

    void setScale(double scale) {
      model_->setScale(scale);
      updateMatrices();
    } 

    void enableEigenDecomposition(bool yn) { eigenDecompose_ = yn; }

    bool enableEigenDecomposition() { return eigenDecompose_; }
	
    /**
     * @brief Tells the model that a parameter value has changed.
     *
     * This updates the matrices consequently.
     */
    virtual void fireParameterChanged(const ParameterList& parameters)
    {
      AbstractParameterAliasable::fireParameterChanged(parameters);
      model_->matchParametersValues(parameters);
      updateRatesModel();
      updateMatrices();
    }
   
    void setNamespace(const std::string& prefix);

  protected:
    
    virtual void updateMatrices();

    /**
     * @brief Update the rates vector, generator and equilibrium frequencies.
     *
     * This method must be implemented by the derived class.
     * It is called by the fireParameterChanged() method.
     */
    virtual void updateRatesModel() = 0;

  };

} //end of namespace bpp.

#endif //_MARKOVMODULATEDSUBSTITUTIONMODEL_H_

