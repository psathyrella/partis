//
// File: SubstitutionModel.h
// Created by: Julien Dutheil
// Created on: Mon May 26 14:52:34 2003
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

#ifndef _SUBSTITUTIONMODEL_H_
#define _SUBSTITUTIONMODEL_H_

#include <cstdlib>
#include <map>
#include <string>

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

//From Seqlib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include "FrequenciesSet/FrequenciesSet.h"

namespace bpp
{
class SubstitutionModel;

/**
 * @brief Exception that may be thrown by susbstitution models.
 *
 * @see SubstitutionModel
 */
class SubstitutionModelException :
  public Exception
{
protected:
  const SubstitutionModel* model_;

public:
  SubstitutionModelException(const std::string& text, const SubstitutionModel* sm = 0);

  SubstitutionModelException(const SubstitutionModelException& sme) :
    Exception(sme), model_(sme.model_) {}

  SubstitutionModelException& operator=(const SubstitutionModelException& sme)
  {
    Exception::operator=(sme);
    model_ = sme.model_;
    return *this;
  }

  ~SubstitutionModelException() throw ();

public:
  /**
   * @brief Get the model that throws the exception.
   *
   * @return The model that throws the exception.
   */
  virtual const SubstitutionModel* getSubstitutionModel() const { return model_; }
};

/**
 * @brief Interface for all substitution models.
 *
 * A substitution model is based on a Markov generator \f$Q\f$, the size of
 * which depends on the alphabet used (4 for nucleotides, 20 for proteins, etc.).
 * Each SubstitutionModel object hence includes a pointer toward an alphabet,
 * and provides a method to retrieve the alphabet used (getAlphabet() method).
 *
 * What we want from a substitution model is to compute the probabilities of state
 * j at time t geven state j at time 0 (\f$P_{i,j}(t)\f$).
 * Typically, this is computed using the formula
 * \f[
 * P(t) = e^{r \times t \times Q},
 * \f]
 * where \f$P(t)\f$ is the matrix with all probabilities \f$P_{i,j}(t)\f$, and
 * \f$ r \f$ the rate.
 * For some models, the \f$P_{i,j}(t)\f$'s can be computed analytically.
 *
 * For more complex models, we need to use a eigen-decomposition of \f$Q\f$:
 * \f[ Q = U^{-1} . D . U, \f]
 * where \f$D = diag(\lambda_i)\f$ is a diagonal matrix.
 * Hence
 * \f[
 * P(t) = e^{r \times t \times Q} = U^{-1} . e^{r \times D \times t} . U,
 * \f]
 * where \f$e^{r \times D \times t} = diag\left(e^{r \times \lambda_i \times t}\right)\f$ is a
 * diagonal matrix with all terms equal to exp the terms in \f$D\f$ multiplied per \f$ r \times t \f$.
 * \f$U\f$ is the matrix of left eigen vectors (by row), and \f$U^{-1}\f$ is the matrix
 * of right eigen vectors (by column).
 * The values in \f$D\f$ are the eigen values of \f$Q\f$.
 * All \f$Q,U,U^{-1}\f$ and \f$D\f$ (its diagonal) may be retrieved from the
 * class (getEigenValues(), getRowRightEigenVectors() and getColumnLeftEigenVectors()
 * functions).
 *
 * First and second order derivatives of \f$P(t)\f$ with respect to \f$t\f$
 * can also be retrieved.
 * These methods may be useful for optimization processes.
 * Derivatives may be computed analytically, or using the general formulas:
 * \f[
 * \frac{\partial P(t)}{\partial t} =
 * U^{-1} . diag\left(r \times \lambda_i \times e^{r \times \lambda_i \times t}\right) . U
 * \f]
 * and
 * \f[
 * \frac{\partial^2 P(t)}{\partial t^2} =
 * U^{-1} . diag\left(r^2 \times \lambda_i^2 \times e^{r \times \lambda_i \times t}\right) . U
 * \f]
 *
 *
 * If Q is not symmetric, then the eigenvalue matrix D is block diagonal
 * with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
 * a + i*b, in 2-by-2 blocks, [a, b; -b, a].  That is, if the complex
 * eigenvalues look like
 * <pre>
 * 
 *           a + ib   .        .    .
 *           .        a - ib   .    .
 *           .        .        x    .
 *           .        .        .    y
 * </pre>
 * then D looks like
 * <pre>
 * 
 *           a          b      .    .
 *          -b          a      .    .
 *           .          .      x    .
 *           .          .      .    y
 * </pre>
 *
 * and exp(tD) equals
 * <pre>
 * 
 *           exp(ta)cos(tb)   exp(ta)sin(tb)  .        .
 *          -exp(ta)sin(tb)   exp(ta)cos(tb)  .        . 
 *           .                .               exp(tx)  .
 *           .                .               .        exp(ty)
 * </pre>
 *
 *
 *
 * If U is singular, it cannot be inverted. In this case exp(tQ) is
 * approximated using Taylor decomposition:
 *
 * \f[
 * P(t) = Id + tQ + \frac{(tQ)^2}{2!} + ... + \frac{(tQ)^n}{n!} + ... 
 * \f]
 *
 * To prevent approximation issues, if @\f$ max(tQ) @\f$ is too high
 * (currently above 0.5), @\f$ t @\f$ is divided in an ad hoc way
 * (e.g. by @\f$ N @\f$), and we compute @\f$ P(t) = (P(t/N))^N @\f$
 * with a Taylor decomposition for @\f$ P(t/N) @\f$.
 *
 * In this case, derivatives according to @\f$ t @\f$ are computed
 * analytically too.
 *
 */

class SubstitutionModel :
  public virtual ParameterAliasable
{
public:
  SubstitutionModel() {}
  virtual ~SubstitutionModel() {}

#ifndef NO_VIRTUAL_COV
  SubstitutionModel* clone() const = 0;
#endif

public:
  /**
   * @brief Get the name of the model.
   *
   * @return The name of this model.
   */
  virtual std::string getName() const = 0;

  /**
   * @return The supported states of the model, as a vector of int codes.
   *
   * @see Alphabet
   */
  virtual const std::vector<int>& getAlphabetChars() const = 0;

  /**
   * @brief Get the char in the alphabet corresponding to a given state in the model.
   *
   * In most cases, this method will return i.
   * @param i The index of the state.
   * @return The corresponding state in the alphabet.
   * @see MarkovModulatedSubstitutionModel
   * @see getStates()
   */
  virtual int getAlphabetChar(size_t i) const = 0;

  /**
   * @brief Get the state in the model corresponding to a particular char in the alphabet.
   *
   * @param i The alphabet char to check.
   * @return A vector of indices of model states.
   */
  virtual std::vector<size_t> getModelStates(int i) const = 0;

  /**
   * @return Equilibrium frequency associated to character i.
   * @see getFrequencies(), getStates()
   */
  virtual double freq(size_t i) const = 0;

  /**
   * @return The rate in the generator of change from state i to state j.
   *
   * @see getStates();
   */
  virtual double Qij(size_t i, size_t j) const = 0;

  /**
   * @return The probability of change from state i to state j during time t.
   * @see getPij_t(), getStates()
   */
  virtual double Pij_t(size_t i, size_t j, double t) const = 0;

  /**
   * @return The first order derivative of the probability of change from state
   * i to state j with respect to time t, at time t.
   * @see getdPij_dt(), getStates()
   */
  virtual double dPij_dt(size_t i, size_t j, double t) const = 0;

  /**
   * @return The second order derivative of the probability of change from state
   * i to state j with respect to time t, at time t.
   * @see getd2Pij_dt2(), getStates()
   */
  virtual double d2Pij_dt2(size_t i, size_t j, double t) const = 0;

  /**
   * @return A vector of all equilibrium frequencies.
   * @see freq()
   */
  virtual const Vdouble& getFrequencies() const = 0;

  /**
   * @return The normalized Markov generator matrix, i.e. all
   * normalized rates of changes from state i to state j. The
   * generator is normalized so that
   * (i) \f$ \forall i; \sum_j Q_{i,j} = 0 \f$, meaning that
   * $\f$ \forall i; Q_{i,i} = -\sum_{j \neq i}Q_{i,j}\f$, and
   * (ii) \f$ \sum_i Q_{i,i} \times \pi_i = -1\f$.
   * This means that, under normalization, the mean rate of replacement at
   * equilibrium is 1 and that time \f$t\f$ are measured in units of
   * expected number of changes per site. Additionnaly, the rate_ attibute provides
   * the possibility to increase or decrease this mean rate.
   *
   * See Kosiol and Goldman (2005), Molecular Biology And Evolution 22(2) 193-9.
   * @see Qij()
   */
  virtual const Matrix<double>& getGenerator() const = 0;

  /**
   * @return The matrix of exchangeability terms.
   * It is recommended that exchangeability matrix be normalized so that the normalized
   * generator be obtained directly by the dot product \f$S . \pi\f$.
   */
  virtual const Matrix<double>& getExchangeabilityMatrix() const = 0;

  /**
   * @return The exchangeability between state i and state j.
   *
   * By definition Sij(i,j) = Sij(j,i).
   */
  
  virtual double Sij(size_t i, size_t j) const = 0;
  /**
   * @return All probabilities of change from state i to state j during time t.
   * @see Pij_t()
   */
  virtual const Matrix<double>& getPij_t(double t) const = 0;

  /**
   * @return Get all first order derivatives of the probability of change from state
   * i to state j with respect to time t, at time t.
   * @see dPij_dt()
   */
  virtual const Matrix<double>& getdPij_dt(double t) const = 0;

  /**
   * @return All second order derivatives of the probability of change from state
   * i to state j with respect to time t, at time t.
   * @see d2Pij_dt2()
   */
  virtual const Matrix<double>& getd2Pij_dt2(double t) const = 0;

  /**
   * @brief Set if eigenValues and Vectors must be computed
   */
  virtual void enableEigenDecomposition(bool yn) = 0;

  /**
   * @brief Tell if eigenValues and Vectors must be computed
   */
  virtual bool enableEigenDecomposition() = 0;

  /**
   * @return A vector with all real parts of the eigen values of the generator of this model;
   */
  virtual const Vdouble& getEigenValues() const = 0;

  /**
   * @return A vector with all imaginary parts of the eigen values of the generator of this model;
   */
  virtual const Vdouble& getIEigenValues() const = 0;

  /**
   * @return True if the model is diagonalizable in R.
   */
  virtual bool isDiagonalizable() const = 0;
  
  /**
   * @return True is the model is non-singular.
   */
  virtual bool isNonSingular() const = 0;

  /**
   * @return A matrix with left eigen vectors.
   * Each row in the matrix stands for an eigen vector.
   */
  virtual const Matrix<double>& getRowLeftEigenVectors() const = 0;

  /**
   * @return A matrix with right eigen vectors.
   * Each column in the matrix stands for an eigen vector.
   */
  virtual const Matrix<double>& getColumnRightEigenVectors() const = 0;

  /**
   * @return Get the alphabet associated to this model.
   */
  virtual const Alphabet* getAlphabet() const = 0;

  /**
   * @brief Get the number of states.
   *
   * For most models, this equals the size of the alphabet.
   *
   * @return The number of different states in the model.
   */
  virtual size_t getNumberOfStates() const = 0;

  /**
   * This method is used to initialize likelihoods in reccursions.
   * It typically sends 1 if i = state, 0 otherwise, where
   * i is one of the possible states of the alphabet allowed in the model
   * and state is the observed state in the considered sequence/site.
   *
   * @param i the index of the state in the model.
   * @param state An observed state in the sequence/site.
   * @return 1 or 0 depending if the two states are compatible.
   * @throw IndexOutOfBoundsException if array position is out of range.
   * @throw BadIntException if states are not allowed in the associated alphabet.
   * @see getStates();
   */
  virtual double getInitValue(size_t i, int state) const throw (IndexOutOfBoundsException, BadIntException) = 0;

  /**
   * @brief Get the scalar product of diagonal elements of the generator
   * and the frequencies vector.
   * If the generator is normalized, then scale=1. Otherwise each element
   * must be multiplied by 1/scale.
   *
   * @return Minus the scalar product of diagonal elements and the frequencies vector.
   */
  virtual double getScale() const = 0;

  /**
   * 
   * @brief Multiplies the current generator by the given scale.
   *
   * @param scale the scale by which the generator is multiplied.
   *
   */

  virtual void setScale(double scale) = 0;

  /**
   * @brief Get the rate
   */
  virtual double getRate() const = 0;

  /**
   * @brief Set the rate of the model (must be positive).
   * @param rate must be positive.
   */
  
  virtual void setRate(double rate) = 0;

  virtual void addRateParameter() = 0;
  
  /**
   * @brief Set equilibrium frequencies equal to the frequencies estimated
   * from the data.
   *
   * @param data The sequences to use.
   * @param pseudoCount A quantity @f$\psi@f$ to add to adjust the observed
   *   values in order to prevent issues due to missing states on small data set.
   * The corrected frequencies shall be computed as
   * @f[
   * \pi_i = \frac{n_i+\psi}{\sum_j (f_j+\psi)}
   * @f]
   */
  virtual void setFreqFromData(const SequenceContainer& data, double pseudoCount = 0) = 0;

  /**
   * @brief Set equilibrium frequencies
   *
   * @param frequencies The map of the frequencies to use.
   */
  virtual void setFreq(std::map<int, double>& frequencies) {}

  /**
   * @brief If the model owns a FrequenciesSet, returns a pointer to
   * it, otherwise return 0.
   *
   */

  virtual const FrequenciesSet* getFrequenciesSet() const {return NULL;}
};


/**
 * @brief Interface for reversible substitution models.
 *
 * For reversible models,
 * \f[ Q = S . \pi, \f]
 * where \f$S\f$ is a symetric matrix called the exchangeability matrix, and
 * \f$\Pi\f$ the diagonal matrix with all equilibrium frequencies.
 * The frequencies may be retrieved as a vector by the getFrequencies() method
 * or individually by the freq() method.
 * The \f$S\f$ matrix may be obtained by the getExchangeabilityMatrix().
 */
class ReversibleSubstitutionModel :
  public virtual SubstitutionModel
{
public:
  ReversibleSubstitutionModel() {}
  virtual ~ReversibleSubstitutionModel() {}

#ifndef NO_VIRTUAL_COV
  ReversibleSubstitutionModel* clone() const = 0;
#endif
};

} //end of namespace bpp.

#endif  //_SUBSTITUTIONMODEL_H_

