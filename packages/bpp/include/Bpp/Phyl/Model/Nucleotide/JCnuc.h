//
// File: JCnuc.h
// Created by: Julien Dutheil
// Created on: Tue May 27 16:04:36 2003
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

#ifndef _JCNUC_H_
#define _JCNUC_H_

#include "NucleotideSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Jukes-Cantor substitution model for nucleotides.
 *
 * All rates equal:
 * \f[
 * S = \begin{pmatrix}
 * \cdots & r & r & r \\
 * r & \cdots & r & r \\
 * r & r & \cdots & r \\
 * r & r & r & \cdots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = diag\left(\frac{1}{4}, \frac{1}{4}, \frac{1}{4}, \frac{1}{4}\right)
 * \f]
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \begin{pmatrix}
 * -4 & \frac{4}{3} & \frac{4}{3} & \frac{4}{3} \\
 * \frac{4}{3} & -4 & \frac{4}{3} & \frac{4}{3} \\
 * \frac{4}{3} & \frac{4}{3} & -4 & \frac{4}{3} \\
 * \frac{4}{3} & \frac{4}{3} & \frac{4}{3} & -4 \\
 * \end{pmatrix}
 * \f]
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$pi\f$:
 * \f[
 * Q = S . \pi = \begin{pmatrix}
 * -1 & \frac{1}{3} & \frac{1}{3} & \frac{1}{3} \\
 * \frac{1}{3} & -1 & \frac{1}{3} & \frac{1}{3} \\
 * \frac{1}{3} & \frac{1}{3} & -1 & \frac{1}{3} \\
 * \frac{1}{3} & \frac{1}{3} & \frac{1}{3} & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * The eigen values are \f$\left(0, -\frac{4}{3}, -\frac{4}{3}, -\frac{4}{3}\right)\f$,
 * the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 *  \frac{1}{4} &  \frac{1}{4} &  \frac{1}{4} &  \frac{1}{4} \\
 * -\frac{1}{4} & -\frac{1}{4} &  \frac{3}{4} & -\frac{1}{4} \\
 * -\frac{1}{4} &  \frac{3}{4} & -\frac{1}{4} & -\frac{1}{4} \\
 *  \frac{3}{4} & -\frac{1}{4} & -\frac{1}{4} & -\frac{1}{4} \\
 * \end{pmatrix}
 * \f]
 * and the right eigen vectors are, by column:
 * \f[
 * U^-1 = \begin{pmatrix}
 * 1 &  0 &  0 &  1 \\
 * 1 &  0 &  1 &  0 \\
 * 1 &  1 &  0 &  0 \\
 * 1 & -1 & -1 & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * In addition, a rate_ factor defines the mean rate of the model.
 *
 * The probabilities of changes are computed analytically using the formulas:
 * \f[
 * P_{i,j}(t) = \begin{cases}
 * \frac{1}{4} + \frac{3}{4}e^{-rate\_*\frac{4}{3}t} & \text{if $i=j$}, \\
 * \frac{1}{4} - \frac{1}{4}e^{-rate\_*\frac{4}{3}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 * First and second order derivatives are also computed analytically using the formulas:
 * \f[
 * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \begin{cases}
 * -e^{-rate\_*\frac{4}{3}t}           & \text{if $i=j$}, \\
 * \frac{1}{3}e^{-rate*\frac{4}{3}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 * \f[
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} = rate\_^2 *\begin{cases}
 * \frac{4}{3}e^{-rate\_*\frac{4}{3}t}  & \text{if $i=j$}, \\
 * -\frac{4}{9}e^{-rate\_*\frac{4}{3}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 * Reference:
 * - Jukes TH and Cantor CR (1969), Evolution_ of proteins molecules_, 121-123, in Mammalian_ protein metabolism_.
 */
class JCnuc :
  public virtual NucleotideSubstitutionModel,
  public AbstractReversibleSubstitutionModel
{
private:
  mutable double exp_;
  mutable RowMatrix<double> p_;

public:
  JCnuc(const NucleicAlphabet* alpha);

  virtual ~JCnuc() {}

#ifndef NO_VIRTUAL_COV
  JCnuc*
#else
  Clonable*
#endif
  clone() const { return new JCnuc(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const;
  double dPij_dt  (size_t i, size_t j, double d) const;
  double d2Pij_dt2(size_t i, size_t j, double d) const;
  const Matrix<double>& getPij_t    (double d) const;
  const Matrix<double>& getdPij_dt  (double d) const;
  const Matrix<double>& getd2Pij_dt2(double d) const;

  std::string getName() const { return "JC69"; }

  /**
   * @brief This method is disabled in this model since frequencies are not free parameters.
   *
   * Consider using the HKY85 model for instance if you want to set frequencies as parameters.
   *
   * @param data Useless parameter.
   */
  void setFreqFromData(const SequenceContainer& data) {}

protected:
  /**
   * In the case of the model of Jukes & Cantor, this method is not usefull since
   * the generator is fully determined. No matrice can be changed... This method is only
   * used in the constructor of the class.
   */
  void updateMatrices();
};
} // end of namespace bpp.

#endif  // _JCNUC_H_

