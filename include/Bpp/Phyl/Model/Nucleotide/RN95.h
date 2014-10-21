//
// File: RN95.h
// Created by: Laurent Guéguen
// Created on: jeudi 24 février 2011, à 20h 43
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

#ifndef _RN95_H_
#define _RN95_H_

#include "NucleotideSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"

#include <Bpp/Numeric/Constraints.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief The model described by Rhetsky \& Nei, where the only
 * hypothesis is that the transversion rates are only dependent of
 * the target nucleotide. This model is not reversible.
 *
 * This model has been thoroughly studied by Schadt & al, and we
 * follow their notations and formula. The parameters are defined to
 * allow allow the direct definition of the stationnary frequencies.
 *
 * After normalization this model has 7 parameters:
 * \f[
 * Q= \frac 1P
 * \begin{pmatrix}
 * . & \gamma & \alpha & \lambda \\
 * \delta & . &  \kappa & \beta \\
 * \epsilon & \gamma & . & \lambda \\
 * \delta & \sigma & \kappa & .\\
 * \end{pmatrix}\f]
 *
 * so in the parametrization process we set: \f[\gamma+\lambda+\delta+\kappa=1\f]
 *
 * The stationnary distribution
 * \f[
 * \pi = \left(\pi_A, \pi_C, \pi_G, \pi_T\right)
 * \f]
 * can be computed analytically, so we define parameters for it (values between 0 and 1):
 *\f[
 * \begin{cases}
 * \theta_R = \pi_A + \pi_G\\
 * \theta_C = \frac{\pi_C}{1 - \theta_R} = \frac{\pi_C}{\pi_C + \pi_G}\\
 * \theta_G = \frac{\pi_G}{\theta_R} = \frac{\pi_G}{\pi_A + \pi_G}\\
 * \end{cases}
 * \f]
 *
 * with parameters with values in [0;1[:
 *
 *\f[
 * \begin{cases}
 * \kappa'=\frac{\kappa}{\theta_R}\\
 * \gamma'=\frac{\gamma}{1-\theta_R}\\
 * \end{cases}
 * \f]
 *
 * and parameters with values > 1:
 *
 *\f[
 * \begin{cases}
 * \alpha'=\frac{\alpha(1-\theta_G)+min(\theta_G,\kappa')(1-\theta_R)}{\theta_G(1-\theta_R)}\\
 * \sigma'=\frac{\sigma(1-\theta_C)+min(\theta_C,\gamma')\theta_R}{\theta_C\theta_R}\\
 * \end{cases}
 * \f]
 *
 * The generator is then computed as:
 *
 *\f[
 * \begin{cases}
 * \kappa=\kappa' \theta_R\\
 * \gamma=\gamma' (1-\theta_R)\\
 * \delta=\theta_R - \kappa\\
 * \lambda=1-\theta_R-\gamma\\
 * \alpha=\frac{\alpha'(1-\theta_R)\theta_G-min(\theta_G,\kappa')(1-\theta_R)}{1-\theta_G}\\
 * \sigma=\frac{\sigma'\theta_R\theta_C-min(\theta_C,\gamma')\theta_R}{1-\theta_C}\\
 * \beta=\frac{\gamma'*\theta_R+\sigma}{\theta_C}-\sigma-\theta_R\\
 * \epsilon=\frac{\kappa'*(1-\theta_R)+\alpha}{\theta_G}-\alpha-(1-\theta_R)\\
 * \end{cases}
 * \f]
 *
 * and \f[P\f] is set for normalization.
 *
 * The parameters are named \c "thetaR", \c "thetaC", \c "thetaG",
 * \c "kappaP", \c "gammaP", \c "sigmaP", \c "alphaP".
 *
 * References:
 * - Rhetsky A. \& Nei M. (1995) MBE 12(1) 131-151.
 * - Schadt, Sinsheimer \& Lange (1998) Genome Research 8 222-233.
 */

class RN95 :
  public virtual NucleotideSubstitutionModel,
  public AbstractSubstitutionModel
{
private:
  double alpha_, beta_, gamma_, delta_, epsilon_, kappa_, lambda_, sigma_;
  double r_;
  /**
   * For calculation purposes as in Schadt & al. (with c1_=1)
   */
  double c1_, c2_, c3_, c4_, c5_, c6_, c7_, c8_, c9_;
  mutable RowMatrix<double> p_;
  mutable double exp1_, exp3_, exp6_, l_;

public:
  RN95(
    const NucleicAlphabet* alphabet,
    double alpha = 1,
    double beta = 1,
    double gamma = 1,
    double delta = 1,
    double epsilon = 1,
    double kappa = 1,
    double lambda = 1,
    double sigma = 1);

  virtual ~RN95() {}

#ifndef NO_VIRTUAL_COV
  RN95*
#else
  Clonable*
#endif
  clone() const { return new RN95(*this); }

public:
  double Pij_t    (int i, int j, double d) const;
  double dPij_dt  (int i, int j, double d) const;
  double d2Pij_dt2(int i, int j, double d) const;
  const Matrix<double>& getPij_t    (double d) const;
  const Matrix<double>& getdPij_dt  (double d) const;
  const Matrix<double>& getd2Pij_dt2(double d) const;
  std::string getName() const { return "RN95"; }

  void updateMatrices();

  void setFreq(std::map<int, double>&);
};
} // end of namespace bpp.

#endif  // _RN95_H_

