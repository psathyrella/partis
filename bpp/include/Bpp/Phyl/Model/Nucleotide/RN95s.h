//
// File: RN95s.h
// Created by: Laurent Guéguen
// Created on: samedi 12 mars 2011, à 06h 49
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

#ifndef _RN95s_H_
#define _RN95s_H_

#include "NucleotideSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"

#include <Bpp/Numeric/Constraints.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief Intersection of models RN95 and L95.
 *
 * The two hypotheses are that the transversion rates are only
 * dependent of the target nucleotide, and strand symmetry.
 *
 * After normalization this model has 3 parameters:
 * \f[
 * Q= \frac 1P
 * \begin{pmatrix}
 * . & \gamma & \alpha & \delta \\
 * \delta & . &  \gamma & \beta \\
 * \beta & \gamma & . & \delta \\
 * \delta & \alpha & \gamma & .\\
 * \end{pmatrix}\f]
 *
 * so in the parametrization process we set: \f[\gamma+\delta=\frac 12\f]
 *
 * The stationnary distribution
 * \f[
 * \pi = \left(\pi_A, \pi_C, \pi_G, \pi_T\right)
 * \f]
 * is such as \f[\pi_A=\pi_T\f] and \f[\pi_C=\pi_G\f]
 * It can be computed analytically.
 *
 * We use as  parameters:
 *
 *\f[
 * \begin{cases}
 * \theta_A (=\pi_A) \in ]0;1/2[\\
 * \gamma \in ]0;1/2[\\
 * \alpha' > 1
 * \end{cases}
 * \f]
 *
 * where \f[\alpha'=\frac{2\alpha\pi_A+min(0.5-\pi_A,\gamma)}{0.5-\pi_A}\f].
 *
 * The generator is then computed as:
 *
 *\f[
 * \begin{cases}
 * \delta=\frac 12-\gamma\\
 * \alpha=\frac{\alpha'(0.5-\pi_A)-min(0.5-\pi_A,\gamma)}{2\pi_A}\\
 * \beta=\frac{2*\pi_A*(\alpha+\frac 12)-\delta}{1-2*\pi_A}\\
 * \end{cases}
 * \f]
 *
 * and @f$P@f$ is set for normalization.
 *
 * The parameters are named \c "thetaA", \c "gamma", \c "alphaP".
 *
 * References:
 * - Rhetsky A. \& Ney M. (1995) MBE 12(1) 131-151.
 * - Lobry J R (1995), Journal_ Of Molecular Evolution_ 40 326-330.
 * - Schadt, Sinsheimer \& Lange (1998) Genome Research 8 222-233.
 */

class RN95s :
  public virtual NucleotideSubstitutionModel,
  public AbstractSubstitutionModel
{
private:
  double alpha_, beta_, gamma_, delta_;
  double r_;
  /**
   * For calculation purposes as in Schadt & al. (with c1_=1)
   */
  double c3_, c4_, c8_;
  mutable RowMatrix<double> p_;
  mutable double exp1_, exp3_, l_;

public:
  RN95s(const NucleicAlphabet* alphabet,
        double alpha = 1,
        double beta = 1,
        double gamma = 1,
        double delta = 1);

  virtual ~RN95s() {}

#ifndef NO_VIRTUAL_COV
  RN95s*
#else
  Clonable*
#endif
  clone() const { return new RN95s(*this); }

public:
  double Pij_t    (int i, int j, double d) const;
  double dPij_dt  (int i, int j, double d) const;
  double d2Pij_dt2(int i, int j, double d) const;
  const Matrix<double>& getPij_t    (double d) const;
  const Matrix<double>& getdPij_dt  (double d) const;
  const Matrix<double>& getd2Pij_dt2(double d) const;

  std::string getName() const { return "RN95s"; }

  void updateMatrices();

  /**
   * @brief This method takes the average value between observed @f$\pi_A@f$ and @f$\pi_T@f$.
   */

  void setFreq(std::map<int, double>&);
};
} // end of namespace bpp.

#endif  // _RN95s_H_

