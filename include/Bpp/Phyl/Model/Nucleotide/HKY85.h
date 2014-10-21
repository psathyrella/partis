//
// File: HKY85.h
// Created by: Julien Dutheil
// Created on: Thu Jan 22 16:17:39 2004
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

#ifndef _HKY85_H_
#define _HKY85_H_

#include "NucleotideSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"

#include <Bpp/Numeric/Constraints.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{

/**
 * @brief The Hasegawa M, Kishino H and Yano T (1985) substitution model for nucleotides.
 *
 * This model is similar to K80 model,
 * but allows distinct equilibrium frequencies between A, C, G and T.
 * \f[
 * \begin{pmatrix}
 * \cdots & r & \kappa r & r \\ 
 * r & \cdots & r & \kappa r \\ 
 * \kappa r & r & \cdots & r \\ 
 * r & \kappa r & r & \cdots \\ 
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\pi_A, \pi_C, \pi_G, \pi_T\right)
 * \f]
 * This models hence includes five parameters, the transition / transversion
 * relative rate \f$\kappa\f$ and four frequencies \f$\pi_A, \pi_C, \pi_G, \pi_T\f$.
 * These four frequencies are not independent parameters, since they have the constraint to
 * sum to 1.
 * We use instead a different parametrization to remove this constraint:
 * \f[
 * \begin{cases}
 * \theta = \pi_C + \pi_G\\
 * \theta_1 = \frac{\pi_A}{1 - \theta} = \frac{\pi_A}{\pi_A + \pi_T}\\
 * \theta_2 = \frac{\pi_G}{\theta} = \frac{\pi_G}{\pi_C + \pi_G}\\
 * \end{cases}
 * \Longleftrightarrow
 * \begin{cases}
 * \pi_A = \theta_1 (1 - \theta)\\
 * \pi_C = (1 - \theta_2) \theta\\
 * \pi_G = \theta_2 \theta\\
 * \pi_T = (1 - \theta_1)(1 - \theta).
 * \end{cases}
 * \f]
 * These parameters can also be measured from the data and not optimized.
 *
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \frac{1}{P}\begin{pmatrix}
 * \frac{-\pi_T-\kappa\pi_G-\pi_C}{\pi_A} & 1 & \kappa & 1 \\ 
 * 1 & \frac{-\kappa\pi_T-\pi_G-\pi_A}{\pi_C} & 1 & \kappa \\ 
 * \kappa & 1 & \frac{-\pi_T-\pi_C-\kappa\pi_A}{\pi_G} & 1 \\ 
 * 1 & \kappa & 1 & \frac{-\pi_G-\kappa\pi_C-\pi_A}{\pi_T} \\ 
 * \end{pmatrix}
 * \f]
 * with \f$P=2\left(\pi_A \pi_C + \pi_C \pi_G + \pi_A \pi_T + \pi_G \pi_T + \kappa \left(\pi_C \pi_T + \pi_A \pi_G\right)\right)\f$.
 *
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -\pi_T-\kappa\pi_G-\pi_C & \pi_C & \kappa\pi_G & \pi_T \\ 
 * \pi_A & -\kappa\pi_T-\pi_G-\pi_A & \pi_G & \kappa\pi_T \\ 
 * \kappa\pi_A & \pi_C & -\pi_T-\pi_C-\kappa\pi_A & \pi_T \\ 
 * \pi_A & \kappa\pi_C & \pi_G & -\pi_G-\kappa\pi_C-\pi_A \\ 
 * \end{pmatrix}
 * \f]
 *  
 * The eigen values are \f$\left(0, -\frac{\pi_R + \kappa\pi_Y}{P}, -\frac{\pi_Y + \kappa\pi_R}{P}, -\frac{1}{P}\right)\f$,
 * with \f$\pi_R=\pi_A+\pi_G\f$ and \f$\pi_Y=\pi_C+\pi_T\f$.
 * The left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 *                    \pi_A &               \pi_C &                    \pi_G & \pi_T \\
 *                        0 & \frac{\pi_T}{\pi_Y} &                        0 & -\frac{\pi_T}{\pi_Y} \\
 *      \frac{\pi_G}{\pi_R} &                   0 &     -\frac{\pi_G}{\pi_R} & 0 \\
 * \frac{\pi_A\pi_Y}{\pi_R} &              -\pi_C & \frac{\pi_G\pi_Y}{\pi_R} & -\pi_T \\
 * \end{pmatrix}
 * \f]
 * and the right eigen vectors are, by column:
 * \f[
 * U^-1 = \begin{pmatrix}
 * 1 &  0 &  1 &  1 \\
 * 1 &  1 &  0 & -\frac{\pi_R}{\pi_Y} \\
 * 1 &  0 & \frac{\pi_A}{\pi_G} &  1 \\
 * 1 & -\frac{\pi_C}{\pi_T} &  0 & -\frac{\pi_R}{\pi_Y} \\
 * \end{pmatrix}
 * \f]
 *
 * In addition, a rate_ factor defines the mean rate of the model.
 *
 * The probabilities of changes are computed analytically using the formulas:
 * \f{multline*}
 * P_{i,j}(t) = \\
 * \begin{pmatrix}
 * \frac{\pi_G}{\pi_R}A + \frac{\pi_A\pi_Y}{\pi_R}B + \pi_A & \pi_C - \pi_CB & -\frac{\pi_G}{\pi_R}A + \frac{\pi_G\pi_Y}{\pi_R}B + \pi_G & \pi_T - \pi_TB \\
 * \pi_A - \pi_AB & \frac{\pi_T}{\pi_Y}A + \frac{\pi_C\pi_R}{\pi_Y}B + \pi_C & \pi_G - \pi_GB & -\frac{\pi_T}{\pi_Y}A + \frac{\pi_T\pi_R}{\pi_Y}B + \pi_T \\
 * -\frac{\pi_A}{\pi_R}A + \frac{\pi_A\pi_Y}{\pi_R}B + \pi_A & \pi_C - \pi_CB & \frac{\pi_A}{\pi_R}A + \frac{\pi_G\pi_Y}{\pi_R}B + \pi_G & \pi_T - \pi_TB \\
 * \pi_A - \pi_AB & -\frac{\pi_C}{\pi_Y}A + \frac{\pi_C\pi_R}{\pi_Y}B + \pi_C & \pi_G - \pi_GB & \frac{\pi_C}{\pi_Y}A + \frac{\pi_R\pi_T}{\pi_Y}B + \pi_T \\
 * \end{pmatrix}
 * \f}
 * with \f$A=e^{-\frac{rate\_*(\pi_Y+\kappa\pi_R)t}{P}}\f$ and \f$B = e^{-\frac{rate\_*t}{P}}\f$. 
 *
 * First and second order derivatives are also computed analytically using the formulas:
 * \f{multline*}
 * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \\
 * \frac{1}{P}
 * \begin{pmatrix}
 * -\frac{\pi_G(\pi_Y+\kappa\pi_R)}{\pi_R}A - \frac{\pi_A\pi_Y}{\pi_R}B & \pi_CB & \frac{\pi_G(\pi_Y+\kappa\pi_R)}{\pi_R}A - \frac{\pi_G\pi_Y}{\pi_R}B & \pi_TB \\
 * \pi_AB & -\frac{\pi_T(\pi_R+\kappa\pi_Y)}{\pi_Y}A - \frac{\pi_C\pi_R}{\pi_Y}B & \pi_GB & \frac{\pi_T(\pi_R+\kappa\pi_Y)}{\pi_Y}A - \frac{\pi_T\pi_R}{\pi_Y}B \\
 * \frac{\pi_A(\pi_Y+\kappa\pi_R)}{\pi_R}A - \frac{\pi_A\pi_Y}{\pi_R}B & \pi_CB & -\frac{\pi_A(\pi_Y+\kappa\pi_R)}{\pi_R}A - \frac{\pi_G\pi_Y}{\pi_R}B & \pi_TB \\
 * \pi_AB & \frac{\pi_C(\pi_R+\kappa\pi_Y)}{\pi_Y}A - \frac{\pi_C\pi_R}{\pi_Y}B & \pi_GB & -\frac{\pi_C(\pi_R+\kappa\pi_Y)}{\pi_Y}A - \frac{\pi_R\pi_T}{\pi_Y}B \\
 * \end{pmatrix}
 * \f}
 * \f{multline*}
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} = rate\_^2 * \\
 * \frac{1}{P^2}
 * \begin{pmatrix}
 * \frac{\pi_G{(\pi_Y+\kappa\pi_R)}^2}{\pi_R}A + \frac{\pi_A\pi_Y}{\pi_R}B & -\pi_CB & -\frac{\pi_G{(\pi_Y+\kappa\pi_R)}^2}{\pi_R}A + \frac{\pi_G\pi_Y}{\pi_R}B & -\pi_TB \\
 * -\pi_AB & \frac{\pi_T{(\pi_R+\kappa\pi_Y)}^2}{\pi_Y}A + \frac{\pi_C\pi_R}{\pi_Y}B & -\pi_GB & -\frac{\pi_T{(\pi_R+\kappa\pi_Y)}^2}{\pi_Y}A + \frac{\pi_T\pi_R}{\pi_Y}B \\
 * -\frac{\pi_A{(\pi_Y+\kappa\pi_R)}^2}{\pi_R}A + \frac{\pi_A\pi_Y}{\pi_R}B & -\pi_CB & \frac{\pi_A{(\pi_Y+\kappa\pi_R)}^2}{\pi_R}A + \frac{\pi_G\pi_Y}{\pi_R}B & -\pi_TB \\
 * -\pi_AB & -\frac{\pi_C{(\pi_R+\kappa\pi_Y)}^2}{\pi_Y}A + \frac{\pi_C\pi_R}{\pi_Y}B & -\pi_GB & \frac{\pi_C{(\pi_R+\kappa\pi_Y)}^2}{\pi_Y}A + \frac{\pi_R\pi_T}{\pi_Y}B \\
  * \end{pmatrix}
 * \f}
 *
 * The parameters are named \c "kappa", \c "theta", \c "theta1", and \c "theta2"
 * and their values may be retrieve with the command 
 * \code
 * getParameterValue("kappa")
 * \endcode
 * for instance.
 *
 * Reference:
 * - Hasegawa M, Kishino H and Yano T (1985), Molecular_ Biology And Evolution_ 22(2) 160-74. 
 */
class HKY85:
  public virtual NucleotideSubstitutionModel,
  public AbstractReversibleSubstitutionModel
{
	private:
		double kappa_, k1_, k2_, r_, piA_, piC_, piG_, piT_, piY_, piR_, theta_, theta1_, theta2_;
    mutable double exp1_, exp21_, exp22_, l_;
    mutable RowMatrix<double> p_;

	public:
		HKY85(
			const NucleicAlphabet * alpha,
			double kappa = 1.,
			double piA = 0.25,
			double piC = 0.25,
			double piG = 0.25,
			double piT = 0.25);
	
		virtual ~HKY85() {}

#ifndef NO_VIRTUAL_COV
    HKY85*
#else
    Clonable*
#endif
    clone() const { return new HKY85(*this); }

  public:
		double Pij_t    (int i, int j, double d) const;
		double dPij_dt  (int i, int j, double d) const;
		double d2Pij_dt2(int i, int j, double d) const;
		const Matrix<double> & getPij_t    (double d) const;
		const Matrix<double> & getdPij_dt  (double d) const;
		const Matrix<double> & getd2Pij_dt2(double d) const;

    std::string getName() const { return "HKY85"; }

  /**
   * @brief This method is redefined to actualize the corresponding parameters piA, piT, piG and piC too.
   */
  void setFreq(std::map<int, double>& freqs);
	
  void updateMatrices();
};

} //end of namespace bpp.

#endif	//_HKY85_H_

