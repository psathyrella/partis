//
// File: K80.h
// Created by: Julien Dutheil
// Created on: Tue May 27 15:24:30 2003
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

#ifndef _K80_H_
#define _K80_H_


#include "NucleotideSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{

  /**
   * @brief The Kimura 2-rates substitution model for nucleotides.
   *
   * All two rates: one for transitions and one for transversions.
   * This models include one parameter, the transition / transversion
   * relative rate \f$\kappa\f$.
   * \f[
   * S = \begin{pmatrix}
   * \cdots & r & \kappa r & r \\ 
   * r & \cdots & r & \kappa r \\ 
   * \kappa r & r & \cdots & r \\ 
   * r & \kappa r & r & \cdots \\ 
   * \end{pmatrix}
   * \f]
   * \f[
   * \pi = \left(\frac{1}{4}, \frac{1}{4}, \frac{1}{4}, \frac{1}{4}\right)
   * \f]
   * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
   * \f[
   * S = \begin{pmatrix}
   * -4 & \frac{4}{\kappa+2} & \frac{4\kappa}{\kappa+2} & \frac{4}{\kappa+2} \\ 
   * \frac{4}{\kappa+2} & -4 & \frac{4}{\kappa+2} & \frac{4\kappa}{\kappa+2} \\ 
   * \frac{4\kappa}{\kappa+2} & \frac{4}{\kappa+2} & -4 & \frac{4}{\kappa+2} \\ 
   * \frac{4}{\kappa+2} & \frac{4\kappa}{\kappa+2} & \frac{4}{\kappa+2} & -4 \\ 
   * \end{pmatrix}
   * \f]
   * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
   * \f[
   * Q = S . \pi = \begin{pmatrix}
   * -1 & \frac{1}{\kappa+2} & \frac{\kappa}{\kappa+2} & \frac{1}{\kappa+2} \\ 
   * \frac{1}{\kappa+2} & -1 & \frac{1}{\kappa+2} & \frac{\kappa}{\kappa+2} \\ 
   * \frac{\kappa}{\kappa+2} & \frac{1}{\kappa+2} & -1 & \frac{1}{\kappa+2} \\ 
   * \frac{1}{\kappa+2} & \frac{\kappa}{\kappa+2} & \frac{1}{\kappa+2} & -1 \\ 
   * \end{pmatrix}
   * \f]
   *  
   * The eigen values are \f$\left(0, -\frac{2\kappa+2}{\kappa+2}, -\frac{2\kappa+2}{\kappa+2}, -\frac{4}{\kappa+2}\right)\f$, 
   * the left eigen vectors are, by row:
   * \f[
   * U = \begin{pmatrix}
   * \frac{1}{4} &  \frac{1}{4} &  \frac{1}{4} &  \frac{1}{4} \\
   *           0 &  \frac{1}{2} &            0 & -\frac{1}{2} \\
   * \frac{1}{2} &            0 & -\frac{1}{2} &            0 \\
   * \frac{1}{4} & -\frac{1}{4} &  \frac{1}{4} & -\frac{1}{4} \\
   * \end{pmatrix}
   * \f]
   * and the right eigen vectors are, by column:
   * \f[
   * U^-1 = \begin{pmatrix}
   * 1 &  0 &  1 &  1 \\
   * 1 &  1 &  0 & -1 \\
   * 1 &  0 & -1 &  1 \\
   * 1 & -1 &  0 & -1 \\
   * \end{pmatrix}
   * \f]
   *
   * In addition, a rate_ factor defines the mean rate of the model.
   *
   * The probabilities of changes are computed analytically using the formulas:
   * \f[
   * P_{i,j}(t) = \begin{pmatrix}
   * \frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B & -\frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B \\
   * \frac{1}{4} - \frac{1}{4}B & \frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B & -\frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} \\
   * -\frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B & \frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B \\
   * \frac{1}{4} - \frac{1}{4}B & -\frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B & \frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} \\
   * \end{pmatrix}
   * \f]
   * with \f$A=e^{-\frac{rate\_ * (2\kappa+2)t}{\kappa+2}}\f$ and \f$B = e^{-\frac{rate\_ * 4t}{\kappa+2}}\f$. 
   *
   * First and second order derivatives are also computed analytically using the formulas:
   * \f[
   * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \begin{pmatrix}
   * -\frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B & \frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B \\
   * \frac{1}{\kappa+2}B & -\frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B & \frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B \\
   * \frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B & -\frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B \\
   * \frac{1}{\kappa+2}B & \frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B & -\frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B \\
   * \end{pmatrix}
   * \f]
   * \f{multline*}
   * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} = rate\_^2 * \\
   * \begin{pmatrix}
   * \frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A - \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B & -\frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B \\
   * -\frac{4}{{(\kappa+2)}^2}B & \frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B & -\frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B \\
   * -\frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B & \frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B \\
   * -\frac{4}{{(\kappa+2)}^2}B & -\frac{2{(\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B & \frac{2{(\kappa+2)}^2}{{(2\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B \\
   * \end{pmatrix}
   * \f}
   *
   * The parameter is named \c "kappa"
   * and its value may be retrieve with the command 
   * \code
   * getParameterValue("kappa")
   * \endcode
   * 
   * Reference:
   * - Kimura M (1980), Journal_ Of Molecular Evolution_ 16(2) 111-20. 
   */
  class K80:
    public virtual NucleotideSubstitutionModel,
    public AbstractReversibleSubstitutionModel
  {
  private:
    double kappa_, r_;
    mutable double l_, k_, exp1_, exp2_;
    mutable RowMatrix<double> p_;

  public:
    K80(const NucleicAlphabet* alpha, double kappa = 1.);

    virtual ~K80() {}

#ifndef NO_VIRTUAL_COV
    K80*
#else
    Clonable*
#endif
    clone() const { return new K80(*this); }

  public:
    double Pij_t    (int i, int j, double d) const;
    double dPij_dt  (int i, int j, double d) const;
    double d2Pij_dt2(int i, int j, double d) const;
    const Matrix<double> & getPij_t    (double d) const;
    const Matrix<double> & getdPij_dt  (double d) const;
    const Matrix<double> & getd2Pij_dt2(double d) const;

    std::string getName() const { return "K80"; }
	   
    /**
     * @brief This method is disabled in this model since frequencies are not free parameters.
     *
     * Consider using the HKY85 model for instance if you want to set frequencies as parameters.
     *
     * @param data Useless parameter.
     */
    void setFreqFromData(const SequenceContainer & data) {}
	
  protected:
    void updateMatrices();

  };

} //end of namespace bpp.

#endif	//_K80_H_

