//
// File: SSR.h
// Created by: Julien Dutheil
// Created on: Tue Nov 4 11:46 2008
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

#ifndef _SSR_H_
#define _SSR_H_

#include "NucleotideSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"

#include <Bpp/Numeric/Constraints.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{

/**
 * @brief The Strand Symmetric Reversible substitution model for
 * nucleotides.
 *
 * We use a parametrization derived from Hobolth et al 2007
 * \f[
 * S = \begin{pmatrix}
 * \cdots & \beta & 1 & \gamma \\ 
 * \beta & \cdots & \delta & 1 \\ 
 * 1 & \delta & \cdots & \beta \\ 
 * \gamma & 1 & \beta & \cdots \\ 
 * \end{pmatrix}
 * \f]
 * The equilibrium frequencies 
 * \f[
 * \pi = \left(1-\frac{\theta}{2}, \frac{\theta}{2}, \frac{\theta}{2}, 1-\frac{\theta}{2}\right)
 * \f]
 * This models hence includes four parameters, three relative rates \f$\beta, \gamma, \delta\f$ and the GC content \f$\theta\f$.
 *
 * Normalization: we set \f$f\f$ to 1, and scale the matrix so that \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -\gamma\pi_T-\pi_G-\beta\pi_C & \beta\pi_C & \pi_G & \gamma\pi_T \\ 
 * \beta\pi_A & -\pi_T-\delta\pi_G-\beta\pi_A & \delta\pi_G & \pi_T \\ 
 * \pi_A & \delta\pi_C & -\beta\pi_T-\delta\pi_C-\pi_A & \beta\pi_T \\ 
 * \gamma\pi_A & \pi_C & \beta\pi_G & -\beta\pi_G-\pi_C-\gamma\pi_A \\ 
 * \end{pmatrix}
 * \f]
 * where P is the normalization constant.
 * For now, the generator of this model is diagonalized numericaly.
 * See AbstractSubstitutionModel for details of how the probabilities are computed.
 *
 * The parameters are named \c "beta", \c "gamma", \c "delta", and \c "theta"
 * and their values may be retrieved with the command 
 * \code
 * getParameterValue("beta")
 * \endcode
 * for instance.
 * 
 * Reference:
 * - Hobolth A, Christensen O Fm Mailund T, Schierup M H (2007), PLoS_ Genetics_ 3(2) e7.
 * - Yap VB, Speed TP (1995), Journal_ Of Molecular Evolution_ 58(1) 12-18
 */
  class SSR:
    public virtual NucleotideSubstitutionModel,
    public AbstractReversibleSubstitutionModel
  {
  private:
    double beta_, gamma_, delta_, theta_, piA_, piC_, piG_, piT_;
  
  public:
    SSR( const NucleicAlphabet* alpha,
         double beta = 1.,
         double gamma = 1.,
         double delta = 1.,
         double theta = 0.5);
  
    virtual ~SSR() {}
  
#ifndef NO_VIRTUAL_COV
    SSR*
#else
    Clonable*
#endif
    clone() const { return new SSR(*this); }
  
  public:
    std::string getName() const { return "Strand Symmetric Reversible"; }
  
    void updateMatrices();
  
    /**
     * @brief This method is redefined to actualize the corresponding parameters theta too.
     */
    void setFreq(std::map<int, double>&);
  };

} //end of namespace bpp.

#endif	//_SSR_H_

