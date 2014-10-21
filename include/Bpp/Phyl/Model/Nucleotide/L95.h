//
// File: L95.h
// Created by: Laurent Guéguen
// Created on: lundi 18 octobre 2010, à 22h 21
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

#ifndef _L95_H_
#define _L95_H_

#include "NucleotideSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"

#include <Bpp/Numeric/Constraints.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{

/**
 * @brief The no-strand bias substitution model for nucleotides, from
 * Lobry 1995. The point of this model is that the substitution rate
 * from a nucleotide N towards another M is the same as the rate from
 * the complement of N towards the complement of M. Note that this
 * model is not reversible.
 *
 * After normalization, this model contains 5 parameters:
 * \f[
 * Q = \frac 1{2*\kappa*\theta*(1-\theta)+\gamma+\theta-2*\theta*\gamma} \begin{pmatrix}
 * \cdots & \kappa.\beta.\theta & \kappa.(1-\beta).\theta & \gamma \\ 
 * \kappa.\alpha.(1-\theta) & \cdots & 1-\gamma & \kappa.(1-\alpha).(1-\theta) \\ 
 * \kappa.(1-\alpha).(1-\theta) & 1-\gamma & \cdots & \kappa.\alpha.(1-\theta) \\ 
 * \gamma & \kappa.(1-\beta).\theta & \kappa.\beta.\theta & \cdots \\ 
 * \end{pmatrix}
 * \f]
 * The equilibrium frequencies are
 * \f[
 * \pi = \left(\frac{1-\theta}{2}, \frac{\theta}{2}, \frac{\theta}{2}, \frac{1-\theta}{2}\right)
 * \f]
 *
 * and then \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 *
 * The generator of this model is diagonalized numerically.
 * See AbstractSubstitutionModel for details of how the probabilities are computed.
 *
 * The parameters are named \c "alpha", \c "beta", \c "gamma", \c
 * "kappa" and \c "theta". The values of \c "alpha", \c "beta", \c
 * "gamma" are between 0 and 1, The values of \c "gamma" are between 0
 * and 1 excluded, the values of "kappa" which are positive. Their
 * values may be retrieved with the command:
 *
 * \code
 * getParameterValue("alpha")
 * \endcode for instance.
 * 
 * Reference:
 * - Lobry J R (1995), Journal_ Of Molecular Evolution_ 40 326-330.
 */
class L95:
  public virtual NucleotideSubstitutionModel,
  public AbstractSubstitutionModel
{
private:
  double alpha_, beta_, gamma_, kappa_, theta_;
  
public:
  L95(
      const NucleicAlphabet* alphabet,
      double alpha = 0.5,
      double beta = 0.5,
      double gamma = 0.5,
      double kappa = 1.,
      double theta = 0.5);
  
  virtual ~L95() {}
  
#ifndef NO_VIRTUAL_COV
  L95*
#else
  Clonable*
#endif
  clone() const { return new L95(*this); }
  
public:
  std::string getName() const { return "L95"; }
  
  void updateMatrices();
  
  /**
   * @brief This method is redefined to actualize the corresponding parameters theta too.
   */
  void setFreq(std::map<int, double>&);
};

} //end of namespace bpp.

#endif	//_L95_H_

