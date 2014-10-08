//
// File: Simplex.h
// Created by: Laurent Guéguen
// Created on: mardi 31 mai 2011, à 11h 02
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for numerical calculus.

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

#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_


// From the STL:
#include <vector>

#include "../AbstractParameterAliasable.h"

namespace bpp
{
/**
 * @brief A Simplex object, used to define sets of probabilities that
 * sum 1.
 *
 * The probabilities are denoted @f$p_i@f$ for @f$i \in <1,n>@f$.
 *
 * If they are parametrized, the parameters are called \c "theta1",
 * ..., \c "theta(n-1)".
 *
 * Two ways of parametrization are available:
 *
 * Global ratio:
 *
 * @f$\forall i<n, \theta_i=\frac{p_i}{1-(p_1+...+p_{i-1})}@f$.
 *
 * In the reverse,
 * @f$\forall i<n, p_i= (1-\theta_1).(1-\theta_2)...\theta_{i}@f$
 * and @f$p_n=(1-\theta_1).(1-\theta_2)...(1-\theta_{n-1})@f$.
 *
 *
 * Local ratio:
 *
 * @f$\theta_i = \frac{p_i}{p_i+p_{i+1}} \forall i \in 1..\textrm{n-1}@f$.
 *
 * In the reverse if we denote @f$\alpha_i=\frac{1-\theta_i}{\theta_i}@f$,
 * @f$p_i=\frac{\alpha_1...\alpha_{i-1}}{1+\sum_{k=1}^{n-1}\alpha_1...\alpha_k}@f$.
 *
 * Binary:
 *
 * This parametrization is based on the binary coding.
 *
 * Given @f$a_b...a_1@f$ the writing of i in binary, we denote
 * @f$i_k=a_k...a_1@f$.
 *
 * Given @f$a_b...a_1@f$ the writing of i in binary where @f$a_b=1@f$,
 * we denote @f$1_i=\sum\{p_{j+1} \text{ such that } j_b=i_b=1i_{b-1}\}@f$ and
 * @f$0_i=\sum\{p_{j+1} \text{ such that }  j_b=0i_{b-1}\}@f$, and then we define:
 *
 *
 * @f$\theta_i=\frac{1_i}{1_i+0_i}@f$
 *
 * and on the reverse, we denote @f$\theta'_{0i_{b-1}}=1-\theta_i@f$
 * and @f$\theta'_{1i_{b-1}}=\theta_i@f$.
 *
 * Then, if @f$c=ceil(log_2(n))@f$, for @f$i \in <0,n-1>@f$.
 *
 * @f$p_{i+1}=\theta'_{i_c}....\theta'_{i_1} @f$
 *
 */
class Simplex :
  public AbstractParameterAliasable
{
private:

  /**
   * @brief The dimension+1 of the space simplex (ie the number of probabilities).
   *
   */
  size_t dim_;

  /**
   * @brief the method of parametrization.
   *
   * 0: No parametrization
   * 1: Global ratio
   * 2: Local ratio
   * 3: Binary
   */
  unsigned short method_;

  std::vector<double> vProb_;

  /**
   * @brief just used with local ratio (method 2)
   */
  std::vector<double> valpha_;

public:
  /**
   * @brief Builds a new Simplex object from a number of
   * probabilities. They are initialized equal.
   *
   * @param dim The number of probabilities.
   * @param method  tells the method of parametrization (default 0)
   *    0: No parametrization
   *    1: Global ratio
   *    2: Local ratio
   *    3: Binary
   * @param allowNull if null probabilites are allowed (default: false)
   * @param name The name passed to AbstractParameterAliasable constructor.
   *
   */
  Simplex(size_t dim, unsigned short method = 0, bool allowNull = false, const std::string& name = "Simplex.");

  /**
   * @brief Builds a new Simplex object from a vector of probabilities
   *
   * @param probas The vector of probabilities.
   * @param method  tells the method of parametrization (default 0)
   *    0: No parametrization
   *    1: Global ratio
   *    2: Local ratio
   *    3: Binary
   * @param allowNull if null probabilites are allowed (default: false)
   * @param name The name passed to AbstractParameterAliasable constructor.
   *
   */
  Simplex(const std::vector<double>& probas, unsigned short method = 0, bool allowNull = false, const std::string& name = "Simplex.");

  virtual ~Simplex() {}

  Simplex* clone() const { return new Simplex(*this); }

public:
  void fireParameterChanged(const ParameterList& parameters);

  size_t dimension() const { return dim_; }

  void setFrequencies(const std::vector<double>&);

  const std::vector<double>& getFrequencies() const { return vProb_;}

  std::vector<double>& getFrequencies() { return vProb_;}

  double prob(size_t i) const { return vProb_[i]; }

  unsigned short getMethod() const { return method_; }
};
} // end of namespace bpp.

#endif  // _SIMPLEDISCRETEDISTRIBUTION_H_

