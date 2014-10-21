//
// File ContingencyTableGenerator.h
// Author: Julien Dutheil
// Created on: Fri Dec 10 2010 16:19
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

#ifndef _CONTINGENCYTABLEGENERATOR_H_
#define _CONTINGENCYTABLEGENERATOR_H_

#include "RandomFactory.h"
#include "RandomTools.h"
#include "../Matrix/Matrix.h"

// From the STL:
#include <cmath>
#include <vector>

namespace bpp
{

  /**
 * @brief Generate a random contingency matrix with given marginal counts.
 *
 * This procedure was adapted from the original fortran code described in:
 * Patefield, W. M. (1981) Algorithm AS159.  An efficient method of
 * generating r x c tables with given row and column totals.
 * _Applied Statistics_ *30*, 91-97.
 * This algorithm is the one also used in R function chisq.test for instance.
 *
 * The code was adpated from file rcont.c, edited by Martin Maechler, Dec 2003,
 * available in the R software source distribution.
 *
 * @param nrowt Marginal counts.
 * @param ncolt Marginal counts.
 * @return A random matrix of counts with the same marginals as specified.
 */
class ContingencyTableGenerator
{
  private:
    std::vector<size_t> nrowt_;
    std::vector<size_t> ncolt_;
    size_t nrow_;
    size_t ncol_;
    size_t nrowm_;
    size_t ncolm_;
    std::vector<size_t> jwork_; //workspace
    size_t ntot_; //total number of observations
    std::vector<double> fact_; //log factorial

  public:
    ContingencyTableGenerator(const std::vector<size_t>& nrowt, const std::vector<size_t>& ncolt);

  public:
    RowMatrix<size_t> rcont2(const RandomFactory& generator = *RandomTools::DEFAULT_GENERATOR); 
};

} //end of namespace bpp.

#endif  //_CONTINGENCYTABLEGENERATOR_H_

