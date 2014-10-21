//
// File: PrincipalComponentAnalysis.h
// Created by: Mathieu Groussin
// Created on: Thu Mar 03 10:13 2011
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide basal and
   utilitary classes. This file belongs to the Bio++ Project.

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


#ifndef _PRINCIPALCOMPONENTANALYSIS_H_
#define _PRINCIPALCOMPONENTANALYSIS_H_

#include "../../Matrix/Matrix.h"
#include "DualityDiagram.h"

namespace bpp
{
/**
 * @brief This class allows to perform a principal component analysis.
 *
 * Two constructors are available. The first one allows the user to specify the row and column weights. The second one specify default weights:
 * uniform weights unit weights are created for rows and columns respectively.
 *
 * The code of this class is deeply inspired from the R code of the dudi.pca function available in the ade4 package.
 */
class PrincipalComponentAnalysis:
  public DualityDiagram
{
private:
  std::vector<double> columnMeans_;
  std::vector<double> columnSd_;

public:
  /**
   * @brief Build a new PrincipalComponentAnalysis object.
   *
   * @param data The input data (a RowMatrix) to analyse.
   * @param nbAxes The number of kept axes during the analysis.
   * @param rowW A vector of values specifying the weights of rows.
   * @param colW A vector of values specifying the weights of columns.
   * @param centered If true the input matrix is centered according to the column means.
   * @param scaled If true the input matrix is normalized according to the standard deviations of columns.
   * @param tol Tolerance threshold for null eigenvalues (a value less than tol times the first one is considered as null)
   * @param verbose Should warnings be dispayed.
   * @throw Exception if an error occured.
   */
  PrincipalComponentAnalysis(
    const Matrix<double>& data,
    unsigned int nbAxes,
    const std::vector<double>& rowW,
    const std::vector<double>& colW,
    bool centered = true,
    bool scaled = true,
    double tol = 0.0000001,
    bool verbose = true)
  throw (Exception);

  /**
   * @brief Build a new PrincipalComponentAnalysis object and specify default row and column weights.
   *
   * @param data The input data (a RowMatrix) to analyse.
   * @param nbAxes The number of kept axes during the analysis.
   * @param centered If true the input matrix is centered according to the column means.
   * @param scaled If true the input matrix is normalized according to the standard deviations of columns.
   * @param tol Tolerance threshold for null eigenvalues (a value less than tol times the first one is considered as null)
   * @param verbose Should warnings be dispayed.
   * @throw Exception if an error occured.
   */
  PrincipalComponentAnalysis(
    const Matrix<double>& data,
    unsigned int nbAxes,
    bool centered = true,
    bool scaled = true,
    double tol = 0.0000001,
    bool verbose = true)
  throw (Exception);

  virtual ~PrincipalComponentAnalysis() {}

  PrincipalComponentAnalysis* clone() const { return new PrincipalComponentAnalysis(*this); }

public:
  /**
   * @brief This function allows to center an input matrix from its column means.
   *
   * @param matrix The input data (a Matrix) to center.
   * @param rowW A vector with row weights.
   */
  static void center(Matrix<double>& matrix, const std::vector<double>& rowW) throw (Exception);

  /**
   * @brief This function allows to center an input matrix from its column means.
   *
   * @param matrix The input data (a Matrix) to center.
   * @param rowW A vector with row weights.
   */
  static void scale(Matrix<double>& matrix, const std::vector<double>& rowW) throw (Exception);

public:
  const std::vector<double>& getColumnMeans() const throw (Exception) { return columnMeans_; }
  const std::vector<double>& getColumnSd() const throw (Exception) { return columnSd_; }
};

} // end of namespace bpp.

#endif  // _PRINCIPALCOMPONENTANALYSIS_H_

