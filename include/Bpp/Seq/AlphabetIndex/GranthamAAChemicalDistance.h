//
// File: GranthamAAChemicalDistance.h
// Created by: Julien Dutheil
// Created on: Mon Feb 21 17:42 2005
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for sequences analysis.

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

#ifndef _GRANTHAMAACHEMICALDISTANCE_H_
#define _GRANTHAMAACHEMICALDISTANCE_H_

// from the STL:
#include <string>

#include "AlphabetIndex2.h"
#include "../Alphabet/ProteicAlphabet.h"
#include "../Alphabet/AlphabetExceptions.h"
#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{
/**
 * @brief Grantham (1974) Amino-Acid chemical distance.
 *
 * Two kinds of matrix can be built:
 * - a symmetric one, with \f$I_{i,j} = I_{i,j}\f$,
 * - or a non-symmetric one, with \f$I_{i,j} = -I_{i,j}\f$.
 * In the second case, which one of the two entries between \f$I_{i,j}\f$ and \f$I_{i,j}\f$ is positive is arbitrary by default.
 * It is also possible to use the coordinate on a principal component analysis between the elementary propoerties of the distance instead (setPC1Sign(true)).
 * The following R code was use in order to get those signs:
 * @code
 * library(seqinr)
 * data(aaindex)
 * data<-data.frame(composition=aaindex[["GRAR740101"]]$I,
 *                     polarity=aaindex[["GRAR740102"]]$I,
 *                       volume=aaindex[["GRAR740103"]]$I)
 * library(ade4)
 * pca<-dudi.pca(data)
 *
 * plot(pca$li[, 1:2], type="n")
 * text(pca$li[, 1:2], rownames(data))
 *
 * s.corcircle(pca$co)
 * layout(matrix(1:3,nrow=1))
 * a1<-pca$li[,1]; names(a1)<-rownames(data); dotchart(sort(a1))
 * a2<-pca$li[,2]; names(a2)<-rownames(data); dotchart(sort(a2))
 * a3<-pca$li[,3]; names(a3)<-rownames(data); dotchart(sort(a3))
 *
 * x<-pca$li[,1] #Contains the coordinates on the first axis.
 * m<-matrix(nrow=20, ncol=20)
 * for(i in 1:length(x))
 *   for(j in 1:length(x))
 *     m[i,j]<-sign(x[j] - x[i])
 *
 * @endcode
 *
 * Reference:
 * Grantham, R.
 * Amino acid difference formula to help explain protein evolution
 * Science 185, 862-864 (1974)
 *
 * Data from AAIndex2 database, Accession Number GRAR740104.
 */
class GranthamAAChemicalDistance :
  public virtual AlphabetIndex2
{
private:
  LinearMatrix<double> distanceMatrix_;
  LinearMatrix<double> signMatrix_;
  const ProteicAlphabet* alpha_;
  short int sign_;

public:
  GranthamAAChemicalDistance();

  GranthamAAChemicalDistance(const GranthamAAChemicalDistance& gd) :
    distanceMatrix_(gd.distanceMatrix_),
    signMatrix_(gd.signMatrix_),
    alpha_(gd.alpha_),
    sign_(gd.sign_)
  {}

  GranthamAAChemicalDistance& operator=(const GranthamAAChemicalDistance& gd)
  {
    distanceMatrix_ = gd.distanceMatrix_;
    signMatrix_ = gd.signMatrix_;
    alpha_ = gd.alpha_;
    sign_ = gd.sign_;
    return *this;
  }

  virtual ~GranthamAAChemicalDistance();

public:
  /**
   * @name Methods from the AlphabetIndex2 interface.
   *
   * @{
   */
  double getIndex(int state1, int state2) const throw (BadIntException);
  double getIndex(const std::string& state1, const std::string& state2) const throw (BadCharException);
  const Alphabet* getAlphabet() const { return alpha_; }
  GranthamAAChemicalDistance* clone() const { return new GranthamAAChemicalDistance(); }
  Matrix<double>* getIndexMatrix() const;
  /** @} */

public:
  void setSymmetric(bool yn) { sign_ = (yn ? SIGN_NONE : SIGN_ARBITRARY); }
  bool isSymmetric() const { return sign_ == SIGN_NONE; }
  /**
   * @brief The sign of the distance is computed using the coordinate on the first axis
   * of a principal component analysis with the 3 elementary properties (Volume, Polarity, Composition).
   * Otherwise, use the default arbitrary sign. Using this option will lead isSymmetric to return false.
   *
   * @param yn Tell is the PC1-based sign should be used instead of the arbitrary one.
   */
  void setPC1Sign(bool yn) { sign_ = (yn ? SIGN_PC1 : SIGN_ARBITRARY); }

  static short int SIGN_ARBITRARY;
  static short int SIGN_PC1;
  static short int SIGN_NONE;
};
} // end of namespace bpp.

#endif // _GRANTHAMAACHEMICALDISTANCE_H_

