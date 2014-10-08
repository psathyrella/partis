//
// File: MiyataAAChemicalDistance.h
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

#ifndef _MIYATAAACHEMICALDISTANCE_H_
#define _MIYATAAACHEMICALDISTANCE_H_

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
 * @brief Miyata et al. (1979) Amino-Acid chemical distance.
 *
 * Two kinds of matrix can be built:
 * - a symmetric one, with \f$I_{i,j} = I_{i,j}\f$,
 * - or a non-symmetric one, with \f$I_{i,j} = -I_{i,j}\f$.
 *
 * Reference:
 * Miyata, T., Miyazawa, S. and Yasunaga, T.
 * Two types of amino acid substitutions in protein evolution
 * J. Mol. Evol. 12, 219-236 (1979)
 *
 * Data from AAIndex2 database, Accession Number MIYT790101.
 */
class MiyataAAChemicalDistance :
  public virtual AlphabetIndex2
{
private:
  LinearMatrix<double> distanceMatrix_;
  const ProteicAlphabet* alpha_;
  bool sym_;

public:
  MiyataAAChemicalDistance();

  MiyataAAChemicalDistance(const MiyataAAChemicalDistance& md) :
    distanceMatrix_(md.distanceMatrix_),
    alpha_(md.alpha_),
    sym_(md.sym_)
  {}

  MiyataAAChemicalDistance& operator=(const MiyataAAChemicalDistance& md)
  {
    distanceMatrix_ = md.distanceMatrix_;
    alpha_ = md.alpha_;
    sym_ = md.sym_;
    return *this;
  }

  virtual ~MiyataAAChemicalDistance() {}

  MiyataAAChemicalDistance* clone() const { return new MiyataAAChemicalDistance(); }

public:
  /**
   * @name Methods from the AlphabetIndex2 interface.
   *
   * @{
   */
  double getIndex(int state1, int state2) const throw (BadIntException);
  double getIndex(const std::string& state1, const std::string& state2) const throw (BadCharException);
  const Alphabet* getAlphabet() const { return alpha_; }
  Matrix<double>* getIndexMatrix() const;
  /** @} */

public:
  void setSymmetric(bool yn) { sym_ = yn; }
  bool isSymmetric() const { return sym_; }
};
} // end of namespace bpp.

#endif // _MIYATAAACHEMICALDISTANCE_H_

