//
// File: CoalaCore.h
// Created by: Mathieu Groussin
// Modified on: Sun Mar 13 12:00:00 2011
//

/*
   Copyright or � or Copr. Bio++ Development Team, (November 16, 2004)
   PCA
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


#ifndef _COALACORE_H_
#define _COALACORE_H_

// From Core:
#include <Bpp/Numeric/Matrix/Matrix.h>

#include <Bpp/Numeric/ParameterList.h>
// From SeqLib:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Container/SequenceContainer.h>


using namespace std;

namespace bpp
{
/**
 * @brief This class is the core class inherited by the Coala class. COaLA is a branch-heterogeneous amino-acid substitution model.
 *
 * This class allows to compute the COA from the alignment, to define the parameter (axis positions), and implements a function used to compute the equilibrium frequencies from a set of
 * coordinates along the principal axes of the COA.
 *
 * @author Mathieu Groussin
 * @param nbAxes The number of principal axes of the COA that have to be taken into account to optimize the 20 branch-specific equilibrium frequencies. This number is common to all branches, as
 * well as on the root, where frequencies are optimized with a MVAprotein object (See the ProteinFrequenciesSet class).
 * @param exch The exchangeability matrix. The matrices currently available are DSO78, JTT92, WAG01 or LG08. A user-defined matrix can be specified with the 'file' argument.
 */

class CoalaCore
{
protected:
  bool init_;
  size_t nbrOfAxes_;
  string exch_;
  RowMatrix<double> P_;
  RowMatrix<double> R_;
  vector<double> colWeights_;
  map<string, string> paramValues_;

public:
  CoalaCore(size_t nbAxes = 0, const string& exch = "LG08");

  virtual ~CoalaCore() {}

  CoalaCore* clone() const { return new CoalaCore(*this); }

public:
  size_t getNbrOfAxes() const { return nbrOfAxes_; }
  const RowMatrix<double>& getTppalAxesMatrix() const { return P_; }
  const RowMatrix<double>& getRowCoordinates() const { return R_; }
  const vector<double>& getColumnWeights() const { return colWeights_; }
  void setParamValues(const map<string, string>& valuesSettings) { paramValues_ = valuesSettings; }

protected:
  ParameterList computeCOA(const SequenceContainer& data, bool param = true);
  
  vector<double> prodMatrixVector(RowMatrix<double>& P, vector<double>& V);
};

} // end of namespace bpp.

#endif  // _COALACORE_H_

