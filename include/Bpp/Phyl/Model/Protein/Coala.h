//
// File: Coala.h
// Created by: Bastien Boussau
// Created on: Tue May 18 15:23:20 2010
// Modified by: Mathieu Groussin
// Modified on: Sun Mar 13 12:00:00 2011
//

/*
   Copyright or � or Copr. CNRS, (November 16, 2004)
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


#ifndef _COALA_H_
#define _COALA_H_

#include "ProteinSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "CoalaCore.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

using namespace std;


namespace bpp
{
/**
 * @brief The Coala branch-heterogeneous amino-acid substitution model.
 *
 * This branch-heterogeneous model allows each branch to have its own set of amino acid equilibrium frequencies. It makes use of a Correspondence Analysis to reduce the number of parameters to be
 * optimized through Maximum Likelihood, focusing on most of the compositional variation observed in the data. The same COA is used for all branches.
 * An empirical exchangeability matrix is used, common to all branches. The choice of this matrix is up to the user. A user-defined empirical matrix can also be employed.
 * The model may also be used as a branch-homogeneous but non-stationary model, where the same estimated equilibrium frequencies are used for all branches, but different equilibrium frequencies
 * (that are optimized) are used on the root. Finally, the model may be employed as a branch-homogeneous and stationary model, where the frequencies at the root are similar to the ones on branches.
 *
 * @author Mathieu Groussin
 * @param alpha The alphabet (Protein)
 * @param nbAxes The number of principal axes of the COA that have to be taken into account to optimize the 20 branch-specific equilibrium frequencies. This number is common to all branches, as
 * well as on the root, where frequencies are optimized with a MVAprotein object (See the ProteinFrequenciesSet class).
 * @param exch The exchangeability matrix. The matrices currently available are DSO78, JTT92, WAG01 or LG08. A user-defined matrix can be specified with the 'file' argument.
 * @param file [optional] Used only to specify the file containing the user-defined exchangeabilities, written in PAML format.
 */

class Coala :
  public ProteinSubstitutionModel,
  public AbstractReversibleSubstitutionModel,
  public CoalaCore
{
protected:
  bool init_;
  unsigned int nbrOfAxes_;
  string exch_;
  string file_;

public:
  Coala(const ProteicAlphabet* alpha,
        unsigned int nbAxes = 0,
        const string exch = "LG08",
        string file = "");

  virtual ~Coala() {}

#ifndef NO_VIRTUAL_COV
  Coala*
#else
  Clonable*
#endif
  clone() const { return new Coala(*this); }

public:
  string getName() const {return "Coala"; }
  string getExch() const {return exch_; }
  void setFreqFromData(const SequenceContainer& data, bool param = true);
  string getEmpiricalMatrixFile() const {return file_; }

protected:
  void readFromFile(string& file);
  void computeEquilibriumFrequencies();
  void updateMatrices();
};
} // end of namespace bpp.
#endif  // _COALA_H_
