//
// File: ProteinFrequenciesSet.h
// Created by: Bastien Boussau
//             Julien Dutheil
// Created on: Tue Aug 21 2007
//

/*
   Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _PROTEINFREQUENCIESSET_H_
#define _PROTEINFREQUENCIESSET_H_

#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include "FrequenciesSet.h"

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies for proteins.
 */
class ProteinFrequenciesSet :
  public virtual FrequenciesSet
{
public:
#ifndef NO_VIRTUAL_COV
  ProteinFrequenciesSet* clone() const = 0;

  const ProteicAlphabet* getAlphabet() const = 0;
#endif
};

/**
 * @brief Protein FrequenciesSet using 19 independent parameters to
 * modelize the 20 frequencies.
 *
 * The parameters are called @f$ \theta_{i \in 1..19} @f$, and are
 * initialized so that all frequencies are equal to 0.005. The
 * parametrization depends on the method used. Default
 * method is 1 (ie global ratio).
 *
 * @see Simplex
 *
 **/
  
class FullProteinFrequenciesSet :
  public virtual ProteinFrequenciesSet,
  public FullFrequenciesSet
{
public:
  FullProteinFrequenciesSet(const ProteicAlphabet* alphabet, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full") :
    FullFrequenciesSet(alphabet, allowNullFreqs, method, name) {}
  FullProteinFrequenciesSet(const ProteicAlphabet* alphabet, const std::vector<double>& initFreqs, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full") :
    FullFrequenciesSet(alphabet, initFreqs, allowNullFreqs, method, name) {}

#ifndef NO_VIRTUAL_COV
  FullProteinFrequenciesSet*
#else
  Clonable*
#endif
  clone() const { return new FullProteinFrequenciesSet(*this); }

public:
#ifndef NO_VIRTUAL_COV
  const ProteicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const ProteicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif
};

/**
 * @brief FrequenciesSet useful for homogeneous and stationary models, protein implementation
 *
 * This set contains no parameter.
 */
class FixedProteinFrequenciesSet :
  public virtual ProteinFrequenciesSet,
  public FixedFrequenciesSet
{
public:
  FixedProteinFrequenciesSet(const ProteicAlphabet* alphabet, const std::vector<double>& initFreqs, const std::string& name = "Fixed") :
    FixedFrequenciesSet(alphabet, initFreqs, name) {}

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet.
   */
  FixedProteinFrequenciesSet(const ProteicAlphabet* alphabet, const std::string& name = "Fixed") :
    FixedFrequenciesSet(alphabet, alphabet->getSize(), name) {}

#ifndef NO_VIRTUAL_COV
  FixedProteinFrequenciesSet*
#else
  FixedFrequenciesSet*
#endif
  clone() const { return new FixedProteinFrequenciesSet(*this); }

#ifndef NO_VIRTUAL_COV
  const ProteicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const ProteicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif
};


} // end of namespace bpp.

#endif // _PROTEINFREQUENCIESSET_H_


