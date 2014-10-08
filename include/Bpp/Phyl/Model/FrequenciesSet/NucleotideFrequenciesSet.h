//
// File: NucleotideFrequenciesSet.h
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

#ifndef _NUCLEOTIDEFREQUENCIESSET_H_
#define _NUCLEOTIDEFREQUENCIESSET_H_

#include "FrequenciesSet.h"
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies for nucleotides.
 */
class NucleotideFrequenciesSet :
  public virtual FrequenciesSet
{
public:
#ifndef NO_VIRTUAL_COV
  NucleotideFrequenciesSet* clone() const = 0;

  const NucleicAlphabet* getAlphabet() const = 0;
#endif
};

/**
 * @brief Nucleotide FrequenciesSet using only one parameter, the GC content.
 */
class GCFrequenciesSet :
  public virtual NucleotideFrequenciesSet,
  public AbstractFrequenciesSet
{
public:
  GCFrequenciesSet(const NucleicAlphabet* alphabet) :
    AbstractFrequenciesSet(4, alphabet, "GC.", "GC")
  {
    addParameter_(new Parameter("GC.theta", 0.5, &Parameter::PROP_CONSTRAINT_IN));
    getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
  }

  GCFrequenciesSet(const NucleicAlphabet* alphabet, double theta) :
    AbstractFrequenciesSet(4, alphabet, "GC.", "GC")
  {
    addParameter_(new Parameter("GC.theta", theta, &Parameter::PROP_CONSTRAINT_IN));
    getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
    getFreq_(1) = getFreq_(2) = theta / 2.;
  }

#ifndef NO_VIRTUAL_COV
  GCFrequenciesSet*
#else
  Clonable*
#endif
  clone() const { return new GCFrequenciesSet(*this); }

public:
#ifndef NO_VIRTUAL_COV
  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif

  void setFrequencies(const std::vector<double>& frequencies);

protected:
  void fireParameterChanged(const ParameterList& parameters);
};

/**
 * @brief Nucleotide FrequenciesSet using three independent parameters
 * (theta, theta1, theta2) to modelize the four frequencies:
 *
 * \f[
 * \begin{cases}
 * \theta = \pi_C + \pi_G\\
 * \theta_1 = \frac{\pi_A}{1 - \theta} = \frac{\pi_A}{\pi_A + \pi_T}\\
 * \theta_2 = \frac{\pi_G}{\theta} = \frac{\pi_G}{\pi_C + \pi_G}\\
 * \end{cases}
 * \Longleftrightarrow
 * \begin{cases}
 * \pi_A = \theta_1 (1 - \theta)\\
 * \pi_C = (1 - \theta_2) \theta\\
 * \pi_G = \theta_2 \theta\\
 * \pi_T = (1 - \theta_1)(1 - \theta).
 * \end{cases}
 * \f]
 *
 * with \f$\pi_x\f$ the frequency of nucleotide \f$x\f$.
 *
 */
  
class FullNucleotideFrequenciesSet :
  public virtual NucleotideFrequenciesSet,
  public AbstractFrequenciesSet
{
public:
  FullNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, bool allowNullFreqs = false, const std::string& name = "Full");

  FullNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, double theta, double theta1, double theta2, bool allowNullFreqs = false, const std::string& name = "Full");

#ifndef NO_VIRTUAL_COV
  FullNucleotideFrequenciesSet*
#else
  Clonable*
#endif
  clone() const { return new FullNucleotideFrequenciesSet(*this); }

public:
#ifndef NO_VIRTUAL_COV
  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif

  void setFrequencies(const std::vector<double>& frequencies);

protected:
  void fireParameterChanged(const ParameterList& parameters);
};


/**
 * @brief FrequenciesSet useful for homogeneous and stationary models, nucleotide implementation
 *
 * This set contains no parameter.
 */
class FixedNucleotideFrequenciesSet :
  public virtual NucleotideFrequenciesSet,
  public FixedFrequenciesSet
{
public:
  FixedNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, const std::vector<double>& initFreqs, const std::string& name = "Fixed") :
    FixedFrequenciesSet(alphabet, initFreqs, name) {}

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet.
   */
  FixedNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, const std::string& name = "Fixed") :
    FixedFrequenciesSet(alphabet, alphabet->getSize(), name) {}

#ifndef NO_VIRTUAL_COV
  FixedNucleotideFrequenciesSet*
#else
  NucleotideFrequenciesSet*
#endif
  clone() const { return new FixedNucleotideFrequenciesSet(*this); }

#ifndef NO_VIRTUAL_COV
  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif
};


} // end of namespace bpp.

#endif // _NUCLEOTIDEFREQUENCIESSET_H_


