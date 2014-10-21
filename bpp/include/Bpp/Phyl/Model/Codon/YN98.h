//
// File: YN98.h
// Created by: Laurent Gueguen
// Created on: July 2009
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _YN98_H_
#define _YN98_H_

#include "../AbstractBiblioSubstitutionModel.h"
#include "CodonDistanceFrequenciesSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Yang and Nielsen (1998) substitution model for codons.
 * @author Laurent Guéguen
 *
 * This model has one rate of transitions and one rate of
 * transversion. It also allows distinct equilibrium frequencies
 * between codons. A multiplicative factor accounts for the selective
 * restraints at the amino acid level, depending on the synonymy of
 * the amino acids.
 *
 * For codons @f$i=i_1i_2i_3@f$ and @f$j=j_1j_2j_3@f$, the generator
 * term @f$Q_{ij} (i \neq j)@f$ is:
 *
 * 0 if 2 or 3 of the pair @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ are different.
 *
 * @f$\mu \pi_j \omega@f$ if exactly 1 of the pairs @f$(i_1,j_1)
 * (i_2,j_2) (i_3,j_3) @f$ is different, that difference is a
 * transversion and amino acids coded by i and j are different.
 *
 * @f$\mu \pi_j @f$ if exactly 1 of the pairs @f$(i_1,j_1) (i_2,j_2)
 * (i_3,j_3) @f$ is different, that difference is a transversion and
 * amino acids coded by i and j are the same.
 *
 * @f$\mu \kappa \pi_j \omega@f$ if exactly 1 of the pairs
 * @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different, that difference
 * is a transition and amino acids coded by i and j are different.
 *
 * @f$\mu \kappa \pi_j @f$ if exactly 1 of the pairs @f$(i_1,j_1)
 * (i_2,j_2) (i_3,j_3) @f$ is different, that difference is a
 * transition and amino acids coded by @f$i@f$ and @f$j@f$ are the
 * same.
 *
 * @f$\mu@f$ is a normalization factor.
 *
 * This model includes 2 parameters (@f$\kappa@f$ and @f$\omega@f$).
 * The codon frequencies @f$\pi_j@f$ are either observed or infered.
 *
 * Reference:
 * -  Yang Z. and Nielsen R. (1998), _Journal of Molecular Evolution_ 46:409--418.
 */
class YN98 :
    public AbstractBiblioSubstitutionModel,
    public virtual ReversibleSubstitutionModel
{
private:
  std::auto_ptr<CodonDistanceFrequenciesSubstitutionModel> pmodel_;

public:
  YN98(const GeneticCode* gc, FrequenciesSet* codonFreqs);

  YN98(const YN98& yn98);

  YN98& operator=(const YN98&);

  virtual ~YN98() {}

  YN98* clone() const { return new YN98(*this); }

public:
  std::string getName() const { return "YN98"; }

  const SubstitutionModel& getModel() const { return *pmodel_.get(); }

private:
  SubstitutionModel& getModel() { return *pmodel_.get(); }

};

} // end of namespace bpp.

#endif  // _YN98_H_

