//
// File: YNGKP_M1.h
// Created by: Laurent Gueguen
// Created on: May 2010
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

#ifndef _YNGKP_M1_H_
#define _YNGKP_M1_H_

#include "../AbstractBiblioMixedSubstitutionModel.h"
#include "../MixtureOfASubstitutionModel.h"
#include "../FrequenciesSet/CodonFrequenciesSet.h"

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{

/**
 * @brief The Yang et al (2000) M1 substitution model for codons, with
 * the more realistic modification in Wong & al (2004).
 * @author Laurent Guéguen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter to allow it to
 * vary among sites. A site is either negatively selected @f$ 0 <
 * \omega_0 < 1 @f$ (with probability @f$p_0 @f$), or neutral (@f$
 * \omega_1 = 1 @f$) with probability @f$1-p_0 @f$.
 *
 * The synonymous rates must be the same between both models, so the
 * overall rates of the models are modified to respect this constraint
 * and such that the mean rate of the mixed model equals one.
 *
 * This model includes 3 parameters (@f$\kappa@f$, @f$ p0 @f$ and
 * @f$\omega@f$). The codon frequencies @f$\pi_j@f$ are either
 * observed or infered.
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 * 
 * Wong, W. S. W., Z. Yang, N. Goldman, and R. Nielsen. (2004)
 * Genetics 168:1041--1051.
 */
class YNGKP_M1:
    public AbstractBiblioMixedSubstitutionModel,
    virtual public ReversibleSubstitutionModel
{
private:
  std::auto_ptr<MixtureOfASubstitutionModel> pmixmodel_;


  /*
   *@brief indexes of 2 codons between which the substitution is
   * synonymous, to set a basis to the homogeneization of the rates.
   *
   */

  int synfrom_, synto_;
  
public:
  YNGKP_M1(const GeneticCode* gc, FrequenciesSet* codonFreqs);

  virtual ~YNGKP_M1();
  
  YNGKP_M1* clone() const { return new YNGKP_M1(*this); }

  YNGKP_M1(const YNGKP_M1&);

  YNGKP_M1& operator=(const YNGKP_M1&);

protected:
  void updateMatrices();

public:
  const SubstitutionModel& getModel() const { return *pmixmodel_.get(); }

  const MixedSubstitutionModel& getMixedModel() const { return *pmixmodel_.get(); }

  std::string getName() const { return "YNGKP_M1"; }

private:
  SubstitutionModel& getModel() { return *pmixmodel_.get(); }
  
  MixedSubstitutionModel& getMixedModel() { return *pmixmodel_.get(); }

  const FrequenciesSet* getFrequenciesSet() const {return pmixmodel_->getNModel(1)->getFrequenciesSet();}
};

} //end of namespace bpp.

#endif	//_YNGKP_M1_H_

