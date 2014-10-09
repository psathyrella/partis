//
// File: LG10_EX_EHO.h
// Created by: Mathieu Groussin
// Created on: Thursday 28 Mar 2013
//

/*
 Copyright or © or Copr. CNRS, (November 16, 2004)
 
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

#ifndef _LG10_EX_EHO_H_
#define _LG10_EX_EHO_H_

#include "../MixtureOfSubstitutionModels.h"
#include "ProteinSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "../AbstractBiblioMixedSubstitutionModel.h"

using namespace std;

namespace bpp
{
/**
 * @brief The Le and Gascuel (2010) EX_EHO substitution model for proteins.
 * @author Mathieu Groussin
 *
 * This model is a mixture of six models. It combines two previously published models, EX2 and EHO.
 * Sites are classified as: exposed & extended, buried & extended, exposed & alpha-helix
 * buried & alpha-helix, exposed & other, or buried & other.
 * The submodels are called BUR_EXT, BUR_HEL, BUR_OTH, EXP_EXT, EXP_HEL and EXP_OTH.
 * 
 * The model includes 10 parameters :
 *
 * - relrate1 to relrate5 are the relative rates of the submodels;
 * - relproba1 to relproba5 are the relative proportions of the submodels;
 *
 * Important: See the relation between these parameters and the
 * rates and probabilities of the models in the description of
 * MixtureOfSubstitutionModels class.
 *
 * Reference:
 *
 * Le S.Q., Gascuel O. (2010) Syst. Biol. 59(3):277–287
 */

class LG10_EX_EHO :
  public AbstractBiblioMixedSubstitutionModel
{
public:
  class EmbeddedModel :
    public virtual ProteinSubstitutionModel,
    public AbstractReversibleSubstitutionModel
  {
private:
    double proportion_;
    string name_;

public:
    EmbeddedModel(const ProteicAlphabet* alpha, string name);
    ~EmbeddedModel(){}
    EmbeddedModel* clone() const { return new EmbeddedModel(*this); }
    string getName() const { return name_;}
    double getProportion() const { return proportion_;}
  };

private:
  std::auto_ptr<MixtureOfSubstitutionModels> pmixmodel_;

public:
  /**
   * @brief Build a EX_EHO model, with original equilibrium frequencies, probabilities and rates.
   *
   * @param alpha A proteic alphabet.
   *
   */

  LG10_EX_EHO(const ProteicAlphabet* alpha);

  ~LG10_EX_EHO();

  LG10_EX_EHO* clone() const { return new LG10_EX_EHO(*this); }

  LG10_EX_EHO(const LG10_EX_EHO&);

  LG10_EX_EHO& operator=(const LG10_EX_EHO&);

  const SubstitutionModel& getModel() const { return *pmixmodel_.get(); }

  const MixedSubstitutionModel& getMixedModel() const { return *pmixmodel_.get(); }
  
  std::string getName() const { return "LG10_EX_EHO"; }

private:
  SubstitutionModel& getModel() { return *pmixmodel_.get(); }
  
  MixedSubstitutionModel& getMixedModel() { return *pmixmodel_.get(); }
  
};
} // end of namespace bpp.

#endif  // _LG10_EX_EHO_H_

