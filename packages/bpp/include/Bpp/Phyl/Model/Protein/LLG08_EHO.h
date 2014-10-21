//
// File: LLG08_EHO.h
// Created by: Laurent Gueguen
// Created on: jeudi 7 octobre 2010, à 21h 47
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

#ifndef _LLG08_EHO_H_
#define _LLG08_EHO_H_

#include "../MixtureOfSubstitutionModels.h"
#include "ProteinSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "../AbstractBiblioMixedSubstitutionModel.h"

using namespace std;

namespace bpp
{
/**
 * @brief The Le et al  (2008) EH0 substitution model for proteins.
 * @author Laurent Guéguen
 *
 * This model is a mixture of three models corresponding to
 * extended/helix/other sites in proteins. The models are considered
 * in this order.
 *
 *
 * This model includes 4 parameters :
 *
 * - relrate1 is the relative rate of model of extended sites;
 * - relrate2 is the relative rate of helix sites;
 * - relproba1 is the proportion  of extended sites;
 * - relproba2 is the ratio of the proportions of helix sites over the sum of
 * the proportion of helix sites plus the proportion of other sites.
 *
 * Important: See the relation between these parameters and the rates
 * and probabilities of the models in the description of
 * MixtureOfSubstitutionModels class.
 *
 * Reference:
 *
 * Le S.Q., Lartillot N., Gascuel O. (2008) Phil. Trans. R. Soc. B 363:3965--3976.
 */

class LLG08_EHO :
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
  auto_ptr<MixtureOfSubstitutionModels> pmixmodel_;

public:
  /**
   * @brief Build a  EH0 model, with original equilibrium frequencies, probabilities and rates.
   *
   * @param alpha A proteic alphabet.
   *
   */

  LLG08_EHO(const ProteicAlphabet* alpha);

  ~LLG08_EHO();

  LLG08_EHO* clone() const { return new LLG08_EHO(*this); }

  LLG08_EHO(const LLG08_EHO&);

  LLG08_EHO& operator=(const LLG08_EHO&);

  const SubstitutionModel& getModel() const { return *pmixmodel_.get(); }

  const MixedSubstitutionModel& getMixedModel() const { return *pmixmodel_.get(); }

  std::string getName() const { return "LLG08_EHO"; }
  
private:
  SubstitutionModel& getModel() { return *pmixmodel_.get(); }

  MixedSubstitutionModel& getMixedModel() { return *pmixmodel_.get(); }

};
} // end of namespace bpp.

#endif  // _LLG08_EHO_H_

