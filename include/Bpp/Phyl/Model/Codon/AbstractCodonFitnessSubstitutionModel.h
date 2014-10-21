//
// File: AbstractCodonFitnessSubstitutionModel.h
// Created by:  Fanny Pouyet
// Created on: mars 2012
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

# ifndef _ABSTRACTCODONFITNESSSUBSTITUTIONMODEL_H_
# define _ABSTRACTCODONFITNESSSUBSTITUTIONMODEL_H_

# include "CodonSubstitutionModel.h"
#include "../FrequenciesSet/CodonFrequenciesSet.h"
namespace bpp
{
  /**
   * @brief Abstract class for modelling of ratios of substitution
   * rates between codons, whatever they are synonymous or not.
   *
   * @author Fanny Pouyet, Laurent Guéguen
   *
   * The fitness of a codon is a value between 0 and 1 defining the
   * relative advantage of a codon, compared to others. If a codon
   * @f$i@f$ has a fitness @f$\phi_i@f$ and another one (@f$j@f$) has
   * a fitness @f$\phi_j@f$, the substitution rate from codon @f$i@f$
   * to codon @f$j@f$ is multiplied by
   * \f[-\frac{\log(\frac{\phi_i}{\phi_j})}{1-\frac{\phi_i}{\phi_j}}\f]
   *
   * The set of fitnesses is implemented through a Codon
   * FrequenciesSet object. The parameters are named \c
   * "fit_NameOfTheParameterInTheFrequenciesSet".
   */


  class AbstractCodonFitnessSubstitutionModel :
    public virtual CodonSubstitutionModel,
    public virtual AbstractParameterAliasable
  {
  private:
    FrequenciesSet* pfitset_;
    std::string fitName_;
  public:
    AbstractCodonFitnessSubstitutionModel(FrequenciesSet* pfitset, const std::string& prefix);
    AbstractCodonFitnessSubstitutionModel(const AbstractCodonFitnessSubstitutionModel& model):
      AbstractParameterAliasable(model),
      pfitset_(model.pfitset_->clone()),
      fitName_(model.fitName_)
    {}

    AbstractCodonFitnessSubstitutionModel& operator=(const AbstractCodonFitnessSubstitutionModel& model){
      AbstractParameterAliasable::operator=(model);
      if (pfitset_) delete pfitset_;
      pfitset_ = model.pfitset_->clone();
      fitName_ = model.fitName_ ;
      return *this;
    }

    virtual ~AbstractCodonFitnessSubstitutionModel();

  public:
    void fireParameterChanged (const ParameterList& parameters);
    void setFreq(std::map<int, double>& frequencies);
    const FrequenciesSet& getFreq() const { return *pfitset_; }
    void setNamespace (const std::string& prefix){
      pfitset_->setNamespace(prefix + fitName_);
    }

    double getCodonsMulRate(size_t i, size_t j) const;

    const FrequenciesSet* getFitness() const { return pfitset_;}

  };
} // end of namespace bpp
# endif
