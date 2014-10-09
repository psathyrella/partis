//
// File: AbstractCodonFrequenciesSubstitutionModel.h
// Created by: jeudi 15 septembre 2011, à 15h 02
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

#ifndef _ABSTRACTCODONFREQUENCIESSUBSTITUTIONMODEL_H_
#define _ABSTRACTCODONFREQUENCIESSUBSTITUTIONMODEL_H_

#include "AbstractCodonSubstitutionModel.h"
#include "../FrequenciesSet/CodonFrequenciesSet.h"

namespace bpp
{
/**
 * @brief Abstract Class for substitution models on codons
 *  parametrized by frequencies.
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, ans does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * @author Laurent Guéguen
 *
 * If we denote @f$F@f$ the equilibrium frequency, the generator term
 * defined from inherited and inheriting classes, @f$Q_{ij})@f$, is
 * multiplied by @f$F_{j}@f$.
 *
 */

class AbstractCodonFrequenciesSubstitutionModel :
  virtual public CodonSubstitutionModel,
  virtual public AbstractParameterAliasable
{
private:
  FrequenciesSet* pfreqset_;
  std::string freqName_;

public:
  /**
   *@brief Build a AbstractCodonFrequenciesSubstitutionModel instance
   *
   *@param pfreq pointer to the AbstractFrequenciesSet equilibrium frequencies.
   *        It is owned by the instance.
   *@param prefix the Namespace
   */

  AbstractCodonFrequenciesSubstitutionModel(FrequenciesSet* pfreq,
                                            const std::string& prefix);

  AbstractCodonFrequenciesSubstitutionModel(const AbstractCodonFrequenciesSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pfreqset_(model.pfreqset_->clone()),
    freqName_(model.freqName_)
  {}

  AbstractCodonFrequenciesSubstitutionModel& operator=(const AbstractCodonFrequenciesSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    if (pfreqset_) delete pfreqset_;
    pfreqset_   = model.pfreqset_->clone();
    freqName_   = model.freqName_;
    return *this;
  }

  virtual ~AbstractCodonFrequenciesSubstitutionModel();

  void fireParameterChanged(const ParameterList& parameters);

  void setFreq(std::map<int, double>& frequencies);

  const FrequenciesSet* getFrequenciesSet() const { return pfreqset_; }

  void setNamespace(const std::string& prefix)
  {
    pfreqset_->setNamespace(prefix + freqName_);
  }

  double getCodonsMulRate(size_t, size_t) const;
};
} // end of namespace bpp.

#endif

