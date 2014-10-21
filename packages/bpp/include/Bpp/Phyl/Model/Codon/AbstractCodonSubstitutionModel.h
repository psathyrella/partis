//
// File: AbstractCodonSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: Tue Dec 24 11:03:53 2003
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

#ifndef _ABSTRACTCODONSUBSTITUTIONMODEL_H_
#define _ABSTRACTCODONSUBSTITUTIONMODEL_H_

#include "../AbstractWordSubstitutionModel.h"
#include "../Nucleotide/NucleotideSubstitutionModel.h"
#include "CodonSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

// From the STL:
#include <memory>

namespace bpp
{
/**
 * @brief Abstract class for substitution models on codons.
 * @author Laurent Guéguen
 *
 * Objects of this class are built from either one (repeated three
 * times) or three different substitution models of NucleicAlphabets.
 * No model is directly accessible. </p>
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 * There is one substitution per codon per unit of time
 * on the equilibrium frequency, and each position has its specific rate.
 *
 * The parameters of this codon are the same as the ones of the models
 * used. Their names have a new prefix, "i_" where i stands for the
 * the phase (1,2 or 3) in the codon.
 */
class AbstractCodonSubstitutionModel :
  public virtual CodonSubstitutionModel,
  public virtual AbstractWordSubstitutionModel
{
private:
  /**
   * @brief boolean for the parametrization of the position relative
   * rates. Default : false.
   */
  bool hasParametrizedRates_;
  const GeneticCode* gCode_;

public:
  const GeneticCode* getGeneticCode() const { return gCode_; }
  
  /**
   * @brief Build a new AbstractCodonSubstitutionModel object from
   * a pointer to a NucleotideSubstitutionModel.
   *
   * @param gCode a pointer toward a genetic code. The codon alphabet from the genetic code will be used by the model class.
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param st string of the Namespace
   * @param paramRates boolean concerning the presence of position
   * relative rates (default: false)
   */
  AbstractCodonSubstitutionModel(
      const GeneticCode* gCode,
      NucleotideSubstitutionModel* pmod,
      const std::string& st,
      bool paramRates = false);

  /**
   * @brief Build a new AbstractCodonSubstitutionModel object
   * from three pointers to NucleotideSubstitutionModels.
   *
   * @param gCode a pointer toward a genetic code. This model instance will own the underlying GeneticCode object and delete it when required.
   *   The codon alphabet from the genetic code will be used by the model class.
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param st string of the Namespace
   * @param paramRates boolean concerning the presence of position
   * relative rates (default: false)
   */
  AbstractCodonSubstitutionModel(
      const GeneticCode* gCode,
      NucleotideSubstitutionModel* pmod1,
      NucleotideSubstitutionModel* pmod2,
      NucleotideSubstitutionModel* pmod3,
      const std::string& st,
      bool paramRates = false);

  virtual ~AbstractCodonSubstitutionModel() {}

  AbstractCodonSubstitutionModel(const AbstractCodonSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractSubstitutionModel(model),
    AbstractWordSubstitutionModel(model),
    hasParametrizedRates_(model.hasParametrizedRates_),
    gCode_(model.gCode_)
  {}

  AbstractCodonSubstitutionModel& operator=(const AbstractCodonSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractSubstitutionModel::operator=(model);
    AbstractWordSubstitutionModel::operator=(model);
    hasParametrizedRates_ = model.hasParametrizedRates_;
    gCode_ = model.gCode_;
    return *this;
  }

#ifndef NO_VIRTUAL_COV
  AbstractCodonSubstitutionModel*
#else
  Clonable*
#endif
  clone() const = 0;

protected:
  /**
   * @brief Method inherited from AbstractWordSubstitutionModel
   *
   * This method sets the rates to/from stop codons to zero and
   * performs the multiplication by the specific codon-codon rate.
   */
  void completeMatrices();

public:
  void updateMatrices();

  /**
   * @brief Method inherited from CodonSubstitutionModel
   *
   * Here this methods returns 1;
   *
   **/
  virtual double getCodonsMulRate(size_t i, size_t j) const { return 1.; }
};
} // end of namespace bpp.

#endif

