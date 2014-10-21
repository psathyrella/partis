//
// File: AbstractCodonDistanceSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: jeudi 15 septembre 2011, à 21h 11
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

#ifndef _ABSTRACTCODONDISTANCESUBSTITUTIONMODEL_H_
#define _ABSTRACTCODONDISTANCESUBSTITUTIONMODEL_H_

#include "CodonSubstitutionModel.h"
#include <Bpp/Numeric/AbstractParameterAliasable.h>


// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Abstract class for modelling of non-synonymous abd
 *  synonymous substitution rates in codon models.
 *
 * @author Laurent Guéguen
 *
 * If a distance @f$d@f$ between amino-acids is defined, the
 *  non-synonymous rate is multiplied with, if the coded amino-acids
 *  are @f$x@f$ and @f$y@f$, @f$\beta*\exp(-\alpha.d(x,y))@f$ with
 *  non-negative parameter \c "alpha" and positive parameter \c
 *  "beta".
 *
 * If such a distance is not defined, the non-synonymous substitution
 *  rate is multiplied with @f$\beta@f$ with positive parameter \c
 *  "beta" (ie @f$d=0@f$).
 *
 * If paramSynRate is true, the synonymous substitution rate is
 *  multiplied with @f$\gamma@f$ (with optional positive parameter \c
 *  "gamma"), else it is multiplied with 1.
 *
 * References:
 * - Goldman N. and Yang Z. (1994), _Molecular Biology And Evolution_ 11(5) 725--736. 
 * - Kosakovsky Pond, S. and Muse, S.V. (2005), _Molecular Biology And Evolution_,
 *   22(12), 2375--2385.
 * - Mayrose, I. and Doron-Faigenboim, A. and Bacharach, E. and Pupko T.
 *   (2007), Bioinformatics, 23, i319--i327.
 */

class AbstractCodonDistanceSubstitutionModel :
  public virtual CodonSubstitutionModel,
  public virtual AbstractParameterAliasable
{
private:
  const AlphabetIndex2* pdistance_;

  double alpha_, beta_;

  double gamma_;
public:
  /**
   * @brief Build a new AbstractCodonDistanceSubstitutionModel object from
   *  a pointer to NucleotideSubstitutionModel.
   *
   * @param pdist optional pointer to a distance between amino-acids
   * @param prefix the Namespace
   * @param paramSynRate is true iff synonymous rate is parametrised
   *       (default=false).
   */
  AbstractCodonDistanceSubstitutionModel(
    const AlphabetIndex2* pdist,
    const std::string& prefix,
    bool paramSynRate = false);

  AbstractCodonDistanceSubstitutionModel(const AbstractCodonDistanceSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pdistance_(model.pdistance_),
    alpha_(model.alpha_),
    beta_(model.beta_),
    gamma_(model.gamma_)
  {}

  AbstractCodonDistanceSubstitutionModel& operator=(
    const AbstractCodonDistanceSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pdistance_ = model.pdistance_;
    alpha_ = model.alpha_;
    beta_ = model.beta_;
    gamma_ = model.gamma_;
    return *this;
  }

  virtual ~AbstractCodonDistanceSubstitutionModel() {}

public:
  void fireParameterChanged(const ParameterList& parameters);

public:
  double getCodonsMulRate(size_t i, size_t j) const;
};

} // end of namespace bpp.

#endif

