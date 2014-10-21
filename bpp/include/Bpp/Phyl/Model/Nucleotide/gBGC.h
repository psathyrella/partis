//
// File: gBGC.h
// Created by: Laurent Gueguen
// Created on: lundi 13 février 2012, à 09h 43
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

#ifndef _GBGC_H_
#define _GBGC_H_

#include "NucleotideSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"

#include <Bpp/Numeric/Constraints.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

// From Utils:

using namespace std;

namespace bpp
{
/**
 * @brief gBGC  model.
 * @author Laurent Guéguen
 *
 * modelling of GC biased gene-conversion.
 *
 * This model adds strand symetric GC biased gene conversion to a
 * given nucleotidic substitution model.
 *
 * In addition to the parameters of the basic nucleic model, the
 * biased gene conversion effect is parametrized by @f$ \gamma @f$
 * (denoted \c "gamma"), and the mutation rates from A<->T nucleotides
 * to C<->G nucleotides are multiplied by @f$ exp(\gamma)@f$.
 *
 * @see AbstractSubstitutionModel
 *
 */

  class gBGC :
    public virtual NucleotideSubstitutionModel,
    public AbstractSubstitutionModel
  {
  private:
    std::auto_ptr<NucleotideSubstitutionModel>  model_;
    std::string nestedPrefix_;

    /*
     * @brief the value of the bias.
     *
     */
  
    double gamma_;
  
  public:
    /*
     * @brief Build a new gBGC substitution model.
     *
     */

    gBGC(const NucleicAlphabet*, NucleotideSubstitutionModel* const, double gamma=0);

    gBGC(const gBGC&);

    gBGC& operator=(const gBGC& gbgc);

#ifndef NOVIRTUAL_COV_
    gBGC*
#else
    Clonable*
#endif
    clone() const { return new gBGC(*this); }

    ~gBGC()  {}

  public:
    std::string getName() const;

    size_t getNumberOfStates() const { return model_->getNumberOfStates(); }

    void fireParameterChanged(const ParameterList&);

    const SubstitutionModel* getNestedModel() const {return model_.get();}
    
    void updateMatrices();

    void setNamespace(const std::string&);
  };
}

#endif // _GBGC_H


