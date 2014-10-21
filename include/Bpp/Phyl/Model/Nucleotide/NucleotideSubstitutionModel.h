//
// File: NucleicSubstitutionModel.h
// Created by: Julien Dutheil
// Created on: Tue May 27 11:03:53 2003
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

#ifndef _NUCLEICSUBSTITUTIONMODEL_H_
#define _NUCLEICSUBSTITUTIONMODEL_H_

#include "../SubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{

/**
 * @brief Specialisation interface for nucleotide substitution model.
 */
  class NucleotideSubstitutionModel :
    public virtual SubstitutionModel
  {
  public:
    virtual ~NucleotideSubstitutionModel() {}

#ifndef NO_VIRTUAL_COV
    NucleotideSubstitutionModel*
#else
    Clonable*
#endif
    clone() const = 0;

  public:
    virtual size_t getNumberOfStates() const { return 4; }

  };

} //end of namespace bpp.

#endif	//_NUCLEICSUBSTITUTIONMODEL_H_

