//
// File: TripletSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: Tue Dec 24 11:03:53 2003
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

#ifndef _TRIPLETSUBSTITUTIONMODEL_H_
#define _TRIPLETSUBSTITUTIONMODEL_H_

#include "../WordSubstitutionModel.h"
#include "../Nucleotide/NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>

namespace bpp
{

/**
 * @brief Class for neutral substitution models on triplets,
 * which correspond to codons that do not have any significance
 * (whether they are STOP or functional).
 * @author Laurent Guéguen
 *
 * Objects of this class are built from three substitution
 * models of NucleicAlphabets. No model is directly accessible. </p>
 *
 */
  
class TripletSubstitutionModel :
  public WordSubstitutionModel
{
public:

  /**
   *@brief Build a new TripletSubstitutionModel object from
   *a pointer to a NucleotideSubstitutionModel. 
   *
   *@param palph pointer to a CodonAlphabet
   *@param pmod  pointer to the NucleotideSubstitutionModel to be used
   *       in the three positions. It is owned by the instance.
   */
  
  TripletSubstitutionModel(const CodonAlphabet* palph,
                                     NucleotideSubstitutionModel* pmod);
  
  /**
   *@brief Build a new TripletSubstitutionModel object
   *from three pointers to NucleotideSubstitutionModels. 
   *
   *@param palph pointer to a CodonAlphabet
   *@param pmod1, pmod2, pmod3 pointers to the
   *   NucleotideSubstitutionModels to use in the three positions. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used
   *   models are owned by the instance.
   */

  TripletSubstitutionModel(const CodonAlphabet* palph,
                                     NucleotideSubstitutionModel* pmod1,
                                     NucleotideSubstitutionModel* pmod2, 
                                     NucleotideSubstitutionModel* pmod3);

  ~TripletSubstitutionModel() {};
  
#ifndef NO_VIRTUAL_COV
  TripletSubstitutionModel*
#else
  Clonable*
#endif
  clone() const { return new TripletSubstitutionModel(*this);}
  
public:
  std::string getName() const;
};

} //end of namespace bpp.

#endif	

