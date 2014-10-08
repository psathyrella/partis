//
// File: DNAToRNA.h
// Created by: Julien Dutheil
// Created on: Sun Oct 12 14:39:29 2003
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#ifndef _DNATORNA_H_
#define _DNATORNA_H_

#include "Transliterator.h"
#include "Alphabet/DNA.h"
#include "Alphabet/RNA.h"
#include "Alphabet/AlphabetTools.h"

namespace bpp
{

/**
 * @brief Reverse iterator to convert RNA to DNA.
 *
 * This only changes U to T characters.
 *
 * @see NucleicAcidsReplication
 */
class DNAToRNA:
  public AbstractReverseTransliterator
{
	private:
		const Alphabet* dna_, * rna_;
	
	public:
		DNAToRNA(): AbstractReverseTransliterator(), dna_(&AlphabetTools::DNA_ALPHABET), rna_(&AlphabetTools::RNA_ALPHABET) {}

    DNAToRNA(const DNAToRNA& d2r): AbstractReverseTransliterator(d2r), dna_(d2r.dna_), rna_(d2r.rna_) {}
    
    DNAToRNA& operator=(const DNAToRNA& d2r)
    {
      AbstractReverseTransliterator::operator=(d2r);
      dna_ = d2r.dna_;
      rna_ = d2r.rna_;
      return *this;
    }
		
    virtual ~DNAToRNA() {}
	
	public:
		virtual const Alphabet* getSourceAlphabet() const { return dna_; }
		virtual const Alphabet* getTargetAlphabet() const { return rna_; }
		int translate(int state) const throw (BadIntException);		
    std::string translate(const std::string& state) const throw (BadCharException);		
		Sequence* translate(const Sequence& sequence) const throw (AlphabetMismatchException, Exception)
    {
      return AbstractReverseTransliterator::translate(sequence);
    }
		int reverse(int state) const throw (BadIntException);
    std::string reverse(const std::string& state) const throw (BadCharException);
		Sequence* reverse(const Sequence& sequence) const throw (AlphabetMismatchException, Exception)
    {
      return AbstractReverseTransliterator::reverse(sequence);
    }

};

} //end of namespace bpp.

#endif	//_DNATORNA_H_

