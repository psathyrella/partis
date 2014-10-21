//
// File: RNA.h
// Created by: Guillaume Deuchst
// Created on: Tue Jul 22 2003
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

#ifndef _RNA_H_
#define _RNA_H_

#include "NucleicAlphabet.h"

//From the STL:
#include <string>

namespace bpp
{

/**
 * @brief This alphabet is used to deal with RNA sequences.
 *
 * It supports all 4 nucleotides (A, U, G and C) with their standard denomination.
 * Gaps are coded by '-', unresolved characters are coded by 'X, N, O, 0 or ?'.
 * Extensive support for generic characters (e.g. 'P', 'Y', etc.) is provided.
 */
class RNA:
  public NucleicAlphabet
{
	public:
    /**
     * @param exclamationMarkCountsAsGap If yes, '!' characters are replaced by gaps.
     * Otherwise, they are counted as unknown characters.
     */
		RNA(bool exclamationMarkCountsAsGap = false);
		virtual ~RNA() {}

	public:
    std::vector<int> getAlias(int state) const throw (BadIntException);
    std::vector<std::string> getAlias(const std::string & state) const throw (BadCharException);
    int getGeneric(const std::vector<int> & states) const throw (BadIntException);
    std::string getGeneric(const std::vector<std::string> & states) const throw (BadCharException);
    std::string getAlphabetType() const { return "RNA alphabet"; }
};

} //end of namespace bpp.

#endif // _RNA_H_
