//
// File: StandardGeneticCode.h
// Created by: Julien Dutheil
// Created on: Mon Oct 13 15:39:17 2003
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

#ifndef _STANDARDGENETICCODE_H_
#define _STANDARDGENETICCODE_H_

#include "GeneticCode.h"
#include "../Alphabet/NucleicAlphabet.h"

namespace bpp
{

/**
 * @brief This class implements the standard genetic code as describe on the NCBI 
 *        web site: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=t#SG1
 */
class StandardGeneticCode:
  public virtual GeneticCode
{
	public:
		StandardGeneticCode(const NucleicAlphabet* alphabet);
		
    virtual ~StandardGeneticCode() {}
	
    virtual StandardGeneticCode* clone() const {
      return new StandardGeneticCode(*this);
    }

  public:
    size_t getNumberOfStopCodons() const { return 3.; }
    
    std::vector<int> getStopCodonsAsInt() const {
      std::vector<int> v(3);
      v[0] = 48;
      v[1] = 50;
      v[02] = 56;
      return v;
    }

    std::vector<std::string> getStopCodonsAsChar() const {
      std::vector<std::string> v(3);
      v[0] = "TAA";
      v[1] = "TAG";
      v[2] = "TGA";
      return v;
    }
    
    bool isStop(int state) const throw (BadIntException) {
      //Test:
      codonAlphabet_.intToChar(state); //throw exception if invalid state!
      return (state == 48 || state == 50 || state == 56);
    }
    
    bool isStop(const std::string& state) const throw (BadCharException) {
      int i = codonAlphabet_.charToInt(state);
      return (i == 48 || i == 50 || i == 56);
    }
 
    bool isAltStart(int state) const throw (BadIntException) {
      //Test:
      codonAlphabet_.intToChar(state); //throw exception if invalid state!
      return (state == 62 || state == 30);
    }
    
    bool isAltStart(const std::string& state) const throw (BadCharException) {
      int i = codonAlphabet_.charToInt(state);
      return (i == 62 || i == 30);
    }

};

} //end of namespace bpp.

#endif	//_STANDARDGENETICCODE_H_

