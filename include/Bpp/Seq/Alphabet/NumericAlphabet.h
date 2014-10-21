//
// File: NumericAlphabet.h
// Created by: Laurent Gueguen
// Created on: March 2010
//

/*
   Copyright or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for sequences analysis.

   This software is governed by the CeCILL license under French law and
   abiding by the rules of distribution of free software. You can use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and rights to copy,
   modify and redistribute granted by the license, users are provided
   only with a limited warranty and the software's author, the holder of
   the economic rights, and the successive licensors have only limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading, using, modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean that it is complicated to manipulate, and that also
   therefore means that it is reserved for developers and experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards
   their requirements in conditions enabling the security of their
   systems and/or data to be ensured and, more generally, to use and
   operate it in the same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */


#ifndef _NUMERICALPHABET_H_
#define _NUMERICALPHABET_H_

#include "AbstractAlphabet.h"
#include "AlphabetNumericState.h"

#include <Bpp/Numeric/Prob.all>

#include <string>

/**
 * @brief This alphabet is used to deal NumericAlphabet
 * @author Laurent Guéguen
 */

namespace bpp
{
class NumericAlphabet : public AbstractAlphabet
{
private:
  const UniformDiscreteDistribution* pdd_;

  std::map<double, unsigned int> values_;
  
public:
  /**
   *@ brief Construction from a UniformDiscreteDistribution. This
   * UniformDiscreteDistribution is cloned.
   *
   */
  
  NumericAlphabet(const UniformDiscreteDistribution&);
  ~NumericAlphabet() { delete pdd_;}
  NumericAlphabet(const NumericAlphabet&);
  NumericAlphabet& operator=(const NumericAlphabet&);
  
public:
  void setState(unsigned int pos, const AlphabetNumericState&);
  void registerState(const AlphabetNumericState& ans);
  
  bool containsGap(const std::string& state) const throw (BadCharException);

  unsigned int getSize() const;
  unsigned int getNumberOfTypes() const;
  int getUnknownCharacterCode() const { return -1; }
  bool isGap(int state) const;
  std::vector<int> getAlias(int state) const throw (BadIntException);
  std::vector<std::string> getAlias(const std::string& state) const throw (BadCharException);
  bool isUnresolved(int state) const;
  bool isUnresolved(const std::string& state) const;

  std::string getAlphabetType() const { return "Numeric alphabet"; }

  AlphabetNumericState& getStateAt(unsigned int pos)  throw (IndexOutOfBoundsException);
  const AlphabetNumericState& getStateAt(unsigned int pos) const throw (IndexOutOfBoundsException);
  /**
   *@ brief Specific methods
   *
   */

  /**
   *@brief Returns the difference between successive values
   *
   */

  double getDelta() const;
  /**
   *@brief Returns the value for the character number 
   *
   */
  
  double intToValue(int state) const throw (BadIntException);

  /**
   *@brief Returns the CategoryIndex of the category to which the value belongs.
   *
   */

  unsigned int valueToInt(double value) const;


  /**
   * @brief Re-update the maps.
   */
  void remap();
  
};
}
#endif // _NUMERICALPHABET_H_

