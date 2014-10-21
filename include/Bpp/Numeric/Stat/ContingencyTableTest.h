//
// File: ContingencyTableTest.h
// Created by: Julien Dutheil
// Created on: Thu Dec 09 14:20 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

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

#ifndef _CONTINGENCYTABLETEST_H_
#define _CONTINGENCYTABLETEST_H_

#include "StatTest.h"

//From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief Implements tests on contingency tables.
 *
 * Performs a chi square test on contingency tables.
 */
class ContingencyTableTest:
  public virtual StatTest
{
  private:
    double statistic_;
    double pvalue_;
    double df_;
    std::vector<size_t> margin1_;
    std::vector<size_t> margin2_;

  public:
    /**
     * @brief Build a new test object and perform computations.
     *
     * @param table The input contingency table.
     * @param nbPermutations If greater than 0, performs a randomization test instead of using the chisquare approximation.
     * @param warn Should a warning message be displayed in case of unsufficient observations?
     */
    ContingencyTableTest(const std::vector< std::vector<size_t> >& table, unsigned int nbPermutations = 0, bool warn = true);
    virtual ~ContingencyTableTest() {}

#ifndef NO_VIRTUAL_COV
    ContingencyTableTest*
#else
    Clonable*
#endif
    clone() const { return new ContingencyTableTest(*this); }

  public:
    std::string getName() const { return "Test on contingency table."; }
    double getStatistic() const { return statistic_; }
    double getPValue() const { return pvalue_; }
    double getDegreesOfFreedom() const { return df_; }
    const std::vector<size_t> getMarginRows() const { return margin1_; }
    const std::vector<size_t> getMarginColumns() const { return margin2_; }

};

} //end of namespace bpp.

#endif //_CONTINGENCYTABLETEST_H_


