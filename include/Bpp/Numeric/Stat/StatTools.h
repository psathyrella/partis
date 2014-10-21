//
// File: StatTools.h
// Created by: Julien Dutheil
// Created on: Sun Jan 30 19:10 2011
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

#ifndef _STATTOOLS_H_
#define _STATTOOLS_H_

//From the STL:
#include <vector>
#include <cstddef>

namespace bpp {

/**
 * @brief Statistics tools and utilitary functions.
 */
class StatTools
{
  private:
    struct PValue_ {
      double pvalue_;
      size_t index_;
      PValue_(double pvalue, size_t index):
        pvalue_(pvalue), index_(index) {}

      bool operator<(const PValue_& pvalue) const {
        return pvalue.pvalue_ < pvalue_;
      }
    };

  public:
    /**
     * @brief Compute the false discovery rate for a set of input p-values, using Benjamini and Hochberg's 'FDR' method.
     * 
     * The false discovery rate is computed by sorting all pvalues.
     * The FDR r is calculated with the formula
     * @f$ r = p * n / i@f$
     * where p is the p-value, n is the number of tests (the size of the input vector) and i is the rank of the p-value, that is the index in the sorted array.
     * 
     * References:
     * - Benjamini, Y and Hochberg, Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society, Series B (Methodological) 57(1):289-300.
     * - Verhoeven, KJF; Simonsen, KL; M. McIntyre, LM (2005). Implementing false discovery rate control: increasing your power. Oikos. 108(3):643-647.
     *
     * @author Julien Dutheil
     * @param pvalues The input p-values.
     * @return The corresponding false discovery rates.
     */
    static std::vector<double> computeFdr(const std::vector<double>& pvalues);
};

} //end of namespace bpp.

#endif //_STATTOOLS_H_

