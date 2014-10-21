//
// File: PhredPoly.h
// Created by: Sylvain Gaillard
// Created on: Fri Oct 31 2008
//

/*
Copyright or © or Copr. CNRS, (October 31, 2008)

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

#ifndef _PHREDPOLY_H_
#define _PHREDPOLY_H_

#include "ISequenceStream.h"
#include "../Sequence.h"

namespace bpp {

  /**
   * @brief The poly sequence file format from phred software.
   *
   * This class read DNA sequence from poly files produced by the phred program
   * from the University of Washington.
   * For now, only read raw sequences and do a basic filter on heterozygous site.
   */
  class PhredPoly: public ISequenceStream {
    protected:
      double ratio_;

    public:

      /**
       * @brief Build a new PhredPoly object.
       */
      PhredPoly(double ratio = 0.8);

      virtual ~PhredPoly() {}

    public:
      /**
       * @name The AbstractISequence interface.
       *
       * @{
       */
      bool nextSequence(std::istream& input, Sequence& seq) const throw (Exception);
      /** @} */

      /**
       * @name The IOSequence interface.
       *
       * @{
       */
      const std::string getDataType() const { return "Sequence"; };
      const std::string getFormatName() const { return "poly file"; };
      const std::string getFormatDescription() const {
        return "Sequences following the poly format as describe in the phred documentation.";
      }
      /** @} */
  };
} //end of namespace bpp

#endif // _PHREDPOLY_H_
