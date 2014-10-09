//
// File IODiscreteDistributionFactory.h
// Created by: Laurent Guéguen
// Created on: lundi 3 septembre 2012, à 14h 45
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

#ifndef _IODISCRETEDISTRIBUTIONFACTORY_H_
#define _IODISCRETEDISTRIBUTIONFACTORY_H_

#include "../Numeric/Prob/DiscreteDistribution.h"
#include "../Exceptions.h"
#include "IoDiscreteDistribution.h"

//From the STL:
#include <string>

namespace bpp
{

/**
 * @brief Utilitary class for creating discrete distribution readers and
 * writers.
 */
class IoDiscreteDistributionFactory
{
public:
  static const std::string BPPO_FORMAT;  

public:

  /**
   * @brief Creates a new factory object.
   *
   * Example:
   * @code
   * IDiscreteDistribution* distReader = IoDiscreteDistributionFactory().createReader(IoDiscreteDistributionFactory::BPP_FORMAT);
   * DiscreteDistribution* dist = distReader->read(...);
   * delete distReader;
   * @endcode
   */
  IoDiscreteDistributionFactory() {}
  virtual ~IoDiscreteDistributionFactory() {}
  
  /**
   * @brief Get a new dynamically created IDiscreteDistribution object.
   *
   * @param format The input file format.
   * @return A pointer toward a new IDiscreteDistribution object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual IDiscreteDistribution* createReader(const std::string& format) throw (Exception);
  
  /**
   * @brief Get a new dynamically created ODiscreteDistribution object.
   *
   * @param format The output file format.
   * @return A pointer toward a new ODiscreteDistribution object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual ODiscreteDistribution* createWriter(const std::string& format) throw (Exception);
};

} //end of namespace bpp.

#endif //_IODISCRETEDISTRIBUTIONFACTORY_H_

