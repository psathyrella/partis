//
// File IODistanceMatrixFactory.h
// Created by: Julien Dutheil
// Created on: Tue 18/04/06 10:24
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

#ifndef _IODISTANCEMATRIXFACTORY_H_
#define _IODISTANCEMATRIXFACTORY_H_

#include "../Distance/DistanceEstimation.h"
#include "IoDistanceMatrix.h"
#include <Bpp/Exceptions.h>

//From the STL:
#include <string>

namespace bpp
{

/**
 * @brief Utilitary class for creating distance matrix readers and writers.
 *
 * @see IOSequenceFactory
 * @see IOTreeFactory
 */
class IODistanceMatrixFactory
{
public:
  static const std::string PHYLIP_FORMAT;  

public:

  /**
   * @brief Creates a new factory object.
   *
   * Example:
   * @code
   * IDistanceMatrix * matReader = IODistanceMatrixFactory().createReader(IODistanceMatrixFactory::PHYLIP);
   * DistanceMatrix * matrix = matReader->read("file.ph");
   * delete matReader;
   * @endcode
   */
  IODistanceMatrixFactory() {}
  virtual ~IODistanceMatrixFactory() {}
  
  /**
   * @brief Get a new dynamically created IDistanceMatrix object.
   *
   * @param format The input file format, and whether names should be
   *      only less than 10 characters, or not (false=10 characters max).
   * @param extended format (default false).
   * @return A pointer toward a new IDistanceMatrix object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual IDistanceMatrix* createReader(const std::string& format, bool extended=false) throw (Exception);
  
  /**
   * @brief Get a new dynamically created ODistanceMatrix object.
   *
   * @param format The output file format, and whether names should be
   *        only less than 10 characters, or not (false=10 characters max).
   * @param extended format (default false).
   * @return A pointer toward a new ODistanceMatrix object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual ODistanceMatrix * createWriter(const std::string& format, bool extended=false) throw (Exception);
};

} //end of namespace bpp.

#endif //_IODISTANCEMATRIXFACTORY_H_

