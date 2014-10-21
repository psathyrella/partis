//
// File IOTreeFactory.h
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

#ifndef _IOTREEFACTORY_H_
#define _IOTREEFACTORY_H_

#include "IoTree.h"

namespace bpp
{

/**
 * @brief Utilitary class for creating tree readers and writers.
 *
 * @see IOSequenceFactory
 * @see IODistanceMatrixFactory
 */
class IOTreeFactory
{
public:
  static const std::string NEWICK_FORMAT;  
  static const std::string NEXUS_FORMAT;  
  static const std::string NHX_FORMAT;  

public:

  /**
   * @brief Creates a new factory object.
   *
   * Example:
   * @code
   * ITree * treeReader = IOTreeFactory().createReader(IOTreeFactory::NEWICK);
   * Tree * tree = treeReader->read("file.dnd");
   * delete treeReader;
   * @endcode
   */
  IOTreeFactory() {}
  virtual ~IOTreeFactory() {}
  
  /**
   * @brief Get a new dynamically created ITree object.
   *
   * @param format The input file format.
   * @return A pointer toward a new ITree object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual ITree* createReader(const std::string& format) throw (Exception);
  
  /**
   * @brief Get a new dynamically created OTree object.
   *
   * @param format The output file format.
   * @return A pointer toward a new OTree object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual OTree* createWriter(const std::string& format) throw (Exception);
};

} //end of namespace bpp.

#endif //_IOTREEFACTORY_H_

