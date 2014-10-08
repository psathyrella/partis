//
// File: ISequence.h
// Created by: Guillaume Deuchst
//             Julien Dutheil
// Created on: Wed Jul 30 2003
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

#ifndef _ISEQUENCE_H_
#define _ISEQUENCE_H_

#include "IoSequence.h"
#include "../Sequence.h"
#include "../Container/SequenceContainer.h"
#include "../Container/SiteContainer.h"
#include <Bpp/Exceptions.h>

//From the STL:
#include <iostream>
#include <string>

namespace bpp
{

/**
 * @brief The ISequence interface.
 *
 * This interface defines the basic methods for reading sequences from a file.
 * NB: This interface is effective only if the VIRTUAL_COV option is enabled (default behavior).
 */
class ISequence :
  public virtual IOSequence
{
  public:
    ISequence() {}
    virtual ~ISequence() {}

  public:
  
    /**
     * @brief Create a new container from a stream.
     *
     * @param input  The input stream to read.
     * @param alpha The alphabet to be associated to the container.
     * @return A new SequenceContainer object.
     * @throw Exception If the file is not in the specified format.
     */
    virtual SequenceContainer* readSequences(std::istream& input, const Alphabet* alpha) const throw (Exception) = 0;
    /**
     * @brief Create a new container from a file.
     *
     * @param path  The path to the file to read.
     * @param alpha The alphabet to be associated to the container.
     * @return A new SequenceContainer object.
     * @throw Exception If the file is not in the specified format.
     */
    virtual SequenceContainer* readSequences(const std::string& path, const Alphabet* alpha) const throw (Exception) = 0;

};

/**
 * @brief The IAlignment interface.
 *
 * This interface defines the basic methods for reading aligned sequences from a file.
 */
class IAlignment:
  public virtual IOSequence
{
  public:
    IAlignment() {}
    virtual ~IAlignment() {}

  public:
  
    /**
     * @brief Create a new container from a stream.
     *
     * @param input  The input stream to read.
     * @param alpha The alphabet to be associated to the container.
     * @return A new SiteContainer object.
     * @throw Exception If the file is not in the specified format.
     */
    virtual SiteContainer* readAlignment(std::istream& input, const Alphabet* alpha) const throw (Exception) = 0;
    /**
     * @brief Create a new container from a file.
     *
     * @param path  The path to the file to read.
     * @param alpha The alphabet to be associated to the container.
     * @return A new SiteContainer object.
     * @throw Exception If the file is not in the specified format.
     */
    virtual SiteContainer* readAlignment(const std::string& path, const Alphabet* alpha) const throw (Exception) = 0;

};

} //end of namespace bpp.

#endif  // _ISEQUENCE_H_

