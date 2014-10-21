//
// File: OSequence.h
// Created by: Guillaume Deuchst
//             Julien Dutheil
// Created on: Tue Aug 21 2003
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

#ifndef _OSEQUENCE_H_
#define _OSEQUENCE_H_

#include "../Container/SequenceContainer.h"
#include "../Container/SiteContainer.h"
#include "IoSequence.h"

#include <Bpp/Exceptions.h>

namespace bpp
{

  /**
   * @brief The OSequence interface.
   * 
   * This interface defines the basic methods for writing sequences to a file.
   */
  class OSequence:
    public virtual IOSequence
  {
    public:
      OSequence() {}
      virtual ~OSequence() {}

    public:

      /**
       * @brief Write a container to a stream.
       *
       * @param output The output stream where to write.
       * @param sc        The container to write.
       * @throw Exception If the file is not in the specified format.
       */
      virtual void writeSequences(std::ostream& output, const SequenceContainer& sc) const throw (Exception) = 0;

      /**
       * @brief Write a container to a file.
       *
       * @param path      The path to the file to write.
       * @param sc        The container to write.
       * @param overwrite If true the sequences are written at the beginning of the file instead of being appended.
       *                  Any previous content will be lost.
       * @throw Exception If the file is not in the specified format.
       */
      virtual void writeSequences(const std::string& path, const SequenceContainer & sc, bool overwrite) const throw (Exception) = 0;

  };


  /**
   * @brief The OAlignment interface.
   * 
   * This interface defines the basic methods for writing alignments to a file.
   */
  class OAlignment:
    public virtual IOSequence
  {
    public:
      OAlignment() {}
      virtual ~OAlignment() {}

    public:

      /**
       * @brief Write a container to a stream.
       *
       * @param output The output stream where to write.
       * @param sc        The container to write.
       * @throw Exception If the file is not in the specified format.
       */
      virtual void writeAlignment(std::ostream& output, const SiteContainer& sc) const throw (Exception) = 0;

      /**
       * @brief Write a container to a file.
       *
       * @param path      The path to the file to write.
       * @param sc        The container to write.
       * @param overwrite If true the sequences are written at the beginning of the file instead of being appended.
       *                  Any previous content will be lost.
       * @throw Exception If the file is not in the specified format.
       */
      virtual void writeAlignment(const std::string& path, const SiteContainer& sc, bool overwrite) const throw (Exception) = 0;

  };
} //end of namespace bpp.

#endif	// _OSEQUENCE_H_

