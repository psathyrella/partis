//
// File: BppOAlphabetIndex1Format.h
// Created by: Julien Dutheil
// Created on: Thursday Februar 07th, 16:30
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

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

#ifndef _BPPOALPHABETINDEX1FORMAT_H_
#define _BPPOALPHABETINDEX1FORMAT_H_

#include <Bpp/Io/IoFormat.h>
#include "../AlphabetIndex/AlphabetIndex1.h"

// From the STL:
#include <string>

namespace bpp
{

  /**
   * @brief AlphabetIndex1 I/O in BppO format.
   *
   * Enables the instanciation of AlphabetIndex1 objects according to
   * the BppO syntax (see the Bio++ Program Suite
   * manual for a detailed description of this syntax).
   *
   */
  class BppOAlphabetIndex1Format:
    public virtual IOFormat
  {
  private:
    const Alphabet* alphabet_;
    std::string message_;
    bool verbose_;

  public:
    /**
      * @param alphabet The alphabet for which indices should be built.
      * The alphabet will be used to check that the instanciated index is compatible.
      * @param message Some text describing what the index is intended for.
      * @param verbose Tell if some messages should be printed while parsing.
      */
    BppOAlphabetIndex1Format(const Alphabet* alphabet, const std::string& message, bool verbose = true):
      alphabet_(alphabet), message_(message), verbose_(verbose) {}

    BppOAlphabetIndex1Format(const BppOAlphabetIndex1Format& format):
      alphabet_(format.alphabet_),
      message_(format.message_),
      verbose_(format.verbose_) {}

    BppOAlphabetIndex1Format& operator=(const BppOAlphabetIndex1Format& format)
    {
      alphabet_ = format.alphabet_;
      message_  = format.message_;
      verbose_  = format.verbose_;
      return *this;
    }

    virtual ~BppOAlphabetIndex1Format() {}

  public:
    const std::string getFormatName() const { return "BppO"; }

    const std::string getFormatDescription() const { return "Bpp Options format."; }

		const std::string getDataType() const { return "AlphabetIndex1"; }

    /**
     * @brief Read a AlphabetIndex1 object from a string.
     *
     * @param description A string describing the index in the keyval syntax.
     * @return A new AlphabetIndex1 object according to options specified.
     * @throw Exception if an error occured.
     */
    AlphabetIndex1* read(const std::string& description) throw (Exception);

  };

} //end of namespace bpp.

#endif //_BPPOALPHABETINDEX1FORMAT_H_

