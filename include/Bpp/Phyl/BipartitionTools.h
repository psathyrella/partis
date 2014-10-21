//
// File: BipartitionTools.h
// Created by: Nicolas Galtier & Julien Dutheil
// Created on: Tue Apr 13 15:09 2007
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

#ifndef _BIPARTITIONTOOLS_H_
#define _BIPARTITIONTOOLS_H_

#include "BipartitionList.h"

// From bpp-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

namespace bpp
{
/**
 * @brief This class provides tools related to the BipartitionList class
 *
 * BipartitionTools includes functions dealing with dynamic arrays of bits
 * and tools for dealing with bipartitions from distinct lists.
 *
 * @see Tree
 * @see BipartitionList
 * @see TreeTools
 */
class BipartitionTools
{
public:
  /**
   * @brief
   *
   * Unit length (in bits) of arrays of bits. Must be a multiple of CHAR_BIT*sizeof(int).
   * Default value is 64.
   */
  static size_t LWORD;

public:
  BipartitionTools() {}
  virtual ~BipartitionTools() {}

public:
  /**
   * @brief Sets bit number num of bit array list to one
   *
   * Note that no control of memory allocation is made
   */
  static void bit1(int* list, int num);

  /**
   * @brief Sets bit number num of bit array plist to zero
   *
   * Note that no control of memory allocation is made
   */
  static void bit0(int* list, int num);

  /**
   * @brief bit-wise logical AND between two arrays of bits
   *
   * (1 AND 1 = 1; 1 AND 0 = 0 AND 1 = 0 AND 0 = 0)
   *
   * Note that no control of memory allocation is made
   *
   * param list1 first array of bit
   * param list2 second array of bit
   * param listet resulting array of bit
   * param len number of int over which the operation is performed
   */
  static void bitAnd(int* listet, int* list1, int* list2, size_t len);

  /**
   * @brief bit-wise logical OR between two arrays of bits
   *
   * (1 OR 1 = 1 OR 0 = 0 OR 1 = 1; 0 OR 0 = 0)
   *
   * Note that no control of memory allocation is made
   *
   * param list1 first array of bit
   * param list2 second array of bit
   * param listou resulting array of bit
   * param len number of int over which the operation is performed
   */
  static void bitOr(int* listou, int* list1, int* list2, size_t len);

  /**
   * @brief bit-wise logical NOT
   *
   * (NOT 1 = 0; NOT 0 = 1)
   *
   * Note that no control of memory allocation is made
   *
   * param list input array of bit
   * param listnot resulting array of bit
   * param len number of int over which the operation is performed
   */
  static void bitNot(int* listnon, int* list, size_t len);

  /**
   * @brief Tells whether bit number num in bit array list is one
   */
  static bool testBit(int* list, int num);

  /**
   * @brief Makes one BipartitionList out of several
   *
   * The input BipartitionList objects must share the same set of elements. This will be checked or not
   * depending on checkElements
   */
  static BipartitionList* mergeBipartitionLists(const std::vector<BipartitionList*>& vecBipartL, bool checkElements = true) throw (Exception);

  /**
   * @brief Construct a BipartitionList containing two bipartitions taken from distinct input lists
   */
  static BipartitionList* buildBipartitionPair(const BipartitionList& bipartL1, size_t i1, const BipartitionList& bipartL2, size_t i2, bool checkElements = true) throw (Exception);

  /**
   * @brief Tells whether two bipartitions from distinct lists are identical
   */
  static bool areIdentical(const BipartitionList& bipart1, size_t i1, const BipartitionList& bipart2, size_t i2, bool checkElements = true);

  /**
   * @brief Tells whether two bipartitions from distinct lists are compatible
   */
  static bool areCompatible(const BipartitionList& bipart1, size_t i1, const BipartitionList& bipart2, size_t i2, bool checkElements = true);

  /**
   * @brief Create a sequence data set corresponding to the Matrix Representation of the input BipartitionList objects
   *
   * The input BipartitionList objects can have distinct sets of elements - missing data will be represented as 'N'.
   * The output alignment (DNA sequences including only A, C and N)) is ready for maximum parsimony analysis
   * according to the MRP supertree method.
   */
  static VectorSiteContainer* MRPEncode(const std::vector<BipartitionList*>& vecBipartL) throw (Exception);
};
} // end of namespace bpp.

#endif // _BIPARTITIONTOOLS_H_

