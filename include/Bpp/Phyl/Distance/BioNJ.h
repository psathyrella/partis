//
// File: BioNJ.cpp
// Created by: Vincent Ranwez
// Created on: Tue Apr 11 14:23 2006
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#ifndef _BIONJ_H_
#define _BIONJ_H_

#include "NeighborJoining.h"

namespace bpp
{
/**
 * @brief The BioNJ distance method.
 *
 * Reference:
 * Gascuel O.
 * BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data.
 * Mol Biol Evol. 1997 Jul;14(7):685-95.
 */
class BioNJ :
  public NeighborJoining
{
private:
  DistanceMatrix variance_;
  double lambda_;

public:
  /**
   * @brief Create a new BioNJ object instance and compute a tree from a distance matrix.
   *
   * @param rooted Tell if the output tree should be rooted.
   * @param positiveLengths Tell if negative lengths should be avoided.
   * @param verbose Allow to display extra information, like progress bars.
   */
  BioNJ(bool rooted = false, bool positiveLengths = false, bool verbose = true) :
    NeighborJoining(rooted, positiveLengths, verbose),
    variance_(0),
    lambda_(0) {}

  /**
   * @brief Create a new BioNJ object instance and compute a tree from a distance matrix.
   *
   * @param matrix Input distance matrix.
   * @param rooted Tell if the output tree should be rooted.
   * @param positiveLengths Tell if negative lengths should be avoided.
   * @param verbose Allow to display extra information, like progress bars.
   */
  BioNJ(const DistanceMatrix& matrix, bool rooted = false, bool positiveLengths = false, bool verbose = true) throw (Exception) :
    NeighborJoining(rooted, positiveLengths, verbose),
    // Use the default constructor, because the other one call computeTree.
    variance_(matrix),
    lambda_(0)
  {
    setDistanceMatrix(matrix);
    outputPositiveLengths(positiveLengths);
    computeTree();
  }

  BioNJ* clone() const { return new BioNJ(*this); }

  virtual ~BioNJ() {}

public:
  std::string getName() const { return "BioNJ"; }

  void setDistanceMatrix(const DistanceMatrix& matrix)
  {
    NeighborJoining::setDistanceMatrix(matrix);
    variance_ = matrix;
  }
  void computeTree() throw (Exception);
  double computeDistancesFromPair(const std::vector<size_t>& pair, const std::vector<double>& branchLengths, size_t pos);
};
} // end of namespace bpp.

#endif // _BIONJ_H_

