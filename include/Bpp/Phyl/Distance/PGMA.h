//
// File: PGMA.h
// Created by: Julien Dutheil
// Created on: Mon jul 11 11:41 2005
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

#ifndef _PGMA_H_
#define _PGMA_H_

#include "AbstractAgglomerativeDistanceMethod.h"
#include "../Tree.h"
#include "../TreeTemplate.h"

namespace bpp
{
/**
 * @brief Inner data structure for WPGMA and UPGMA distance methods.
 */
struct PGMAInfos
{
  size_t numberOfLeaves;
  double time;
};

/**
 * @brief Compute WPGMA and UPGMA trees from a distance matrix.
 *
 * WPGMA = Weighted pair group method using arithmetic averaging,
 * is equivalent to the average linkage hierarchical clustering method.
 * The distance between two taxa is the average distance between all individuals in each taxa.
 * The unweighted version (named UPGMA), uses a weighted average, with the number of individuals in a group as a weight.
 */
class PGMA :
  public AbstractAgglomerativeDistanceMethod
{
protected:
  bool weighted_;

public:
  PGMA(bool weighted = true) :
    AbstractAgglomerativeDistanceMethod(true, true),
    weighted_(weighted) {}

  /**
   * @brief Create a (U/W)PGMA object instance.
   *
   * @param matrix Input distance matrix.
   * @param weighted Tell if we should perform Weighted or Unweighted pair group method.
   * @param verbose Allow to display extra information, like progress bars.
   */
  PGMA(const DistanceMatrix& matrix, bool weighted = true, bool verbose = true) throw (Exception) :
    AbstractAgglomerativeDistanceMethod(matrix, verbose, true),
    weighted_(weighted)
  {
    computeTree();
  }
  virtual ~PGMA() {}

  PGMA* clone() const { return new PGMA(*this); }

public:
  std::string getName() const { return std::string(weighted_ ? "W" : "U") + "PGMA"; }

  void setDistanceMatrix(const DistanceMatrix& matrix)
  {
    AbstractAgglomerativeDistanceMethod::setDistanceMatrix(matrix);
  }

  TreeTemplate<Node>* getTree() const;

  void setWeighted(bool weighted) { weighted_ = weighted; }
  bool isWeighted() const { return weighted_; }

protected:
  std::vector<size_t> getBestPair() throw (Exception);
  std::vector<double> computeBranchLengthsForPair(const std::vector<size_t>& pair);
  double computeDistancesFromPair(const std::vector<size_t>& pair, const std::vector<double>& branchLengths, size_t pos);
  void finalStep(int idRoot);
  virtual Node* getLeafNode(int id, const std::string& name);
  virtual Node* getParentNode(int id, Node* son1, Node* son2);
};
} // end of namespace bpp.

#endif // _PGMA_H_

