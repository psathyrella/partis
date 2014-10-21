//
// File: HierarchicalClustering.h
// From file Cluster.h in CoMap package.
// Created by: Julien Dutheil
// Created on: Tue Aug 30 17:19 2005
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

   This software is a computer program whose purpose is to map substitutions
   on a tree and to detect co-evolving positions in a dataset.

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

#ifndef _HIERARCHICALCLUSTERING_H_
#define _HIERARCHICALCLUSTERING_H_

#include "AbstractAgglomerativeDistanceMethod.h"

namespace bpp
{
class ClusterInfos
{
public:
  size_t numberOfLeaves;
  double length;

public:
  ClusterInfos() : numberOfLeaves(0),
    length(0) {}
};

/**
 * @brief Hierarchical clustering.
 *
 * This class implements the complete, single, average (= UPGMA), median, ward and centroid linkage methods.
 */
class HierarchicalClustering :
  public AbstractAgglomerativeDistanceMethod
{
public:
  static const std::string COMPLETE;
  static const std::string SINGLE;
  static const std::string AVERAGE;
  static const std::string MEDIAN;
  static const std::string WARD;
  static const std::string CENTROID;

protected:
  std::string method_;

public:
  /**
   * @brief Builds a new clustering object.
   *
   * @param method The linkage method to use. should be one of COMPLETE, SINGLE, AVERAGE, MEDIAN, WARD, CENTROID.
   * @param verbose Tell if some progress information should be displayed.
   */
  HierarchicalClustering(const std::string& method, bool verbose = false) :
    AbstractAgglomerativeDistanceMethod(verbose),
    method_(method) {}
  HierarchicalClustering(const std::string& method, const DistanceMatrix& matrix, bool verbose = false) throw (Exception) :
    AbstractAgglomerativeDistanceMethod(matrix, verbose, true),
    method_(method)
  {
    computeTree();
  }

  virtual ~HierarchicalClustering() {}

  HierarchicalClustering* clone() const { return new HierarchicalClustering(*this); }

public:
  std::string getName() const { return "Hierarchical clustering: " + method_; }

  TreeTemplate<Node>* getTree() const;

protected:
  std::vector<size_t> getBestPair() throw (Exception);
  std::vector<double> computeBranchLengthsForPair(const std::vector<size_t>& pair);
  double computeDistancesFromPair(const std::vector<size_t>& pair, const std::vector<double>& branchLengths, size_t pos);
  void finalStep(int idRoot);
  virtual Node* getLeafNode(int id, const std::string& name);
  virtual Node* getParentNode(int id, Node* son1, Node* son2);
};
} // end of namespace bpp.

#endif // _HIERARCHICALCLUSTERING_H_

