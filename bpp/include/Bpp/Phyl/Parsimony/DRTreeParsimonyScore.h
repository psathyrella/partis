//
// File: DRTreeParsimonyScore.h
// Created by: Julien Dutheil
// Created on: Thu Jul 28 18:31 2005
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

#ifndef _DRTREEPARSIMONYSCORE_H_
#define _DRTREEPARSIMONYSCORE_H_

#include "AbstractTreeParsimonyScore.h"
#include "DRTreeParsimonyData.h"
#include "../NNISearchable.h"
#include "../TreeTools.h"

namespace bpp
{
/**
 * @brief Double recursive implementation of interface TreeParsimonyScore.
 *
 * Uses a DRTreeParsimonyData object for data storage.
 */
class DRTreeParsimonyScore :
  public AbstractTreeParsimonyScore,
  public virtual NNISearchable
{
private:
  DRTreeParsimonyData* parsimonyData_;
  size_t nbDistinctSites_;

public:
  DRTreeParsimonyScore(
    const Tree& tree,
    const SiteContainer& data,
    bool verbose = true,
    bool includeGaps = false)
  throw (Exception);

  DRTreeParsimonyScore(
    const Tree& tree,
    const SiteContainer& data,
    const StateMap* statesMap,
    bool verbose = true)
  throw (Exception);

  DRTreeParsimonyScore(const DRTreeParsimonyScore& tp);

  DRTreeParsimonyScore& operator=(const DRTreeParsimonyScore& tp);

  virtual ~DRTreeParsimonyScore();

#ifndef NO_VIRTUAL_COV
  DRTreeParsimonyScore*
#else
  Clonable*
#endif
  clone() const { return new DRTreeParsimonyScore(*this); }

private:
  void init_(const SiteContainer& data, bool verbose);

protected:
  /**
   * @brief Compute all scores.
   *
   * Call the computeScoresPreorder and computeScoresPostorder methods, and then initialize rootBitsets_ and rootScores_.
   */
  virtual void computeScores();
  /**
   * @brief Compute scores (preorder algorithm).
   */
  virtual void computeScoresPreorder(const Node*);
  /**
   * @brief Compute scores (postorder algorithm).
   */
  virtual void computeScoresPostorder(const Node*);

public:
  unsigned int getScore() const;
  unsigned int getScoreForSite(size_t site) const;

  /**
   * @brief Compute bitsets and scores for each site for a node, in postorder.
   *
   * @param pData    The node data to use.
   * @param rBitsets The bitset array where to store the resulting bitsets.
   * @param rScores  The score array where to write the resulting scores.
   */
  static void computeScoresPostorderForNode(
    const DRTreeParsimonyNodeData& pData,
    std::vector<Bitset>& rBitsets,
    std::vector<unsigned int>& rScores);

  /**
   * @brief Compute bitsets and scores for each site for a node, in preorder.
   *
   * @param pData    The node data to use.
   * @param source   The node where we are coming from.
   * @param rBitsets The bitset array where to store the resulting bitsets.
   * @param rScores  The score array where to write the resulting scores.
   */
  static void computeScoresPreorderForNode(
    const DRTreeParsimonyNodeData& pData,
    const Node* source,
    std::vector<Bitset>& rBitsets,
    std::vector<unsigned int>& rScores);

  /**
   * @brief Compute bitsets and scores for each site for a node, in all directions.
   *
   * @param pData    The node data to use.
   * @param rBitsets The bitset array where to store the resulting bitsets.
   * @param rScores  The score array where to write the resulting scores.
   */
  static void computeScoresForNode(
    const DRTreeParsimonyNodeData& pData, std::vector<Bitset>& rBitsets,
    std::vector<unsigned int>& rScores);

  /**
   * @brief Compute bitsets and scores from an array of arrays.
   *
   * This method is the more general score computation.
   * Depending on what is passed as input, it may computes scroes fo a subtree
   * or the whole tree.
   *
   * @param iBitsets The vector of bitset arrays to use.
   * @param iScores  The vector of score arrays to use.
   * @param oBitsets The bitset array where to store the resulting bitsets.
   * @param oScores  The score array where to write the resulting scores.
   */
  static void computeScoresFromArrays(
    const std::vector<const std::vector<Bitset>*>& iBitsets,
    const std::vector<const std::vector<unsigned int>*>& iScores,
    std::vector<Bitset>& oBitsets,
    std::vector<unsigned int>& oScores);

  /**
   * @name Thee NNISearchable interface.
   *
   * @{
   */
  double getTopologyValue() const throw (Exception) { return getScore(); }

  double testNNI(int nodeId) const throw (NodeException);

  void doNNI(int nodeId) throw (NodeException);

  // Tree& getTopology() { return getTree(); } do we realy need this one?
  const Tree& getTopology() const { return getTree(); }

  void topologyChangeTested(const TopologyChangeEvent& event)
  {
    parsimonyData_->reInit();
    computeScores();
  }

  void topologyChangeSuccessful(const TopologyChangeEvent& event) {}
  /**@} */
};
} // end of namespace bpp.

#endif // _DRTREEPARSIMONYSCORE_H_

