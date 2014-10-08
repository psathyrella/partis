//
// File: DRTreeParsimonyData.h
// Created by: Julien Dutheil
// Created on: Tue Jan 09 17:15 2007
// From file DRTreeParsimonyScore.h
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

#ifndef _DRTREEPARSIMONYDATA_H_
#define _DRTREEPARSIMONYDATA_H_

#include "AbstractTreeParsimonyData.h"
#include "../Model/StateMap.h"

// From SeqLib
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <bitset>

namespace bpp
{
typedef std::bitset<21> Bitset; // 20AA + gaps, codon not lalowed so far :s

/**
 * @brief Parsimony data structure for a node.
 *
 * This class is for use with the DRTreeParsimonyData class.
 *
 * Store for each neighbor node
 * - a vector of bitsets,
 * - a vector of score for the corresponding subtree.
 *
 * @see DRTreeParsimonyData
 */
class DRTreeParsimonyNodeData :
  public TreeParsimonyNodeData
{
private:
  mutable std::map<int, std::vector<Bitset> > nodeBitsets_;
  mutable std::map<int, std::vector<unsigned int> > nodeScores_;
  const Node* node_;

public:
  DRTreeParsimonyNodeData() :
    nodeBitsets_(),
    nodeScores_(),
    node_(0)
  {}

  DRTreeParsimonyNodeData(const DRTreeParsimonyNodeData& tpnd) :
    nodeBitsets_(tpnd.nodeBitsets_),
    nodeScores_(tpnd.nodeScores_),
    node_(tpnd.node_)
  {}

  DRTreeParsimonyNodeData& operator=(const DRTreeParsimonyNodeData& tpnd)
  {
    nodeBitsets_ = tpnd.nodeBitsets_;
    nodeScores_  = tpnd.nodeScores_;
    node_        = tpnd.node_;
    return *this;
  }

  DRTreeParsimonyNodeData* clone() const { return new DRTreeParsimonyNodeData(*this); }

public:
  const Node* getNode() const { return node_; }

  void setNode(const Node* node) { node_ = node; }

  std::vector<Bitset>& getBitsetsArrayForNeighbor(int neighborId)
  {
    return nodeBitsets_[neighborId];
  }
  const std::vector<Bitset>& getBitsetsArrayForNeighbor(int neighborId) const
  {
    return nodeBitsets_[neighborId];
  }
  std::vector<unsigned int>& getScoresArrayForNeighbor(int neighborId)
  {
    return nodeScores_[neighborId];
  }
  const std::vector<unsigned int>& getScoresArrayForNeighbor(int neighborId) const
  {
    return nodeScores_[neighborId];
  }

  bool isNeighbor(int neighborId) const
  {
    return nodeBitsets_.find(neighborId) != nodeBitsets_.end();
  }

  void eraseNeighborArrays()
  {
    nodeBitsets_.erase(nodeBitsets_.begin(), nodeBitsets_.end());
    nodeScores_.erase(nodeScores_.begin(), nodeScores_.end());
  }
};

/**
 * @brief Parsimony data structure for a leaf.
 *
 * This class is for use with the DRTreeParsimonyData class.
 *
 * Store the vector of bitsets associated to a leaf.
 *
 * @see DRTreeParsimonyData
 */
class DRTreeParsimonyLeafData :
  public TreeParsimonyNodeData
{
private:
  mutable std::vector<Bitset> leafBitsets_;
  const Node* leaf_;

public:
  DRTreeParsimonyLeafData() :
    leafBitsets_(),
    leaf_(0)
  {}

  DRTreeParsimonyLeafData(const DRTreeParsimonyLeafData& tpld) :
    leafBitsets_(tpld.leafBitsets_),
    leaf_(tpld.leaf_)
  {}

  DRTreeParsimonyLeafData& operator=(const DRTreeParsimonyLeafData& tpld)
  {
    leafBitsets_ = tpld.leafBitsets_;
    leaf_        = tpld.leaf_;
    return *this;
  }


  DRTreeParsimonyLeafData* clone() const { return new DRTreeParsimonyLeafData(*this); }

public:
  const Node* getNode() const { return leaf_; }
  void setNode(const Node* node) { leaf_ = node; }

  std::vector<Bitset>& getBitsetsArray()
  {
    return leafBitsets_;
  }
  const std::vector<Bitset>& getBitsetsArray() const
  {
    return leafBitsets_;
  }
};

/**
 * @brief Parsimony data structure for double-recursive (DR) algorithm.
 *
 * States are coded using bitsets for faster computing (@see AbstractTreeParsimonyData).
 * For each inner node in the tree, we store a DRTreeParsimonyNodeData object in nodeData_.
 * For each leaf node in the tree, we store a DRTreeParsimonyLeafData object in leafData_.
 *
 * The dataset is first compressed, removing all identical sites.
 * The resulting dataset is stored in shrunkData_.
 * The corresponding positions are stored in rootPatternLinks_, inherited from AbstractTreeParsimonyData.
 */
class DRTreeParsimonyData :
  public AbstractTreeParsimonyData
{
private:
  mutable std::map<int, DRTreeParsimonyNodeData> nodeData_;
  mutable std::map<int, DRTreeParsimonyLeafData> leafData_;
  mutable std::vector<Bitset> rootBitsets_;
  mutable std::vector<unsigned int> rootScores_;
  SiteContainer* shrunkData_;
  size_t nbSites_;
  size_t nbStates_;
  size_t nbDistinctSites_;

public:
  DRTreeParsimonyData(const TreeTemplate<Node>* tree) :
    AbstractTreeParsimonyData(tree),
    nodeData_(),
    leafData_(),
    rootBitsets_(),
    rootScores_(),
    shrunkData_(0),
    nbSites_(0),
    nbStates_(0),
    nbDistinctSites_(0)
  {}

  DRTreeParsimonyData(const DRTreeParsimonyData& data);

  DRTreeParsimonyData& operator=(const DRTreeParsimonyData& data);

  virtual ~DRTreeParsimonyData() { delete shrunkData_; }

  DRTreeParsimonyData* clone() const { return new DRTreeParsimonyData(*this); }

public:
  /**
   * @brief Set the tree associated to the data.
   *
   * All node data will be actualized accordingly by calling the setNode() method on the corresponding nodes.
   * @warning: the old tree and the new tree must be two clones! And particularly, they have to share the
   * same topology and nodes id.
   *
   * @param tree The tree to be associated to this data.
   */
  void setTree(const TreeTemplate<Node>* tree)
  {
    AbstractTreeParsimonyData::setTreeP_(tree);
    for (std::map<int, DRTreeParsimonyNodeData>::iterator it = nodeData_.begin(); it != nodeData_.end(); it++)
    {
      int id = it->second.getNode()->getId();
      it->second.setNode(tree_->getNode(id));
    }
    for (std::map<int, DRTreeParsimonyLeafData>::iterator it = leafData_.begin(); it != leafData_.end(); it++)
    {
      int id = it->second.getNode()->getId();
      it->second.setNode(tree_->getNode(id));
    }
  }

  DRTreeParsimonyNodeData& getNodeData(int nodeId)
  {
    return nodeData_[nodeId];
  }
  const DRTreeParsimonyNodeData& getNodeData(int nodeId) const
  {
    return nodeData_[nodeId];
  }

  DRTreeParsimonyLeafData& getLeafData(int nodeId)
  {
    return leafData_[nodeId];
  }
  const DRTreeParsimonyLeafData& getLeafData(int nodeId) const
  {
    return leafData_[nodeId];
  }

  std::vector<Bitset>& getBitsetsArray(int nodeId, int neighborId)
  {
    return nodeData_[nodeId].getBitsetsArrayForNeighbor(neighborId);
  }
  const std::vector<Bitset>& getBitsetsArray(int nodeId, int neighborId) const
  {
    return nodeData_[nodeId].getBitsetsArrayForNeighbor(neighborId);
  }

  std::vector<unsigned int>& getScoresArray(int nodeId, int neighborId)
  {
    return nodeData_[nodeId].getScoresArrayForNeighbor(neighborId);
  }
  const std::vector<unsigned int>& getScoresArray(int nodeId, int neighborId) const
  {
    return nodeData_[nodeId].getScoresArrayForNeighbor(neighborId);
  }

  size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const
  {
    return currentPosition;
  }

  std::vector<Bitset>& getRootBitsets() { return rootBitsets_; }
  const std::vector<Bitset>& getRootBitsets() const { return rootBitsets_; }
  const Bitset& getRootBitset(size_t i) const { return rootBitsets_[i]; }

  std::vector<unsigned int>& getRootScores() { return rootScores_; }
  const std::vector<unsigned int>& getRootScores() const { return rootScores_; }
  unsigned int getRootScore(size_t i) const { return rootScores_[i]; }

  size_t getNumberOfDistinctSites() const { return nbDistinctSites_; }
  size_t getNumberOfSites() const { return nbSites_; }
  size_t getNumberOfStates() const { return nbStates_; }

  void init(const SiteContainer& sites, const StateMap& stateMap) throw (Exception);
  void reInit() throw (Exception);

protected:
  void init(const Node* node, const SiteContainer& sites, const StateMap& stateMap) throw (Exception);
  void reInit(const Node* node) throw (Exception);
};
} // end of namespace bpp.

#endif // _DRTREEPARSIMONYDATA_H_

