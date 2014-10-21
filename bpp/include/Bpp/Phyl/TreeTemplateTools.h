//
// File: TreeTemplateTools.h
// Created by:  Julien Dutheil
// Created on: Fri Oct  13 13:00 2006
// From file TreeTools.h
// Created on: Wed Aug  6 13:45:28 2003
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

#ifndef _TREETEMPLATETOOLS_H_
#define _TREETEMPLATETOOLS_H_

#include "TreeTools.h"
#include <Bpp/Numeric/Random/RandomTools.h>

// From the STL:
#include <string>
#include <vector>

namespace bpp
{
template<class N> class TreeTemplate;


/**
 * @brief Utilitary methods working with TreeTemplate and Node objects.
 *
 * @see TreeTools for more generic methods.
 */
class TreeTemplateTools
{
public:
  TreeTemplateTools() {}
  virtual ~TreeTemplateTools() {}

public:
  /**
   * @name Retrieve topology information
   *
   * @{
   */

  /**
   * @brief Retrieve all leaves from a subtree.
   *
   * @param node The node that defines the subtree.
   * @return A vector of pointers toward each leaf in the subtree.
   */
  template<class N>
  static std::vector<N*> getLeaves(N& node)
  {
    std::vector<N*> leaves;
    getLeaves<N>(node, leaves);
    return leaves;
  }

  /**
   * @brief Retrieve all leaves from a subtree.
   *
   * @param node The node that defines the subtree.
   * @param leaves A vector of pointers toward each leaf in the subtree.
   */
  template<class N>
  static void getLeaves(N& node, std::vector<N*>& leaves)
  {
    if (node.isLeaf())
    {
      leaves.push_back(&node);
    }
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      getLeaves<N>(*node.getSon(i), leaves);
    }
  }

  /**
   * @brief Retrieve all leaves ids from a subtree.
   *
   * @param node The node that defines the subtree.
   * @return A vector of ids.
   */
  static std::vector<int> getLeavesId(const Node& node)
  {
    std::vector<int> ids;
    getLeavesId(node, ids);
    return ids;
  }

  /**
   * @brief Retrieve all leaves ids from a subtree.
   *
   * @param node The node that defines the subtree.
   * @param ids A vector of ids.
   */
  static void getLeavesId(const Node& node, std::vector<int>& ids)
  {
    if (node.isLeaf())
    {
      ids.push_back(node.getId());
    }
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      getLeavesId(*node.getSon(i), ids);
    }
  }

  /**
   * @brief Retrieve all nodes ids that are ancestors of a node.
   *
   * @param node The node
   * @return A vector of ids.
   */
  static std::vector<int> getAncestorsId(const Node& node)
  {
    std::vector<int> ids;
    const Node* n = &node;
    while (n->hasFather())
    {
      n = n->getFather();
      ids.push_back(n->getId());
    }
    return ids;
  }

  /**
   * @brief Get the id of a leaf given its name in a subtree.
   *
   * @param node The node defining the subtree to search.
   * @param name The name of the node.
   * @return The id of the node.
   * @throw NodeNotFoundException If the node is not found.
   */
  static int getLeafId(const Node& node, const std::string& name) throw (NodeNotFoundException)
  {
    int* id = 0;
    searchLeaf(node, name, id);
    if (id == 0) throw NodeNotFoundException("TreeTemplateTools::getLeafId().", name);
    else
    {
      int i = *id;
      delete id;
      return i;
    }
  }

  /**
   * @brief Get the id of a leaf given its name in a subtree.
   *
   * @param node The node defining the subtree to search.
   * @param name The name of the node.
   * @param id The id of the node.
   * @throw NodeNotFoundException If the node is not found.
   */
  static void searchLeaf(const Node& node, const std::string& name, int*& id) throw (NodeNotFoundException)
  {
    if (node.isLeaf())
    {
      if (node.getName() == name)
      {
        id = new int(node.getId());
        return;
      }
    }
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      searchLeaf(*node.getSon(i), name, id);
    }
  }

  /**
   * @brief Remove a leaf node and its parent node, while correcting for branch lengths.
   *
   * @param tree The tree to edit.
   * @param leafName The name of the leaf node.
   * @throw NodeNotFoundException If the node is not found.
   */
  template<class N>
  static void dropLeaf(TreeTemplate<N>& tree, const std::string& leafName) throw (NodeNotFoundException, Exception)
  {
    N* leaf = tree.getNode(leafName);
    if (!leaf->hasFather())
      throw Exception("TreeTemplateTools::dropLeaf(). Leaf is the only node in the tree, can't remove it.");
    N* parent = leaf->getFather();
    if (parent->getNumberOfSons() > 2)
    {
      // The easy case:
      parent->removeSon(leaf);
      delete leaf;
    }
    else if (parent->getNumberOfSons() == 2)
    {
      // We have to delete the parent node as well:
      N* brother = parent->getSon(0);
      if (brother == leaf) brother = parent->getSon(1);
      if (!parent->hasFather())
      {
        // The brother becomes the root:
        if (leaf->hasDistanceToFather() && brother->hasDistanceToFather())
        {
          brother->setDistanceToFather(brother->getDistanceToFather() + leaf->getDistanceToFather());
        }
        brother->removeFather();
        tree.setRootNode(brother);
        delete parent;
        delete leaf;
      }
      else
      {
        N* gParent = parent->getFather();
        if (brother->hasDistanceToFather() && parent->hasDistanceToFather())
        {
          brother->setDistanceToFather(brother->getDistanceToFather() + parent->getDistanceToFather());
        }
        size_t pos = gParent->getSonPosition(parent);
        gParent->setSon(pos, brother);
        delete parent;
        delete leaf;
      }
    }
    else
    {
      // Dunno what to do in that case :(
      throw Exception("TreeTemplateTools::dropLeaf. Parent node as only one child, I don't know what to do in that case :(");
    }
  }

  /**
   * @brief Remove a subtree defined by its root node and its parent node, while correcting for branch lengths.
   *
   * @param tree The tree to edit.
   * @param subtree The subtree to remove, defined by its root node.
   * @throw Exception If something unexpected happens :s
   */
  template<class N>
  static void dropSubtree(TreeTemplate<N>& tree, Node* subtree) throw (Exception)
  {
    if (!subtree->hasFather())
      throw Exception("TreeTemplateTools::dropSubtree(). Trying to remove the full tree!");
    N* parent = subtree->getFather();
    if (parent->getNumberOfSons() > 2)
    {
      // The easy case:
      parent->removeSon(subtree);
      deleteSubtree(subtree);
    }
    else if (parent->getNumberOfSons() == 2)
    {
      // We have to delete the parent node as well:
      N* brother = parent->getSon(0);
      if (brother == subtree) brother = parent->getSon(1);
      if (!parent->hasFather())
      {
        // The brother becomes the root:
        if (subtree->hasDistanceToFather() && brother->hasDistanceToFather())
        {
          brother->setDistanceToFather(brother->getDistanceToFather() + subtree->getDistanceToFather());
        }
        tree.setRootNode(brother);
        delete parent;
        deleteSubtree(subtree);
      }
      else
      {
        N* gParent = parent->getFather();
        if (brother->hasDistanceToFather() && parent->hasDistanceToFather())
        {
          brother->setDistanceToFather(brother->getDistanceToFather() + parent->getDistanceToFather());
        }
        size_t pos = gParent->getSonPosition(parent);
        gParent->setSon(pos, brother);
        delete parent;
        deleteSubtree(subtree);
      }
    }
    else
    {
      // Dunno what to do in that case :(
      throw Exception("TreeTemplateTools::dropSubtree. Parent node as only one child, I don't know what to do in that case :(");
    }
  }

  /**
   * @brief Sample a subtree by removing leaves randomly.
   *
   * @param tree The tree to edit.
   * @param leaves The leafs names that should be sampled. They must be found in the tree otherwise an exception will be thrown.
   * @param size The number of leaves in the final sample. If greater or equal to the number of leaf names, the function returns without doing anything.
   */
  template<class N>
  static void sampleSubtree(TreeTemplate<N>& tree, const std::vector<std::string>& leaves, size_t size)
  {
    std::vector<std::string> names = leaves;
    for (size_t n = names.size(); n > size; --n)
    {
      size_t i = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(n);
      dropLeaf(tree, names[i]);
      names.erase(names.begin() + i);
    }
  }

  /**
   * @brief Retrieve all son nodes from a subtree.
   *
   * @param node The node that defines the subtree.
   * @return A vector of pointers toward each son node in the subtree.
   */
  template<class N>
  static std::vector<N*> getNodes(N& node)
  {
    std::vector<N*> nodes;
    getNodes<N>(node, nodes);
    return nodes;
  }

  /**
   * @brief Retrieve all son nodes from a subtree.
   *
   * @param node The node that defines the subtree.
   * @param nodes A vector of pointers toward each son node in the subtree.
   */
  template<class N>
  static void getNodes(N& node, std::vector<N*>& nodes)
  {
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      getNodes<N>(*node.getSon(i), nodes);
    }
    nodes.push_back(&node);
  }

  /**
   * @brief Retrieve all nodes ids from a subtree.
   *
   * @param node The node that defines the subtree.
   * @return A vector of ids.
   */
  static std::vector<int> getNodesId(const Node& node)
  {
    std::vector<int> ids;
    getNodesId(node, ids);
    return ids;
  }

  /**
   * @brief Retrieve all branches ids from a subtree.
   *
   * @param node The node that defines the subtree.
   * @return A vector of ids.
   */
  static std::vector<int> getBranchesId(const Node& node)
  {
    std::vector<int> ids;
    getBranchesId(node, ids);
    return ids;
  }

  /**
   * @brief Retrieve all nodes ids from a subtree.
   *
   * @param node The node that defines the subtree.
   * @param ids A vector of ids.
   */
  static void getNodesId(const Node& node, std::vector<int>& ids)
  {
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      getNodesId(*node.getSon(i), ids);
    }
    ids.push_back(node.getId());
  }

  /**
   * @brief Retrieve all branches ids from a subtree.
   *
   * @param node The node that defines the subtree.
   * @param ids A vector of ids.
   */
  static void getBranchesId(const Node& node, std::vector<int>& ids)
  {
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      getNodesId(*node.getSon(i), ids);
    }
  }

  /**
   * @brief Retrieve all inner nodes from a subtree.
   *
   * @param node The node that defines the subtree.
   * @return A vector of pointers toward each inner node in the subtree.
   */
  template<class N>
  static std::vector<N*> getInnerNodes(N& node)
  {
    std::vector<N*> nodes;
    getInnerNodes<N>(node, nodes);
    return nodes;
  }

  /**
   * @brief Retrieve all inner nodes from a subtree.
   *
   * A inner node is a node with degree > 1, that is, all nodes but the leaves, be they terminal or not.
   *
   * @param node The node that defines the subtree.
   * @param nodes A vector to be filled with pointers toward each inner node in the subtree.
   */
  template<class N>
  static void getInnerNodes(N& node, std::vector<N*>& nodes)
  {
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      getInnerNodes<N>(*node.getSon(i), nodes);
    }
    if (!node.isLeaf())
      nodes.push_back(&node);  // Do not add leaves!
  }

  /**
   * @brief Retrieve all inner nodes ids from a subtree.
   *
   * A inner node is a node with degree > 1, that is, all nodes but the leaves, be they terminal or not.
   *
   * @param node The node that defines the subtree.
   * @return A vector of ids.
   */
  static std::vector<int> getInnerNodesId(const Node& node)
  {
    std::vector<int> ids;
    getInnerNodesId(node, ids);
    return ids;
  }

  /**
   * @brief Retrieve all inner nodes ids from a subtree.
   *
   * @param node The node that defines the subtree.
   * @param ids  A vector to be filled with the resulting ids.
   */
  static void getInnerNodesId(const Node& node, std::vector<int>& ids)
  {
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      getInnerNodesId(*node.getSon(i), ids);
    }
    if (!node.isLeaf())
      ids.push_back(node.getId());  // Do not add leaves!
  }

  /**
   * @param node The node defining the subtree to be searched.
   * @param id   The id to search for.
   * @return     Nodes with the specified id.
   */
  template<class N>
  static std::vector<N*> searchNodeWithId(N& node, int id)
  {
    std::vector<N*> nodes;
    searchNodeWithId<N>(node, id, nodes);
    return nodes;
  }

  /**
   * @param node  The node defining the subtree to be searched.
   * @param id    The id to search for.
   * @param nodes A vector to be filled with the matching nodes.
   */
  template<class N>
  static void searchNodeWithId(N& node, int id, std::vector<N*>& nodes)
  {
    for (size_t i = 0; i < node.getNumberOfSons(); ++i)
    {
      searchNodeWithId<N>(*node.getSon(i), id, nodes);
    }
    if (node.getId() == id) nodes.push_back(&node);
  }

  /**
   * @param node  The node defining the subtree to be searched.
   * @param id    The id to search for.
   * @return The first node encountered with the given id, or 0 if no node with the given id is found.
   */
  static Node* searchFirstNodeWithId(Node& node, int id)
  {
    if (node.getId() == id)
      return &node;
    else
    {
      for (size_t i = 0; i < node.getNumberOfSons(); ++i)
      {
        Node* result = searchFirstNodeWithId(*node.getSon(i), id);
        if (result)
          return result;
      }
    }
    return 0;
  }

  /**
   * @param node  The node defining the subtree to be searched.
   * @param id    The id to search for.
   * @return The first node encountered with the given id, or 0 if no node with the given id is found.
   */
  static const Node* searchFirstNodeWithId(const Node& node, int id)
  {
    if (node.getId() == id)
      return &node;
    else
    {
      for (size_t i = 0; i < node.getNumberOfSons(); ++i)
      {
        const Node* result = searchFirstNodeWithId(*node.getSon(i), id);
        if (result)
          return result;
      }
    }
    return 0;
  }

  /**
   * @param node The node defining the subtree to be searched.
   * @param id   The id to search for.
   * @return     True if the subtree contains a node with the specified id.
   */
  template<class N>
  static bool hasNodeWithId(const N& node, int id)
  {
    if (node.getId() == id) return true;
    else
    {
      for (size_t i = 0; i < node.getNumberOfSons(); i++)
      {
        if (hasNodeWithId(*node.getSon(i), id)) return true;
      }
      return false;
    }
  }

  /**
   * @param node The node defining the subtree to be searched.
   * @param name The name to search for.
   * @return     Nodes with the specified name.
   */
  template<class N>
  static std::vector<N*> searchNodeWithName(N& node, const std::string& name)
  {
    std::vector<N*> nodes;
    searchNodeWithId<N>(node, name, nodes);
    return nodes;
  }

  /**
   * @param node  The node defining the subtree to be searched.
   * @param name  The name to search for.
   * @param nodes A vector to be filled with the matching nodes.
   */
  template<class N>
  static void searchNodeWithName(N& node, const std::string& name, std::vector<N*>& nodes)
  {
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      searchNodeWithName<N>(*node.getSon(i), name, nodes);
    }
    if (node.hasName() && node.getName() == name) nodes.push_back(&node);
  }

  /**
   * @param node The node defining the subtree to be searched.
   * @param name The name to search for.
   * @return     True if the subtree contains a node with the specified name.
   */
  template<class N>
  static bool hasNodeWithName(const N& node, const std::string& name)
  {
    if (node.hasName() & node.getName() == name) return true;
    else
    {
      for (size_t i = 0; i < node.getNumberOfSons(); i++)
      {
        if (hasNodeWithName(*node.getSon(i), name)) return true;
      }
      return false;
    }
  }

  /**
   * @brief Tell if a particular node is the root of a tree
   * i.e. if it has a father node.
   *
   * @param node The node to check.
   * @return True if node does not have a father.
   */
  static bool isRoot(const Node& node) { return !node.hasFather(); }

  /**
   * @brief Get the number of leaves of a subtree defined by a particular node.
   *
   * @param node The node defining the subtree to check.
   * @return The number of leaves.
   */
  static unsigned int getNumberOfLeaves(const Node& node);

  /**
   * @brief Get the number of nodes of a subtree defined by a particular node.
   *
   * @param node The node defining the subtree to check.
   * @return The number of nodes.
   */
  static unsigned int getNumberOfNodes(const Node& node);

  /**
   * @brief Get the leaves names of a subtree defined by a particular node.
   *
   * @param node The node defining the subtree to check.
   * @return The list of all leaves names.
   */
  static std::vector<std::string> getLeavesNames(const Node& node);

  /**
   * @brief Get the depth of the subtree defined by node 'node', i.e. the maximum
   * number of sons 'generations'.
   *
   * ex:
   * @verbatim
   *    +----------A
   *    |
   * ---+ N1     +-------B
   *    |        |
   *    +--------+ N2
   *             |
   *             +------C
   * @endverbatim
   * Depth of node 'N1' id 2, depth of node 'N2' is 1, depth of leaves is 0.
   *
   * @param node The node defining the subtree to check.
   * @return The depth of the subtree.
   */
  static unsigned int getDepth(const Node& node);

  /**
   * @brief Get the depths for all nodes of the subtree defined by node 'node', i.e. the maximum
   * number of sons 'generations'.
   *
   * ex:
   * @verbatim
   *    +----------A
   *    |
   * ---+ N1     +-------B
   *    |        |
   *    +--------+ N2
   *             |
   *             +------C
   * @endverbatim
   * Depth of node 'N1' id 2, depth of node 'N2' is 1, depth of leaves is 0.
   *
   * @param node The node defining the subtree to check.
   * @param depths The map that will contain all the depths of the nodes, with node pointers as keys.
   * @return The depth of the subtree.
   */
  static unsigned int getDepths(const Node& node, std::map<const Node*, unsigned int>& depths);

  /**
   * @brief Get the height of the subtree defined by node 'node', i.e. the maximum
   * distance between leaves and the root of the subtree.
   *
   * The distance do not include the branch length of the subtree root node.
   * The height of a leaf is hence 0.
   *
   * @param node The node defining the subtree to check.
   * @return The height of the subtree.
   * @throw NodePException If a branch length is lacking.
   */
  static double getHeight(const Node& node);

  /**
   * @brief Get the heights of all nodes within a subtree defined by node 'node', i.e. the maximum
   * distance between leaves and the root of the subtree.
   *
   * The height of a leaf is 0.
   *
   * @param node The node defining the subtree to check.
   * @param heights The map that will contain all the heights of the nodes, with node pointers as keys.
   * @return The height of the subtree.
   * @throw NodePException If a branch length is lacking.
   */
  static double getHeights(const Node& node, std::map<const Node*, double>& heights);

  /**
   * @brief Tell is a subtree is multifurcating.
   *
   * @param node The root node of the subtree.
   * @return True is the subtree contains at least one multifurcating
   * node (including the root node).
   */
  static bool isMultifurcating(const Node& node);

  /**
   * @brief Tells if two subtrees have the same topology.
   *
   * The comparison is based on parental relationships and leaf names only, node ids and all branch/node properties are ignored.
   * The ordering of son nodes is taken into account so that ((A,B),C) will be considered different from ((B,A),C). Considerer
   * ordering the trees first if you want to perform a strict topological comparison.
   *
   * @param n1 Root node of the first subtree.
   * @param n2 Root node of the second subtree.
   * @return true if the two subtrees have the same topology.
   */
  static bool haveSameOrderedTopology(const Node& n1, const Node& n2);

  static std::vector<Node*> getPathBetweenAnyTwoNodes(Node& node1, Node& node2, bool includeAncestor = true);

  static std::vector<const Node*> getPathBetweenAnyTwoNodes(const Node& node1, const Node& node2, bool includeAncestor = true);

  /**
   * @brief Recursively clone a subtree structure.
   *
   * This is a template function allowing to specify the class of the copy.
   * The template class has to have a constructor accepting const Node& as single argument.
   *
   * @param node The basal node of the subtree.
   * @return The basal node of the new copy.
   */
  template<class N>
  static N* cloneSubtree(const Node& node)
  {
    // First we copy this node using default copy constuctor:
    N* clone = new N(node);
    // We remove the link toward the father:
    // clone->removeFather();

    // Now we perform a hard copy:
    for (int i = 0; i < static_cast<int>(node.getNumberOfSons()); i++)
    {
      clone->addSon(cloneSubtree<N>(*node[i]));
    }
    return clone;
  }

  /**
   * @brief Recursively delete a subtree structure.
   *
   * @param node The basal node of the subtree.
   */
  template<class N>
  static void deleteSubtree(N* node)
  {
    for (size_t i = 0; i < node->getNumberOfSons(); ++i)
    {
      N* son = node->getSon(i);
      deleteSubtree(son);
      delete son;
    }
  }


  template<class N>
  static N* cloneSubtree(const Tree& tree, int nodeId)
  {
    // First we copy this node using default copy constuctor:
    N* clone = tree.hasNodeName(nodeId) ? new N(nodeId, tree.getNodeName(nodeId)) : new N(nodeId);
    // Then we set the length:
    if (tree.hasDistanceToFather(nodeId))
      clone->setDistanceToFather(tree.getDistanceToFather(nodeId));
    // Now we copy all sons:
    std::vector<int> sonsId = tree.getSonsId(nodeId);
    for (size_t i = 0; i < sonsId.size(); i++)
    {
      clone->addSon(cloneSubtree<N>(tree, sonsId[i]));
    }
    // Must copy all properties too:
    std::vector<std::string> names;
    names = tree.getNodePropertyNames(nodeId);
    for (size_t i = 0; i < names.size(); i++)
    {
      clone->setNodeProperty(names[i], *tree.getNodeProperty(nodeId, names[i]));
    }
    names = tree.getBranchPropertyNames(nodeId);
    for (size_t i = 0; i < names.size(); i++)
    {
      clone->setBranchProperty(names[i], *tree.getBranchProperty(nodeId, names[i]));
    }

    return clone;
  }
  /** @} */

  /**
   * @name Act on branch lengths.
   *
   * @{
   */

  /**
   * @brief Get all the branch lengths of a subtree.
   *
   * @param node The root node of the subtree.
   * @return A vector with all branch lengths.
   * @throw NodePException If a branch length is lacking.
   */
  static Vdouble getBranchLengths(const Node& node) throw (NodePException);

  /**
   * @brief Get the total length (sum of all branch lengths) of a subtree.
   *
   * @param node The root node of the subtree.
   * @param includeAncestor Tell if the branch length of the most ancient node should be included in the counting.
   * (this should be set to false if this node is the root of the tree for instance).
   * @return The total length of the subtree.
   * @throw NodePException If a branch length is lacking.
   */
  static double getTotalLength(const Node& node, bool includeAncestor = true) throw (NodePException);

  /**
   * @brief Set all the branch lengths of a subtree.
   *
   * @param node  The root node of the subtree.
   * @param brLen The branch length to apply.
   */
  static void setBranchLengths(Node& node, double brLen);

  /**
   * @brief Remove all the branch lengths of a subtree.
   *
   * @param node  The root node of the subtree.
   */
  static void deleteBranchLengths(Node& node);

  /**
   * @brief Give a length to branches that don't have one in a subtree.
   *
   * @param node  The root node of the subtree.
   * @param brLen The branch length to apply.
   */
  static void setVoidBranchLengths(Node& node, double brLen);

  /**
   * @brief Scale a given tree.
   *
   * Multiply all branch lengths by a given factor.
   *
   * @param node   The root node of the subtree to scale.
   * @param factor The factor to multiply all branch lengths with.
   * @throw NodePException If a branch length is lacking.
   */
  static void scaleTree(Node& node, double factor) throw (NodePException);

  /**
   * @brief Get the total distance between to nodes.
   *
   * Sum all branch lengths between two nodes.
   *
   * @param node1 The first node.
   * @param node2 The second node.
   * @return The sum of all branch lengths between the two nodes.
   */
  static double getDistanceBetweenAnyTwoNodes(const Node& node1, const Node& node2);

  /**
   * @brief Compute a distance matrix from a tree.
   *
   * Compute all distances between each leaves and store them in a matrix.
   * A new DistanceMatrix object is created, and a pointer toward it is returned.
   * The destruction of this matrix is left up to the user.
   *
   * From version 1.9 of Bio++, this function has been rewritten in a more efficient way
   * and does not use getDistanceBetweenAnyTwoNodes anymore, but makes use of a more clever
   * pass on the tree. The new function now works well on trees with thousands of leaves.
   *
   * @see getDistanceBetweenAnyTwoNodes
   *
   * @author Nicolas Rochette
   *
   * @param tree The tree to use.
   * @return The distance matrix computed from tree.
   */
  static DistanceMatrix* getDistanceMatrix(const TreeTemplate<Node>& tree);

private:
  /**
   * @brief Inner function used by getDistanceMatrix.
   *
   * (1) Retrieves leaf-leaf distances in node's subtree and
   *  writes them in the distance matrix.
   * (2) Returns distances from node's father to those leaves.
   *
   * @param node The current node in the recursion.
   * @param matrix The output matrix which will be filled.
   * @param distsToNodeFather Intermediate computations contianing the distances of the node to the leaves.
   */
  static void processDistsInSubtree_(const Node* node, DistanceMatrix& matrix, std::vector< std::pair<std::string, double> >& distsToNodeFather);

public:
  /** @} */

  /**
   * @name Conversion tools.
   *
   * Convert from Newick standard tree description.
   * The description is for a node, and hence is to be surrounded with
   * parenthesis. ex: (A:0.001, (B:0.001, C:0.02)90:0.005)50:0.0005
   *
   * @{
   */

  struct Element
  {
public:
    std::string content;
    std::string length;
    std::string bootstrap;
    bool isLeaf;

public:
    Element() : content(),
      length(),
      bootstrap(),
      isLeaf(false) {}
  };

  static Element getElement(const std::string& elt) throw (IOException);

  /**
   * @brief Parse a string in the parenthesis format and convert it to
   * a subtree.
   *
   * @param description the string to parse;
   * @param bootstrap Tell is real bootstrap values are expected. If so, a property with name TreeTools::BOOTSTRAP will be created and stored at the corresponding node.
   * The property value will be of type Number<double>. Otherwise, an object of type String will be created and stored with the property name propertyName.
   * @param propertyName The name of the property to store. Only used if bootstrap = false.
   * @param withId Tells if node ids have been stored in the tree. If set at "true", no bootstrap or property values can be read. Node ids are positioned as bootstrap values for internal nodes, and are concatenated to leaf names after a "_" sign.
   * @return A pointer toward a dynamically created subtree.
   */
  static Node* parenthesisToNode(const std::string& description, bool bootstrap = true, const std::string& propertyName = TreeTools::BOOTSTRAP, bool withId = false);

  /**
   * @brief Parse a string in the parenthesis format and convert it to
   * a tree.
   *
   * @param description the string to parse;
   * @param bootstrap Tells if real bootstrap values are expected. If so, a property with name TreeTools::BOOTSTRAP will be created and stored at the corresponding node.
   * The property value will be of type Number<double>. Otherwise, an object of type String will be created and stored with the property name propertyName.
   * @param propertyName The name of the property to store. Only used if bootstrap = false.
   * @param withId Tells if node ids have been stored in the tree. If set at "true", no bootstrap or property values can be read. Node ids are positioned as bootstrap values for internal nodes, and are concatenated to leaf names after a "_" sign.
   * @return A pointer toward a dynamically created tree.
   * @throw Exception in case of bad format.
   */
  static TreeTemplate<Node>* parenthesisToTree(const std::string& description, bool bootstrap = true, const std::string& propertyName = TreeTools::BOOTSTRAP, bool withId = false) throw (Exception);

  /**
   * @brief Get the parenthesis description of a subtree.
   *
   * @param node The node defining the subtree.
   * @param writeId Tells if node ids must be printed.
   *                This will overwrite bootstrap values if there are ones.
   *                Leaves id will be added to the leave names, separated by a '_' character.
   * @return A string in the parenthesis format.
   */
  static std::string nodeToParenthesis(const Node& node, bool writeId = false);

  /**
   * @brief Get the parenthesis description of a subtree.
   *
   * @param node The node defining the subtree.
   * @param bootstrap Tell is bootstrap values must be writen.
   * If so, the content of the property with name TreeTools::BOOTSTRAP will be written as bootstrap value.
   * The property should be a Number<double> object.
   * Otherwise, the content of the property with name 'propertyName' will be written.
   * In this later case, the property should be a String object.
   * @param propertyName The name of the property to use. Only used if bootstrap = false.
   * @return A string in the parenthesis format.
   */
  static std::string nodeToParenthesis(const Node& node, bool bootstrap, const std::string& propertyName);

  /**
   * @brief Get the parenthesis description of a tree.
   *
   * @param tree The tree to convert.
   * @param writeId Tells if node ids must be printed.
   *                This will overwrite bootstrap values if there are ones.
   *                Leaves id will be added to the leave names, separated by a '_' character.
   * @return A string in the parenthesis format.
   */
  static std::string treeToParenthesis(const TreeTemplate<Node>& tree, bool writeId = false);

  /**
   * @brief Get the parenthesis description of a tree.
   *
   * @param tree The tree to convert.
   * @param bootstrap Tell is bootstrap values must be writen.
   * If so, the content of the property with name TreeTools::BOOTSTRAP will be written as bootstrap value.
   * The property should be a Number<double> object.
   * Otherwise, the content of the property with name 'propertyName' will be written.
   * In this later case, the property should be a String object.
   * @param propertyName The name of the property to use. Only used if bootstrap = false.
   * @return A string in the parenthesis format.
   */
  static std::string treeToParenthesis(const TreeTemplate<Node>& tree, bool bootstrap, const std::string& propertyName);

  /** @} */

  /**
   * @name Random trees
   *
   * @{
   */

  /**
   * @brief Draw a random tree from a list of taxa, using a Yule process.
   *
   * @param leavesNames A list of taxa.
   * @param rooted Tell is the output tree should be rooted.
   * @return A random tree with all corresponding taxa.
   */
  static TreeTemplate<Node>* getRandomTree(std::vector<std::string>& leavesNames, bool rooted = true);

  /** @} */

  /**
   * @brief Get a subset of node neighbors.
   *
   * Get all neighbors of node node1 that are neither node1 nor node2.
   * This method is useful for topology manipulations, like NNI.
   *
   * @param node1 The node whose neighbors must be retrieved.
   * @param node2 One neighbor to exclude.
   * @param node3 Another neighbor to exclude.
   * @return A vector of neighbors.
   */
  static std::vector<const Node*> getRemainingNeighbors(const Node* node1, const Node* node2, const Node* node3);

  /**
   * @brief This method will add a given value (possibly negative) to all identifiers in a (sub)tree.
   *
   * @param node The root node of the (sub)tree to use.
   * @param increment The value to add.
   */
  static void incrementAllIds(Node* node, int increment);

  /**
   * @name Retrieve properties from a (sub)tree.
   *
   * @{
   */

  /**
   * @brief Retrieve the names of all available node properties in the tree.
   *
   * @param node [in] The root node of the (sub)tree to use.
   * @param propertyNames [out] a vector where names will be added.
   */
  static void getNodePropertyNames(const Node& node, std::vector<std::string>& propertyNames);

  /**
   * @brief Retrieve all node property objects with a given name over a (sub) tree (const version).
   *
   * @param node [in] The root node of the (sub)tree to use.
   * @param propertyName [in] The name of the property to retrieve.
   * @param properties [out] A map with pointers toward the properties as values, and node ids as key.
   * If a node does not contain the given property, then no entry in the map is created.
   * If an entry already exists in the map, it will be replaced, but the underlying property will not be destroyed.
   * Property objects are not cloned when added to the map, but passed as pointers.
   */
  static void getNodeProperties(const Node& node, const std::string& propertyName, std::map<int, const Clonable*>& properties);

  /**
   * @brief Retrieve all node property objects with a given name over a (sub) tree.
   *
   * @param node [in] The root node of the (sub)tree to use.
   * @param propertyName [in] The name of the property to retrieve.
   * @param properties [out] A map with pointers toward the properties as values, and node ids as key.
   * If a node does not contain the given property, then no entry in the map is created.
   * If an entry already exists in the map, it will be replaced, but the underlying property will not be destroyed.
   * Property objects are not cloned when added to the map, but passed as pointers.
   */
  static void getNodeProperties(Node& node, const std::string& propertyName, std::map<int, Clonable*>& properties);

  /**
   * @brief Retrieve the names of all available branch properties in the tree.
   *
   * @param node [in] The root node of the (sub)tree to use.
   * @param propertyNames [out] a vector where names will be added.
   */
  static void getBranchPropertyNames(const Node& node, std::vector<std::string>& propertyNames);

  /**
   * @brief Retrieve all branch property objects with a given name over a (sub) tree (const version).
   *
   * @param node [in] The root node of the (sub)tree to use.
   * @param propertyName [in] The name of the property to retrieve.
   * @param properties [out] A map with pointers toward the properties as values, and node ids as key.
   * If a node does not contain the given property, then no entry in the map is created.
   * If an entry already exists in the map, it will be replaced, but the underlying property will not be destroyed.
   * Property objects are not cloned when added to the map, but passed as pointers.
   */
  static void getBranchProperties(const Node& node, const std::string& propertyName, std::map<int, const Clonable*>& properties);

  /**
   * @brief Retrieve all branch property objects with a given name over a (sub) tree.
   *
   * @param node [in] The root node of the (sub)tree to use.
   * @param propertyName [in] The name of the property to retrieve.
   * @param properties [out] A map with pointers toward the properties as values, and node ids as key.
   * If a node does not contain the given property, then no entry in the map is created.
   * If an entry already exists in the map, it will be replaced, but the underlying property will not be destroyed.
   * Property objects are not cloned when added to the map, but passed as pointers.
   */
  static void getBranchProperties(Node& node, const std::string& propertyName, std::map<int, Clonable*>& properties);

  /**
   * @brief Swap nodes in the subtree so that they are ordered according to the underlying number of leaves.
   *
   * @param node The root node of the (sub)tree to use.
   * @param downward If yes, biggest subtrees (in terms of number of leaves) will come first. Otherwise, the smallest subtrees will come first.
   * @param orderLeaves Tell if leaves have to be ordered alphabetically. This ensures that two identical topology will always have the same ordered tree, whatever the initial ordering of nodes.
   */
  static void orderTree(Node& node, bool downward = true, bool orderLeaves = false)
  {
    orderTree_(node, downward, orderLeaves);
  }
  /** @} */

  /**
   * @brief Midroot the tree by minimizing a given criterion ("variance" or "sum of squares")
   *
   * @details
   * For each branch, the best root position, according to the given criterion, is computed analytically.
   *
   * For the 'variance' criterion :
   * \f[
   *  (n_1+n_2)^2 V(x)
   *   = (n_1+n_2) \left[ \sum_{F_1} (d_i + x \delta )^2 + \sum_{F_2} (d_i + (1-x) \delta )^2 \right]
   *     - \left[ \sum_{F_1} (d_i + x \delta) + \sum_{F_2} (d_i + (1-x) \delta) \right]^2
   *   = A x^2 + B x + C
   * \f]
   * With
   * \f[ \begin{array}{rcl}
   * A &=& 4 n_1 n_2 \delta^2 \\
   * B &=& 4 \delta ( n_2 S_1 - n_1 S_2 - n_1 n_2 \delta ) \\
   * C &=& (n_1+n_2) (C_1+C_2) + n_1 n_2 \delta^2 + 2 n_1 S_2 \delta - 2 n_2 S_1 \delta - (S_1+S_2)^2 \\
   * \end{array} \f]
   *
   * Where \f$F_1\f$ and \f$F_2\f$ are the sets of leaves on either side of
   * the root branch,
   * \f$d_i\f$ is the distance of leaf \f$i\f$ to the nearest end of the root branch,
   * \f$\delta\f$ is the length of the root branch, and \f$S_k\f$ and \f$C_k\f$ are respectively
   * \f$\sum_{F_k} d_i\f$ and \f$\sum_{F_k} d_i^2\f$
   *
   * ~
   *
   * If force_branch_root==true, then the function will always root the tree on a branch.
   * To do so, in cases where the root is placed on a node, a new node new_root is created between the root and its nearest child.
   * If force_branch_root==false, it may be placed on a node.
   *
   *
   * @param tree
   * @param criterion The criterion upon which to reroot. Legal values : TreeTemplateTools::MIDROOT_VARIANCE
   *   to minimize root-leaf distance variance (molecular clock assumption) or
   *   TreeTemplateTools::MIDROOT_SUM_OF_SQUARES to minimize the sum of root-leaf distance squares.
   * @param force_branch_root If true, the root must be placed on a branch, otherwise it may also be placed on a node. 
   *
   * @author Nicolas Rochette
   */
  static void midRoot(TreeTemplate<Node>& tree, short criterion, const bool force_branch_root);

  /**
   * @brief Get the caracteristic radius of a tree (average distance to the root minimizing the sum of squared distances).
   *
   * @param tree The tree (which is rerooted in the process).
   */
  static double getRadius(TreeTemplate<Node>& tree);


  /**
   * @brief Unresolve nodes with low confidence value.
   *
   * The underlying branches will be removed, resulting in a multifurcation.
   * the branch length of the removed node is added to the length of its son nodes,
   * so that pairwise phylogenetic distances are conserved along the tree.
   * Leaves are not checked. Node with missing values are ignored.
   *
   * @author Julien Dutheil.
   *
   * @param subtree   The node defining the subtree where nodes should be collapsed.
   * @param threshold The minimum value for which a node is considered to be confident.
   * @param property  The branch property to be considered as a confidence value (bootstrap vlaues by default).
   */
  static void unresolveUncertainNodes(Node& subtree, double threshold, const std::string& property = TreeTools::BOOTSTRAP);

private:
  struct OrderTreeData_
  {
    size_t size;
    std::string firstLeaf;
    OrderTreeData_() : size(0),
      firstLeaf("") {}
  };

  static OrderTreeData_ orderTree_(Node& node, bool downward, bool orderLeaves);

  /**
   * @brief
   * A <i>structure</i> recording, for a subtree, the sum of root-leaf distances, the sum of their squares,
   * and the number of elements in these sums (ie. the number of leaves).
   *
   * @details
   * The branch at the base of the subtree should never be included,
   * as the subtree of the root does not have one.
   *
   */
  struct Moments_
  {
    double sum;
    double squaresSum;
    int numberOfLeaves;
  };

  /**
   * @brief
   * Computes the moment of a subtree
   *
   * @param node The root of the subtree
   * @return A Moments_ structure
   */
  static Moments_ getSubtreeMoments_(const Node* node);

  /**
   * @brief Find, in the branches of a subtree, the root that minimizes a criterion over the tree.
   *
   * @details
   * The branches are explored recursively. For each branch leaving the input node, the method
   * computes the best root position, possibly updates the bestRoot parameter, then recurses.
   *
   * @param tree The tree to which the subtree belongs. (The root is moved.)
   * @param criterion The criterion to minimize. Legal values are TreeTemplateTools::MIDROOT_VARIANCE and TreeTemplateTools::MIDROOT_SUM_OF_SQUARES.
   * @param node The root of the subtree.
   * @param bestRoot The object storing the best root found, if it is better than the initial one, or otherwise left unchanged.
   *
   * @author Nicolas Rochette, Manolo Gouy
   */
  static void getBestRootInSubtree_(bpp::TreeTemplate<bpp::Node>& tree, short criterion,  bpp::Node* node, std::pair<bpp::Node*, std::map<std::string, double> >& bestRoot);

public:
  static const short MIDROOT_VARIANCE;
  static const short MIDROOT_SUM_OF_SQUARES;
};
} // end of namespace bpp.

#endif // _TREETEMPLATETOOLS_H_

