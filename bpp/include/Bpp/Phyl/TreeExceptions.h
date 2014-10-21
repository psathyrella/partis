//
// File: TreeExceptions.h
// Created by: Julien Dutheil
// Created on: Mon Nov  3 17:04:46 2003
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _TREEEXCEPTIONS_H_
#define _TREEEXCEPTIONS_H_

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

// From the STL:
#include <string>

namespace bpp
{
class Node;
class Tree;

/**
 * @brief General exception thrown when something is wrong with a particular node.
 */
class NodeException :
  public Exception
{
protected:
  int nodeId_;

public:
  /**
   * @brief Build a new NodePException.
   * @param text A message to be passed to the exception hierarchy.
   * @param nodeId The id of the node that threw the exception.
   */
  NodeException(const std::string& text, int nodeId) :
    Exception("NodeException: " + text + "(id:" + TextTools::toString(nodeId) + ")"),
    nodeId_(nodeId) {}

  virtual ~NodeException() throw () {}

public:
  /**
   * @brief Get the id of node that threw the exception.
   *
   * @return The id of the faulty node.
   */
  virtual int getNodeId() const { return nodeId_; }
};


/**
 * @brief General exception thrown when something is wrong with a particular node.
 */
class NodePException :
  public NodeException
{
private:
  const Node* node_;

public:
  /**
   * @brief Build a new NodePException.
   * @param text A message to be passed to the exception hierarchy.
   * @param node A const pointer toward the node that threw the exception.
   */
  NodePException(const std::string& text, const Node* node = 0);

  /**
   * @brief Build a new NodePException.
   * @param text A message to be passed to the exception hierarchy.
   * @param nodeId The id of the node that threw the exception.
   */
  NodePException(const std::string& text, int nodeId) :
    NodeException(text, nodeId), node_(0) {}

  NodePException(const NodePException& nex) :
    NodeException(nex),
    node_(nex.node_)
  {}

  NodePException& operator=(const NodePException& nex)
  {
    NodeException::operator=(nex);
    node_ = nex.node_;
    return *this;
  }

  virtual ~NodePException() throw () {}

public:
  /**
   * @brief Get the node that threw the exception.
   *
   * @return A pointer toward the faulty node.
   */
  virtual const Node* getNode() const { return node_; };
  /**
   * @brief Get the id of node that threw the exception.
   *
   * @return The id of the faulty node.
   */
  virtual int getNodeId() const { return nodeId_; }
};

/**
 * @brief General exception thrown if a property could not be found.
 */
class PropertyNotFoundException :
  public NodePException
{
private:
  std::string propertyName_;

public:
  /**
   * @brief Build a new PropertyNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param node A const pointer toward the node that threw the exception.
   */
  PropertyNotFoundException(const std::string& text, const std::string& propertyName, const Node* node = 0) :
    NodePException("Property not found: " + propertyName + ". " + text, node),
    propertyName_(propertyName) {}

  /**
   * @brief Build a new PropertyNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param nodeId The id of the node that threw the exception.
   */
  PropertyNotFoundException(const std::string& text, const std::string& propertyName, int nodeId) :
    NodePException("Property not found: " + propertyName + ". " + text, nodeId),
    propertyName_(propertyName) {}

  virtual ~PropertyNotFoundException() throw () {}

public:
  /**
   * @brief Get the name of the property that could not be found.
   *
   * @return The name of the missing property.
   */
  virtual const std::string& getPropertyName() const { return propertyName_; }
};

/**
 * @brief Exception thrown when something is wrong with a particular node.
 */
class NodeNotFoundException :
  public Exception
{
private:
  std::string id_;

public:
  /**
   * @brief Build a new NodeNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param id   A string describing the node.
   */
  NodeNotFoundException(const std::string& text, const std::string& id);

  /**
   * @brief Build a new NodeNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param id   A node identifier.
   */
  NodeNotFoundException(const std::string& text, int id);

  virtual ~NodeNotFoundException() throw () {}

public:
  /**
   * @brief Get the node id that threw the exception.
   *
   * @return The id of the node.
   */
  virtual std::string getId() const { return id_; }
};

/**
 * @brief General exception thrown when something wrong happened in a tree.
 */
class TreeException :
  public Exception
{
private:
  const Tree* tree_;

public:
  /**
   * @brief Build a new TreeException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param tree A const pointer toward the tree that threw the exception.
   */
  TreeException(const std::string& text, const Tree* tree = 0);

  TreeException(const TreeException& tex) :
    Exception(tex),
    tree_(tex.tree_)
  {}

  TreeException& operator=(const TreeException& tex)
  {
    Exception::operator=(tex);
    tree_ = tex.tree_;
    return *this;
  }

  virtual ~TreeException() throw () {}

public:
  /**
   * @brief Get the tree that threw the exception.
   *
   * @return The faulty tree
   */
  virtual const Tree* getTree() const { return tree_; }
};

/**
 * @brief Exception thrown when a tree is expected to be rooted.
 */
class UnrootedTreeException :
  public TreeException
{
public:
  /**
   * @brief Build a new UnrootedTreeException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param tree A const pointer toward the tree that threw the exception.
   */
  UnrootedTreeException(const std::string& text, const Tree* tree = 0);

  virtual ~UnrootedTreeException() throw () {}
};

} // end of namespace bpp.

#endif  // _TREEEXCEPTIONS_H_

