//
// File: Node.h
// Created by: Julien Dutheil
// Created on: Thu Mar 13 12:03:18 2003
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

#ifndef _NODE_H_
#define _NODE_H_

#include "TreeExceptions.h"

#include <Bpp/Clonable.h>
#include <Bpp/Utils/MapTools.h>
#include <Bpp/BppString.h>
#include <Bpp/Numeric/Number.h>

// From the STL:
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>

namespace bpp
{
/**
 * @brief The phylogenetic node class.
 *
 * This class is for use with the TreeTemplate class, an implementation of the Tree interface.
 * TreeTemplates are made made of nodes, instances of this class.
 * Since trees are oriented (rooted), each node has one <i>father node</i> and possibly
 * many <i>son nodes</i>. Leaves are nodes without descendant and root is defined has the without
 * father. Inner nodes will generally contain two descendants (the tree is then called
 * <i>bifurcating</i>), but mutlifurcating trees are also allowed with this kind of description.
 * In the rooted case, each inner node also defines a <i>subtree</i>.
 * This allows to work recursively on trees, which is very convenient in most cases.</p>
 *
 * This class is made the more general as possible, while keeping it very simple. It contains:</p>
 * - An identity tag, to identity it in the tree;
 * - A name, necessary for leaf nodes, optionnal else;
 * - A pointer toward the father node;
 * - A std::vector of pointer toward son nodes;
 * - The distance from the father node:
 * - A property map, that may contain any information to link to each node, e.g. bootstrap
 * value or GC content.
 *
 * Methods are provided to help the building of trees from scratch.
 * Trees are more easily built from root to leaves:
 * The addSon(Node) method adds a node to the list of direct descendants of a
 * given node. The son node will also have its father set to the current node.
 * It is also possible to build a tree starting from the leaves using the setFather method.
 * Changing the parent node will automatically append the current node to the son nodes of the new father.
 *
 * @see Tree, TreeTemplate
 */
class Node
{

protected:
  int id_;
  std::string* name_;
  std::vector<Node*> sons_;
  Node* father_;
  double* distanceToFather_;
  mutable std::map<std::string, Clonable*> nodeProperties_;
  mutable std::map<std::string, Clonable*> branchProperties_;

public:
  /**
   * @brief Build a new void Node object.
   */
  Node() :
    id_(0),
    name_(0),
    sons_(),
    father_(0),
    distanceToFather_(0),
    nodeProperties_(),
    branchProperties_()
  {}

  /**
   * @brief Build a new Node with specified id.
   */
  Node(int id) :
    id_(id),
    name_(0),
    sons_(),
    father_(0),
    distanceToFather_(0),
    nodeProperties_(),
    branchProperties_()
  {}

  /**
   * @brief Build a new Node with specified name.
   */
  Node(const std::string& name) :
    id_(0),
    name_(new std::string(name)),
    sons_(),
    father_(0),
    distanceToFather_(0),
    nodeProperties_(),
    branchProperties_()
  {}

  /**
   * @brief Build a new Node with specified id and name.
   */
  Node(int id, const std::string& name) :
    id_(id),
    name_(new std::string(name)),
    sons_(),
    father_(0),
    distanceToFather_(0),
    nodeProperties_(),
    branchProperties_()
  {}

  /**
   * @brief Copy constructor.
   *
   * @warning This constructor copies all fields, excepted father and son node pointers.
   *
   * @param node The node to copy.
   */
  Node(const Node& node);

  /**
   * @brief Assignation operator.
   *
   * @warning This operator copies all fields, excepted father and son node pointers.
   *
   * @param node the node to copy.
   * @return A reference toward this node.
   */
  Node& operator=(const Node& node);

  Node* clone() const { return new Node(*this); }

public:
  virtual ~Node()
  {
    if (name_) delete name_;
    if (distanceToFather_) delete distanceToFather_;
    for (std::map<std::string, Clonable*>::iterator i = nodeProperties_.begin(); i != nodeProperties_.end(); i++)
    {
      delete i->second;
    }
    for (std::map<std::string, Clonable*>::iterator i = branchProperties_.begin(); i != branchProperties_.end(); i++)
    {
      delete i->second;
    }
  }

public:
  /**
   * @name Identity
   *
   * @{
   */

  /**
   * @brief Get the node's id.
   *
   * @return The identity tag of this node.
   */
  virtual int getId() const { return id_; }

  /**
   * @brief Set this node's id.
   *
   * @param id The new identity tag.
   */
  virtual void setId(int id) { id_ = id; }

  virtual std::vector<int> getSonsId() const
  {
    std::vector<int> sonsId(sons_.size());
    for (size_t i = 0; i < sons_.size(); i++)
    {
      sonsId[i] = sons_[i]->getId();
    }
    return sonsId;
  }

  /** @} */

  /**
   * @name Name:
   *
   * @{
   */

  /**
   * @brief Get the name associated to this node, if there is one,
   * otherwise throw a NodeException.
   *
   * @return The name associated to this node.
   */
  virtual std::string getName() const throw (NodePException)
  {
    if (!hasName()) throw NodePException("Node::getName: no name associated to this node.", this);
    return *name_;
  }

  /**
   * @brief Give a name or update the name associated to the node.
   *
   * @param name The name to give to the node.
   */
  virtual void setName(const std::string& name)
  {
    if (name_) delete name_;
    name_ = new std::string(name);
  }

  /**
   * @brief Delete the name associated to this node (do nothing if there is no name).
   */
  virtual void deleteName()
  {
    if (name_) delete name_;
    name_ = 0;
  }

  /**
   * @brief Tell is this node has a name.
   *
   * @return True if name != 0.
   */
  virtual bool hasName() const { return name_ != 0; }

  /** @} */

  /**
   * @name Distances:
   *
   * @{
   */

  /**
   * @brief Get the distance to the father node is there is one,
   * otherwise throw a NodeException.
   *
   * @return The distance to the father node.
   */
  virtual double getDistanceToFather() const
  {
    if (!hasDistanceToFather())
      throw NodePException("Node::getDistanceToFather: Node has no distance.", this);
    return *distanceToFather_;
  }

  /**
   * @brief Set or update the distance toward the father node.
   *
   * Warning: a distance to the father node may be set even if no father node is specified.
   * This is used by several tree reconstruction methods.
   * It may also be useful for manipulating subtrees.
   *
   * @param distance The new distance to the father node.
   */
  virtual void setDistanceToFather(double distance)
  {
    if (distanceToFather_)
      delete distanceToFather_;
    distanceToFather_ = new double(distance);
  }

  /**
   * @brief Delete the distance to the father node.
   */
  virtual void deleteDistanceToFather()
  {
    if (distanceToFather_)
      delete distanceToFather_;
    distanceToFather_ = 0;
  }

  /**
   * @brief Tell is this node has a distance to the father.
   *
   * @return True if distanceToFather != 0.
   */
  virtual bool hasDistanceToFather() const
  {
    return distanceToFather_ != 0;
  }

  /** @} */

  /**
   * @name Father:
   *
   * @{
   */

  /**
   * @brief Get the father of this node is there is one.
   *
   * @return A pointer toward the father node, 0 if there is not.
   */
  virtual const Node* getFather() const { return father_; }

  /**
   * @brief Get the father of this node is there is one.
   *
   * @return A pointer toward the father node, 0 if there is not.
   */
  virtual Node* getFather() { return father_; }

  virtual int getFatherId() const { return father_->getId(); }

  /**
   * @brief Set the father node of this node.
   *
   * @param node The father node.
   */
  virtual void setFather(Node* node) throw (NullPointerException)
  {
    if (!node)
      throw NullPointerException("Node::setFather(). Empty node given as input.");
    father_ = node;
    if (find(node->sons_.begin(), node->sons_.end(), this) == node->sons_.end())
      node->sons_.push_back(this);
    else // Otherwise node is already present.
      std::cerr << "DEVEL warning: Node::setFather. Son node already registered! No pb here, but could be a bug in your implementation..." << std::endl;
  }

  /**
   * @brief Remove the father of this node.
   */
  virtual Node* removeFather()
  {
    Node* f = father_;
    father_ = 0;
    return f;
  }

  /**
   * @brief Tell if this node has a father node.
   */
  virtual bool hasFather() const { return father_ != 0; }

  /** @} */

  /**
   * @name Sons:
   *
   * @{
   */
  virtual size_t getNumberOfSons() const { return sons_.size(); }

  virtual std::vector<Node*>& getSons()
  {
    return sons_;
  }

  virtual const Node* getSon(size_t pos) const throw (IndexOutOfBoundsException)
  {
    if (pos >= sons_.size()) throw IndexOutOfBoundsException("Node::getSon().", pos, 0, sons_.size() - 1);
    return sons_[pos];
  }

  virtual Node* getSon(size_t pos) throw (IndexOutOfBoundsException)
  {
    if (pos >= sons_.size()) throw IndexOutOfBoundsException("Node::getSon().", pos, 0, sons_.size() - 1);
    return sons_[pos];
  }

  virtual void addSon(size_t pos, Node* node) throw (NullPointerException, NodePException)
  {
    if (!node)
      throw NullPointerException("Node::addSon(). Empty node given as input.");
    if (find(sons_.begin(), sons_.end(), node) == sons_.end())
      sons_.insert(sons_.begin() + pos, node);
    else // Otherwise node is already present.
      std::cerr << "DEVEL warning: Node::addSon. Son node already registered! No pb here, but could be a bug in your implementation..." << std::endl;

    node->father_ = this;
  }

  virtual void addSon(Node* node) throw (NullPointerException, NodePException)
  {
    if (!node)
      throw NullPointerException("Node::addSon(). Empty node given as input.");
    if (find(sons_.begin(), sons_.end(), node) == sons_.end())
      sons_.push_back(node);
    else // Otherwise node is already present.
      throw NodePException("Node::addSon. Trying to add a node which is already present.");
    node->father_ = this;
  }

  virtual void setSon(size_t pos, Node* node) throw (IndexOutOfBoundsException, NullPointerException, NodePException)
  {
    if (!node)
      throw NullPointerException("Node::setSon(). Empty node given as input.");
    if (pos >= sons_.size())
      throw IndexOutOfBoundsException("Node::setSon(). Invalid node position.", pos, 0, sons_.size() - 1);
    std::vector<Node*>::iterator search = find(sons_.begin(), sons_.end(), node);
    if (search == sons_.end() || search == sons_.begin() + pos)
      sons_[pos] = node;
    else
      throw NodePException("Node::setSon. Trying to set a node which is already present.");
    node->father_ = this;
  }

  virtual Node* removeSon(size_t pos) throw (IndexOutOfBoundsException)
  {
    if (pos >= sons_.size())
      throw IndexOutOfBoundsException("Node::removeSon(). Invalid node position.", pos, 0, sons_.size() - 1);
    Node* node = sons_[pos];
    sons_.erase(sons_.begin() + pos);
    node->removeFather();
    return node;
  }

  virtual void removeSon(Node* node) throw (NodeNotFoundException, NullPointerException)
  {
    if (!node)
      throw NullPointerException("Node::removeSon(). Empty node given as input.");
    for (size_t i = 0; i < sons_.size(); i++)
    {
      if (sons_[i] == node)
      {
        sons_.erase(sons_.begin() + i);
        node->removeFather();
        return;
      }
    }
    throw NodeNotFoundException("Node::removeSon.", node->getId());
  }

  virtual void removeSons()
  {
    while (sons_.size() != 0)
      removeSon(static_cast<size_t>(0));
  }

  virtual void swap(size_t branch1, size_t branch2) throw (IndexOutOfBoundsException);

  virtual size_t getSonPosition(const Node* son) const throw (NodeNotFoundException, NullPointerException);

  /** @} */

  // These functions must not be declared as virtual!!

  std::vector<const Node*> getNeighbors() const;

  std::vector<Node*> getNeighbors();

  virtual size_t degree() const { return getNumberOfSons() + (hasFather() ? 1 : 0); }

  /**
   * @name Operators:
   *
   * - a positive value send the corresponding son;
   * - a negative value send the father.
   *
   * @{
   */
  Node* operator[](int i) { return (i < 0) ? father_ : sons_[i]; }

  const Node* operator[](int i) const { return (i < 0) ? father_ : sons_[i]; }

  /** @} */

  /**
   * @name Node properties:
   *
   * @{
   */

  /**
   * @brief Set/add a node property.
   *
   * If no property with the same name is found, the new property will be added to the list.
   * Conversely, the property will be deleted and replaced by the new one.
   * If you want to keep a copy of the old property, consider using the removeNodeProperty function before.
   *
   * @param name The name of the property to set.
   * @param property The property object (will be cloned).
   */
  virtual void setNodeProperty(const std::string& name, const Clonable& property)
  {
    if (hasNodeProperty(name))
      delete nodeProperties_[name];
    nodeProperties_[name] = property.clone();
  }

  virtual Clonable* getNodeProperty(const std::string& name) throw (PropertyNotFoundException)
  {
    if (hasNodeProperty(name))
      return nodeProperties_[name];
    else
      throw PropertyNotFoundException("", name, this);
  }

  virtual const Clonable* getNodeProperty(const std::string& name) const throw (PropertyNotFoundException)
  {
    if (hasNodeProperty(name))
      return const_cast<const Clonable*>(nodeProperties_[name]);
    else
      throw PropertyNotFoundException("", name, this);
  }

  virtual Clonable* removeNodeProperty(const std::string& name) throw (PropertyNotFoundException)
  {
    if (hasNodeProperty(name))
    {
      Clonable* removed = nodeProperties_[name];
      nodeProperties_.erase(name);
      return removed;
    }
    else
      throw PropertyNotFoundException("", name, this);
  }

  virtual void deleteNodeProperty(const std::string& name) throw (PropertyNotFoundException)
  {
    if (hasNodeProperty(name))
    {
      delete nodeProperties_[name];
      nodeProperties_.erase(name);
    }
    else
      throw PropertyNotFoundException("", name, this);
  }

  /**
   * @brief Remove all node properties.
   *
   * Attached objects will not be deleted.
   */
  virtual void removeNodeProperties()
  {
    nodeProperties_.clear();
  }

  /**
   * @brief Delete all node properties.
   */
  virtual void deleteNodeProperties()
  {
    for (std::map<std::string, Clonable*>::iterator i = nodeProperties_.begin(); i != nodeProperties_.end(); i++)
    {
      delete i->second;
    }
    nodeProperties_.clear();
  }

  virtual bool hasNodeProperty(const std::string& name) const { return nodeProperties_.find(name) != nodeProperties_.end(); }

  virtual std::vector<std::string> getNodePropertyNames() const { return MapTools::getKeys(nodeProperties_); }

  /** @} */

  /**
   * @name Branch properties:
   *
   * @{
   */

  /**
   * @brief Set/add a branch property.
   *
   * If no property with the same name is found, the new property will be added to the list.
   * Conversely, the property will be deleted and replaced by the new one.
   * If you want to keep a copy of the old property, consider using the removeBranchProperty function before.
   *
   * @param name The name of the property to set.
   * @param property The property object (will be cloned).
   */
  virtual void setBranchProperty(const std::string& name, const Clonable& property)
  {
    if (hasBranchProperty(name))
      delete branchProperties_[name];
    branchProperties_[name] = property.clone();
  }

  virtual Clonable* getBranchProperty(const std::string& name) throw (PropertyNotFoundException)
  {
    if (hasBranchProperty(name))
      return branchProperties_[name];
    else
      throw PropertyNotFoundException("", name, this);
  }

  virtual const Clonable* getBranchProperty(const std::string& name) const throw (PropertyNotFoundException)
  {
    if (hasBranchProperty(name))
      return const_cast<const Clonable*>(branchProperties_[name]);
    else
      throw PropertyNotFoundException("", name, this);
  }

  virtual Clonable* removeBranchProperty(const std::string& name) throw (PropertyNotFoundException)
  {
    if (hasBranchProperty(name))
    {
      Clonable* removed = branchProperties_[name];
      branchProperties_.erase(name);
      return removed;
    }
    else
      throw PropertyNotFoundException("", name, this);
  }

  virtual void deleteBranchProperty(const std::string& name) throw (PropertyNotFoundException)
  {
    if (hasBranchProperty(name))
    {
      delete branchProperties_[name];
      branchProperties_.erase(name);
    }
    else
      throw PropertyNotFoundException("", name, this);
  }

  /**
   * @brief Remove all branch properties.
   *
   * Attached objects will not be deleted.
   */
  virtual void removeBranchProperties()
  {
    branchProperties_.clear();
  }

  /**
   * @brief Delete all branch properties.
   */
  virtual void deleteBranchProperties()
  {
    for (std::map<std::string, Clonable*>::iterator i = branchProperties_.begin(); i != branchProperties_.end(); i++)
    {
      delete i->second;
    }
    branchProperties_.clear();
  }

  virtual bool hasBranchProperty(const std::string& name) const { return branchProperties_.find(name) != branchProperties_.end(); }

  virtual std::vector<std::string> getBranchPropertyNames() const { return MapTools::getKeys(branchProperties_); }

  virtual bool hasBootstrapValue() const;

  virtual double getBootstrapValue() const throw (PropertyNotFoundException);
  /** @} */
  // Equality operator:

  virtual bool operator==(const Node& node) const { return id_ == node.id_; }

  // Tests:

  virtual bool isLeaf() const { return degree() <= 1; }

};
} // end of namespace bpp.

#endif  // _NODE_H_

