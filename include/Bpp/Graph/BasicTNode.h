// 
// File:    BasicTNode.h
// Author:  Sylvain Gaillard
// Created: 13/01/2011 16:39:23
// 

/*
Copyright or © or Copr. CNRS, (January 12, 2011)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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

#ifndef _BPP_GRAPH_BASICTNODE_H_
#define _BPP_GRAPH_BASICTNODE_H_

#include "TNode.h"
#include "../Exceptions.h"

namespace bpp {
  /**
   * @brief Simple implementation of TNode
   *
   * Contains only methods for node manipulation.
   *
   * @author Sylvain Gaillard
   */
  class BasicTNode: public TNode {
    private:
      std::vector< BasicTNode * > sons_;
      BasicTNode * father_;

    public:
      /**
       * @brief Simple constructor.
       */
      BasicTNode(): sons_(), father_(0) {};

      /**
       * @brief Destructor.
       *
       * When destroyed, the node remove himself as son of its father and as
       * father of its sons.
       */
      virtual ~BasicTNode();

      /**
       * @brief Copy constructor.
       */
      BasicTNode(const BasicTNode& node);

      /**
       * @brief Assignation operator.
       */
      BasicTNode& operator=(const BasicTNode& node);

      BasicTNode* clone() const {
        return new BasicTNode(* this);
      }

      // Neighbors

      const BasicTNode* getNeighbor(int pos) const;
      BasicTNode* getNeighbor(int pos);

      int degree() const { return sons_.size() - 1 + father_ ? 1 : 0; }

      const BasicTNode* operator[](int i) const;
      BasicTNode* operator[](int i);

      // Fathers

      /**
       * @brief Tell if the node has a father.
       */
      bool hasFather() const { return father_ ? true : false; }
      bool hasFathers() const { return father_ ? true : false; }
      int getNumberOfFathers() const { return father_ ? 1 : 0; }

      const BasicTNode* getFather(int pos) const;
      BasicTNode* getFather(int pos);
      const BasicTNode* getFather() const;
      BasicTNode* getFather();

      /**
       * @brief Tell if the node is a father of this node.
       */
      virtual bool isFather(const BasicTNode* node) const;

      /**
       * @brief Add a father to this node.
       */
      virtual void addFather(BasicTNode* node);

      /**
       * @brief Remove the father of this node.
       *
       * @return A pointer to the removed father node.
       */
      virtual BasicTNode* removeFather();

      // Sons

      bool hasSons() const { return !sons_.empty(); }
      int getNumberOfSons() const { return static_cast<int>(sons_.size()); }
      const BasicTNode* getSon(int pos) const;
      BasicTNode* getSon(int pos);

      /**
       * @brief Tell if a node is son of this node.
       */
      virtual bool isSon(const BasicTNode* node) const;

      /**
       * @brief Add a son to this node.
       */
      virtual void addSon(BasicTNode* node);

      /**
       * @brief Remove a son of this node.
       */
      virtual void removeSon(BasicTNode* son);

      /**
       * @brief Remove a son of this node.
       *
       * @return A pointer to the removed son node or a Null pointer if son is
       * not found.
       */
      virtual BasicTNode* removeSon(int pos);
  };
}

#endif //_BPP_GRAPH_BASICTNODE_H_
