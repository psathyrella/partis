// 
// File:    ONode.h
// Author:  Sylvain Gaillard
// Created: 12/01/2011 08:36:47
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

#ifndef _BPP_GRAPH_ONODE_H_
#define _BPP_GRAPH_ONODE_H_

#include "UNode.h"

namespace bpp {
  /**
   * @brief Oriented Node interface
   *
   * ONode is an interface for oriented nodes aimed to build oriented graphs.
   *
   * @author Sylvain Gaillard
   */
  class ONode: public virtual UNode {
    public:
      /**
       * @name Neighbors
       *
       * @{
       */

      virtual const ONode * getNeighbor(int pos) const = 0;
      virtual ONode * getNeighbor(int pos) = 0;

      /** @} */

      /**
       * @name The Clonable interface.
       *
       * @{
       */
#ifndef NO_VIRTUAL_COV      
      ONode * clone() const = 0;
#endif
      /** @} */

      /**
       * @name Fathers
       *
       * @{
       */

      /**
       * @brief Get a particular father in const environment.
       */
      virtual const ONode * getFather(int pos) const = 0;

      /**
       * @brief Get a particular father.
       */
      virtual ONode * getFather(int pos) = 0;

      /**
       * @brief Tell if this node has one or more father nodes.
       */
      virtual bool hasFathers() const = 0;

      /**
       * @brief Give the number of father nodes for this node.
       */
      virtual int getNumberOfFathers() const = 0;

      /** @} */

      /**
       * @name Sons
       *
       * @{
       */

      /**
       * @brief Get a particular son in const environment.
       */
      virtual const ONode * getSon(int pos) const = 0;

      /**
       * @brief Get a particular son.
       */
      virtual ONode * getSon(int pos) = 0;

      /**
       * @brief Tell if this node has one or more son nodes.
       */
      virtual bool hasSons() const = 0;

      /**
       * @brief Give the number of son nodes for this node.
       */
      virtual int getNumberOfSons() const = 0;

      /** @} */

      /**
       * @name Operators
       *
       * @{
       */

      /**
       * @brief Direct access to a neighbor in const context.
       *
       * - a positive i gives access to sons (from 0 to n - 1)
       * - a negative i gives access to fathers (from 1 to m)
       *
       * No check is done, you have to ensure that you query an existing
       * neighbor.
       *
       * @param i the position of the neighbor
       * @return A pointer toward the neighbor
       */
      virtual const ONode * operator[] (int i) const = 0;

      /**
       * @brief Direct access to a neighbor.
       *
       * - a positive i gives access to sons (from 0 to n - 1)
       * - a negative i gives access to fathers (from 1 to m)
       *
       * No check is done, you have to ensure that you query an existing
       * neighbor.
       *
       * @param i the position of the neighbor
       * @return A pointer toward the neighbor
       */
      virtual ONode * operator[] (int i) = 0;

      /** @} */
  };
}

#endif //_BPP_GRAPH_ONODE_H_
