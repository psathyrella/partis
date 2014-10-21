// 
// File:    UNode.h
// Author:  Sylvain Gaillard
// Created: 12/01/2011 08:33:55
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

#ifndef _BPP_GRAPH_UNODE_H_
#define _BPP_GRAPH_UNODE_H_

#include "../Clonable.h"

// From the STL
#include <vector>

namespace bpp {
  /**
   * @brief Unoriented node interface
   *
   * UNode is an interface for unoriented nodes aimed to build unoriented
   * graphs.
   *
   * In these classes we choose to use int for positions rathed than size_t
   * because negative positions are used in some implementations to distinguish
   * between two types of neighbors in operator[]. @see bpp::ONode
   *
   * @author Sylvain Gaillard
   */
  class UNode: public virtual Clonable {
    public:
      /**
       * @name Neighbors
       *
       * @{
       */

      /**
       * @brief Get a neighbor of this node in const context.
       *
       * @param pos the position of the neighbor to get.
       * @return A pointer toward the neighbor node.
       */
      virtual const UNode * getNeighbor(int pos) const = 0;

      /**
       * @brief Get a neighbor of this node.
       *
       * @param pos the position of the neighbor to get.
       * @return A pointer toward the neighbor node.
       */
      virtual UNode * getNeighbor(int pos) = 0;

      /**
       * @brief Get the degree i.e. the number of neighbors of this node.
       */
      virtual int degree() const = 0;

      /** @} */

      /**
       * @name The Clonable interface.
       *
       * @{
       */
#ifndef NO_VIRTUAL_COV      
      UNode * clone() const = 0;
#endif
      /** @} */

      /**
       * @name Operators
       *
       * @{
       */

      /**
       * @brief Direct access to a neighbor in const context.
       */
      virtual const UNode * operator[] (int i) const = 0;

      /**
       * @brief Direct access to a neighbor.
       */
      virtual UNode * operator[] (int i) = 0;

      /** @} */

  };
}

#endif //_BPP_GRAPH_UNODE_H_
