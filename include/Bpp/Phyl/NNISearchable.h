//
// File: NNISearchable.h
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Thu Jul 28 16:32 2005
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

#ifndef _NNISEARCHABLE_H_
#define _NNISEARCHABLE_H_

#include "Node.h"
#include "TreeTemplate.h"
#include "TopologySearch.h"

namespace bpp
{

/**
 * @brief Interface for Nearest Neighbor Interchanges algorithms.
 *
 * This interface defines the methods to work with NNI algorithms.
 * 
 * NNISearchable objects are supposed to work with TreeTemplate objects.
 * NNI are defined as follow:
 * <pre>
 * ------------->
 *     +------- C
 *     |
 * D --+ X  +-- B
 *     |    |
 *     +----+ F
 *          |
 *          +-- A
 * </pre>
 * Where:
 * -F is the focus (parent) node,
 * -A and B are the son of F
 * -X is the parent of F and so on.
 * Two NNI's are possible for branch (XF):
 * - swaping B and C, which is the same as D and A
 * - swaping A and C, which is the same as D and B
 * Because of the rooted representation, we'll consider B \f$\leftrightarrow\f$ C and A \f$\leftrightarrow\f$ C,
 * which are simpler to perform.
 *
 * For unrooted tree, we have at the 'root' node:
 * <pre>
 * ------------->
 *   +------- D
 *   |
 * X +------- C
 *   |
 *   |    +-- B
 *   |    |
 *   +----+ F
 *        |
 *        +-- A
 * </pre>
 * In this case, we swap A or B with one of C or D.
 * Which one of C or D depends on the implementation, but will always be the same,
 * so that swapping A or B involve 2 distinct NNI.
 * 
 */
class NNISearchable:
  public TopologyListener,
  public virtual Clonable
{
	public:
		NNISearchable() {}
		virtual ~NNISearchable() {}

#ifndef NO_VIRTUAL_COV
    virtual NNISearchable* clone() const = 0;
#endif

	public:

		/**
		 * @brief Send the score of a NNI movement, without performing it.
		 *
		 * This methods sends the score variation.
		 * This variation must be negative if the new point is better,
		 * i.e. the object is to be used with a minimizing optimization
		 * (for consistence with Optimizer objects).
		 *
		 * @param nodeId The id of the node defining the NNI movement.
		 * @return The score variation of the NNI.
		 * @throw NodeException If the node does not define a valid NNI.
		 */
		virtual double testNNI(int nodeId) const throw (NodeException) = 0;

		/**
		 * @brief Perform a NNI movement.
		 *
		 * @param nodeId The id of the node defining the NNI movement.
		 * @throw NodeException If the node does not define a valid NNI.
		 */
		virtual void doNNI(int nodeId) throw (NodeException) = 0;

		/**
		 * @brief Get the tree associated to this NNISearchable object.
		 *
		 * @return The tree associated to this instance.
		 */
		virtual const Tree& getTopology() const = 0;
    
    /**
     * @brief Get the current score of this NNISearchable object.
     *
     * @return The current score of this instance.
     */
    virtual double getTopologyValue() const throw (Exception) = 0;

};

} //end of namespace bpp.

#endif //_NNISEARCHABLE_H_

