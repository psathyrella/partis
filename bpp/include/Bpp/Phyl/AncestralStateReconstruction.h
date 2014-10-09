//
// File: AncestralStateReconstruction.h
// Created by: Julien Dutheil
// Created on: Fri Jul 08 13:32 2005
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

#ifndef _ANCESTRALSTATESRECONSTRUCTION_H_
#define _ANCESTRALSTATESRECONSTRUCTION_H_

// From SeqLib:
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <vector>
#include <map>

namespace bpp
{

class Node;

/**
 * @brief Interface for ancestral states reconstruction methods.
 */
class AncestralStateReconstruction
{
  public:
    AncestralStateReconstruction() {}
    virtual ~AncestralStateReconstruction() {}

  public:
    /**
     * @brief Get ancestral states for a given node as a vector of int.
     *
     * The size of the vector depends on the implementation.
     * This method is mainly for efficient internal use in other classes.
     * Consider using the getAncestralSequenceForNode() method for a more
     * general output.
     *
     * @param nodeId the id of the node at which the states must be reconstructed.
     * @return A vector of states indices.
     * @see getAncestralSequenceForNode
     */ 
    virtual std::vector<size_t> getAncestralStatesForNode(int nodeId) const = 0;

    /**
     * @brief Get all ancestral states for all nodes.
     *
     * Call the getAncestralSequenceForNode() method on each node in the tree.
     *
     * @return A map with nodes id as key, and a vector of states indices as value.
     * @see getAncestralSequenceForNode
     */
    virtual std::map<int, std::vector<size_t> > getAllAncestralStates() const = 0;

    /**
     * @brief Get the ancestral sequence for a given node.
     *
     * @param nodeId The id of the node at which the sequence must be reconstructed.
     * @return A sequence object.
     */ 
    virtual Sequence* getAncestralSequenceForNode(int nodeId) const = 0;

    /**
     * @brief Get all the ancestral sequences for all nodes.
     *
     * @return A new SiteContainer object.
     */ 
    virtual SiteContainer* getAncestralSequences() const = 0;
    
};

} //end of namespace bpp.

#endif // _ANCESTRALSTATESRECONSTRUCTION_H_

