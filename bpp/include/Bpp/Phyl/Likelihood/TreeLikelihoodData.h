//
// File: TreeLikelihoodData.h
// Created by: Julien Dutheil
// Created on: Sat Dec 30 12:48 2006
// From file AbstractTreeLikelihood.h
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

#ifndef _TREELIKELIHOODDATA_H_
#define _TREELIKELIHOODDATA_H_

#include "../Node.h"
#include "../TreeTemplate.h"

//From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{

/**
 * @brief TreeLikelihood partial data structure.
 *
 * Stores inner computation for a given node.
 *
 * @see TreeLikelihoodData
 */
class TreeLikelihoodNodeData:
  public virtual Clonable
{
  public:
    TreeLikelihoodNodeData() {}
    virtual ~TreeLikelihoodNodeData() {}

#ifndef NO_VIRTUAL_COV
    TreeLikelihoodNodeData*
#else
    Clonable*
#endif
    clone() const = 0;

  public:
    /**
     * @brief Get the node associated to this data structure.
     *
     * @return The node associated to this structure.
     */
    virtual const Node* getNode() const = 0;

    /**
     * @brief Set the node associated to this data
     *
     * A pointer toward this node will be created and associated to this data.
     *
     * @param node The node to be associated to this data.
     */
    virtual void setNode(const Node* node) = 0;
};

/**
 * @brief TreeLikelihood data structure.
 *
 * Stores all the inner computations:
 * - conditionnal likelihoods for each node,
 * - correspondance between sites in the dataset and array indices.
 * 
 * @see TreeLikelihoodNodeData
 */
class TreeLikelihoodData:
  public virtual Clonable
{
  public:
    TreeLikelihoodData() {}
    virtual ~TreeLikelihoodData() {}
    
#ifndef NO_VIRTUAL_COV
    TreeLikelihoodData* clone() const = 0;
#endif

  public:
    virtual const Alphabet* getAlphabet() const = 0;
    virtual const TreeTemplate<Node>* getTree() const = 0;
    virtual size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const = 0;
    virtual size_t getRootArrayPosition(size_t site) const = 0;
    virtual       TreeLikelihoodNodeData& getNodeData(int nodeId) = 0;
    virtual const TreeLikelihoodNodeData& getNodeData(int nodeId) const = 0;

    /**
     * @return The number of non redundant patterns.
     */
    virtual size_t getNumberOfDistinctSites() const = 0;
    
    /**
     * @return The total number of sites.
     */
    virtual size_t getNumberOfSites() const = 0;
    
    /**
     * @return Get the number of states used in the model.
     */
    virtual size_t getNumberOfStates() const = 0;

    /**
     * @return The frequency of a given pattern.
     */
    virtual unsigned int getWeight(size_t pos) const = 0;
    
    /**
     * @return Frequencies for each pattern.
     */
    virtual const std::vector<unsigned int>& getWeights() const = 0;

};

} //end of namespace bpp.

#endif //_TREELIKELIHOODDATA_H_

