//
// File: DRASRTreeLikelihoodData.h
// Created by: Julien Dutheil
// Created on: Sat Dec 30 14:20 2006
// From file HomogeneousTreeLikelihood.h
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

#ifndef _DRASRHOMOGENEOUSTREELIKELIHOODDATA_H_
#define _DRASRHOMOGENEOUSTREELIKELIHOODDATA_H_

#include "AbstractTreeLikelihoodData.h"
#include "../Model/SubstitutionModel.h"
#include "../SitePatterns.h"

#include <Bpp/Numeric/VectorTools.h>

// From the STL:
#include <map>

namespace bpp
{

/**
 * @brief Likelihood data structure for a node.
 * 
 * This class is for use with the DRASRTreeParsimonyData class.
 * 
 * Store all conditionnal likelihoods:
 * <pre>
 * x[i][c][s]
 *   |---------> Site i
 *      |------> Rate class c
 *         |---> Ancestral state s
 * </pre> 
 * We call this the <i>likelihood array</i> for each node.
 * In the same way, we store first and second order derivatives.
 *
 * @see DRASRTreeLikelihoodData
 */
class DRASRTreeLikelihoodNodeData :
  public virtual TreeLikelihoodNodeData
{
  private:
    mutable VVVdouble nodeLikelihoods_;
    mutable VVVdouble nodeDLikelihoods_;
    mutable VVVdouble nodeD2Likelihoods_;
    const Node* node_;

  public:
    DRASRTreeLikelihoodNodeData() : nodeLikelihoods_(), nodeDLikelihoods_(), nodeD2Likelihoods_(), node_(0) {}
    
    DRASRTreeLikelihoodNodeData(const DRASRTreeLikelihoodNodeData& data) :
      nodeLikelihoods_(data.nodeLikelihoods_),
      nodeDLikelihoods_(data.nodeDLikelihoods_),
      nodeD2Likelihoods_(data.nodeD2Likelihoods_),
      node_(data.node_)
    {}
    
    DRASRTreeLikelihoodNodeData& operator=(const DRASRTreeLikelihoodNodeData& data)
    {
      nodeLikelihoods_   = data.nodeLikelihoods_;
      nodeDLikelihoods_  = data.nodeDLikelihoods_;
      nodeD2Likelihoods_ = data.nodeD2Likelihoods_;
      node_              = data.node_;
      return *this;
    }
 
#ifndef NO_VIRTUAL_COV
    DRASRTreeLikelihoodNodeData*
#else
    Clonable*
#endif
    clone() const
    {
      return new DRASRTreeLikelihoodNodeData(*this);
    }

  public:
    const Node* getNode() const { return node_; }
    void setNode(const Node* node) { node_ = node; }

    VVVdouble& getLikelihoodArray() { return nodeLikelihoods_; }
    const VVVdouble& getLikelihoodArray() const { return nodeLikelihoods_; }
    
    VVVdouble& getDLikelihoodArray() { return nodeDLikelihoods_; }
    const VVVdouble& getDLikelihoodArray() const { return nodeDLikelihoods_; }

    VVVdouble& getD2LikelihoodArray() { return nodeD2Likelihoods_; }
    const VVVdouble& getD2LikelihoodArray() const { return nodeD2Likelihoods_; }
};

/**
 * @brief discrete Rate Across Sites, (simple) Recursive likelihood data structure.
 */
class DRASRTreeLikelihoodData :
  public virtual AbstractTreeLikelihoodData
{
  private:
    /**
     * @brief This contains all likelihood values used for computation.
     *
     */
    mutable std::map<int, DRASRTreeLikelihoodNodeData> nodeData_;
      
    /**
     * @brief This map defines the pattern network.
     *
     * Let n1 be the id of a node in the tree, and n11 and n12 the ids of its sons.
     * Providing the likelihood array is known for nodes n11 and n12,
     * the likelihood array for node n1 and site <i>i</i> (_likelihood[n1][i]) must be computed  
     * using arrays patternLinks_[n1][n11][i] and patternLinks_[n1][n12][i].
     * This network is intialized once for all in the constructor of this class.
     *
     * The double map contains the position of the site to use (second dimension)
     * of the likelihoods array.
     */
    mutable std::map<int, std::map<int, std::vector<size_t> > > patternLinks_;
    SiteContainer* shrunkData_;
    size_t nbSites_; 
    size_t nbStates_;
    size_t nbClasses_;
    size_t nbDistinctSites_; 
    bool usePatterns_;

  public:
    DRASRTreeLikelihoodData(const TreeTemplate<Node>* tree, size_t nbClasses, bool usePatterns = true) :
      AbstractTreeLikelihoodData(tree),
      nodeData_(), patternLinks_(), shrunkData_(0), nbSites_(0), nbStates_(0),
      nbClasses_(nbClasses), nbDistinctSites_(0), usePatterns_(usePatterns)
    {}

    DRASRTreeLikelihoodData(const DRASRTreeLikelihoodData& data):
      AbstractTreeLikelihoodData(data),
      nodeData_(data.nodeData_),
      patternLinks_(data.patternLinks_),
      shrunkData_(0),
      nbSites_(data.nbSites_), nbStates_(data.nbStates_),
      nbClasses_(data.nbClasses_), nbDistinctSites_(data.nbDistinctSites_),
      usePatterns_(data.usePatterns_)
    {
      if (data.shrunkData_)
        shrunkData_      = dynamic_cast<SiteContainer *>(data.shrunkData_->clone());
    }

    DRASRTreeLikelihoodData& operator=(const DRASRTreeLikelihoodData & data)
    {
      AbstractTreeLikelihoodData::operator=(data);
      nodeData_          = data.nodeData_;
      patternLinks_      = data.patternLinks_;
      nbSites_           = data.nbSites_;
      nbStates_          = data.nbStates_;
      nbClasses_         = data.nbClasses_;
      nbDistinctSites_   = data.nbDistinctSites_;
      if (shrunkData_) delete shrunkData_;
      if (data.shrunkData_)
        shrunkData_      = dynamic_cast<SiteContainer*>(data.shrunkData_->clone());
      else
        shrunkData_      = 0;
      usePatterns_       = data.usePatterns_;
      return *this;
    }

    virtual ~DRASRTreeLikelihoodData() { delete shrunkData_; }

    DRASRTreeLikelihoodData* clone() const { return new DRASRTreeLikelihoodData(*this); }

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
      tree_ = tree;
      for (std::map<int, DRASRTreeLikelihoodNodeData>::iterator it = nodeData_.begin(); it != nodeData_.end(); it++)
      {
        int id = it->second.getNode()->getId();
        it->second.setNode(tree_->getNode(id));
      }
    }
 
    DRASRTreeLikelihoodNodeData& getNodeData(int nodeId)
    { 
      return nodeData_[nodeId];
    }
    const DRASRTreeLikelihoodNodeData& getNodeData(int nodeId) const
    { 
      return nodeData_[nodeId];
    }
    size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const
    {
      return patternLinks_[parentId][sonId][currentPosition];
    }
    size_t getRootArrayPosition(size_t currentPosition) const
    {
      return rootPatternLinks_[currentPosition];
    }
    const std::vector<size_t>& getArrayPositions(int parentId, int sonId) const
    {
      return patternLinks_[parentId][sonId];
    }
    std::vector<size_t>& getArrayPositions(int parentId, int sonId)
    {
      return patternLinks_[parentId][sonId];
    }
    size_t getArrayPosition(int parentId, int sonId, size_t currentPosition)
    {
      return patternLinks_[parentId][sonId][currentPosition];
    }

    VVVdouble& getLikelihoodArray(int nodeId)
    {
      return nodeData_[nodeId].getLikelihoodArray();
    }
    
    VVVdouble& getDLikelihoodArray(int nodeId)
    {
      return nodeData_[nodeId].getDLikelihoodArray();
    }
    
    VVVdouble& getD2LikelihoodArray(int nodeId)
    {
      return nodeData_[nodeId].getD2LikelihoodArray();
    }

    size_t getNumberOfDistinctSites() const { return nbDistinctSites_; }
    size_t getNumberOfSites() const { return nbSites_; }
    size_t getNumberOfStates() const { return nbStates_; }
    size_t getNumberOfClasses() const { return nbClasses_; }
    
    void initLikelihoods(const SiteContainer& sites, const SubstitutionModel& model) throw (Exception);

  protected:
    /**
     * @brief This method initializes the leaves according to a sequence file.
     * likelihood is set to 1 for the state corresponding to the sequence site,
     * otherwise it is set to 0.
     *
     * All likelihood arrays at each nodes are initialized according to alphabet
     * size and sequences length, and filled with 1.
     *
     * NB: This method is recursive.
     *
     * @param node      The node defining the subtree to analyse.
     * @param sequences The data to be used for initialization.
     * @param model     The model to use.
     */
    virtual void initLikelihoods(const Node* node, const SiteContainer& sequences, const SubstitutionModel& model) throw (Exception);

    /**
     * @brief This method initializes the leaves according to a sequence file.
     *
     * likelihood is set to 1 for the state corresponding to the sequence site,
     * otherwise it is set to 0.
     *
     * All likelihood arrays at each nodes are initialized according to alphabet
     * size and sequences length, and filled with 1.
     *
     * NB: This method is recursive.
     *
     * @param node      The node defining the subtree to analyse.
     * @param sequences The data to be used for initialization.
     * @param model     The model to use.
     * @return The shrunk sub-dataset + indices for the subtree defined by <i>node</i>.
     */
    virtual SitePatterns* initLikelihoodsWithPatterns(const Node* node, const SiteContainer& sequences, const SubstitutionModel& model) throw (Exception);
  
};

} //end of namespace bpp.

#endif //_DRASRHOMOGENEOUSTREELIKELIHOODDATA_H_

