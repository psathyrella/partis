//
// File: Tree.h
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

#ifndef _TREE_H_
#define _TREE_H_

#include "TreeExceptions.h"

// From Utils:
#include <Bpp/Clonable.h>

// From the STL:
#include <string>
#include <vector>
#include <map>

/**
 * @mainpage
 *
 * @section tree Tree data storage and manipulation
 *
 * @par
 * The Bio++ Phylogenetics Library (PhylLib) provides classes and methods for phylogenetics and molecular evolution.
 * The bpp::Tree interface provides general methods to store and manipulate phylogenetic trees.
 * Several utilitary methods can also be found in the bpp::TreeTools static class.
 * The only implementation for now of the bpp::Tree interface is the bpp::TreeTemplate class.
 * It uses a recursive storage of bpp::Node objects (or any class inheriting from it).
 * The bpp::Node object contains methods to access the father and son nodes in the hierarchy, and several fields like a name,
 * an id or the length of the branch connected to its father. It also includes support for node and branch properties (like
 * bootstrap values) that can be attached to and manipulated together with the tree.
 * The bpp::NodeTemplate class can be used in order to extend the tree structure and add more complex
 * data to the tree. The corresponding bpp::TreeTemplateTools provide specific methods, in most cases more efficient
 * than their equivalent in the bpp::TreeTools.
 *
 * @par
 * Trees can also be read and written from/to files, using the bpp::Newick class.
 *
 *
 * @section reconstruction Phylogenetic reconstruction methods
 *
 * @par
 * PhylLib provides tools to reconstruct phylogenies from sequence data, using maximum parsimony, distance-based methods 
 * and maximum likelihood, all of them implemented in an object-oriented way, and hence involving several classes.
 *
 * @par Maximum parcimony
 * See bpp::TreeParsimonyScore for parsimony score computation. Only a Nearest Neighbor Interchange (NNI) algorithm for
 * topology estimation is provided for now, see bpp::NNISearchable, bpp::NNITopologySearch and bpp::OptimizationTools for
 * more user-friendly methods.
 *
 * @par Distance methods
 * The bpp::DistanceEstimation class allows you to compute pairwise distances from a large set of models (see next section),
 * and store them as a bpp::DistanceMatrix. This matrix is the input of any distance-based method.
 * The (U/W)PGMA (bpp::PGMA), neighbor-joining (bpp::NeighborJoining) and BioNJ (bpp::BioNJ) methods are implemented.
 *
 * @par Maximum likelihood methods
 * Use a model to describe the evolutionary process, among many available (see next section).
 * Support for homogeneous (reversible or not) and non-homogeneous models is provided. Several likelihood computation
 * algorithms are provided, depending on the final usage. All classes are instances of the bpp::TreeLikelihood interface.
 * The bpp::DiscreteRatesAcrossSitesTreeLikelihood interface adds support for rate heterogeneity across sites.
 * - The bpp::RHomogeneousTreeLikelihood class is the most simple implementation. It uses Felsenstein's recursion.
 *   The possibility of compressing sites save comptation time and memory.
 * - The bpp::DRHomogeneousTreeLikelihood class is similar to the previous one, but uses a double-recursive algorithm.
 *   It is more CPU expensive when computing likelihoods, and uses more memory. Computation of branch length analytical 
 *   derivatives is nonetheless faster, since they do not involve any additionnal recursion.
 *   You also have to use this class in order to perform substitution mapping (bpp::SubstitutionMappingTools) or reconstruct
 *   ancestral sequences (bpp::AncestralStateReconstruction).
 * - The bpp::NNIHomogeneousTreeLikelihood class inherits from bpp::DRHomogeneousTreeLikelihood, and implements the bpp::NNISearchable
 *   interface. This class should hence be used in order to optimize the tree topology.
 * - The bpp::RNonHomogeneousTreeLikelihood and bpp::DRNonHomogeneousTreeLikelihood are similar to their homogeneous homologues,
 *   but are designed for non-reversible or non-homogeneous models of substitution.
 * - Finally, the bpp::ClockTreeLikelihood interface uses a different parametrization by assuming a global molecular clock.
 *   It is implemented in the bpp::RHomogeneousClockTreeLikelihood class, which inherits from bpp::RHomogeneousTreeLikelihood.
 *
 * @par 
 * The bpp::TreeLikelihood class inherits from the bpp::Function interface, which means that any optimization method from
 * the NumCalc library can be used to estimate numerical parameters. The bpp::OptimizationTools static class provides
 * general methods with predefined options, including for topology estimation.
 *
 *
 * @section models Evolutionary models
 *
 * @par
 * The Bio++ phylogenetic library provides different kinds of models.
 * Substitution models are provided via the bpp::SubstitutionModel interface.
 * All commonly used models for nucleotides and proteins are provided (see for instance bpp::JCnuc, bpp::K80, bpp::GTR, bpp::JTT92, etc.).
 * You can add your own model by implementing the bpp::SubstitutionModel interface.
 * Rate across sites (RAS) models are integrated thanks to the bpp::DiscreteDistribution interface, providing support for the gamma
 * (bpp::GammaDiscreteDistribution) and gamma+invariant (bpp::InvariantMixedDiscreteDistribution) rate distributions.
 * Here again, this is very easy to add support for new rate distributions by implementing the corresponding interface.
 * 
 * @par
 * Markov-modulated Markov models (of which the covarion model is a particular case) are included via the bpp::MarkovModulatedSubstitutionModel interface,
 * and its implementation bpp::G2001 and bpp::TS98.
 *
 * @par
 * Finally from version 1.5, it is possible to build virtually any kind of non-homogeneous model thanks to the bpp::SubstitutionModelSet class.
 *
 *
 * @section more And more...
 *
 * @par
 * PhylLib allows you to perform a lot of analysis, like evolutionary rate estimation, tree consensus, etc.
 *
 */

namespace bpp
{

  /**
   * @brief Interface for phylogenetic tree objects. 
   */
  class Tree:
    public virtual Clonable
  {

  public: // Constructors and destructor:
		
    Tree() {}
    virtual ~Tree() {}

    Tree* clone() const = 0;

    /**
     * @brief clones a Subtree rooted at given node Id
     *
     **/
    
    virtual Tree* cloneSubtree(int newRootId) const = 0;
  public:
		
    /**
     * @brief Tree name.
     *
     * @{
     */
    virtual std::string getName() const = 0;
		
    virtual void setName(const std::string& name) = 0;
    /** @} */
		
    virtual size_t getNumberOfLeaves() const = 0;
		
    virtual size_t getNumberOfNodes() const = 0;

    virtual std::vector<double> getBranchLengths() const = 0;

    virtual std::vector<std::string> getLeavesNames() const = 0;

    /**
     * @name Retrieving ids.
     *
     * @{
     */
    virtual int getRootId() const = 0;

    virtual int getLeafId(const std::string& name) const throw (NodeNotFoundException)= 0;
	
    virtual std::vector<int> getLeavesId() const = 0;

    virtual std::vector<int> getNodesId() const = 0;
		
    virtual std::vector<int> getInnerNodesId() const = 0;
		
    virtual std::vector<int> getBranchesId() const = 0;

    virtual std::vector<int> getSonsId(int parentId) const throw (NodeNotFoundException) = 0;

    virtual std::vector<int> getAncestorsId(int nodeId) const throw (NodeNotFoundException) = 0;

    virtual int getFatherId(int parentId) const throw (NodeNotFoundException) = 0;

    virtual bool hasFather(int nodeId) const throw (NodeNotFoundException) = 0;

/** @} */

    /**
     * @name Dealing with node names.
     *
     * @{
     */
    virtual std::string getNodeName(int nodeId) const throw (NodeNotFoundException) = 0;
		
    virtual void setNodeName(int nodeId, const std::string& name) throw (NodeNotFoundException) = 0;
		
    virtual void deleteNodeName(int nodeId) throw (NodeNotFoundException) = 0;
		
    virtual bool hasNodeName(int nodeId) const throw (NodeNotFoundException) = 0;
    /** @} */
		
    /**
     * @name Several tests.
     *
     * @{
     */
    virtual bool hasNode(int nodeId) const = 0;

    virtual bool isLeaf(int nodeId) const throw (NodeNotFoundException) = 0;

    virtual bool isRoot(int nodeId) const throw (NodeNotFoundException) = 0;
    /** @} */

    /**
     * @name Acting on topology.
     *
     * @{
     */

    /**
     * @brief Swap two son nodes.
     *
     * @param tree The tree.
     * @param nodeId The node.
     * @param i1 First son node index.
     * @param i2 Second son node index.
     * @throw NodeNotFoundException If the node is not found.
     * @throw IndexOutOfBoundsException If one node index is not valid, or if the node
     */
    void swapNodes(const Tree& tree, int nodeId, size_t i1 = 0, size_t i2 = 1) throw (NodeNotFoundException,IndexOutOfBoundsException);
  
    /** @} */

    /**
     * @name Dealing with branch lengths.
     *
     * @{
     */
    virtual double getDistanceToFather(int nodeId) const = 0;
		
    virtual void setDistanceToFather(int nodeId, double length) = 0;
		
    virtual void deleteDistanceToFather(int nodeId) = 0;
		
    virtual bool hasDistanceToFather(int nodeId) const = 0;
    /** @} */

    /**
     * @name Node properties.
     *
     * @{
     */
    virtual bool hasNodeProperty(int nodeId, const std::string& name) const throw (NodeNotFoundException) = 0;
		
    virtual void setNodeProperty(int nodeId, const std::string& name, const Clonable& property) throw (NodeNotFoundException) = 0;
				
    virtual Clonable* getNodeProperty(int nodeId, const std::string& name) throw (NodeNotFoundException) = 0;
				
    virtual const Clonable* getNodeProperty(int nodeId, const std::string& name) const throw (NodeNotFoundException) = 0;
				
    virtual Clonable* removeNodeProperty(int nodeId, const std::string& name) throw (NodeNotFoundException) = 0;

    virtual std::vector<std::string> getNodePropertyNames(int nodeId) const throw (NodeNotFoundException) = 0;
    /** @} */
		
    /**
     * @name Branch properties.
     *
     * @{
     */
    virtual bool hasBranchProperty(int nodeId, const std::string& name) const throw (NodeNotFoundException) = 0;
		
    virtual void setBranchProperty(int nodeId, const std::string& name, const Clonable & property) throw (NodeNotFoundException) = 0;
				
    virtual Clonable* getBranchProperty(int nodeId, const std::string& name) throw (NodeNotFoundException) = 0;
				
    virtual const Clonable* getBranchProperty(int nodeId, const std::string& name) const throw (NodeNotFoundException) = 0;
				
    virtual Clonable* removeBranchProperty(int nodeId, const std::string& name) throw (NodeNotFoundException) = 0;

    virtual std::vector<std::string> getBranchPropertyNames(int nodeId) const throw (NodeNotFoundException) = 0;
    /** @} */

    /**
     * @brief Change the root node.
     *
     * Works on unrooted tree.
     * If the tree is rooted, the method unroots it first.
     *
     * @param nodeId The id of the node that will be the new root.
     */
    virtual void rootAt(int nodeId) throw (NodeNotFoundException) = 0;

    /**
     * @brief Root a tree by specifying an outgroup.
     *
     * If the tree is rooted, unroot it first, change the root node and then
     * reroot the tree using the previous root id.
     * If the tree is unrooted, change the root node and then create a new root node.
     *
     * @param nodeId The id of the node that will be the new root.
     */
    virtual void newOutGroup(int nodeId) throw (NodeNotFoundException) = 0;
		
    /**
     * @brief Tell if the tree is rooted.
     * 
     * @return True if the tree is rooted.
     */
    virtual bool isRooted() const = 0;
		
    /**
     * @brief Unroot a rooted tree.
     *
     * @return True if the tree has been unrooted.
     * @throw UnrootedTreeException If the tree is already rooted.
     */
    virtual bool unroot() throw (UnrootedTreeException) = 0;

    /**
     * @brief Number nodes.
     */
    virtual void resetNodesId() = 0;
		
    // Works on (multi)furcations:
		
    /**
     * @brief Tell if the tree is multifurcating.
     * 
     * @return True if the tree is multifurcating.
     */
    virtual bool isMultifurcating() const = 0;
		
    /**
     * @brief Get all the branch lengths of a tree.
     *
     * @return A vector with all branch lengths.
     * @throw NodeException If a branch length is lacking.
     */
    virtual std::vector<double> getBranchLengths() throw (NodeException) = 0;

    /**
     * @brief Get the total length (sum of all branch lengths) of a tree.
     *
     * @return The total length of the subtree.
     * @throw NodeException If a branch length is lacking.
     */
    virtual double getTotalLength() throw (NodeException) = 0;

    /**
     * @brief Set all the branch lengths of a tree.
     *
     * @param brLen The branch length to apply.
     */
    virtual void setBranchLengths(double brLen) = 0;
		
    /**
     * @brief Give a length to branches that don't have one in a tree.
     *
     * @param brLen The branch length to apply.
     */
    virtual void setVoidBranchLengths(double brLen) = 0;
	
    /**
     * @brief Scale a given tree.
     *
     * Multiply all branch lengths by a given factor.
     *
     * @param factor The factor to multiply all branch lengths with.
     * @throw NodeException If a branch length is lacking.
     */
    virtual void scaleTree(double factor) throw (NodeException) = 0;

    /**
     * @brief Get an id.
     *
     * @return an unused node id.
     */
    virtual int getNextId() = 0;

  };

} //end of namespace bpp.

#endif	//_TREE_H_

