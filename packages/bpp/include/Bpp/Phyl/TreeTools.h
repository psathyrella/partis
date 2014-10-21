//
// File: TreeTools.h
// Created by:  Julien Dutheil
// Created on: Wed Aug  6 13:45:28 2003
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

#ifndef _TREETOOLS_H_
#define _TREETOOLS_H_

#include "TreeExceptions.h"
#include "Node.h"
#include "Tree.h"
#include "BipartitionList.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/VectorTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/DistanceMatrix.h>

namespace bpp
{

/**
 * @brief Generic utilitary methods dealing with trees.
 *
 * These methods work with all Tree object.
 * However, depending on the tree implementation, they may not be the most efficient.
 *
 * @see TreeTemplateTools
 */
class TreeTools
{
  public:
    TreeTools() {}
    virtual ~TreeTools() {}
  
  public:

    /**
     * @name Retrieve topology information
     *
     * @{
     */
    
    /**
     * @brief Retrieve all leaves from a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @return A vector with the ids of all leaves in the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static std::vector<int> getLeavesId(const Tree& tree, int nodeId) throw (NodeNotFoundException);

    /**
     * @brief Retrieve all leaves from a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @param leaves A vector with the ids of all leaves in the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void getLeavesId(const Tree& tree, int nodeId, std::vector<int>& leaves) throw (NodeNotFoundException);
 
    /**
     * @brief Count the number of leaves from a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static size_t getNumberOfLeaves(const Tree& tree, int nodeId) throw (NodeNotFoundException);
 
    /**
     * @brief Get the id of a leaf given its name in a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @param name The name of the node.
     * @return The id of the node.
     * @throw NodeNotFoundException If the node is not found.
     */
    static int getLeafId(const Tree& tree, int nodeId, const std::string& name) throw (NodeNotFoundException);

    /**
     * @brief Get the id of a leaf given its name in a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @param name The name of the node.
     * @param id The id of the node.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void searchLeaf(const Tree& tree, int nodeId, const std::string& name, int*& id) throw (NodeNotFoundException);

    /**
     * @brief Get a vector of ancestor nodes between to nodes.
     *
     * @param tree The tree to use.
     * @param nodeId1 Id of first node.
     * @param nodeId2 Id of second node.
     * @param includeAncestor Tell if the common ancestor must be included in the vector.
     * @return A vector of ancestor nodes ids.
     * @throw NodeNotFoundException If the node is not found.
     */
    static std::vector<int> getPathBetweenAnyTwoNodes(const Tree& tree, int nodeId1, int nodeId2, bool includeAncestor = true) throw (NodeNotFoundException);
 
    /**
     * @brief Get a list of all ids of parents nodes, from the current node (not included) to the root of the tree.
     *
     * @param tree The tree to use.
     * @param nodeId The id of node defining the subtree.
     * @return The list of ancestors ids.
     * @throw NodeNotFoundException If the node is not found.
     */
    static std::vector<int> getAncestors(const Tree& tree, int nodeId) throw (NodeNotFoundException);

    /**
     * @brief Get the id of the last common ancestors of all specified nodes.
     *
     * Nodes id need not correspond to leaves.
     *
     * @author Simon Carrignon
     * @param tree The tree to use.
     * @param nodeIds The ids of the input nodes.
     * @throw NodeNotFoundException If at least of of input node is not found.
     */
    static int getLastCommonAncestor(const Tree& tree, const std::vector<int>& nodeIds) throw (NodeNotFoundException, Exception);

    /**
     * @brief Get the depth of the subtree defined by node 'node', i.e. the maximum
     * number of sons 'generations'.
     *
     * ex:
     * @code
     *    +----------A
     *    |
     * ---+ N1     +-------B
     *    |        |
     *    +--------+ N2
     *             |
     *             +------C
     * @endcode
     * Depth of node 'N1' id 2, depth of node 'N2' is 1, depth of leaves is 0.
     *
     * @param tree The tree.
     * @param nodeId The id of node defining the subtree.
     * @return The depth of the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static size_t getDepth(const Tree& tree, int nodeId) throw (NodeNotFoundException);

    /**
     * @brief Get the depths for all nodes of the subtree defined by node 'node', i.e. the maximum
     * number of sons 'generations'.
     *
     * ex:
     * @verbatim
     *    +----------A
     *    |
     * ---+ N1     +-------B
     *    |        |
     *    +--------+ N2
     *             |
     *             +------C
     * @endverbatim
     * Depth of node 'N1' id 2, depth of node 'N2' is 1, depth of leaves is 0.
     *
     * @param tree The tree.
     * @param nodeId The id of node defining the subtree.
     * @param depths The map that will contain all the depths of the nodes, with node ids as keys.
     * @return The depth of the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static size_t getDepths(const Tree& tree, int nodeId, std::map<int, size_t>& depths) throw (NodeNotFoundException);

    /**
     * @brief Get the height of the subtree defined by node 'node', i.e. the maximum
     * distance between leaves and the root of the subtree.
     *
     * The distance do not include the branch length of the subtree root node.
     * The height of a leaf is hence 0.
     *
     * @param tree The tree.
     * @param nodeId The id of node defining the subtree.
     * @return The height of the subtree.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If a branch length is lacking.
     */ 
    static double getHeight(const Tree& tree, int nodeId) throw (NodeNotFoundException,NodeException);

    /**
     * @brief Get the heights of all nodes within a subtree defined by node 'node', i.e. the maximum
     * distance between leaves and the root of the subtree.
     *
     * The height of a leaf is 0.
     *
     * @param tree The tree.
     * @param nodeId The id of node defining the subtree.
     * @param heights The map that will contain all the heights of the nodes, with node ids as keys.
     * @return The height of the subtree.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If a branch length is lacking.
     */ 
    static double getHeights(const Tree& tree, int nodeId, std::map<int, double>& heights) throw (NodeNotFoundException,NodeException);

    /** @} */

    /**
     * @name Act on branch lengths.
     *
     * @{
     */
    
     /**
     * @brief Get all the branch lengths of a subtree.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @return A vector with all branch lengths.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If a branch length is lacking.
     */
    static Vdouble getBranchLengths(const Tree& tree, int nodeId) throw (NodeNotFoundException,NodeException);
    
    /**
     * @brief Get the total length (sum of all branch lengths) of a subtree.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @param includeAncestor Tell if the branch length of the most ancient node should be included in the counting.
     * (this should be set to false if this node is the root of the tree for instance).
     * @return The total length of the subtree.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If a branch length is lacking.
     */
    static double getTotalLength(const Tree& tree, int nodeId, bool includeAncestor = true) throw (NodeNotFoundException,NodeException);

    /**
     * @brief Set all the branch lengths of a subtree.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @param brLen The branch length to apply.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void setBranchLengths(Tree& tree, int nodeId, double brLen) throw (NodeNotFoundException);
          
    /**
     * @brief Give a length to branches that don't have one in a subtree.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @param brLen The branch length to apply.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void setVoidBranchLengths(Tree& tree, int nodeId, double brLen) throw (NodeNotFoundException);
        
    /**
     * @brief Scale a given tree.
     *
     * Multiply all branch lengths by a given factor.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @param factor The factor to multiply all branch lengths with.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If a branch length is lacking.
     */
    static void scaleTree(Tree& tree, int nodeId, double factor) throw (NodeNotFoundException,NodeException);

    /**
     * @brief Grafen's method to initialize branch lengths.
     *
     * Each height of the node (total distance from the leaves) is set equal to the number of
     * leaf nodes for the corresponding subtrees - 1 for inner nodes, 0 for leaves.
     * 
     * If the tree already has branch lengths, they will be ignored.
     * 
     * Reference:
     * Grafen A. The phylogenetic regression. Philos Trans R Soc Lond B Biol Sci. 1989; 326(1233):119-57
     * 
     * @param tree The tree.
     */ 
    static void initBranchLengthsGrafen(Tree& tree);
    
    /**
     * @brief Compute branch lengths using Grafen's method.
     *
     * The 'height' of each node is devided by the total height of the tree, and the ratio is raised at power 'rho'.
     * A value of rho=0 hence returns a star tree.
     *
     * Reference:
     * Grafen A. The phylogenetic regression. Philos Trans R Soc Lond B Biol Sci. 1989; 326(1233):119-57
     * 
     * @param tree The tree to use.
     * @param power The rho parameter.
     * @param init Tell if the height must be initialized by calling the initBranchLengthsGrafen() method.
     *             Otherwise use branch lengths.
     * @throw NodeException If init=false and one branch length is lacking.
     */
    static void computeBranchLengthsGrafen(Tree& tree, double power = 1, bool init = true) throw (NodeException);
   
  private:
    static size_t initBranchLengthsGrafen(Tree& tree, int nodeId) throw (NodeNotFoundException);
    static void computeBranchLengthsGrafen(Tree& tree, int nodeId, double power, double total, double& height, double& heightRaised) throw (NodeNotFoundException,NodeException);

  public:
    /**
     * @brief Modify a tree's branch lengths to make a clock tree, by rebalancing branch lengths.
     *
     * The height of each node is set to the mean height of all son nodes.
     * This may however lead to negative branch lengths, since the mean heigth
     * may be inferior to one of the son heights, due to short branch lengths.
     * If the 'noneg' is set to yes, the mean height is checked against all son
     * heights. If it is inferior to one of the son heights, the maximum son
     * height is used instead. This results in a multifurcation.
     * 
     * This method is recursive and will be applied on all sons nodes.
     * 
     * @param tree The tree to use.
     * @param nodeId The node defining the subtree.
     * @param noneg Tell if the correction for non negative branch lengths must be used.
     * @return The modified height of the node.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If one branch length is lacking.
     */
    static double convertToClockTree(Tree& tree, int nodeId, bool noneg = false);
    
    /**
     * @brief Modify a tree's branch lengths to make a clock tree, by rescaling subtrees.
     *
     * The height of each node is set to the mean height of all son nodes.
     * All branch lengths of the corresponding subtrees are updated proportionally.
     * This algorithm is smaller than the convertToClockTree method, but may be more accurate.
     * 
     * This method is recursive and will be applied on all sons nodes.
     * 
     * @param tree The tree to use.
     * @param nodeId The node defining the subtree.
     * @return The modified height of the node.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If one branch length is lacking.
     */
    static double convertToClockTree2(Tree& tree, int nodeId);
 
    /**
     * @brief Get the total distance between two nodes.
     *
     * Sum all branch lengths between two nodes.
     *
     * @param tree The tree to consider.
     * @param nodeId1 First node id.
     * @param nodeId2 Second node id.
     * @return The sum of all branch lengths between the two nodes.
     * @throw NodeNotFoundException If the node is not found.
     */
    static double getDistanceBetweenAnyTwoNodes(const Tree& tree, int nodeId1, int nodeId2);
    
    /**
     * @brief Compute a distance matrix from a tree.
     *
     * Compute all distances between each leaves and store them in a matrix.
     * A new DistanceMatrix object is created, and a pointer toward it is returned.
     * The destruction of this matrix is left up to the user.
     *
     * @see getDistanceBetweenAnyTwoNodes
     *
     * @param tree The tree to use.
     * @return The distance matrix computed from tree.
     */
    static DistanceMatrix* getDistanceMatrix(const Tree& tree);

    /**
     * @brief (Re)root the tree using the midpoint method.
     *
     * This methods compute the pairwise distance matrix from the tree and get the maximum distance.
     * The root is then set on the branch located at half this distance.
     *
     * @param tree The tree to (re)root.
     * @deprecated Use TreeTemplateTools::midRoot instead!
     */
    static void midpointRooting(Tree& tree);
    /** @} */


    /**
     * @name Conversion tools.
     *
     * Convert to Newick standard tree description.
     * The description is for a node, and hence is to be surrounded with
     * parenthesis. ex: (A:0.001, (B:0.001, C:0.02)90:0.005)50:0.0005
     *
     * @{
     */
   
    /**
     * @brief Get the parenthesis description of a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @param writeId Tells if node ids must be printed.
     *                This will overwrite bootstrap values if there are ones.
     *                Leaves id will be added to the leave names, separated by a '_' character.
     * @return A string in the parenthesis format.
     * @throw NodeNotFoundException If the node is not found.
     */
    static std::string nodeToParenthesis(const Tree& tree, int nodeId, bool writeId = false) throw (NodeNotFoundException);

    /**
     * @brief Get the parenthesis description of a subtree.
     *
     * @param tree The tree
     * @param nodeId The node defining the subtree.
     * @param bootstrap Tell is bootstrap values must be writen.
     * If so, the content of the property with name TreeTools::BOOTSTRAP will be written as bootstrap value.
     * The property should be a Number<double> object.
     * Otherwise, the content of the property with name 'propertyName' will be written.
     * In this later case, the property should be a String object.
     * @param propertyName The name of the property to use. Only used if bootstrap = false.
     * @return A string in the parenthesis format.
     * @throw NodeNotFoundException If the node is not found.
     */
    static std::string nodeToParenthesis(const Tree& tree, int nodeId, bool bootstrap, const std::string& propertyName) throw (NodeNotFoundException);

    /**
     * @brief Get the parenthesis description of a tree.
     *
     * @param tree The tree to convert.
     * @param writeId Tells if node ids must be printed.
     *                This will overwrite bootstrap values if there are ones.
     *                Leaves id will be added to the leave names, separated by a '_' character.
     * @return A string in the parenthesis format.
     */
    static std::string treeToParenthesis(const Tree& tree, bool writeId = false);
    
    /**
     * @brief Get the parenthesis description of a tree.
     *
     * @param tree The tree to convert.
     * @param bootstrap Tell is bootstrap values must be writen.
     * If so, the content of the property with name TreeTools::BOOTSTRAP will be written as bootstrap value.
     * The property should be a Number<double> object.
     * Otherwise, the content of the property with name 'propertyName' will be written.
     * In this later case, the property should be a String object.
     * @param propertyName The name of the property to use. Only used if bootstrap = false.
     * @return A string in the parenthesis format.
     */
    static std::string treeToParenthesis(const Tree& tree, bool bootstrap, const std::string& propertyName);
    
    /** @} */

    /**
     * @name Deal with identifiers
     *
     * @{
     */

    /**
     * @brief Retrieve all nodes ids nodes from a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of the node that defines the subtree.
     * @return A vector of ids of each node in the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static std::vector<int> getNodesId(const Tree& tree, int nodeId) throw (NodeNotFoundException);

    /**
     * @brief Retrieve all nodes ids from a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of the node that defines the subtree.
     * @param nodes A vector of ids of each node in the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void getNodesId(const Tree& tree, int nodeId, std::vector<int>& nodes) throw (NodeNotFoundException);

    /**
     * @brief Get the maximum identifier used in a (sub)tree.
     *
     * This is a recursive method.
     *
     * @param tree The tree to check.
     * @param id The identifier of the subtree from which the recursion will be performed.
     * Use id=tree.getRootId() to search for the whole tree.
     * @return The identifier number with maximum value.
     */
    static int getMaxId(const Tree& tree, int id);

    /**
     * @brief Get the minimum positive non-used identifier in a (sub)tree.
     *
     * This method uses the recursive method getNodesId, and then sort the ids.
     *
     * @param tree The tree to check.
     * @param id The identifier of the subtree from which the recursion will be performed.
     * Use id=tree.getRootId() to search for the whole tree.
     * @return A non-used identifier number.
     */
    static int getMPNUId(const Tree& tree, int id);

    /**
     * @brief Check if the ids are uniques.
     *
     * @param tree The tree to check.
     * @param throwException If set to true, the function throws qn exception if a duplicated is found.
     * @return true if the tree has uniqe ids.
     */
    static bool checkIds(const Tree& tree, bool throwException) throw (Exception);

    /** @} */


    /**
     * @name Topology methods
     *
     * @{
     */

    /**
     * @brief Creates a sequence data set corresponding to the Matrix Representation of the input trees
     *
     * @author Nicolas Galtier
     * Trees can have distinct sets of elements - missing data will be represented as 'N'.
     * The output alignment (DNA sequences including only A, C and N)) is ready for maximum parsimony analysis
     * according to the MRP supertree method.
     */
    static VectorSiteContainer* MRPEncode(const std::vector<Tree*>& vecTr);

    /**
     * @brief Tells whether two trees have the same unrooted topology
     *
     * @author Nicolas Galtier
     * Note that the location of the root, if any, is ignored.
     */
    static bool haveSameTopology(const Tree& tr1, const Tree& tr2);

    /**
     * @brief Calculates the Robinson-Foulds topological distance between two trees
     *
     * The two trees must share a common set of leaves (checked if checkNames is true)
     * Three numbers are calculated:
     *
     * @author Nicolas Galtier
     * @param tr1 First input tree.
     * @param tr2 Second input tree.
     * @param missing_in_tr2 Output as the number of bipartitions occurring in the first tree but not the second
     * @param missing_in_tr1 Output as the number of bipartitions occurring in the second tree but not the first
     * @param checkNames Tell whether we should check the trees first.
     * @return Robinson-Foulds distance = *missing_in_tr1 + *missing_in_tr2
     * @throw Exception If checkNames is set to true and trees do not share the same leaves names.
     */
    static int robinsonFouldsDistance(const Tree& tr1, const Tree& tr2, bool checkNames = true, int* missing_in_tr2 = NULL, int* missing_in_tr1 = NULL) throw (Exception);

    /**
     * @brief Counts the total number of occurrences of every bipartition from the input trees
     *
     * Returns the list of distinct bipartitions found at least once in the set of input trees,
     * and writes the number of occurrence of each of these bipartitions in vector bipScore.
     *
     * @author Nicolas Galtier
     * @param vecTr Vector of input trees (must share a common set of leaves - not checked in this function)
     * @param bipScore Output as the numbers of occurrences of the returned distinct bipartitions
     * @return A BipartitionList object including only distinct bipartitions
     */
    static BipartitionList* bipartitionOccurrences(const std::vector<Tree*>& vecTr, std::vector<size_t>& bipScore);

    /**
     * @brief General greedy consensus tree method
     *
     * Calculates the consensus tree of a set of trees defined from the number of occurrences
     * of bipartitions. Bipartitions are considered in decreasing score order.
     * A bipartition is included if it is compatible with all previously included bipartitions, and if its score
     * is higher than a threshold.
     *
     * @author Nicolas Galtier
     * @param vecTr Vector of input trees (must share a common set of leaves - checked if checkNames is true)
     * @param threshold Minimal acceptable score =number of occurrence of a bipartition/number of trees (0.<=threshold<=1.)
     * @param checkNames Tell whether we should check the trees first.
     */
    static TreeTemplate<Node>* thresholdConsensus(const std::vector<Tree*>& vecTr, double threshold, bool checkNames = true) throw (Exception);

    /**
     * @brief Fully-resolved greedy consensus tree method
     *
     * Calls thresholdConsensus with threshold=0, i.e. no constraint on the number of occurrence of bipartitions.
     * The resulting tree is fully resolved.
     *
     * @author Nicolas Galtier
     * @param vecTr Vector of input trees (must share a common set of leaves - checked if checkNames is true)
     * @param checkNames Tell whether we should check the trees first.
     */
    static TreeTemplate<Node>* fullyResolvedConsensus(const std::vector<Tree*>& vecTr, bool checkNames = true);

    /**
     * @brief Majority consensus tree method
     *
     * Calls thresholdConsensus with threshold=0.5: internal branches present in a majority of trees are kept.
     *
     * @author Nicolas Galtier
     * @param vecTr Vector of input trees (must share a common set of leaves - checked if checkNames is true)
     * @param checkNames Tell whether we should check the trees first.
     */
    static TreeTemplate<Node>* majorityConsensus(const std::vector<Tree*>& vecTr, bool checkNames = true);

    /**
     * @brief Strict consensus tree method
     *
     * Calls thresholdConsensus with threshold=1: only internal branches present in all trees are kept.
     *
     * @author Nicolas Galtier
     * @param vecTr Vector of input trees (must share a common set of leaves - checked if checkNames is true)
     * @param checkNames Tell whether we should check the trees first.
     */
    static TreeTemplate<Node>* strictConsensus(const std::vector<Tree*>& vecTr, bool checkNames = true);

    /** @} */

    /**
     * @brief Matrix Representation Parsimony supertree method
     *
     * This implementation of the MRP method takes a BIONJ tree (Jukes-Cantor distances)
     * as the starting tree and optimizes the parsimony score using only NNI (in a
     * PHYML-like way).
     *
     * @author Nicolas Galtier
     * @param vecTr A vector of trees.
     * @return The MRP super tree.
     */
    static Tree* MRP(const std::vector<Tree*>& vecTr);

    /**
     * @brief Compute bootstrap values.
     *
     * @param tree Input tree. the BOOTSTRAP banch property of the tree will be modified if it already exists.
     * @param vecTr A list of trees to compare to 'tree'.
     * @param verbose Tell if a progress bar should be displayed.
     */
    static void computeBootstrapValues(Tree& tree, const std::vector<Tree*>& vecTr, bool verbose = true);
	
    /**
     * @brief Determine the mid-point position of the root along the branch that already contains the root. Consequently, the topology of the rooted tree remains identical.
     * 
     * This code uses two inner functions to compute the mid-point position: statFromNode_ and bestRootPosition_.
     * This code is inspired by a code performing a similar calculation in Seaview (Guindon et al., 2010, Mol. Biol. Evol. 27(2):221-4).
     * 
     * @param tree The rooted tree for which the root has to be moved to its mid-point position, along the branch where it already stands.
     */    
    static void constrainedMidPointRooting(Tree& tree);
	
    /**
     * @name Some properties.
     *
     * @{
     */
     
    /**
     * @brief Bootstrap tag.
     */
    static const std::string BOOTSTRAP;

  private:
	  struct Moments_ {
	    double N;
	    double sum, squaredSum;
	    Moments_(): N(0), sum(0), squaredSum(0) {}
	  };	  

	  static Moments_ statFromNode_(Tree& tree, int rootId);
	  static double bestRootPosition_(Tree& tree, int nodeId1, int nodeId2, double length);


    /** @} */

};

} //end of namespace bpp.

#endif  //_TREETOOLS_H_

