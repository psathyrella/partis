//
// File: PHyloStatistics.h
// Created by: Julien Dutheil
// Created on: Sat Aug 08 07:29 2009
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

#ifndef _PHYLOSTATISTICS_H_
#define _PHYLOSTATISTICS_H_

#include "Tree.h"

#include <Bpp/Clonable.h>

//From the STL:
#include <vector>

namespace bpp
{
  /**
   *  @brief Compute several quantities on a tree simulateously, optimizing the recursions on the tree.
   *  
   *  This class uses a TreeTemplate. If the input tree is not a TreeTemplate, then a copy is performed
   *  before any computation.
   *
   *  @see TreeTools, TreeTemplateTools.
   */
  class PhyloStatistics:
    public virtual Clonable
  {
    private:
      size_t numberOfLeaves_;
      size_t numberOfAncestors_;
      std::vector<double> branchLengths_;
      std::vector<double> nodeHeights_;
      std::vector<size_t> nodeDepths_;
      std::vector<size_t> nodeNumberOfSons_;
      std::vector<int> nodeIds_;

    public:
      PhyloStatistics() : 
        numberOfLeaves_(0), numberOfAncestors_(0),
        branchLengths_(), nodeHeights_(), nodeDepths_(), nodeNumberOfSons_(), nodeIds_()
      {}
      virtual ~PhyloStatistics() {}

#ifndef NO_VIRTUAL_COV
      Clonable*
#else
      PhyloStatistics*
#endif
      clone() const { return new PhyloStatistics(*this); }

      /**
       * @brief Compute statistics for a given input tree.
       *
       * @param tree The tree for which the statistics should be computed.
       */
      void setTree(const Tree& tree);

      size_t getNumberOfLeaves() const { return numberOfLeaves_; }
      size_t getNumberOfAncestors() const { return numberOfAncestors_; }
      const std::vector<double>& getBranchLengths() const { return branchLengths_; }
      const std::vector<double>& getNodeHeights() const { return nodeHeights_; }
      const std::vector<size_t>& getNodeDepths() const { return nodeDepths_; }
      const std::vector<size_t>& getNodeNumberOfSons() const { return nodeNumberOfSons_; }
      const std::vector<int>& getNodeIds() const { return nodeIds_; }

    private:
      void computeForSubtree_(const Node* node, double& height, size_t& depth);

  };

} //end of namespace bpp.

#endif //_PHYLOSTATISTICS_H_

