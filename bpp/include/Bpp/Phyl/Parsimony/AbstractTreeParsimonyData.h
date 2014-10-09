//
// File: AbstractTreeParsimonyData.h
// Created by: Julien Dutheil
// Created on: Tue Jan 09 17:15 2007
// From file AbstractTreeParsimonyScore.h
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

#ifndef _ABSTRACTTREEPARSIMONYDATA_H_
#define _ABSTRACTTREEPARSIMONYDATA_H_

#include "TreeParsimonyData.h"

namespace bpp
{
/**
 * @brief Partial implementation of the TreeParsimonyData interface.
 *
 * This data structure provides a simple compression, by performing and storing computations
 * only one time per identical sites.
 *
 * The compression is achieved by the TreeParsimonyScore object.
 * The correspondance between sites in the dataset and the arrays in the structures is given
 * by the rootPatternLinks_ array: the array indice for site @f$i@f$ if given by:
 * @code
 * rootPatternLinks_[i]
 * @endcode
 *
 * Finally, the rootWeights_ array gives for each array position, the number of sites with this
 * pattern.
 * The global parsimony score is then given by the sum of all scores for each array position,
 * weighted by the corresponding number of sites.
 */
class AbstractTreeParsimonyData :
  public TreeParsimonyData
{
protected:
  std::vector<size_t> rootPatternLinks_;
  std::vector<unsigned int> rootWeights_;
  const TreeTemplate<Node>* tree_;

public:
  AbstractTreeParsimonyData(const TreeTemplate<Node>* tree) :
    rootPatternLinks_(),
    rootWeights_(),
    tree_(tree)
  {}

  AbstractTreeParsimonyData(const AbstractTreeParsimonyData& atpd) :
    rootPatternLinks_(atpd.rootPatternLinks_),
    rootWeights_(atpd.rootWeights_),
    tree_(atpd.tree_)
  {}

  AbstractTreeParsimonyData& operator=(const AbstractTreeParsimonyData& atpd)
  {
    rootPatternLinks_ = atpd.rootPatternLinks_;
    rootWeights_      = atpd.rootWeights_;
    tree_             = atpd.tree_;
    return *this;
  }


  virtual ~AbstractTreeParsimonyData() {}

public:
  size_t getRootArrayPosition(size_t site) const
  {
    return rootPatternLinks_[site];
  }

  unsigned int getWeight(size_t pos) const
  {
    return rootWeights_[pos];
  }

  const TreeTemplate<Node>* getTree() const { return tree_; }

protected:
  void setTreeP_(const TreeTemplate<Node>* tree) { tree_ = tree; }
  const TreeTemplate<Node>* getTreeP_() const { return tree_; }
};
} // end of namespace bpp.

#endif // _ABSTRACTTREEPARSIMONYDATA_H_

