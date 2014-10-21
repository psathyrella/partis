//
// File: AbstractTreeLikelihoodData.h
// Created by: Julien Dutheil
// Created on: Sat Dec 30 12:48 2006
// From file AbstractTreeLikelihood.h
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

#ifndef _ABSTRACTTREELIKELIHOODDATA_H_
#define _ABSTRACTTREELIKELIHOODDATA_H_

#include "TreeLikelihoodData.h"

//From the STL:
#include <vector>
#include <map>

namespace bpp
{

/**
 * @brief Partial implementation of the TreeLikelihoodData interface.
 *
 * This data structure provides a simple compression, by performing and storing computations
 * only one time per identical sites.
 *
 * The compression is achieved by the TreeLikelihood object.
 * The correspondance between sites in the dataset and the arrays in the structures is given
 * by the rootPatternLinks_ array: the array indice for site @f$i@f$ if given by:
 * @code
 * rootPatternLinks_[i]
 * @endcode
 *
 * Finally, the rootWeights_ array gives for each array position, the number of sites with this
 * pattern.
 * The global likelihood is then given by the product of all likelihoods for each array position,
 * weighted by the corresponding number of sites.
 */
class AbstractTreeLikelihoodData :
	public TreeLikelihoodData
{
	protected:
		/**
		 * @brief Links between sites and patterns.
		 * 
		 * The size of this vector is equal to the number of sites in the container,
		 * each element corresponds to a site in the container and points to the
		 * corresponding column in the likelihood array of the root node.
		 * If the container contains no repeated site, there will be a strict
		 * equivalence between each site and the likelihood array of the root node.
		 * However, if this is not the case, some pointers may point toward the same
		 * element in the likelihood array.
		 */
    std::vector<size_t> rootPatternLinks_;

		/**
		 * @brief The frequency of each site.
		 */
    std::vector<unsigned int> rootWeights_;

		const TreeTemplate<Node>* tree_;

		const Alphabet* alphabet_;

  public:
		AbstractTreeLikelihoodData(const TreeTemplate<Node>* tree):
      rootPatternLinks_(), rootWeights_(), tree_(tree), alphabet_(0) {}

		AbstractTreeLikelihoodData(const AbstractTreeLikelihoodData& atd) :
      rootPatternLinks_(atd.rootPatternLinks_),
      rootWeights_(atd.rootWeights_),
      tree_(atd.tree_),
      alphabet_(atd.alphabet_)
    {}

    AbstractTreeLikelihoodData& operator=(const AbstractTreeLikelihoodData& atd)
    {
      rootPatternLinks_ = atd.rootPatternLinks_;
      rootWeights_      = atd.rootWeights_;
      tree_             = atd.tree_;
      alphabet_         = atd.alphabet_;
      return *this;
    }


		virtual ~AbstractTreeLikelihoodData() {}

	public:
    std::vector<size_t>& getRootArrayPositions() { return rootPatternLinks_; }
		const std::vector<size_t>& getRootArrayPositions() const { return rootPatternLinks_; }
		size_t getRootArrayPosition(const size_t site) const
		{
			return rootPatternLinks_[site];
		}
		unsigned int getWeight(size_t pos) const
		{
			return rootWeights_[pos];
		}
		const std::vector<unsigned int>& getWeights() const
		{ 
			return rootWeights_;
		}

		const Alphabet* getAlphabet() const { return alphabet_; }

		const TreeTemplate<Node>* getTree() const { return tree_; }  

};

} //end of namespace bpp.

#endif //_ABSTRACTTREELIKELIHOODDATA_H_

