//
// File: MarginalAncestralStateReconstruction.h
// Created by: Julien Dutheil
// Created on: Fri Jul 08 13:32 2005
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

#ifndef _MARGINALANCESTRALSTATESRECONSTRUCTION_H_
#define _MARGINALANCESTRALSTATESRECONSTRUCTION_H_

#include "../AncestralStateReconstruction.h"
#include "DRTreeLikelihood.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Sequence.h>

// From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief Likelihood ancestral states reconstruction: marginal method.
 *
 * Reference:
 * Z Yang, S Kumar and M Nei (1995), _Genetics_ 141(4) 1641-50.
 */
class MarginalAncestralStateReconstruction:
  public virtual AncestralStateReconstruction
{
	private:
		const DRTreeLikelihood* likelihood_;
    TreeTemplate<Node> tree_;
		const Alphabet* alphabet_;
		size_t nbSites_;
		size_t nbDistinctSites_;
		size_t nbClasses_;
		size_t nbStates_;
    std::vector<size_t> rootPatternLinks_;
    std::vector<double> r_;
    std::vector<double> l_;
		
	public:
		MarginalAncestralStateReconstruction(const DRTreeLikelihood* drl) :
      likelihood_      (drl),
      tree_            (drl->getTree()),
			alphabet_        (drl->getAlphabet()),
			nbSites_         (drl->getLikelihoodData()->getNumberOfSites()),
			nbDistinctSites_ (drl->getLikelihoodData()->getNumberOfDistinctSites()),
			nbClasses_       (drl->getLikelihoodData()->getNumberOfClasses()),
			nbStates_        (drl->getLikelihoodData()->getNumberOfStates()),
			rootPatternLinks_(drl->getLikelihoodData()->getRootArrayPositions()),
      r_               (drl->getRateDistribution()->getProbabilities()),
      l_               (drl->getLikelihoodData()->getRootRateSiteLikelihoodArray())
    {}

    MarginalAncestralStateReconstruction(const MarginalAncestralStateReconstruction& masr) :
      likelihood_      (masr.likelihood_),
      tree_            (masr.tree_),
      alphabet_        (masr.alphabet_),
      nbSites_         (masr.nbSites_),
      nbDistinctSites_ (masr.nbDistinctSites_),
      nbClasses_       (masr.nbClasses_),
      nbStates_        (masr.nbStates_),
      rootPatternLinks_(masr.rootPatternLinks_),
      r_               (masr.r_),
      l_               (masr.l_)
    {}

    MarginalAncestralStateReconstruction& operator=(const MarginalAncestralStateReconstruction& masr)
    {
      likelihood_       = masr.likelihood_;
      tree_             = masr.tree_;
      alphabet_         = masr.alphabet_;
      nbSites_          = masr.nbSites_;
      nbDistinctSites_  = masr.nbDistinctSites_;
      nbClasses_        = masr.nbClasses_;
      nbStates_         = masr.nbStates_;
      rootPatternLinks_ = masr.rootPatternLinks_;
      r_                = masr.r_;
      l_                = masr.l_;
      return *this;
    }


#ifndef NO_VIRTUAL_COV
    MarginalAncestralStateReconstruction*
#else
    Clonable*
#endif
    clone() const { return new MarginalAncestralStateReconstruction(*this); }

		virtual ~MarginalAncestralStateReconstruction() {}

	public:

		/**
		 * @brief Get ancestral states for a given node as a vector of int.
		 *
		 * The size of the vector is the number of distinct sites in the container
		 * associated to the likelihood object.
		 * This method is mainly for efficient internal use in other classes.
		 * Consider using the getAncestralSequenceForNode() method for a more
		 * general output.
		 *
		 * @param nodeId The id of the node at which the states must be reconstructed.
     * @param probs  A vector to be filled with the probability for each state at each position (will be the same size as the returned vector for states).
     * @param sample Tell if the sequence should be sample from the posterior distribution instead of taking the one with maximum probability.
		 * @return A vector of states indices.
		 * @see getAncestralSequenceForNode
		 */ 
    std::vector<size_t> getAncestralStatesForNode(int nodeId, VVdouble& probs, bool sample) const;
		
    std::vector<size_t> getAncestralStatesForNode(int nodeId) const
    {
      VVdouble probs(nbSites_);
      return getAncestralStatesForNode(nodeId, probs, false);
    }
		
    std::map<int, std::vector<size_t> > getAllAncestralStates() const;

		/**
		 * @brief Get the ancestral sequence for a given node.
		 *
		 * The name of the sequence will be the name of the node if there is one, its id otherwise.
		 * A new sequence object is created, whose destruction is up to the user.
		 *
		 * @param nodeId The id of the node at which the sequence must be reconstructed.
     * @param probs  A pointer toward a vector to be filled with the probability for each state at each site (set to NULL if you don't want these probabilities).
     * @param sample Tell if the sequence should be sample from the posterior distribution instead of taking the one with maximum probability.
		 * @return A sequence object.
		 */ 
		Sequence* getAncestralSequenceForNode(int nodeId, VVdouble* probs, bool sample) const;
		
    Sequence* getAncestralSequenceForNode(int nodeId) const
    {
      return getAncestralSequenceForNode(nodeId, 0, false);
    }

    AlignedSequenceContainer* getAncestralSequences() const
    {
      return getAncestralSequences(false);
    }

#ifndef NO_VIRTUAL_COV
    AlignedSequenceContainer *
#else
    SequenceContainer *
#endif
    getAncestralSequences(bool sample) const;
	
  private:
		void recursiveMarginalAncestralStates(
			const Node* node,
			std::map<int, std::vector<size_t> >& ancestors,
			AlignedSequenceContainer& data) const;

		
};

} //end of namespace bpp.

#endif // _MARGINALANCESTRALSTATESRECONSTRUCTION_H_

