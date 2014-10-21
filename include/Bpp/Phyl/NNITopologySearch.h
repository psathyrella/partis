//
// File: NNITopologySearch.h
// Created by: Julien Dutheil
// Created on: Wed Oct 12 10:52 2005
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

#ifndef _NNITOPOLOGYSEARCH_H_
#define _NNITOPOLOGYSEARCH_H_

#include "TopologySearch.h"
#include "NNISearchable.h"

namespace bpp
{

/**
 * @brief NNI topology search method.
 *
 * Several algorithm are implemented:
 * - Fast algorithm: loop over all nodes, check all NNIs and perform the corresponding change if it improve the score.
 *   When a NNI is done, reloop from the first node without checking the remaining ones.
 * - Better algorithm: loop over all nodes, check all NNIS.
 *   Then choose the NNI corresponding to the best improvement and perform it.
 *   Then re-loop over all nodes.
 * - PhyML algorithm (not fully tested, use with care): as the previous one, but perform all NNI improving the score at the same time.
 *   Leads to faster convergence.
 */
class NNITopologySearch :
  public virtual TopologySearch
{
	public:
		const static std::string FAST;
		const static std::string BETTER;
		const static std::string PHYML;
		
	private:
		NNISearchable* searchableTree_;
    std::string algorithm_;
		unsigned int verbose_;
    std::vector<TopologyListener*> topoListeners_;
		
	public:
		NNITopologySearch(
        NNISearchable& tree,
        const std::string& algorithm = FAST,
        unsigned int verbose = 2) :
      searchableTree_(&tree), algorithm_(algorithm), verbose_(verbose), topoListeners_()
    {}

    NNITopologySearch(const NNITopologySearch& ts) :
      searchableTree_(ts.searchableTree_),
      algorithm_(ts.algorithm_),
      verbose_(ts.verbose_),
      topoListeners_(ts.topoListeners_)
    {
      //Hard-copy all listeners:
      for (unsigned int i = 0; i < topoListeners_.size(); i++)
        topoListeners_[i] = dynamic_cast<TopologyListener*>(ts.topoListeners_[i]->clone());
    }
	
    NNITopologySearch& operator=(const NNITopologySearch& ts)
    {
      searchableTree_ = ts.searchableTree_;
      algorithm_      = ts.algorithm_;
      verbose_        = ts.verbose_;
      topoListeners_  = ts.topoListeners_;
      //Hard-copy all listeners:
      for (unsigned int i = 0; i < topoListeners_.size(); i++)
        topoListeners_[i] = dynamic_cast<TopologyListener*>(ts.topoListeners_[i]->clone());
      return *this;
    }
	
	
		virtual ~NNITopologySearch()
    {
      for (std::vector <TopologyListener*>::iterator it = topoListeners_.begin();
           it != topoListeners_.end();
           it++)
        delete *it;
    }

	public:
		void search() throw (Exception);
    
    /**
     * @brief Add a listener to the list.
     *
     * All listeners will be notified in the order of the list.
     * The first listener to be notified is the NNISearchable object itself.
     *
     * The listener will be owned by this instance, and copied when needed.
     */
    void addTopologyListener(TopologyListener* listener)
    {
      if (listener)
        topoListeners_.push_back(listener);
    }

	public:
		/**
		 * @brief Retrieve the tree.
		 *
		 * @return The tree associated to this instance.
		 */
		const Tree& getTopology() const { return searchableTree_->getTopology(); }
		
    /**
     * @return The NNISearchable object associated to this instance.
     */
    NNISearchable* getSearchableObject() { return searchableTree_; }
    /**
     * @return The NNISearchable object associated to this instance.
     */
    const NNISearchable* getSearchableObject() const { return searchableTree_; }

	protected:
		void searchFast()   throw (Exception);
		void searchBetter() throw (Exception);
		void searchPhyML()  throw (Exception);

    /**
     * @brief Process a TopologyChangeEvent to all listeners.
     */
    void notifyAllPerformed(const TopologyChangeEvent& event);
    /**
     * @brief Process a TopologyChangeEvent to all listeners.
     */
    void notifyAllTested(const TopologyChangeEvent& event);
    /**
     * @brief Process a TopologyChangeEvent to all listeners.
     */
    void notifyAllSuccessful(const TopologyChangeEvent& event);
		
};

} //end of namespace bpp.

#endif //_NNITOPOLOGYSEARCH_H_

