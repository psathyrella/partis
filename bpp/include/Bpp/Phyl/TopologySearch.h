//
// File: TopologySearch.h
// Created by: Julien Dutheil
// Created on: Wed Oct 12 10:18 2005
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

#ifndef _TOPOLOGYSEARCH_H_
#define _TOPOLOGYSEARCH_H_

#include <Bpp/Exceptions.h>
#include <Bpp/Clonable.h>

// From the STL:
#include <string>
#include <vector>

namespace bpp
{

/**
 * @brief Class for notifying new toplogy change events.
 */
class TopologyChangeEvent
{
	protected:
    std::string message_;
		
	public:
		TopologyChangeEvent(): message_("") {}
		TopologyChangeEvent(const std::string& message): message_(message) {}
		virtual ~TopologyChangeEvent() {}

	public:
		/**
		 * @brief Get the message associated to this event.
		 *
		 * @return The message associated to this event.
		 */
		virtual const std::string& getMessage() const { return message_; }

};

class TopologySearch;

/**
 * @brief Implement this interface to be notified when the topology of a tree
 * has changed during topology search.
 */
class TopologyListener :
  public virtual Clonable
{
  public:
    TopologyListener() {}
    virtual ~TopologyListener() {}

#ifndef NO_VIRTUAL_COV
    TopologyListener* clone() const = 0;
#endif

  public:
		/**
		 * @brief Notify a topology change event.
		 *
		 * This method is to be invoked after one or several NNI are performed.
		 * It allows appropriate recomputations.
     *
     * In most case, this is the same as
     * topologyChangeTested() + topologyChangeSuccessful().
		 *
		 * @param event The topology change event.
		 */
    virtual void topologyChangePerformed(const TopologyChangeEvent& event)
    {
      topologyChangeTested(event);
      topologyChangeSuccessful(event);
    }
		/**
		 * @brief Notify a topology change event.
		 *
		 * @param event The topology change event.
		 */
    virtual void topologyChangeTested(const TopologyChangeEvent& event) = 0;

    /**
     * @brief Tell that a topology change is definitive.
     *
     * This method is called after the topologyChangeTested() method.
     *
		 * @param event The topology change event.
     */
    virtual void topologyChangeSuccessful(const TopologyChangeEvent& event) = 0;

};



/**
 * @brief Interface for topology search methods.
 */
class TopologySearch
{
	public:
		TopologySearch() {}
		virtual ~TopologySearch() {}

	public:

		/**
		 * @brief Performs the search.
		 */
		virtual void search() throw (Exception) = 0;

    /**
     * @brief Add a topology listener to this class.
     *
     * TopologyListeners will be notified when the topology of the tree is modified. 
     */
    virtual void addTopologyListener(TopologyListener* listener) = 0;			
};

} //end of namespace bpp.

#endif //_TOPOLOGYSEARCH_H_

