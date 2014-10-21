//
// File: SiteContainerIterator.h
// Created by: Julien Dutheil
// Created on: Sun Oct 19 12:47:16 2003
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#ifndef _SITECONTAINERITERATOR_H_
#define _SITECONTAINERITERATOR_H_

#include "../Site.h"
#include "../SiteIterator.h"
#include "SiteContainer.h"

namespace bpp
{

/**
 * @brief Partial implementation of the SiteIterator interface, allowing to loop over a site container.
 */
class AbstractSiteContainerIterator :
  public virtual ConstSiteIterator
{
	protected:
		const SiteContainer* sites_;
		int currentPosition_;
	
	public:
		AbstractSiteContainerIterator(const SiteContainer& sites);
		
    AbstractSiteContainerIterator(const AbstractSiteContainerIterator& asi) :
      sites_(asi.sites_),
      currentPosition_(asi.currentPosition_) {}
    
    AbstractSiteContainerIterator& operator=(const AbstractSiteContainerIterator& asi)
    {
      sites_ = asi.sites_;
      currentPosition_ = asi.currentPosition_;
      return *this;
    }

		virtual ~AbstractSiteContainerIterator() {}
	
};

/**
 * @brief Loop over all sites in a SiteContainer.
 */
class SimpleSiteContainerIterator: public AbstractSiteContainerIterator
{
	public:
		SimpleSiteContainerIterator(const SiteContainer& sites);
		virtual ~SimpleSiteContainerIterator() {}
	
	public:
		const Site* nextSite();
		bool hasMoreSites() const;
};

/**
 * @brief Loop over all sites without gaps in a SiteContainer.
 */
class NoGapSiteContainerIterator: public AbstractSiteContainerIterator
{
	public:
		NoGapSiteContainerIterator(const SiteContainer & sites);
		virtual ~NoGapSiteContainerIterator() {}
	
	public:
		const Site* nextSite();
		bool hasMoreSites() const;
		int nextSiteWithoutGapPosition(int current) const;
		int previousSiteWithoutGapPosition(int current) const;
};

/**
 * @brief Loop over all complete sites in a SiteContainer
 * (i.e. sites without gap and unresolved characters).
 */
class CompleteSiteContainerIterator: public AbstractSiteContainerIterator
{
	public:
		CompleteSiteContainerIterator(const SiteContainer & sites);
		virtual ~CompleteSiteContainerIterator() {}
	
	public:
		const Site * nextSite();
		bool hasMoreSites() const;
		int nextCompleteSitePosition(int current) const;
		int previousCompleteSitePosition(int current) const;
};

} //end of namespace bpp.

#endif	//_SITEITERATOR_H_

