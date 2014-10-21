//
// File: SiteExceptions.h
// Author : Julien Dutheil
// Created on: dim mar 7 2004
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

#ifndef _SITEEXCEPTIONS_H_
#define _SITEEXCEPTIONS_H_

#include <Bpp/Exceptions.h>

namespace bpp
{

class Site;

/**
 * @brief The site exception base class.
 *
 * @see Exception
 */
class SiteException:
  public Exception
{
	private:

		/**
		 * @brief A pointer toward a site object.
		 */
		const Site* site_;

	public:	// Class constructor

		/**
		 * @brief Build a new SiteException object.
		 *
		 * @param text A message to be passed to the exception hierarchy.
		 * @param s    A const pointer toward the site that threw the exception.
		 */
		SiteException(const std::string& text, const Site* s = 0);

    SiteException(const SiteException& se): Exception(se), site_(se.site_) {}
    
    SiteException& operator=(const SiteException& se)
    {
      Exception::operator=(se);
      site_ = se.site_;
      return *this;
    }

		// Class destructor
		virtual ~SiteException() throw() {}

	public:

		/**
		 * @brief Get the site that threw the exception.
		 *
		 * @return A const pointer toward the site.
		 */
		virtual const Site* getSite() const { return site_; };
};

/**
 * @brief Exception sent when a empty site is found.
 */
class EmptySiteException:
  public SiteException
{
	public:
		EmptySiteException(const std::string& text, const Site* s = 0): SiteException(text, s) {}

		virtual ~EmptySiteException() throw() {}
};

/**
 * @brief Exception sent when a site containing gap is found.
 */
class SiteWithGapException:
  public SiteException
{
	public:
		SiteWithGapException(const std::string& text, const Site* s = 0): SiteException(text, s) {}

		virtual ~SiteWithGapException() throw() {}
};

} //end of namespace bpp.

#endif // _SITEEXCEPTIONS_H_

