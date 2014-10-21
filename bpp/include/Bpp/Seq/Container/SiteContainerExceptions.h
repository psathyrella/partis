//
// File SiteContainerExceptions.h
// Author: Julien Dutheil
// Created on: mer mar 31 2004
// 

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#ifndef _SITECONTAINEREXCEPTIONS_H_
#define _SITECONTAINEREXCEPTIONS_H_

// From STL
#include <string>
#include <Bpp/Exceptions.h>

namespace bpp
{

/**
 * @brief The site not found exception base class.
 *
 * @see Exception
 */
class SiteNotFoundException:
  public Exception
{

	protected:

		/**
		 * @brief The id of the site that was to be found.
		 */
		const std::string id;

	public:	// Class constructor
	
		/**
		 * @brief Build a new SiteNotFoundException object.
		 *
		 * @param text A message to be passed to the exception hierarchy.
		 * @param sId  A the id of the site that was to be found.
		 */
		SiteNotFoundException(const char *   text, const char * sId = "");

		/**
		 * @brief Build a new SiteNotFoundException object.
		 *
		 * @param text A message to be passed to the exception hierarchy.
		 * @param sId  A the id of the site that was to be found.
		 */
		SiteNotFoundException(const std::string & text, const std::string & sId = "");

		// Class destructor
		~SiteNotFoundException() throw();

	public:

		/**
		 * @brief Get the id of the site that was to be found.
		 *
		 * @return The id of the site that was to be found.
		 */
		virtual const std::string getSiteId() const;
};

} //end of namespace bpp.

#endif // _SITECONTAINEREXCEPTIONS_H_

