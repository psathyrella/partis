//
// File Site.h
// Created by: Guillaume Deuchst
//             Julien Dutheil
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

#ifndef _SITE_H_
#define _SITE_H_

#include "SymbolList.h"
#include "SiteExceptions.h"

namespace bpp
{

/**
 * @brief The Site class.
 *
 * Define specific attributes and methods for sites manipulation.
 * It is very similar to the Sequence object (a site is a vertical sequence!),
 * and characters at each position are coded as integers.
 * Sites have a 'position' attribute.
 * This attribute stands for an indice in a an alignment, and may be used as a unique identifier,
 * in the same manner that names identify sequence objects.
 * But for now, we do not allow to construct a Site directly from a string.
 * This should not be a constraint, since you never read sites directly from a file.
 */
class Site:
  public BasicSymbolList 
{  
  private:
    /**
     * @brief The position associated to this site.
     */
    int position_;

  public:
    
    /**
     * @brief Build a new void Site object with the specified alphabet.
     *
     * @param alpha The alphabet to use.
     */
    Site(const Alphabet* alpha);

    /**
     * @brief Build a new void Site object with the specified alphabet and position.
     *
     * @param alpha    The alphabet to use.
     * @param position The position attribute for this site.
     */
    Site(const Alphabet* alpha, int position);

    /**
     * @brief Build a new Site object with the specified alphabet.
     * The content of the site is initialized from a vector of characters.
     *
     * @param site     The content of the site.
     * @param alpha    The alphabet to use.
     * @throw BadCharException If the content does not match the specified alphabet.
     */
    Site(const std::vector<std::string>& site, const Alphabet* alpha) throw (BadCharException);

    /**
     * @brief Build a new Site object with the specified alphabet and position.
     * The content of the site is initialized from a vector of characters.
     *
     * @param site     The content of the site.
     * @param alpha    The alphabet to use.
     * @param position The position attribute for this site.
     * @throw BadCharException If the content does not match the specified alphabet.
     */
    Site(const std::vector<std::string>& site, const Alphabet* alpha, int position) throw (BadCharException);

    /**
     * @brief Build a new Site object with the specified alphabet.
     * The content of the site is initialized from a vector of integers.
     *
     * @param site     The content of the site.
     * @param alpha    The alphabet to use.
     * @throw BadIntException If the content does not match the specified alphabet.
     */
    Site(const std::vector<int>& site, const Alphabet* alpha) throw (BadIntException);

    /**
     * @brief Build a new Site object with the specified alphabet and position.
     * The content of the site is initialized from a vector of integers.
     *
     * @param site     The content of the site.
     * @param alpha    The alphabet to use.
     * @param position The position attribute for this site.
     * @throw BadIntException If the content does not match the specified alphabet.
     */
    Site(const std::vector<int>& site, const Alphabet* alpha, int position) throw (BadIntException);

    /**
     * @brief The copy constructor.
     */
    Site(const Site& site);

    /**
     * @brief The assignment operator.
     */
    Site& operator=(const Site& s);

    virtual ~Site() {}
  
  public:

    /**
     * @name The Clonable interface
     *
     * @{
     */
    Site* clone() const { return new Site(*this); }
    /** @} */

    /**
     * @name Setting/getting the position of the site.
     *
     * @{
     */

    /**
     * @brief Get the position of this site.
     *
     * @return This site position.
     */
    virtual int getPosition() const { return position_; }

    /**
     * @brief Set the position of this site.
     *
     * @param position The new position of the site.
     */
    virtual void setPosition(int position) { position_ = position; }
};

// Sites comparison operators overload
bool operator == (const Site& site1, const Site& site2);
bool operator < (const Site& site1, const Site& site2);

} //end of namespace bpp.

#endif  //_SITE_H_

