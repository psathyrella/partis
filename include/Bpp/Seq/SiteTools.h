//
// File SiteTools.h
// Author: Julien Dutheil
//         Guillaume Deuchst
// Last modification : Friday August 8 2003
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

#ifndef _SITETOOLS_H_
#define _SITETOOLS_H_

#include "SymbolListTools.h"
#include "Site.h"
#include <Bpp/Exceptions.h>

// From the STL:
#include <map>

namespace bpp
{
/**
 * @brief Utilitary methods dealing with sites.
 */

class SiteTools :
  public SymbolListTools
{
public:
  SiteTools() {}
  virtual ~SiteTools() {}

public:
  /**
   * @param site A site.
   * @return True if the site contains one or several gap(s).
   */
  static bool hasGap(const Site& site);

  /**
   * @param site A site.
   * @return True if the site contains only gaps.
   */
  static bool isGapOnly(const Site& site);

  /**
   * @param site A site.
   * @return True if the site contains only gaps.
   */
  static bool isGapOrUnresolvedOnly(const Site& site);

  /**
   * @param site A site.
   * @return True if the site contains one or several unknwn characters.
   */
  static bool hasUnknown(const Site& site);

  /**
   * @param site A site.
   * @return True if the site contains no gap and no unknown characters.
   */
  static bool isComplete(const Site& site);

  /**
   * @brief Tell if a site is constant, that is displaying the same state in all sequences that do not present a gap.
   *
   * @param site A site.
   * @param ignoreUnknown If true, positions with unknown positions will be ignored.
   * Otherwise, a site with one single state + any uncertain state will not be considered as constant.
   * @param unresolvedRaisesException In case of ambiguous case (gap only site for instance), throw an exception. Otherwise returns false.
   * @return True if the site is made of only one state.
   * @throw EmptySiteException If the site has size 0 or if the site cannot be resolved (for instance is made of gaps only) and unresolvedRaisesException is set to true.
   */
  static bool isConstant(const Site& site, bool ignoreUnknown = false, bool unresolvedRaisesException = true) throw (EmptySiteException);

  /**
   * @param site1 The first site.
   * @param site2 The second site.
   * @return True if the two states have the same content (and, of course, alphabet).
   */
  static bool areSitesIdentical(const Site& site1, const Site& site2);

  /**
   * @brief Compute the Shannon entropy index of a site.
   *
   * \f[
   * I = - \sum_x f_x\cdot \ln(f_x)
   * \f]
   * where \f$f_x\f$ is the frequency of state \f$x\f$.
   *
   * @author J. Dutheil
   * @param site A site.
   * @param resolveUnknowns Tell is unknown characters must be resolved.
   * @return The Shannon entropy index of this site.
   * @throw EmptySiteException If the site has size 0.
   */
  static double variabilityShannon(const Site& site, bool resolveUnknowns) throw (EmptySiteException);

  /**
   * @brief Compute the factorial diversity index of a site.
   *
   * \f[
   * F = \frac{log\left(\left(\sum_x p_x\right)!\right)}{\sum_x \log(p_x)!}
   * \f]
   * where \f$p_x\f$ is the number of times state \f$x\f$ is observed in the site.
   *
   * @author J. Dutheil
   * @param site A site.
   * @return The factorial diversity index of this site.
   * @throw EmptySiteException If the site has size 0.
   */
  static double variabilityFactorial(const Site& site) throw (EmptySiteException);

  /**
   * @brief Compute the mutual information between two sites.
   *
   * \f[
   * MI = \sum_x \sum_y p_{x,y}\ln\left(\frac{p_{x,y}}{p_x \cdot p_y}\right)
   * \f]
   * where \f$p_x\f$ and \f$p_y\f$ are the frequencies of states \f$x\f$ and \f$y\f$, and
   * \f$p_{x,y}\f$ is the frequency of the pair \f$(x,y)\f$.
   *
   * @author J. Dutheil
   * @param site1 First site
   * @param site2 Second site
   * @param resolveUnknowns Tell is unknown characters must be resolved.
   * @return The mutual information for the pair of sites.
   * @throw DimensionException If the sites do not have the same length.
   * @throw EmptySiteException If the sites have size 0.
   */
  static double mutualInformation(const Site& site1, const Site& site2, bool resolveUnknowns) throw (DimensionException,EmptySiteException);

  /**
   * @brief Compute the entropy of a site. This is an alias of method variabilityShannon.
   *
   * \f[
   * I = - \sum_x f_x\cdot \ln(f_x)
   * \f]
   * where \f$f_x\f$ is the frequency of state \f$x\f$.
   *
   * @author J. Dutheil
   * @param site A site.
   * @param resolveUnknowns Tell is unknown characters must be resolved.
   * @return The Shannon entropy index of this site.
   * @throw EmptySiteException If the site has size 0.
   */
  static double entropy(const Site& site, bool resolveUnknowns) throw (EmptySiteException) {
    return variabilityShannon(site, resolveUnknowns); 
  }


  /**
   * @brief Compute the joint entropy between two sites.
   *
   * \f[
   * H_{i,j} = - \sum_x \sum_y p_{x,y}\ln\left(p_{x,y}\right)
   * \f]
   * where \f$p_{x,y}\f$ is the frequency of the pair \f$(x,y)\f$.
   *
   * @author J. Dutheil
   * @param site1 First site
   * @param site2 Second site
   * @param resolveUnknowns Tell is unknown characters must be resolved.
   * @return The mutual information for the pair of sites.
   * @throw DimensionException If the sites do not have the same length.
   * @throw EmptySiteException If the sites have size 0.
   */
  static double jointEntropy(const Site& site1, const Site& site2, bool resolveUnknowns) throw (DimensionException,EmptySiteException);

  /**
   * @brief Compute the heterozygosity index of a site.
   *
   * \f[
   * H = 1 - \sum_x f_x^2
   * \f]
   * where \f$f_x\f$ is the frequency of state \f$x\f$.
   *
   * @param site A site.
   * @return The heterozygosity index of this site.
   * @throw EmptySiteException If the site has size 0.
   */
  static double heterozygosity(const Site& site) throw (EmptySiteException);

  /**
   * @brief Give the number of distinct characters at a site.
   *
   * @param site a Site
   * @return The number of distinct characters in the given site.
   */
  static size_t getNumberOfDistinctCharacters(const Site& site) throw (EmptySiteException);

  /**
   * @brief Tell if a site has singletons
   *
   *
   * @param site a Site.
   * @return True if the site has singletons.
   */
  static bool hasSingleton(const Site& site) throw (EmptySiteException);

  /**
   * @brief Tell if a site is a parsimony informative site.
   *
   * At least two distinct characters must be present.
   *
   * @param site a Site.
   * @return True if the site is parsimony informative.
   */
  static bool isParsimonyInformativeSite(const Site& site) throw (EmptySiteException);


  /**
   * @brief Tell if a site has more than 2 distinct characters
   *
   * @param site a Site.
   * @return True if the site has more than 2 distinct characters
   */
  static bool isTriplet(const Site& site) throw (EmptySiteException);
};
} // end of namespace bpp.

#endif  // _SITETOOLS_H_

