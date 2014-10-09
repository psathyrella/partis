//
// File: TreeParsimonyScore.h
// Created by: Julien Dutheil
// Created on: Thu Jul 28 17:15 2005
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

#ifndef _TREEPARSIMONYSCORE_H_
#define _TREEPARSIMONYSCORE_H_

#include "../TreeTemplate.h"

#include <Bpp/Clonable.h>

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief Compute a parsimony score.
 */
class TreeParsimonyScore :
  public virtual Clonable
{
public:
  TreeParsimonyScore() {}
  virtual ~TreeParsimonyScore() {}

#if defined(NO_VIRTUAL_COV)
  Clonable* clone() const = 0;
#else
  TreeParsimonyScore* clone() const = 0;
#endif

public:
  /**
   * @brief Get the score for the current tree, i.e. the total minimum number of changes in the tree.
   *
   * @return The minimum total number of changes in the tree.
   */
  virtual unsigned int getScore() const = 0;

  /**
   * @brief Get the score for a given site for the current tree, i.e. the total minimum number of changes in the tree for this site.
   *
   * @param site The corresponding site.
   * @return The minimum total number of changes in the tree for site 'site'.
   */
  virtual unsigned int getScoreForSite(size_t site) const = 0;

  /**
   * @brief Get the score for each site for the current tree, i.e. the total minimum number of changes in the tree for each site.
   *
   * @return The minimum total number of changes in the tree for each site.
   */
  virtual std::vector<unsigned int> getScoreForEachSite() const = 0;

  /**
   * @brief Get the tree for wich scores are computed.
   *
   * @return The tree associated to this object.
   */
  virtual const Tree& getTree() const = 0;
};
} // end of namespace bpp.

#endif // _TREEPARSIMONYSCORE_H_


