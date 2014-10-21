//
// File: RewardMapping.h
// Created by: Laurent Guéguen
// Created on: vendredi 29 mars 2013, à 11h 32
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#ifndef _REWARDMAPPING_H_
#define _REWARDMAPPING_H_

#include "Mapping.h"

#include <Bpp/Clonable.h>

//From the STL:
#include <vector>
#include <memory>

namespace bpp
{

/**
 * @brief General interface for storing reward mapping data.
 *
 * Since only probabilistic reward mapping is implemented for now, the basal 
 * interface only contains a few methods.
 * More methods are expected to be added later.
 */
  
  class RewardMapping:
    virtual public Mapping
  {
    
  public:
    RewardMapping() {}
    virtual ~RewardMapping() {}

#ifndef NO_VIRTUAL_COV
    RewardMapping* clone() const = 0;
#endif

  public:
    
    virtual double& operator()(size_t nodeIndex, size_t siteIndex) = 0;
    virtual const double& operator()(size_t nodeIndex, size_t siteIndex) const = 0;
};



/**
 * @brief Partial implementation of the reward mapping interface.
 *
 * This implementation copies the input tree in a TreeTemplate<Node> object.
 */
class AbstractRewardMapping:
    virtual public RewardMapping,
    virtual public AbstractMapping
{
public:
  //  AbstractRewardMapping() : AbstractMapping(){}

  AbstractRewardMapping(const Tree& tree) : AbstractMapping(tree){}

  AbstractRewardMapping(const AbstractRewardMapping& arm):
    AbstractMapping(arm) {}

#ifndef NO_VIRTUAL_COV
  AbstractRewardMapping* clone() const = 0;
#endif

  AbstractRewardMapping& operator=(const AbstractRewardMapping& arm)
    {
      AbstractMapping::operator=(arm);
      return *this;
    }

  virtual ~AbstractRewardMapping() {}

};

} //end of namespace bpp.

#endif //_REWARDMAPPING_H_

