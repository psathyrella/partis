//
// File: ProbabilisticRewardMapping.h
// Created by: Laurent Guéguen
// Created on: mercredi 27 mars 2013, à 15h 20
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

#ifndef _PROBABILISTICREWARDMAPPING_H_
#define _PROBABILISTICREWARDMAPPING_H_

#include "RewardMapping.h"
#include "Reward.h"
#include "../TreeExceptions.h"

#include <Bpp/Text/TextTools.h>

// From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief Data storage class for probabilistic rewards mappings.
 *
 * A 'probabilistic' mapping contains an single value for each branch and each site.
 * This number is an average reward.
 * Probabilistic was coined there by opposition to the'stochastic'
 * mapping, where a path (sequence of rewards along the branch) is
 * available for each branch and site.
 */
  class ProbabilisticRewardMapping:
    public AbstractRewardMapping
  {
  private:
    const Reward* reward_;
    /**
     * @brief Rewards storage.
     *
     * Rewards are stored by sites.
     */
    std::vector< std::vector<double> > mapping_;
    
  public:
    
    /**
     * @brief Build a new ProbabilisticRewardMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     * @param reward A pointer toward the Reward object that has been used for the mapping, if any.
     * @param numberOfSites The number of sites to map.
     */
    ProbabilisticRewardMapping(const Tree& tree, const Reward* reward, size_t numberOfSites) :
      AbstractMapping(tree), AbstractRewardMapping(tree), reward_(reward), mapping_(0)
    {
      setNumberOfSites(numberOfSites);
    }
    
    /**
     * @brief Build a new ProbabilisticRewardMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     */
    ProbabilisticRewardMapping(const Tree& tree) :
    AbstractMapping(tree), AbstractRewardMapping(tree), reward_(0), mapping_(0)
    {}
    

    ProbabilisticRewardMapping* clone() const { return new ProbabilisticRewardMapping(*this); }

    ProbabilisticRewardMapping(const ProbabilisticRewardMapping& prm):
    AbstractMapping(prm), AbstractRewardMapping(prm), reward_(prm.reward_), mapping_(prm.mapping_)
    {}

    ProbabilisticRewardMapping& operator=(const ProbabilisticRewardMapping& prm)
    {
      AbstractRewardMapping::operator=(prm);
      reward_  = prm.reward_;
      mapping_ = prm.mapping_;
      return *this;
    }

    virtual ~ProbabilisticRewardMapping() {}

  public:

    virtual double getReward(int nodeId, size_t siteIndex) const
    {
      return mapping_[siteIndex][getNodeIndex(nodeId)];
    }
    
    /**
     * @brief (Re)-set the phylogenetic tree associated to this mapping.
     *
     * @param tree The new tree.
     */
    virtual void setTree(const Tree& tree);

    virtual void setNumberOfSites(size_t numberOfSites);
    
    /**
     * @brief Direct access to rewards.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual double& operator()(size_t nodeIndex, size_t siteIndex)
    {
      return mapping_[siteIndex][nodeIndex];
    }

    /**
     * @brief Direct access to rewards.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual const double& operator()(size_t nodeIndex, size_t siteIndex) const
    {
      return mapping_[siteIndex][nodeIndex];
    }
     
    /**
     * @brief Direct access to rewards.
     *
     * @warning No index checking is performed, use with care!
     */
    std::vector<double>& operator[](size_t siteIndex)
    {
      return mapping_[siteIndex];
    }

    /**
     * @brief Direct access to rewards.
     *
     * @warning No index checking is performed, use with care!
     */
    const std::vector<double>& operator[](size_t siteIndex) const
    {
      return mapping_[siteIndex];
    }
};

} //end of namespace bpp.

#endif //_PROBABILISTICREWARDMAPPING_H_

