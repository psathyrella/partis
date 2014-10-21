//
// File: ProbabilisticSubstitutionMapping.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 10:47 2006
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

#ifndef _PROBABILISTICSUBSTITUTIONMAPPING_H_
#define _PROBABILISTICSUBSTITUTIONMAPPING_H_

#include "SubstitutionMapping.h"
#include "SubstitutionCount.h"
#include "../TreeExceptions.h"

#include <Bpp/Text/TextTools.h>

// From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief Data storage class for probabilistic substitution mappings.
 *
 * A 'probabilistic' mapping contains an single value for each branch and each site.
 * This number can be an average number of substitutions, optionally waited, or a probability of observing a certain number of substitutions.
 * Probabilistic was coined there by opposition to the'stochastic' mapping, where a path (number of susbstitutions + there position along the branch)
 * is available for each branch and site. The probabilistic mapping can however be extended to contain a matrix will all types of substitutions, instead of their total number.
 */
class ProbabilisticSubstitutionMapping:
  public AbstractSubstitutionMapping
{
  private:
    const SubstitutionCount* substitutionCount_;
    /**
     * @brief Substitution numbers storage.
     *
     * Numbers are stored by sites.
     */
    std::vector< std::vector< std::vector<double> > > mapping_;
  
  public:
    
    /**
     * @brief Build a new ProbabilisticSubstitutionMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     * @param sc A pointer toward the substitution count object that has been used for the mapping, if any.
     * This object allows to get the substitution types description, if there are several. If set to 0, then
     * the mapping will be considered as having only one type of substitution mapped.
     * @param numberOfSites The number of sites to map.
     */
    ProbabilisticSubstitutionMapping(const Tree& tree, const SubstitutionCount* sc, size_t numberOfSites) :
      AbstractMapping(tree), AbstractSubstitutionMapping(tree), substitutionCount_(sc), mapping_(0)
    {
      setNumberOfSites(numberOfSites);
    }

    /**
     * @brief Build a new ProbabilisticSubstitutionMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     */
    ProbabilisticSubstitutionMapping(const Tree& tree) :
    AbstractMapping(tree), AbstractSubstitutionMapping(tree), substitutionCount_(0), mapping_(0)
    {}
    

    ProbabilisticSubstitutionMapping* clone() const { return new ProbabilisticSubstitutionMapping(*this); }

    ProbabilisticSubstitutionMapping(const ProbabilisticSubstitutionMapping& psm):
    AbstractMapping(psm), AbstractSubstitutionMapping(psm), substitutionCount_(psm.substitutionCount_), mapping_(psm.mapping_)
    {}

    ProbabilisticSubstitutionMapping& operator=(const ProbabilisticSubstitutionMapping& psm)
    {
      AbstractSubstitutionMapping::operator=(psm);
      substitutionCount_ = psm.substitutionCount_;
      mapping_           = psm.mapping_;
      return *this;
    }

    virtual ~ProbabilisticSubstitutionMapping() {}

  public:

    size_t getNumberOfSubstitutionTypes() const
    {
      if (!substitutionCount_) return 1;
      return substitutionCount_->getNumberOfSubstitutionTypes();
    }
     
    virtual double getNumberOfSubstitutions(int nodeId, size_t siteIndex, size_t type) const
    {
      return mapping_[siteIndex][getNodeIndex(nodeId)][type];
    }
    
    virtual std::vector<double> getNumberOfSubstitutions(int nodeId, size_t siteIndex) const
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
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual double& operator()(size_t nodeIndex, size_t siteIndex, size_t type)
    {
      return mapping_[siteIndex][nodeIndex][type];
    }

    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual const double& operator()(size_t nodeIndex, size_t siteIndex, size_t type) const
    {
      return mapping_[siteIndex][nodeIndex][type];
    }
     
    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    std::vector< std::vector<double> >& operator[](size_t siteIndex)
    {
      return mapping_[siteIndex];
    }

    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    const std::vector< std::vector<double> >& operator[](size_t siteIndex) const
    {
      return mapping_[siteIndex];
    }
};

} //end of namespace bpp.

#endif //_PROBABILISTICSUBSTITUTIONMAPPING_H_

