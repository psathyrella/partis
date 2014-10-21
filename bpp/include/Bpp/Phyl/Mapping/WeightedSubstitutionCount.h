//
// File: WeightedSubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 15:12 2006
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

#ifndef _WEIGHTEDSUBSTITUTIONCOUNT_H_
#define _WEIGHTEDSUBSTITUTIONCOUNT_H_

#include "SubstitutionCount.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{

/**
 * @brief Interface allowing for weighting of  substitution counts according to state properties.
 */
class WeightedSubstitutionCount:
  public virtual SubstitutionCount
{
  public:
    virtual void setWeights(const AlphabetIndex2* index, bool ownWeights) = 0;
    virtual bool hasWeights() const = 0;
    virtual const AlphabetIndex2* getWeights() const = 0;
};

/**
 * @brief Partial implementation of the WeightedSubstitutionCount interface.
 */
class AbstractWeightedSubstitutionCount:
  public virtual WeightedSubstitutionCount
{
	protected:
		const AlphabetIndex2* weights_;
		bool ownWeights_;
	
	public:
		AbstractWeightedSubstitutionCount(const AlphabetIndex2* weights, bool ownWeights) :
      weights_(weights),
      ownWeights_(ownWeights)
    {
    }

    AbstractWeightedSubstitutionCount(const AbstractWeightedSubstitutionCount& index) :
      weights_(index.weights_),
      ownWeights_(index.ownWeights_)
    {
      if (ownWeights_)
        weights_ = dynamic_cast<AlphabetIndex2*>(index.weights_->clone());
    }

    AbstractWeightedSubstitutionCount& operator=(const AbstractWeightedSubstitutionCount& index)
    {
      ownWeights_ = index.ownWeights_;
      if (ownWeights_) weights_ = dynamic_cast<AlphabetIndex2*>(index.weights_->clone());
      else weights_ = index.weights_;
      
      return *this;
    }
		
		virtual ~AbstractWeightedSubstitutionCount()
    {
			if (ownWeights_)
        delete weights_;
		}
		
  public:
    void setWeights(const AlphabetIndex2* weights, bool ownWeights);
    bool hasWeights() const { return weights_ != 0; }
    const AlphabetIndex2* getWeights() const { return weights_; }

  protected:
    virtual void weightsHaveChanged() = 0;

};

} //end of namespace bpp.

#endif //_WEIGHTEDSUBSTITUTIONCOUNT_H_

