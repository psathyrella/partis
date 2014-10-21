//
// File: SitePatterns.h
// Created by: Julien Dutheil
// Created on: Tue Nov 29 15:37 2005
//  from file PatternTools.h
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
 
#ifndef _SITEPATTERNS_H_
#define _SITEPATTERNS_H_

#include "Tree.h"

#include <Bpp/Clonable.h>
#include <Bpp/Numeric/VectorTools.h>

//From bpp-seq:
#include <Bpp/Seq/Site.h>
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <map>
#include <vector>
#include <string>

namespace bpp
{

/**
 * @brief Data structure for site patterns.
 * 
 * 'names' are the sequence names
 * 'sites' points toward a unique site
 * 'weights' is the number of sites identical to this sites
 * 'indices' are the positions in the original container
 */
class SitePatterns :
  public virtual Clonable
{
  private:
    /**
     * @brief Class used for site pattern sorting.
     */
    class SortableSite
    {
	    public:
        std::string siteS; 
		    const Site* siteP;
		    size_t originalPosition;
	    
      public:
		    SortableSite() : siteS(), siteP(0), originalPosition(0) {}
        SortableSite(const SortableSite& ss) : siteS(ss.siteS), siteP(ss.siteP), originalPosition(ss.originalPosition) {}
        SortableSite& operator=(const SortableSite& ss)
        {
          siteS = ss.siteS;
          siteP = ss.siteP;
          originalPosition = ss.originalPosition;
          return *this;
        }

        bool operator<(const SortableSite& ss) const { return siteS < ss.siteS; }

		    virtual ~SortableSite() {}
    };

    ///**
    // * @brief Class used for site pattern sorting.
    // */
    //struct SSComparator :
    //  std::binary_function<SortableSite, SortableSite, bool>
    //{
	  //  bool operator()(const SortableSite& ss1, const SortableSite& ss2) const { return ss1.siteS < ss2.siteS; }
    //};

  private: 
    std::vector<std::string> names_;
    std::vector<const Site *> sites_;
    std::vector<unsigned int> weights_;
    std::vector<size_t> indices_;
    const SiteContainer* sequences_;
    const Alphabet* alpha_;
    bool own_;

  public:
   /**
     * @brief Build a new SitePattern object.
     *
     * Look for patterns (unique sites) within a site container.
     *
     * @param sequences The container to look in.
     * @param own       Tel is the class own the sequence container.
     * If yes, the sequences wll be deleted together with this instance.
     */
    SitePatterns(const SiteContainer* sequences, bool own = false);

    virtual ~SitePatterns()
    {
      if(own_) delete sequences_;
    }

    SitePatterns(const SitePatterns& patterns) :
  	  names_(patterns.names_),
	    sites_(patterns.sites_),
	    weights_(patterns.weights_),
	    indices_(patterns.indices_),
      sequences_(0),
      alpha_(patterns.alpha_),
      own_(patterns.own_)
    {
      if(!patterns.own_) sequences_ = patterns.sequences_;
      else               sequences_ = dynamic_cast<SiteContainer*>(patterns.sequences_->clone());
    }

    SitePatterns& operator=(const SitePatterns& patterns)
    {
  	  names_     = patterns.names_;
	    sites_     = patterns.sites_;
	    weights_   = patterns.weights_;
	    indices_   = patterns.indices_;
      if(!patterns.own_) sequences_ = patterns.sequences_;
      else               sequences_ = dynamic_cast<SiteContainer*>(patterns.sequences_->clone());
      alpha_     = patterns.alpha_;
      own_       = patterns.own_;
      return *this;
    }

#ifdef NO_VIRTUAL_COV
    Clonable*
#else
    SitePatterns *
#endif
    clone() const { return new SitePatterns(*this); }

  public:
    /**
     * @return The number of times each unique site was found.
     */
		const std::vector<unsigned int>& getWeights() const { return weights_; }
    /**
     * @return The position of each unique site.
     */
		const std::vector<size_t>& getIndices() const { return indices_; }

    /**
     * @return A new container with each unique site.
     */
		SiteContainer* getSites() const;
    
};

} //end of namespace bpp.

#endif // _SITEPATTERNS_H_

