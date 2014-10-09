//
// File: DetailedSiteSimulator.h
// Created by: Julien Dutheil
// Created on: Tue Mar  14 10:51 2006
// from old file DetailedSequenceSimulator.h
// Created on: Wed Aug  24 15:20 2005
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

#ifndef _DETAILEDSITESIMULATOR_H_
#define _DETAILEDSITESIMULATOR_H_

#include "SiteSimulator.h"
#include "MutationProcess.h"
#include "../TreeTemplate.h"

// From the STL:
#include <map>
#include <vector>

namespace bpp
{

/**
 * @brief Data structure to store the result of a DetailedSiteSimulator.
 *
 * This data structure stores each transitional state, and the time when it occured.
 */
class SiteSimulationResult
{
  private:
    mutable std::map<int, size_t> indexes_;
    size_t currentIndex_;
    std::vector<MutationPath> paths_;
    std::vector<int> ancestralStates_;
    const Tree* tree_;
    std::vector<int> leavesId_;
    const Alphabet* alphabet_;
    
  public:
    SiteSimulationResult(const Tree* tree, const Alphabet* alphabet, int ancestralState) :
      indexes_        (),
      currentIndex_   (0),
      paths_          (),
      ancestralStates_(),
      tree_           (tree),
      leavesId_       (tree->getLeavesId()),
      alphabet_       (alphabet)
    {
      indexes_[tree->getRootId()] = 0;
      //Warning, watch out the indices there!
      ancestralStates_.push_back(ancestralState);
    }

    SiteSimulationResult(const SiteSimulationResult& ssr) :
      indexes_        (ssr.indexes_),
      currentIndex_   (ssr.currentIndex_),
      paths_          (ssr.paths_),
      ancestralStates_(ssr.ancestralStates_),
      tree_           (ssr.tree_),
      leavesId_       (ssr.leavesId_),
      alphabet_       (ssr.alphabet_)
    {}
 
    SiteSimulationResult& operator=(const SiteSimulationResult& ssr)
    {
      indexes_         = ssr.indexes_;
      currentIndex_    = ssr.currentIndex_;
      paths_           = ssr.paths_;
      ancestralStates_ = ssr.ancestralStates_;
      tree_            = ssr.tree_;
      leavesId_        = ssr.leavesId_;
      alphabet_        = ssr.alphabet_;
      return *this;
    }

    virtual ~SiteSimulationResult() {}
  
  public:
    /**
     * @return The alphabet associated to this simulation.
     */
    virtual const Alphabet* getAlphabet() const { return alphabet_; }
    
    virtual void addNode(int nodeId, MutationPath path)
    {
      indexes_[nodeId] = currentIndex_;
      currentIndex_++;
      paths_.push_back(path);
      ancestralStates_.push_back(path.getFinalState());
    }

    virtual int getAncestralState(size_t i)    const { return ancestralStates_[i]; }

    virtual int getAncestralState(int nodeId) const { return ancestralStates_[1 + indexes_[nodeId]]; }

    virtual const MutationPath& getMutationPath(size_t i) const { return paths_[i]; }

    virtual const MutationPath& getMutationPath(int nodeId) const { return paths_[indexes_[nodeId]]; }

    virtual size_t getSubstitutionCount(size_t i) const { return paths_[i].getNumberOfEvents(); }
    
    virtual void getSubstitutionCount(size_t i, const SubstitutionRegister& reg, std::vector<double>& counts) const {
      paths_[i].getEventCounts(counts, reg);
    }
    
    virtual size_t getSubstitutionCount(int nodeId) const { return paths_[indexes_[nodeId]].getNumberOfEvents(); }
    
    virtual void getSubstitutionCount(int nodeId, const SubstitutionRegister& reg, std::vector<double>& counts) const {
      paths_[indexes_[nodeId]].getEventCounts(counts, reg);
    }
    
    virtual VVdouble getSubstitutionVector(const SubstitutionRegister& reg) const
    {
      size_t n = paths_.size();
      VVdouble counts(n);
      for (size_t i = 0; i < n; ++i) {
        counts[i].resize(reg.getNumberOfSubstitutionTypes());
        paths_[i].getEventCounts(counts[i], reg);
      }
      return counts;
    }

    /**
     * @return The states at the leaves.
     */
    virtual std::vector<int> getFinalStates() const
    {
      size_t n = leavesId_.size(); 
      std::vector<int> states(n);
      for (size_t i = 0; i < n; i++)
      {
        states[i] = ancestralStates_[1 + indexes_[leavesId_[i]]];
      }
      return states;
    }

    /**
     * @return The site corresponding to this simulation.
     */
    virtual Site* getSite() const { return new Site(getFinalStates(), alphabet_); }

    /**
     * @return A vector with the leaves names.
     */
    virtual std::vector<std::string> getLeaveNames() const
    {
      size_t n = leavesId_.size(); 
      std::vector<std::string> names(n);
      for (size_t i = 0; i < n; i++)
      {
        names[i] = tree_->getNodeName(leavesId_[i]);
      }
      return names;
    }

};

//---------------------------------------------------------------------------

/**
 * @brief Data structure to store the result of a DetailedSiteSimulator.
 *
 * This sructure inherits from the SequenceSimulationResult class, and add support for
 * rate variation across sites.
 */
class RASiteSimulationResult:
  public SiteSimulationResult
{
  protected:
    double rate_;
    
  public:
    RASiteSimulationResult(const Tree* tree, const Alphabet * alphabet, int ancestralState, double rate):
      SiteSimulationResult(tree, alphabet, ancestralState),
      rate_(rate) {}

    virtual ~RASiteSimulationResult() {}
  
  public:
    /**
     * @return The rate of this simulation.
     */
    virtual double getRate() const { return rate_; }
};

//---------------------------------------------------------------------------

/**
 * @brief This interface adds the dSimulate method to the SiteSimulator interface.
 *
 * Instances of this class should be used when a detailed output of the simulation is needed.
 */
class DetailedSiteSimulator:
  public virtual SiteSimulator
{
  public:
    DetailedSiteSimulator() {}
    virtual ~DetailedSiteSimulator() {}
  
  public:
    /**
     * @brief Get a detailed simulation result for one site.
     *
     * @return A SiteSimulationResult object with all ancestral
     * states for all nodes and branches.
     */
    virtual SiteSimulationResult* dSimulate() const = 0;
    virtual SiteSimulationResult* dSimulate(int ancestralState) const = 0;
    virtual SiteSimulationResult* dSimulate(int ancestralState, double rate) const = 0;
    virtual SiteSimulationResult* dSimulate(double rate) const = 0;
    
};

} //end of namespace bpp.

#endif // _DETAILEDSITESIMULATOR_H_

