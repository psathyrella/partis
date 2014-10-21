//
// File: AbstractTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:57:21 2003
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

#ifndef _ABSTRACTTREELIKELIHOOD_H_
#define _ABSTRACTTREELIKELIHOOD_H_

#include "TreeLikelihood.h"
#include "../Tree.h"
#include "../TreeTemplate.h"

#include <Bpp/Numeric/AbstractParametrizable.h>

//From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{

/**
 * @brief Partial implementation of the TreeLikelihood interface. 
 *
 * This class implements a few methods useful for most of likelihood
 * computation methods.
 *
 * It includes a tree_ and a data_ pointers.
 * This objects are owned by the class, and hence hard copied when cloning, and destroyed by the destructor.
 * 
 * - The Parametrizable interface;
 * - The getTree() method;
 * 
 * It also adds an abstract method for recursive computations.
 */
class AbstractTreeLikelihood :
	public virtual TreeLikelihood,
	public AbstractParametrizable
{
  public:
    /**
     * @brief A very simple branch iterator.
     *
     * The constructor takes a vector of nodes id to iterate over.
     */
    class SimpleBranchIterator :
      public BranchIterator
    {
      private:
        std::vector<int> nodesId_;
        size_t index_;

      public:
        SimpleBranchIterator(const std::vector<int>& nodesId) :
          nodesId_(nodesId), index_(0) {}

      public:
        int next() throw (Exception)
        {
          if (!hasNext())
            throw Exception("AbstractTreeLikelihood::SimpleBranchIterator::next(). No more branch in the set.");
          return nodesId_[index_++];
        }

        bool hasNext() const { return index_ < nodesId_.size(); }
        
    };

    /**
     * @brief A very simple site iterator.
     *
     * This iterator loops over a continuous range of sites.
     * The ocnstructor takes as input the number of sites to iterate over,
     * and optionally an offset argument, specifying the index of the first site.
     */
    class SimpleSiteIterator :
      public SiteIterator
    {
      private:
        size_t maxIndex_;
        size_t index_;
        size_t offset_;

      public:
        SimpleSiteIterator(size_t nbSites, size_t offset = 0) :
          maxIndex_(nbSites), index_(0), offset_(offset) {}

      public:
        size_t next() throw (Exception)
        {
          if (!hasNext())
            throw Exception("AbstractTreeLikelihood::SimpleSiteIterator::next(). No more site in the set.");
          return offset_ + index_++;
        }

        bool hasNext() const { return index_ < maxIndex_; }
        
    };
   
    /**
     * @name Branch iterator for models without site partition.
     *
     * @{
     */
    class ConstNoPartitionBranchModelDescription :
      public ConstBranchModelDescription
    {
      private:
        const SubstitutionModel* model_;
        size_t nbSites_;

      public:
        ConstNoPartitionBranchModelDescription(const SubstitutionModel* model, size_t nbSites) :
          model_(model), nbSites_(nbSites) {}

        ConstNoPartitionBranchModelDescription(const ConstNoPartitionBranchModelDescription& bmd) :
          model_(bmd.model_),
          nbSites_(bmd.nbSites_)
        {}

        ConstNoPartitionBranchModelDescription& operator=(const ConstNoPartitionBranchModelDescription& bmd)
        {
          model_ = bmd.model_;
          nbSites_ = bmd.nbSites_;
          return *this;
        }

      public:
        const SubstitutionModel* getModel() const { return model_; }
        
        SiteIterator* getNewSiteIterator() const { return new SimpleSiteIterator(nbSites_); }
    };

    class ConstNoPartitionBranchModelIterator :
      public ConstBranchModelIterator
    {
      private:
        ConstNoPartitionBranchModelDescription branchModelDescription_;
        size_t index_;

      public:
        ConstNoPartitionBranchModelIterator(const SubstitutionModel* model, size_t nbSites) :
          branchModelDescription_(model, nbSites), index_(0) {}

      public:
        ConstNoPartitionBranchModelDescription* next() throw (Exception)
        {
          if (!hasNext())
            throw Exception("AbstractHomogeneousTreeLikelihood::ConstHomogeneousBranchModelIterator::next(). No more branch in the set.");
          index_++;
          return &branchModelDescription_;
        }

        bool hasNext() const { return index_ == 0; }
    };
 
    class ConstNoPartitionSiteModelDescription :
      public ConstSiteModelDescription
    {
      private:
        const SubstitutionModel* model_;
        std::vector<int> nodesId_;

      public:
        ConstNoPartitionSiteModelDescription(const SubstitutionModel* model, const std::vector<int> nodesId) :
          model_(model), nodesId_(nodesId) {}

        ConstNoPartitionSiteModelDescription(const ConstNoPartitionSiteModelDescription& smd) :
          model_(smd.model_),
          nodesId_(smd.nodesId_)
        {}

        ConstNoPartitionSiteModelDescription& operator=(const ConstNoPartitionSiteModelDescription& smd)
        {
          model_ = smd.model_;
          nodesId_ = smd.nodesId_;
          return *this;
        }

      public:
        const SubstitutionModel* getModel() const { return model_; }
        
        BranchIterator* getNewBranchIterator() const { return new SimpleBranchIterator(nodesId_); }
    };

   /** @} */




	protected:
		const SiteContainer* data_;
		mutable TreeTemplate<Node>* tree_;
		bool computeFirstOrderDerivatives_;
		bool computeSecondOrderDerivatives_;
    bool initialized_;

	public:
		AbstractTreeLikelihood():
      AbstractParametrizable(""),
      data_(0),
      tree_(0),
      computeFirstOrderDerivatives_(true),
      computeSecondOrderDerivatives_(true),
      initialized_(false) {}

    AbstractTreeLikelihood(const AbstractTreeLikelihood & lik):
      AbstractParametrizable(lik),
      data_(0),
      tree_(0),
      computeFirstOrderDerivatives_(lik.computeFirstOrderDerivatives_),
      computeSecondOrderDerivatives_(lik.computeSecondOrderDerivatives_),
      initialized_(lik.initialized_) 
    {
      if (lik.data_) data_ = dynamic_cast<SiteContainer*>(lik.data_->clone());
      if (lik.tree_) tree_ = lik.tree_->clone();
    }

    AbstractTreeLikelihood & operator=(const AbstractTreeLikelihood& lik)
    {
      AbstractParametrizable::operator=(lik);
      if (data_) delete data_;
      if (lik.data_) data_ = dynamic_cast<SiteContainer*>(lik.data_->clone());
      else           data_ = 0;
      if (tree_) delete tree_;
      if (lik.tree_) tree_ = lik.tree_->clone();
      else           tree_ = 0;
      computeFirstOrderDerivatives_ = lik.computeFirstOrderDerivatives_;
      computeSecondOrderDerivatives_ = lik.computeSecondOrderDerivatives_;
      initialized_ = lik.initialized_;
      return *this;
    }

    /**
     * @brief Abstract class destructor
     *
     * This destructor is empty.
     */
		virtual ~AbstractTreeLikelihood()
    {
      if (data_) delete data_;
      if (tree_) delete tree_;
    }
	
	public:
		/**
		 * @name The TreeLikelihood interface.
		 *
		 * @{
		 */
		const SiteContainer* getData() const { return data_; }
		const Alphabet* getAlphabet() const { return data_->getAlphabet(); }	
		Vdouble getLikelihoodForEachSite()                 const;
		Vdouble getLogLikelihoodForEachSite()              const;
		VVdouble getLikelihoodForEachSiteForEachState()    const;
		VVdouble getLogLikelihoodForEachSiteForEachState() const;
		size_t getNumberOfSites() const { return data_->getNumberOfSites(); }
		size_t getNumberOfStates() const { return data_->getAlphabet()->getSize(); }
		const Tree& getTree() const { return *tree_; }
		void enableDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
		void enableFirstOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = yn; }
		void enableSecondOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
		bool enableFirstOrderDerivatives() const { return computeFirstOrderDerivatives_; }
		bool enableSecondOrderDerivatives() const { return computeSecondOrderDerivatives_; }
    bool isInitialized() const { return initialized_; }
    void initialize() throw (Exception) { initialized_ = true; }
		/** @} */

//	protected:
//		
//		/**
//		 * @brief Recompute pxy_, dpxy_ and d2pxy_ arrays, and derivatives if needed.
//		 *
//		 * This method is called when some parameter has changed.
//		 *
//		 * @param params The parameters that changed.
//		 */
//		virtual void fireParameterChanged(const ParameterList & params) = 0;
		
};

} //end of namespace bpp.

#endif	//_ABSTRACTTREELIKELIHOOD_H_

