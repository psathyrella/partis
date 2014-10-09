//
// File: TreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:36:44 2003
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

#ifndef _TREELIKELIHOOD_H_
#define _TREELIKELIHOOD_H_

#include "../Node.h"
#include "../Tree.h"
#include "../Model/SubstitutionModel.h"
#include "TreeLikelihoodData.h"

#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/VectorTools.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{

/**
 * @brief The TreeLikelihood interface.
 *
 * This interface defines the methods needed for computing the likelihood
 * of a phylogenetic tree, given a dataset.
 */ 
class TreeLikelihood:
  public virtual DerivableSecondOrder
{
  public:
    /**
     * @brief An iterator over a set of branches, specified by their node ids.
     */
    class BranchIterator
    {
      public:
        virtual ~BranchIterator() {}

      public:
        /**
         * @return The id of the next node in the set.
         */
        virtual int next() throw (Exception) = 0;
        /**
         * @return True if there is at least another node in the set.
         */
        virtual bool hasNext() const = 0;
    };

    /**
     * @brief An iterator over a set of sites, speicfied by their position.
     *
     * In most cases, the position will reflect the index of an inner array used for likelihood storage.
     */
    class SiteIterator
    {
      public:
        virtual ~SiteIterator() {}

      public:
        /**
         * @return The position of the next site in the set.
         */
        virtual size_t next() throw (Exception) = 0;
        /**
         * @return True is there is at least another site in the set.
         */
        virtual bool hasNext() const = 0;
    };

    /**
     * @brief A pair of SubstitutionModel / SiteIterator.
     */
    class ConstBranchModelDescription
    {
      public:
        virtual ~ConstBranchModelDescription() {}

      public:
        virtual const SubstitutionModel* getModel() const = 0;
        virtual SiteIterator* getNewSiteIterator() const = 0;
    };

    /**
     * @brief Iterates through all models used for all sites on a given branch.
     */
    class ConstBranchModelIterator
    {
      public:
        virtual ~ConstBranchModelIterator() {}

      public:
        virtual ConstBranchModelDescription* next() throw (Exception) = 0;
        virtual bool hasNext() const = 0;
    };

    /**
     * @brief A pair of SubstitutionModel / BranchIterator.
     */
    class ConstSiteModelDescription
    {
      public:
        virtual ~ConstSiteModelDescription() {}

      public:
        virtual const SubstitutionModel* getModel() const = 0;
        virtual BranchIterator* getNewBranchIterator() const = 0;
    };

    /**
     * @brief Iterates through all models used for all branches on a given site.
     */
    class ConstSiteModelIterator
    {
      public:
        virtual ~ConstSiteModelIterator() {}

      public:
        virtual ConstSiteModelDescription* next() throw (Exception) = 0;
        virtual bool hasNext() const = 0;
    };

  public:
    TreeLikelihood() {}
    virtual ~TreeLikelihood() {}

#ifndef NO_VIRTUAL_COV
    TreeLikelihood* clone() const = 0;
#endif

  public:

    /**
     * @brief Set the dataset for which the likelihood must be evaluated.
     *
     * @param sites The data set to use.
     */
    virtual void setData(const SiteContainer& sites) = 0;
    
    /**
     * @brief Get the dataset for which the likelihood must be evaluated.
     *
     * @return A pointer toward the site container where the sequences are stored.
     */
    virtual const SiteContainer* getData() const = 0;

    /**
     * @brief Init the likelihood object.
     *
     * This method is used to initialize all parameters.
     * It is typically called after the constructor and the setData method.
     * It contains virtual methods that can't be called in the constructor.
     * @throw Exception if something bad happened, for instance if no data are associated to the likelihood function.
     */
    virtual void initialize() throw (Exception) = 0;

    /**
     * @return 'true' is the likelihood function has been initialized.
     */
    virtual bool isInitialized() const = 0;

    /**
     * @return The underlying likelihood data structure.
     */
    virtual TreeLikelihoodData* getLikelihoodData() = 0;

    /**
     * @return The underlying likelihood data structure.
     */
    virtual const TreeLikelihoodData* getLikelihoodData() const = 0;

    /**
     * @brief Get the likelihood for a site.
     *
     * @param site The site index to analyse.
     * @return The likelihood for site <i>site</i>.
     */
    virtual double getLikelihoodForASite(size_t site) const = 0;

    /**
     * @brief Get the logarithm of the likelihood for a site.
     *
     * @param site The site index to analyse.
     * @return The logarithm of the likelihood for site <i>site</i>.
     */
    virtual double getLogLikelihoodForASite(size_t site) const = 0;

    /**
     * @brief Get the likelihood for a site and for a state.
     *
     * @param site The site index to analyse.
     * @param state The state to consider.
     * @return The likelihood for site <i>site</i> and state <i>state</i>.
     */
    virtual double getLikelihoodForASiteForAState(size_t site, int state) const = 0;

    /**
     * @brief Get the logarithm of the likelihood for a site and for a state.
     *
     * @param site The site index to analyse.
     * @param state The state to consider.
     * @return The logarithm of the likelihood for site <i>site</i> and state <i>state</i>.
     */
    virtual double getLogLikelihoodForASiteForAState(size_t site, int state) const = 0;

    /**
     * @brief Get the likelihood for each site.
     *
     * @return A vector with all likelihoods for each site.
     */
    virtual Vdouble getLikelihoodForEachSite() const = 0;

    /**
     * @brief Get the logarithm of the likelihood for each site.
     *
     * @return A vector with all log likelihoods for each site.
     */
    virtual Vdouble getLogLikelihoodForEachSite() const = 0;

    /**
     * @brief Get the likelihood for each site and for each state.
     *
     * @return A 2d vector with all likelihoods for each site and for each state.
     */
    virtual VVdouble getLikelihoodForEachSiteForEachState() const = 0;

    /**
     * @brief Get the logarithm of the likelihood for each site and for each state.
     *
     * @return A 2d vector with all log likelihoods for each site and for each state.
     */
    virtual VVdouble getLogLikelihoodForEachSiteForEachState() const = 0;
    
    /**
     * @brief Get the likelihood for the whole dataset.
     *
     * @return The likelihood of the dataset.
     */
    virtual double getLikelihood() const = 0;

    /**
     * @brief Get the logarithm of the likelihood for the whole dataset.
     *
     * @return The logarithm of the likelihood of the dataset.
     */
    virtual double getLogLikelihood() const = 0;
  
    /**
     * @brief Get the tree (topology and branch lengths).
     *
     * @return The tree of this TreeLikelihood object.
      */
    virtual const Tree& getTree() const = 0;

    /**
     * @brief Get the number of sites in the dataset.
     *
     * @return the number of sites in the dataset.
     */
    virtual size_t getNumberOfSites() const = 0;

    /**
     * @brief Get the number of states in the alphabet associated to the dataset.
     *
     * @return the number of states in the alphabet associated to the dataset.
     */    
    virtual size_t getNumberOfStates() const = 0;
    
    /**
     * @brief Get the alphabet associated to the dataset.
     *
     * @return the alphabet associated to the dataset.
     */    
    virtual const Alphabet* getAlphabet() const = 0;
   
    /**
     * @name Retrieve some particular parameters subsets.
     *
     * @{
     */
    
    /**
     * @brief Get the branch lengths parameters.
     *
     * @return A ParameterList with all branch lengths.
     */
    virtual ParameterList getBranchLengthsParameters() const = 0;
    
    /**
     * @brief Get the parameters associated to substitution model(s).
     *
     * @return A ParameterList.
     */
    virtual ParameterList getSubstitutionModelParameters() const = 0;

    /**
     * @brief Get the substitution model associated to a given node and alignment column.
     *
     * @param nodeId The id of the request node.
     * @param siteIndex The index of the alignment position.
     * @see getSiteIndex
     * @return A pointer toward the corresponding model.
     * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
     */
    virtual const SubstitutionModel* getSubstitutionModel(int nodeId, size_t siteIndex) const throw (NodeNotFoundException) = 0;

    /**
     * @brief Get the substitution model associated to a given node and alignment column.
     *
     * @param nodeId The id of the request node.
     * @param siteIndex The index of the alignment position.
     * @see getSiteIndex
     * @return A pointer toward the corresponding model.
     * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
     */
    virtual SubstitutionModel* getSubstitutionModel(int nodeId, size_t siteIndex) throw (NodeNotFoundException) = 0;

    /**
     * @brief Retrieves all Pij(t) for a particular branch, defined by the upper node and site.
     *
     * These intermediate results may be used by other methods.
     *
     * @param nodeId The node defining the branch of interest.
     * @param siteIndex The index of the alignment position.
     * @see getSiteIndex
     * @return An array of dimension 2, where a[x][y] is the probability of substituting from x to y.
     */
    virtual VVdouble getTransitionProbabilities(int nodeId, size_t siteIndex) const = 0;

    virtual ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const = 0;

    virtual ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const = 0;

    /**
     * @brief Get the index (used for inner computations) of a given site (original alignment column).
     *
     * @param site An alignment position.
     * @return The site index corresponding to the given input alignment position.
     */
    virtual size_t getSiteIndex(size_t site) const throw (IndexOutOfBoundsException) = 0;

    /**
     * @brief Get the values of the frequencies for each state in the alphabet at the root node.
     *
     * For reversible models, these are the equilibrium frequencies.
     * For non-reversible models, these usually are distinct parameters.
     *
     * For models without site partitioning, the set of frequencies is the same for all positions.
     * For partition models, the frequencies may differ from one site to another.
     *
     * @param siteIndex The index of the alignment position.
     * @see getSiteIndex
     * @return A vector with ancestral frequencies for each state in the alphabet;
     */
    virtual const std::vector<double>& getRootFrequencies(size_t siteIndex) const = 0;
    
    /** @} */

    /**
     * @brief Tell if derivatives must be computed.
     *
     * This methods calls the enableFirstOrderDerivatives and enableSecondOrderDerivatives.
     *
     * @param yn Yes or no.
     */
    virtual void enableDerivatives(bool yn) = 0;

    /**
     * @brief All derivable parameters.
     *
     * Usually, this contains all branch lengths parameters.
     *
     * @return A ParameterList.
     */
    virtual ParameterList getDerivableParameters() const = 0;

    /**
     * @brief All non derivable parameters.
     *
     * Usually, this contains all substitution model parameters and rate distribution.
     *
     * @return A ParameterList.
     */
    virtual ParameterList getNonDerivableParameters() const = 0;

};

} //end of namespace bpp.

#endif  //_TREELIKELIHOOD_H_

