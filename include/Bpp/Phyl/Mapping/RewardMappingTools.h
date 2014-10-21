//
// File: RewardMappingTools.h
// Created by: Laurent Guéguen
// Created on: vendredi 29 mars 2013, à 14h 08
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

#ifndef _REWARDMAPPINGTOOLS_H_
#define _REWARDMAPPINGTOOLS_H_

#include "ProbabilisticRewardMapping.h"
#include "Reward.h"

#include "../Likelihood/DRTreeLikelihood.h"

namespace bpp
{

  /**
   * @brief Provide methods to compute reward mappings.
   *
   * For now, 4 methods are implemented, and provide reward mappings.
   *
   * See:
   * Minin, V.N. and Suchard, M.A., 
   * Fast, accurate and simulation-free stochastic mapping
   * Philosophical Transactions of the Royal Society B 2008 363:3985-95.
   *
   * @author Laurent Guéguen
   */
  class RewardMappingTools
  {
  public:
    RewardMappingTools() {}
    virtual ~RewardMappingTools() {}
    
  public:
		
    /**
     * @brief Compute the reward vectors for a particular dataset
     * using the double-recursive likelihood computation.
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The Ids of the nodes the reward vectors
     *                                are computed on.
     * @param reward            The Reward to use.
     * @param verbose           Print info to screen.
     * @return A vector of reward vectors (one for each site).
     * @throw Exception If the likelihood object is not initialized.
     */
    static ProbabilisticRewardMapping* computeRewardVectors(
                                                            const DRTreeLikelihood& drtl,
                                                            const std::vector<int>& ids,
                                                            Reward& reward,
                                                            bool verbose = true) throw (Exception);
		
    
    /**
     * @brief Write the reward vectors to a stream.
     *
     * @param rewards The reward vectors to write.
     * @param sites         The dataset associated to the vectors
     * (needed to know the position of each site in the dataset).
     * @param out           The output stream where to write the vectors.
     * @throw IOException If an output error happens.
     */
    static void writeToStream(
                              const ProbabilisticRewardMapping& rewards,
                              const SiteContainer& sites,
                              std::ostream& out)
      throw (IOException);
	

    /**
     * @brief Read the reward vectors from a stream.
     *
     * @param in            The input stream where to read the vectors.
     * @param rewards       The mapping object to fill.
     * @throw IOException If an input error happens.
     */
    static void readFromStream(std::istream& in, ProbabilisticRewardMapping& rewards)
      throw (IOException);


    /**
     * @brief Sum all rewards of a given branch (specified by its index).
     *
     * @param smap The reward map to use.
     * @param branchIndex The index of the reward vector for which the counts should be computed.
     * @return A vector will all rewards summed.
     */
    static double computeSumForBranch(const RewardMapping& smap, size_t branchIndex);
 

    /**
     * @brief Sum all substitutions for each type of a given site (specified by its index).
     *
     * @param smap The substitution map to use.
     * @param siteIndex The index of the substitution vector for which the counts should be computed.
     * @return A vector will all counts summed for each types of substitutions. 
     */
    static double computeSumForSite(const RewardMapping& smap, size_t siteIndex);
  };

} //end of namespace bpp.

#endif //_REWARDMAPPINGTOOLS_H_

