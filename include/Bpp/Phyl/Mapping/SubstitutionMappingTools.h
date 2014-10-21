//
// File: SubstitutionMappingTools.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 13:04 2006
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

#ifndef _SUBSTITUTIONMAPPINGTOOLS_H_
#define _SUBSTITUTIONMAPPINGTOOLS_H_

#include "ProbabilisticSubstitutionMapping.h"
#include "SubstitutionCount.h"
#include "OneJumpSubstitutionCount.h"
#include "../Likelihood/DRTreeLikelihood.h"

namespace bpp
{

  /**
   * @brief Provide methods to compute substitution mappings.
   *
   * For now, only 4 methods are implemented, and provide probabilistic substitution mappings.
   *
   * See:
   * Dutheil J, Pupko T, Jean-Marie A, Galtier N.
   * A model-based approach for detecting coevolving positions in a molecule.
   * Mol Biol Evol. 2005 Sep;22(9):1919-28. Epub 2005 Jun 8.
   *
   * @author Julien Dutheil
   */
  class SubstitutionMappingTools
  {
  public:
    SubstitutionMappingTools() {}
    virtual ~SubstitutionMappingTools() {}
    
  public:
		
    /**
     * @brief Compute the substitutions vectors for a particular dataset using the
     * double-recursive likelihood computation.
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The Ids of the nodes the substitutions
     *                                are counted on.
     * @param substitutionCount The SubstitutionCount to use.
     * @param verbose           Print info to screen.
     * @return A vector of substitutions vectors (one for each site).
     * @throw Exception If the likelihood object is not initialized.
     */
    static ProbabilisticSubstitutionMapping* computeSubstitutionVectors(
        const DRTreeLikelihood& drtl,
        const std::vector<int>& ids,
        SubstitutionCount& substitutionCount,
        bool verbose = true) throw (Exception);
		
    /**
     * @brief Compute the substitutions vectors for a particular dataset using the
     * double-recursive likelihood computation.
     *
     * In this method, substitution counts are computed using the pair of ancestral
     * states with maximum likelihood.
     * This is a kind of joint-pair ancestral reconstruction, as in Galtier and Boursot (1998).
     * This reconstruction possibly takes into account several rate classes, and
     * substitution counts are averaged over all rate classes, weighted by their conditional
     * likelihood.
     *
     * This function is mainly for testing purpose (see Dutheil et al. 2005).
     * For practical use, consider using the 'getSubstitutionVectors' method instead.
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param substitutionCount The substitutionsCount to use.
     * @param verbose           Print info to screen.
     * @return A vector of substitutions vectors (one for each site).
     * @throw Exception If the likelihood object is not initialized.
     */
    static ProbabilisticSubstitutionMapping* computeSubstitutionVectorsNoAveraging(
        const DRTreeLikelihood& drtl,
        SubstitutionCount& substitutionCount,
        bool verbose = true) throw (Exception);

		
    /**
     * @brief Compute the substitutions vectors for a particular dataset using the
     * double-recursive likelihood computation.
     *
     * In this method, all ancestral states are estimated using marginal likelihoods,
     * putatively intregated over several rate classes.
     * For each branch, the number of substitution given marginal states is used.
     * This method, used with a SimpleSubstitutionCount objet is equivalent to
     * Tufféry and Darlu's (2000) computation of substitution vectors.
     *
     * Use with another substitution count objet is in most cases irrelevent.
     * 
     * @param drtl              A DRTreeLikelihood object.
     * @param substitutionCount The substitutionsCount to use.
     * @param verbose           Print info to screen.
     * @return A vector of substitutions vectors (one for each site).
     * @throw Exception If the likelihood object is not initialized.
     */
    static ProbabilisticSubstitutionMapping* computeSubstitutionVectorsNoAveragingMarginal(
        const DRTreeLikelihood& drtl,
        SubstitutionCount& substitutionCount,
        bool verbose = true) throw (Exception);

		
    /**
     * @brief Compute the substitutions vectors for a particular dataset using the
     * double-recursive likelihood computation.
     *
     * The marginal probability is used for weighting, i.e. the product of probabilities for the pair.
     *
     * This function is mainly for testing purpose (see Dutheil et al. 2005).
     * For practical use, consider using the 'getSubstitutionVectors' method instead.
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param substitutionCount The substitutionsCount to use.
     * @param verbose           Print info to screen.
     * @return A vector of substitutions vectors (one for each site).
     * @throw Exception If the likelihood object is not initialized.
     */
    static ProbabilisticSubstitutionMapping* computeSubstitutionVectorsMarginal(
        const DRTreeLikelihood& drtl,
        SubstitutionCount& substitutionCount,
        bool verbose = true) throw (Exception);
	

    /**
     * @brief This method computes for each site and for each branch the probability that
     * at least one jump occurred.
     *
     * Here 'jump' refer to a change in the model state. Depending on the model, this might
     * not be the same as a substitution (an alphabet state change).
     */
    static ProbabilisticSubstitutionMapping* computeOneJumpProbabilityVectors(
                                                                              const DRTreeLikelihood& drtl,
                                                                              bool verbose = true) throw (Exception)
    {
      OneJumpSubstitutionCount ojsm(0);
      return computeSubstitutionVectors(drtl, drtl.getTree().getNodesId(), ojsm, 0);
    }


    /**
     * @brief Write the substitutions vectors to a stream.
     *
     * @param substitutions The substitutions vectors to write.
     * @param sites         The dataset associated to the vectors
     * (needed to know the position of each site in the dataset).
     * @param type          The type of substitutions to be output. See SubstitutionCount class.
     * Only one type of substitution can be output at a time.
     * @param out           The output stream where to write the vectors.
     * @throw IOException If an output error happens.
     */
    static void writeToStream(
                              const ProbabilisticSubstitutionMapping& substitutions,
                              const SiteContainer& sites,
                              size_t type,
                              std::ostream& out)
      throw (IOException);
	

    /**
     * @brief Read the substitutions vectors from a stream.
     *
     * @param in            The input stream where to read the vectors.
     * @param substitutions The mapping object to fill.
     * @param type          The type of substitutions that are read. Should be in supported by the substittuion count obect assiciated to the mapping, if any.
     * @throw IOException If an input error happens.
     */
    static void readFromStream(std::istream& in, ProbabilisticSubstitutionMapping& substitutions, size_t type)
      throw (IOException);


    /**
     * @brief Sum all type of substitutions for each branch of a given position (specified by its index).
     *
     * @param smap The substitution map to use.
     * @param siteIndex The index of the substitution vector for which the counts should be computed.
     * @return A vector will all counts for all types of substitutions summed. 
     */
    static std::vector<double> computeTotalSubstitutionVectorForSite(const SubstitutionMapping& smap, size_t siteIndex);


    /**
     * @brief Compute the norm of a substitution vector for a given position (specified by its index).
     *
     * The norm is computed as:
     * @f$ N_i = \sqrt{\left(\sum_l {\left(\sum_t n_{l, i, t}\right)}^2\right)}@f$,
     * where @f$n_{l, i, t}@f$ is the number of substitutions of type t on site i on branch l, obtained using the () operator for the SubstitutionMapping object.
     *
     * @param smap The substitution map to use.
     * @param siteIndex The index of the substitution vector for which the norm should be computed.
     * @return The norm of the substitution vector.
     */
    static double computeNormForSite(const SubstitutionMapping& smap, size_t siteIndex);
    

    /**
     * @brief Sum all substitutions for each type of a given branch (specified by its index).
     *
     * @param smap The substitution map to use.
     * @param branchIndex The index of the substitution vector for which the counts should be computed.
     * @return A vector will all counts summed for each types of substitutions. 
     */
    static std::vector<double> computeSumForBranch(const SubstitutionMapping& smap, size_t branchIndex);
 

    /**
     * @brief Sum all substitutions for each type of a given site (specified by its index).
     *
     * @param smap The substitution map to use.
     * @param siteIndex The index of the substitution vector for which the counts should be computed.
     * @return A vector will all counts summed for each types of substitutions. 
     */
    static std::vector<double> computeSumForSite(const SubstitutionMapping& smap, size_t siteIndex);


    /**
     * @brief Returns the counts on each branch.
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param model             The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @param threshold         value above which counts are considered saturated
     *                                        (default: -1 means no threshold).
     * @return A vector of substitutions vectors (one per branch per type).
     */
    static std::vector< std::vector<double> > getCountsPerBranch(
        DRTreeLikelihood& drtl,
        const std::vector<int>& ids,
        SubstitutionModel* model,
        const SubstitutionRegister& reg,
        double threshold = -1,
        bool verbose = true);


    /**
     * @brief Returns the normalization factors due to the null model
     * on each branch, for each register
     *
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param nullModel         The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @return A vector of normalization vectors (one per branch per type).
     */
    static std::vector< std::vector<double> > getNormalizationsPerBranch(
        DRTreeLikelihood& drtl,
        const std::vector<int>& ids,
        const SubstitutionModel* nullModel,
        const SubstitutionRegister& reg,
        bool verbose = true);


    /**
     * @brief Returns the normalization factors due to the set of null
     * models on each branch, for each register.
     *
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param nullModelSet      The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @return A vector of normalization vectors (one per branch per type).
     */
    static std::vector< std::vector<double> > getNormalizationsPerBranch(
        DRTreeLikelihood& drtl,
        const std::vector<int>& ids,
        const SubstitutionModelSet* nullModelSet,
        const SubstitutionRegister& reg,
        bool verbose = true);

    
    /**
     * @brief Returns the counts normalized by a null model
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param model             The model on which the SubstitutionCount is built
     * @param nullModel         The null model used for normalization.
     * @param reg               the Substitution Register
     * @param complete          boolean if the counts should be
     *                              completed to sum 1.
     *
     */
    
    static std::vector< std::vector<double> >  getNormalizedCountsPerBranch(
        DRTreeLikelihood& drtl,
        const std::vector<int>& ids,
        SubstitutionModel* model,
        SubstitutionModel* nullModel,
        const SubstitutionRegister& reg,
        bool complete = true,
        bool verbose = true);

    /**
     * @brief Returns the counts normalized by a null model set
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param modelSet          The model set on which the SubstitutionCount is built
     * @param nullModelSet      The null model set used for normalization.
     * @param reg               the Substitution Register
     * @param complete          boolean if the counts should be
     *                              completed to sum 1.
     *
     */
    
    static std::vector< std::vector<double> > getNormalizedCountsPerBranch(
        DRTreeLikelihood& drtl,
        const std::vector<int>& ids,
        SubstitutionModelSet* modelSet,
        SubstitutionModelSet* nullModelSet,
        const SubstitutionRegister& reg,
        bool complete = true,
        bool verbose = true);

    /**
     * @brief Returns the counts relative to the frequency of the
     * states in case of non-stationarity.
     *
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param model             The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     * @param stationarity      if false, a correction is made if the SubstitutionRegister
     *                             is a CategorySubstitutionRegister
     * @param threshold         value above which counts are considered saturated
     *                                        (default: -1 means no threshold).
     *
     */
    static std::vector< std::vector<double> >  getRelativeCountsPerBranch(
        DRTreeLikelihood& drtl,
        const std::vector<int>& ids,
        SubstitutionModel* model,
        const SubstitutionRegister& reg,
        bool stationarity = true,
        double threshold = -1);

    /**
     * @brief Output the sum of the counts par branch per site, in a
     * file.
     *
     * @param filename          The name of the output file
     * @param drtl              A DRTreeLikelihood object.
     * @param ids               The numbers of the nodes of the tree
     * @param model             The model on which the SubstitutionCount is built
     * @param reg               the Substitution Register
     *
     */
    static void outputTotalCountsPerBranchPerSite(
        std::string& filename,
        DRTreeLikelihood& drtl,
        const std::vector<int>& ids,
        SubstitutionModel* model,
        const SubstitutionRegister& reg);

   

  };

} //end of namespace bpp.

#endif //_SUBSTITUTIONMAPPINGTOOLS_H_

