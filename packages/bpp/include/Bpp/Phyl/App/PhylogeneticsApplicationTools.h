//
// File: PhylogeneticsApplicationTools.h
// Created by: Julien Dutheil
// Created on: Fri Oct 21 16:49 2005
// from old file ApplicationTools.h created on Sun Dec 14 09:36:26 2003
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

#ifndef _PHYLOGENETICSAPPLICATIONTOOLS_H_
#define _PHYLOGENETICSAPPLICATIONTOOLS_H_

#include "../Tree.h"
#include "../Model/SubstitutionModel.h"
#include "../Model/SubstitutionModelSet.h"
#include "../Model/MixedSubstitutionModelSet.h"
#include "../Model/MarkovModulatedSubstitutionModel.h"
#include "../Likelihood/HomogeneousTreeLikelihood.h"
#include "../Likelihood/ClockTreeLikelihood.h"
#include "../Mapping/SubstitutionCount.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Io/OutputStream.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/MultipleDiscreteDistribution.h>
#include <Bpp/Numeric/Function/Optimizer.h>

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// From the STL:
#include <string>
#include <map>

namespace bpp
{


  /**
   * @brief This class provides some common tools for applications.
   *
   * The functions parse some option file, create corresponding objects and send
   * a pointer toward it.
   * 
   * The option files are supposed to follow this simple format:
   * @code
   * parameterName = parameterContent
   * @endcode
   * with one parameter per line.
   *
   * @see ApplicationTools
   */
  class PhylogeneticsApplicationTools
  {
  
  public:
    PhylogeneticsApplicationTools();
    virtual ~PhylogeneticsApplicationTools();
  

    /**
     * @brief Build a Tree object according to options.
     *
     * See the Bio++ Program Suite manual for a description of available options.
     *
     * @param params  The attribute map where options may be found.
     * @param prefix  A prefix to be applied to each attribute name.
     * @param suffix  A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @return A new Tree object according to the specified options.
     * @throw Exception if an error occured.
     */
    static Tree* getTree(
                         std::map<std::string, std::string>& params,
                         const std::string& prefix = "input.",
                         const std::string& suffix = "",
                         bool suffixIsOptional = true,
                         bool verbose = true) throw (Exception);
 
    /**
     * @brief Build a list ofTree objects according to options.
     *
     * See the Bio++ Program Suite manual for a description of available options.
     *
     * @param params  The attribute map where options may be found.
     * @param prefix  A prefix to be applied to each attribute name.
     * @param suffix  A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @return A new vector of Tree objects according to the specified options.
     * @throw Exception if an error occured.
     */
    static std::vector<Tree*> getTrees(
                                       std::map<std::string, std::string>& params,
                                       const std::string& prefix = "input.",
                                       const std::string& suffix = "",
                                       bool suffixIsOptional = true,
                                       bool verbose = true) throw (Exception);
  
    /**
     * @brief Build a SubstitutionModel object according to options.
     *
     * Creates a new substitution model object according to model description syntax
     * (see the Bio++ Progam Suite manual for a detailed description of this syntax). The
     * function also parses the parameter values and set them accordingly.
     *
     * @param alphabet The alphabet to use in the model.
     * @param gCode    The genetic code to use (only for codon models, otherwise can be set to 0).
     *                 If set to NULL and a codon model is requested, an Exception will be thrown.
     * @param data     A pointer toward the SiteContainer for which the substitution model is designed.
     *                 The alphabet associated to the data must be of the same type as the one specified for the model.
     *                 May be equal to NULL, but in this case use_observed_freq option will be unavailable.
     * @param params   The attribute map where options may be found.
     * @param suffix   A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @return A new SubstitutionModel object according to options specified.
     * @throw Exception if an error occured.
     */
    static SubstitutionModel* getSubstitutionModel(
        const Alphabet* alphabet,
        const GeneticCode* gCode,
        const SiteContainer* data, 
        std::map<std::string, std::string>& params,
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool verbose = true) throw (Exception);
  

    /**
     * @brief Set parameter initial values of a given model in a set according to options.
     *
     * Parameters actually depends on the model passed as argument.
     * See getSubstitutionModelSet for more information.
     *
     * This function is mainly for internal usage, you're probably looking for the getSubstitutionModel or getSubstitutionModelSet function.
     *
     * @param model                   The model to set.
     * @param unparsedParameterValues A map that contains all the model parameters
     *                                names and their corresponding unparsed value, if they were found.
     * @param modelNumber The number of this model in the SubstitutionModelSet.
     * @param data   A pointer toward the SiteContainer for which the substitution model is designed.
     *               The alphabet associated to the data must be of the same type as the one specified for the model.
     *               May be equal to NULL, but in this case use_observed_freq option will be unavailable.
     * @param existingParams (in/out) A map with already existing value that have been found in previous calls, and may be recalled here.
     *                       New parameters found here will be added.
     * @param sharedParams (out) remote parameters will be recorded here.
     * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
     */
    static void setSubstitutionModelParametersInitialValuesWithAliases(
                                                                       SubstitutionModel& model,
                                                                       std::map<std::string, std::string>& unparsedParameterValues,
                                                                       size_t modelNumber,
                                                                       const SiteContainer* data,
                                                                       std::map<std::string, double>& existingParams,
                                                                       std::map<std::string, std::string>& sharedParams,
                                                                       bool verbose) throw (Exception);

    /**
     * @brief Get A FrequenciesSet object for root frequencies (NH models) according to options.
     *
     * @param alphabet The alpabet to use.
     * @param gCode    The genetic code to use (only for codon alphabets, otherwise can be set to 0).
     *                 If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
     * @param data      A pointer toward the SiteContainer for which the substitution model is designed.
     *                  The alphabet associated to the data must be of the same type as the one specified for the model.
     *                  May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param params    The attribute map where options may be found.
     * @param rateFreqs A vector of rate catégories frequencies in case of a Markov Modulated Markov Model.
     *                  Ignored if a vector with size 0 is passed.
     * @param suffix    A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose   Print some info to the 'message' output stream.
     * @return A new FrequenciesSet object according to options specified.
     * @throw Exception if an error occured.
     */
    static FrequenciesSet* getRootFrequenciesSet(
        const Alphabet* alphabet,
        const GeneticCode* gCode,
        const SiteContainer* data, 
        std::map<std::string, std::string>& params,
        const std::vector<double>& rateFreqs,
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool verbose = true) throw (Exception);

    /**
     * @brief Get A FrequenciesSet object according to options.
     *
     * @param alphabet The alpabet to use.
     * @param gCode    The genetic code to use (only for codon alphabets, otherwise can be set to 0).
     *                 If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
     * @param freqDescription A string in the keyval syntaxe describing the frequency set to use.:if expand("%") == ""|browse confirm w|else|confirm w|endif
     * 
     * @param data      A pointer toward the SiteContainer for which the substitution model is designed.
     *                  The alphabet associated to the data must be of the same type as the one specified for the model.
     *                  May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param rateFreqs A vector of rate catégories frequencies in case of a Markov Modulated Markov Model.
     *                  Ignored if a vector with size 0 is passed.
     * @param verbose   Print some info to the 'message' output stream.
     * @return A new FrequenciesSet object according to options specified.
     * @throw Exception if an error occured.
     */
    static FrequenciesSet* getFrequenciesSet(
        const Alphabet* alphabet,
        const GeneticCode* gCode,
        const std::string& freqDescription,
        const SiteContainer* data, 
        const std::vector<double>& rateFreqs,
        bool verbose = true)
      throw (Exception);

    /**
     * @brief Gets a SubstitutionModelSet object according to options.
     *
     * See setSubstitutionModelSet and setMixedSubstitutionModelSet
     * methods.
     */
     static SubstitutionModelSet* getSubstitutionModelSet(
         const Alphabet* alphabet,
         const GeneticCode* gcode,
         const SiteContainer* data, 
         std::map<std::string, std::string>& params,
         const std::string& suffix = "",
         bool suffixIsOptional = true,
         bool verbose = true);    

     /**
     * @brief Sets a SubstitutionModelSet object according to options.
     *
     * This model set is meant to be used with non-homogeneous substitution models of sequence evolution.
     *
     * Recognized options are:
     * - number_of_models: the number of distinct SubstitutionModel to use.
     *
     * Then, for each of the models, the following information must be provided:
     * - model1='model name(parameters'='value',...)
     * Model names and parameters follow the same syntaxe as for the getSubstitutionModel method.
     * - model1.nodes='list of nodes id, separated by comas'.
     * And then
     * - model2=...
     * etc.
     *
     * All models must be fully specified, and at the end of the description, all nodes must be attributed to a model,
     * otherwise an exception is thrown.
     * 
     * Finally, this is also allowed for models to share one or several parameters.
     * for instance:
     * @code
     * model1=T92(kappa=2.0, theta=0.5)
     * model2=T92(kappa=model1.kappa, theta=0.5)
     * @endcode
     * In this case model1 and model2 with have their own and independent theta parameter, but only one kappa parameter will be used for both models.
     * Note that
     * @code
     * model1=T92(kappa=2.0, theta=0.5)
     * model1.nodes=1,2,3
     * model2=T92(kappa= model1.kappa, theta=model1.theta)
     * model2.nodes=4,5,6
     * @endcode
     * is equivalent to
     * @code
     * model1=T92(kappa=2.0, theta=0.5)
     * model1.nodes=1,2,3,4,5,6
     * @endcode
     * but will require more memory and use more CPU, since some calculations will be performed twice.
     *
     * @param modelSet The modified SubstitutionModelSet object according to options specified.
     * @param alphabet The alpabet to use in all models.
     * @param gcode    The genetic code to use (only for codon models, otherwise can be set to 0).
     *                 If set to NULL and a codon model is requested, an Exception will be thrown.
     * @param data     A pointer toward the SiteContainer for which the substitution model is designed.
     *                  The alphabet associated to the data must be of the same type as the one specified for the model.
     *                 May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param params   The attribute map where options may be found.
     * @param suffix   A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
     */
    static void setSubstitutionModelSet(
        SubstitutionModelSet& modelSet,
        const Alphabet* alphabet,
        const GeneticCode* gcode,
        const SiteContainer* data, 
        std::map<std::string, std::string>& params,
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool verbose = true);
    
    /**
     * @brief Complete a MixedSubstitutionModelSet object according to
     * options, given this model has already been filled through
     * setSubstitutionModelSet method.
     *
     * In addition, this method builds the allowed combinations of
     * submodels of the different mixed models.
     *
     * If none combination is given, then all possible submodels
     * combinations will be considered.
     *
     * The submodels dependencies are given a sets of combinations of
     * the mixed variables of the mixed models. For instance, if we
     * have:
     *
     * @code
     * model1=MixedModel(model=T92(kappa=Gamma(n=4), theta=0.5))
     * model2=MixedModel(model=T92(kappa=Gaussian(n=5), theta=Beta(n=3)))
     * @endcode
     *
     * In this case model1 is a mixture of 4 T92 submodels and model2
     * a mixture of 15 T92 submodels. These submodels are denoted with
     * the parameter name and the class number. For example, the
     * submodels of model1 are denoted model1[kappa_1], ...,
     * model1[kappa_4], and the submodels of model2 are denoted
     * model2[kappa_1,theta_1], ..., model2[kappa_5, theta_3].
     * Additionnaly, for instance, model2[kappa_2] denotes all the
     * submodels whose description has kappa_2.
     *
     * By default, when switching from model1 to model2, a site is
     * allowed to switch between any submodel of model1 and any
     * submodel of model2. If the only allowed combination is that a
     * site follows submodels model1(kappa_1) and
     * model2(kappa_3,theta_2), it is denoted:
     *
     * @code
     * site.allowedpaths= model1[kappa_1] & model2[kappa_3,theta_2]
     * @endcode
     *
     * With additional combination saying that a site can follow
     * submodels model1[kappa_2] and any submodel of model2[kappa_3]
     * is denoted:
     *
     * @code
     * site.allowedpaths= model1[kappa_1] & model2[kappa_3,theta_2] |
     *                    model1[kappa_2] & model2[kappa_3]
     * @endcode
     *
     * See MixedSubstitutionModelSet description for further
     * information.
     *
     * @param mixedModelSet The modified MixedSubstitutionModelSet object according to options specified.
     * @param alphabet The alpabet to use in all models.
     * @param data     A pointer toward the SiteContainer for which the substitution model is designed.
     *                  The alphabet associated to the data must be of the same type as the one specified for the model.
     *                 May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param params   The attribute map where options may be found.
     * @param suffix   A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
     */
    static void completeMixedSubstitutionModelSet(
        MixedSubstitutionModelSet& mixedModelSet,
        const Alphabet* alphabet,
        const SiteContainer* data, 
        std::map<std::string, std::string>& params,
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool verbose = true);

    /**
     * @brief Build a multi-dimension distribution as a
     * MultipleDiscreteDistribution object with default parameter
     * values according to a keyval description.
     *
     * Check the Bio++ Program Suite documentation for a description of the syntax.
     * It is mainly for internal usage, you're probably looking for the getRateDistribution function.
     *
     * @param distDescription         A string describing the model in the keyval syntax.
     * @param unparsedParameterValues [out] a map that will contain all the distribution parameters
     *                                names and their corresponding unparsed value, if they were found.
     * @param verbose                 Print some info to the 'message' output stream.
     * @return A new MultipleDiscreteDistribution object according to options specified.
     * @throw Exception if an error occured.
     */
    
    static MultipleDiscreteDistribution* getMultipleDistributionDefaultInstance(
        const std::string& distDescription,
        std::map<std::string, std::string>& unparsedParameterValues,
        bool verbose = true);

    /**
     * @brief Build a DiscreteDistribution object according to options.
     *
     * Creates a new rate distribution object according to distribution description syntax
     * (see the Bio++ Progam Suite manual for a detailed description of this syntax). The
     * function also parses the parameter values and set them accordingly.
     *
     * @param params  The attribute map where options may be found.
     * @param suffix  A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @return A new DiscreteDistribution object according to options specified.
     * @throw Exception if an error occured.
     */
    static DiscreteDistribution* getRateDistribution(
        std::map<std::string, std::string>& params,
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool verbose = true)
      throw (Exception);
      
    /**
     * @brief Optimize parameters according to options.
     *
     * Options used are:
     * - optimization = Tell if optimization must be performed.
     * - optimization.message_handler = [std, file_path]
     *   A path to a specific path (existing will be overwritten) or std for use
     *   of the standard output.
     * - optimization.profiler = [std, file_path], idem for the profiling (history
     *   of all functions evaluations).
     * - optimization.max_number_f_eval = The maximum number of function evaluation.
     * - optimization.tolerance = The tolerance parameter (when to stop the optimization).
     * - optimization.scale_first = Tell if we must scale the tree first.
     * - optimization.ignore_parameter = A coma-separated list of parameter
     *   names to ignore in the optimizing process.
     * - optimization.method = [DB|fullD] Algorithm to use: Derivatives+Brent or full derivatives, with numerical derivatives when required.
     * - optimization.method.derivatives = [gradient|newton] Use Conjugate Grandient or Newton-Rhaphson algorithm.
     * - optimization.final = [none|simplex|powell] Perform a downhill simplex or a Powell multidimensions optimization
     * - optimization.topology = Tell if we must optimize tree topology. Toplogy estimation uses the DB algorithm with Newton-Raphson during estimation.
     *   The previous options will be used only for final estimation of numerical parameters.
     * Options depending on other options:
     * - If optimization.scale_first is set to true:
     *   - optimization.scale_first.tolerance = The tolerance of the scaling alogrithm.
     *   - optimization.scale_first.max_number_f_eval = the maximum number of function evaluations
     *     for the scaling algorithm.
     * - optimization.method_DB.nstep = number of progressive steps to use in DB algorithm.
     * - optimization.topology.algorithm = [nni] algorithm to use (for now, only Nearest Neighbor Interchanges (NNI) are implemented). 
     * - optimization.topology.algorithm_nni.method = [fast,better,phyml]
     * - optimization.topology.nstep = Estimate numerical parameters every 'n' NNI rounds.
     * - optimization.topology.numfirst = Tell if numerical parameters must be estimated prior to topology search.
     * - optimization.topology.tolerance.before = Numerical parameters estimation prior to topology search.
     * - optimization.topology.tolerance.during = Numerical parameters estimation during topology search.
     *
     * @param tl               The TreeLikelihood function to optimize.
     * @param parameters       The initial list of parameters to optimize.
     *                         Use tl->getIndependentParameters() in order to estimate all parameters.
     * @param params           The attribute map where options may be found.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @throw Exception        Any exception that may happen during the optimization process.
     * @return A pointer toward the final likelihood object.
     * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
     * clone this object. We may change this bahavior in the future...
     * You hence should write something like
     * @code
     * tl = PhylogeneticsApplicationTools::optimizeParameters(tl, ...);
     * @endcode
     */
    static TreeLikelihood* optimizeParameters(
        TreeLikelihood* tl,
        const ParameterList& parameters,
        std::map<std::string, std::string>& params,
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool verbose = true)
      throw (Exception);
    
    /**
     * @brief Optimize parameters according to options, with a molecular clock.
     *
     * Options used are:
     * - optimization = Tell if optimization must be performed.
     * - optimization.message_handler = [std, file_path]
     *   A path to a specific path (existing will be overwritten) or std for use
     *   of the standard output.
     * - optimization.profiler = [std, file_path], idem for the profiling (history
     *   of all functions evaluations).
     * - optimization.max_number_f_eval = The maximum number of function evaluation.
     * - optimization.tolerance = The tolerance parameter (when to stop the optimization).
     * - optimization.ignore_parameter = A coma-separated list of parameter
     *   names to ignore in the optimizing process.
     * - optimization.method = [DB|fullD] Algorithm to use: Derivatives+Brent or full derivatives, with numerical derivatives when required.
     * - optimization.method.derivatives = [gradient|newton] Use Conjugate Grandient or Newton-Rhaphson algorithm.
     * - optimization.final = [none|simplex|powell] Perform a downhill simplex or a Powell multidimensions optimization
     * - optimization.method_DB.nstep = number of progressive steps to use in DB algorithm.
     *
     * @param tl               The ClockTreeLikelihood function to optimize.
     * @param parameters       The initial list of parameters to optimize.
     *                         Use tl->getIndependentParameters() in order to estimate all parameters.
     * @param params           The attribute map where options may be found.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @throw Exception        Any exception that may happen during the optimization process.
     */
    static void optimizeParameters(
        DiscreteRatesAcrossSitesClockTreeLikelihood* tl,
        const ParameterList& parameters,
        std::map<std::string, std::string>& params,
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool verbose = true)
      throw (Exception);
    
    /**
     * @brief Check if parameter values are close to their definition boundary.
     *
     * This allows the detection of potential optimization issues.
     * A warning message will be output for each problematic parameter.
     *
     * @param pl A list of parameters. Parameters without constraint will be ignored.
     */
    static void checkEstimatedParameters(const ParameterList& pl);

    /**
     * @brief Get a SubstitutionCount instance.
     *
     * @param alphabet The alphabet to use.
     * @param model The model to use.
     * @param params The attribute map where options may be found.
     * @param suffix Optional suffix for command name.
     */
    static SubstitutionCount* getSubstitutionCount(
        const Alphabet* alphabet,
        const SubstitutionModel* model,
        map<string, string>& params,
        string suffix = "");
   
    /**
     * @brief Write a tree according to options.
     *
     * See the Bio++ Program Suite manual for a descriptio of all available options.
     *
     * @param tree    The tree to write.
     * @param params  The attribute map where options may be found.
     * @param prefix  A prefix to be applied to each attribute name.
     * @param suffix  A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @param checkOnly If this parameter is set to true, then all options are
     * checked and error messages sent, but no file is written.
     * @throw Exception if an error occured.
     */
    static void writeTree(
        const TreeTemplate<Node>& tree,
        std::map<std::string, std::string>& params,
        const std::string& prefix = "output.",
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool verbose = true,
        bool checkOnly = false) throw (Exception);
    
    /**
     * @brief Write a tree according to options.
     *
     * See the Bio++ Program Suite manual for a descriptio of all available options.
     *
     * @param trees   The trees to write.
     * @param params  The attribute map where options may be found.
     * @param prefix  A prefix to be applied to each attribute name.
     * @param suffix  A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @param checkOnly If this parameter is set to true, then all options are
     * checked and error messages sent, but no file is written.
     * @throw Exception if an error occured.
     */
    static void writeTrees(
        const std::vector<Tree*>& trees,
        std::map<std::string, std::string>& params,
        const std::string& prefix = "output.",
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool verbose = true,
        bool checkOnly = false) throw (Exception);


    
    /**
     * @brief Output a SubstitutionModel description to a file.
     *
     * @param model The model to serialize.
     * @param out   The stream where to print.
     */
    static void printParameters(const SubstitutionModel* model, OutputStream& out);



    /**
     * @brief Output a SubstitutionModelSet description to a file.
     *
     * @param modelSet The model set to serialize.
     * @param out      The stream where to print.
     */
    static void printParameters(const SubstitutionModelSet* modelSet, OutputStream& out);



    /**
     * @brief Output a DiscreteDistribution description to a file.
     *
     * @param rDist The rate distribution to serialize.
     * @param out   The stream where to print.
     */
    static void printParameters(const DiscreteDistribution* rDist, OutputStream& out);

  };

} //end of namespace bpp.

#endif  //_PHYLOGENETICSAPPLICATIONTOOLS_H_

