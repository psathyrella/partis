//
// File: SequenceApplicationTools.h
// Created by: Julien Dutheil
// Created on: Fri Oct 21 13:13
// from file old ApplicationTools.h created on Sun Dec 14 09:36:26 2003
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for sequences analysis.

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

#ifndef _SEQUENCEAPPLICATIONTOOLS_H_
#define _SEQUENCEAPPLICATIONTOOLS_H_

#include "../Alphabet/Alphabet.h"
#include "../GeneticCode/GeneticCode.h"
#include "../AlphabetIndex/AlphabetIndex1.h"
#include "../AlphabetIndex/AlphabetIndex2.h"
#include "../Container/SequenceContainer.h"
#include "../Container/VectorSiteContainer.h"

#include <map>
#include <string>

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
class SequenceApplicationTools
{
public:
  SequenceApplicationTools() {}
  virtual ~SequenceApplicationTools() {}

public:
  /**
   * @brief Build an Alphabet object according to options.
   *
   * Options used are:
   * - alphabet = [DNA|RNA|Protein], the alphabet type to use.
   *            = [DNA|RNA|Protein](length=n) a word-alphabet of
   *                 words with length n
   *            = [EchinodermMitochondrialCodonAlphabet
   *                   | InvertebrateMitochondrialCodonAlphabet
   *                   | InvertebrateMitochondrialCodonAlphabet
   *                   | StandardCodonAlphabet
   *                   | VertebrateMitochondrialCodonAlphabet]([alphn=NA|RNA])
   *                  a codon-alphabet
   *
   * @param params  The attribute map where options may be found.
   * @param suffix  A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose Print some info to the 'message' output stream.
   * @param allowGeneric Tell if generic alphabets can be used.
   * @return A new Alphabet object according to options specified.
   */
  static Alphabet* getAlphabet(
    std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    bool allowGeneric = false) throw (Exception);

  /**
   * @brief Build a GeneticCode object according to options.
   *
   * @param alphabet pointer to the NucleicAlphabet
   * @param description for the name of the GeneticCode:
   *    [EchinodermMitochondrialGeneticCode
   *    | InvertebrateMitochondrialGeneticCode
   *    | InvertebrateMitochondrialGeneticCode
   *    | StandardGeneticCode
   *    | VertebrateMitochondrialGeneticCode]
   * @return A new GeneticCode object
   * @throw Exception in case of bad description.
   */
  static GeneticCode* getGeneticCode(const NucleicAlphabet* alphabet, const std::string& description) throw (Exception);


  /**
   * @brief Build a AlphabetIndex1 object for a given alphabet.
   *
   * @param alphabet The alphabet to use. This is currently only used for assessing the type of distance allowed.
   * @param description Which distance to use. See the Bio++ Program Suite reference manual for a description of the syntax.
   * @param message To be displayed when parsing.
   * @param verbose Tell if some info should be displayed while parsing.
   * @return A new AlphabetIndex1 object.
   * @throw Exception in case of bad description.
   */
  static AlphabetIndex1* getAlphabetIndex1(const Alphabet* alphabet, const std::string& description, const std::string& message = "Alphabet distance:", bool verbose = true) throw (Exception);


  /**
   * @brief Build a AlphabetIndex2 object for a given alphabet.
   *
   * @param alphabet The alphabet to use. This is currently only used for assessing the type of distance allowed.
   * @param description Which distance to use. See the Bio++ Program Suite reference manual for a description of the syntax.
   * @param message To be displayed when parsing.
   * @return A new AlphabetIndex2 object.
   * @param verbose Tell if some info should be displayed while parsing.
   * @throw Exception in case of bad description.
   */
  static AlphabetIndex2* getAlphabetIndex2(const Alphabet* alphabet, const std::string& description, const std::string& message = "Alphabet distance:", bool verbose = true) throw (Exception);


  /**
   * @brief Build a SequenceContainer object according to options.
   *
   * The sequences do not have to be aligned.
   * The supported sequence formats are Fasta, DCSE, Clustal, Mase, Phylip and GenBank.
   *
   * See the Bio++ program suite manual for a full description of the syntax.
   *
   * @param alpha   The alphabet to use in the container.
   * @param params  The attribute map where options may be found.
   * @param suffix  A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose Print some info to the 'message' output stream.
   * @return A new VectorSequenceContainer object according to options specified.
   * @see getSiteContainer to read an alignment.
   */

  static SequenceContainer* getSequenceContainer(
    const Alphabet* alpha,
    std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true);

  /**
   * @brief Build a SiteContainer object according to options.
   *
   * Sequences in file must be aligned.
   * The supported sequence formats are Fasta, DCSE, Clustal, Mase and Phylip.
   *
   * See the Bio++ program suite manual for a full description of the syntax.
   *
   * @param alpha   The alphabet to use in the container.
   * @param params  The attribute map where options may be found.
   * @param suffix  A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose Print some info to the 'message' output stream.
   * @return A new VectorSiteContainer object according to options specified.
   */
  static VectorSiteContainer* getSiteContainer(
    const Alphabet* alpha,
    std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true);

  /**
   * @brief Retrieves sites suitable for the analysis.
   *
   * Options used are:
   * - sequence.sites_to_use = [all|complete|nogap].
   *
   * If the 'complete' option is used, only fully resolve site will be taken
   * into account.
   * If the 'nogap' option is used, only sites without gap will be taken into
   * account.
   * If 'gapAsUnknown' is set to true and the all option is selected, gaps will
   * be changed to 'unknown' character is sequences.
   *
   * - sequence.max_gap_allowed = [57%|30]
   * If a % sign fallow the number, it is taken to be a frequence (in percent).
   * This specify the maximum amount of gaps allowed for each site.
   * Sites not satisfying this amount will be removed.
   * A value of 100% will remove all gap-only sites, a value >100% will keep all sites.
   *
   * @param allSites The site container from which sites must be retrieved.
   * @param params   The attribute map where options may be found.
   * @param suffix   A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param gapAsUnknown Convert gaps to unknown characters.
   * @param verbose Print some info to the 'message' output stream.
   * @return A new VectorSiteContainer object containing sites of interest.
   */
  static VectorSiteContainer* getSitesToAnalyse(
    const SiteContainer& allSites,
    std::map<std::string, std::string>& params,
    std::string suffix = "",
    bool suffixIsOptional = true,
    bool gapAsUnknown = true,
    bool verbose = true);

  /**
   * @brief Write a sequence file according to options.
   *
   * The supported sequence formats are Fasta and Mase.
   *
   * See the Bio++ program suite manual for a full description of the syntax.
   *
   * @see writeSequenceFile(SiteContainer) for writing alignments, with more output formats.
   *
   * @param sequences The sequences to write.
   * @param params  The attribute map where options may be found.
   * @param suffix  A suffix to be applied to each attribute name.
   * @param verbose Print some info to the 'message' output stream.
   */
  static void writeSequenceFile(
    const SequenceContainer& sequences,
    std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool verbose = true);

  /**
   * @brief Write a sequence alignment file according to options.
   *
   * The supported sequence formats are Fasta, Mase and Phylip.
   *
   * See the Bio++ program suite manual for a full description of the syntax.
   *
   * @param sequences The aligned sequences to write.
   * @param params  The attribute map where options may be found.
   * @param suffix  A suffix to be applied to each attribute name.
   * @param verbose Print some info to the 'message' output stream.
   */
  static void writeAlignmentFile(
    const SiteContainer& sequences,
    std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool verbose = true);
};
} // end of namespace bpp.

#endif // _SEQUENCEAPPLICATIONTOOLS_H_

