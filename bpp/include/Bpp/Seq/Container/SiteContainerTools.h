//
// File: SiteContainerTools.h
// Created by: Julien Dutheil
// Created on: Fri Dec 12 18:55:06 2003
//

#ifndef _SITECONTAINERTOOLS_H_
#define _SITECONTAINERTOOLS_H_

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

#include "SiteContainer.h"
#include "VectorSiteContainer.h"
#include "AlignedSequenceContainer.h"
#include "../AlphabetIndex/AlphabetIndex2.h"
#include "../DistanceMatrix.h"
#include "../GeneticCode/GeneticCode.h"
#include <Bpp/Numeric/Matrix/Matrix.h>

//From the STL:
#include <vector>
#include <map>

namespace bpp
{

  typedef std::vector<size_t> SiteSelection;

/**
 * @brief Some utililitary methods to deal with site containers.
 */
  class SiteContainerTools
  {
  public:
    SiteContainerTools() {}
    virtual ~SiteContainerTools() {}

  public:

    /**
     * @brief Retrieves sites without gaps from SiteContainer.
     *
     * This function build a new SiteContainer instance with only sites without gaps.
     * The container passed as input is not modified, all sites are copied.
     *
     * @param sites The container to analyse.
     * @return A pointer toward a new SiteContainer with only sites with no gaps.
     */
    static SiteContainer* getSitesWithoutGaps(const SiteContainer& sites);

    /**
     * @brief Retrieves complete sites from SiteContainer.
     *
     * This function build a new SiteContainer instance with only complete sites,
     * i.e. site with fully resolved states (no gap, no unknown caracters).
     * The container passed as input is not modified, all sites are copied.
     *
     * @param sites The container to analyse.
     * @return A pointer toward a new SiteContainer with only complete sites.
     */
    static SiteContainer* getCompleteSites(const SiteContainer& sites);

    /**
     * @brief Get a site set without gap-only sites.
     *
     * This function build a new SiteContainer instance without sites with only gaps.
     * The container passed as input is not modified, all sites are copied.
     *
     * @see removeGapOnlySites(SiteContainer& sites)
     * @param sites The container to analyse.
     * @return A pointer toward a new SiteContainer.
     */
    static SiteContainer* removeGapOnlySites(const SiteContainer& sites);

    /**
     * @brief Remove gap-only sites from a site set.
     *
     * @param sites The container where the sites have to be removed.
     */
    static void removeGapOnlySites(SiteContainer& sites);

    /**
     * @brief Get a site set without gap/unresolved-only sites.
     *
     * This function build a new SiteContainer instance without sites with only gaps or unresolved characters.
     * The container passed as input is not modified, all sites are copied.
     *
     * @param sites The container to analyse.
     * @return A pointer toward a new SiteContainer.
     */
    static SiteContainer* removeGapOrUnresolvedOnlySites(const SiteContainer& sites);
		
    /**
     * @brief Remove gap/unresolved-only sites from a site set.
     *
     * @param sites The container where the sites have to be removed.
     */
    static void removeGapOrUnresolvedOnlySites(SiteContainer& sites);

    /**
     * @brief Get a siteset with sites with less than a given amount of gaps.
     *
     * @param sites The container from which the sites have to be removed.
     * @param maxFreqGaps The maximum frequency of gaps in each site.
     * @return A pointer toward a new SiteContainer.
     */
    static SiteContainer* removeGapSites(const SiteContainer& sites, double maxFreqGaps);

    /**
     * @brief Remove sites with a given amount of gaps.
     *
     * @param sites The container from which the sites have to be removed.
     * @param maxFreqGaps The maximum frequency of gaps in each site.
     */
    static void removeGapSites(SiteContainer& sites, double maxFreqGaps);

    /**
     * @brief Get a site set without stop codons, if the alphabet is a CodonAlphabet, otherwise throws an Exception.
     *
     * This function build a new SiteContainer instance without sites that have at least a stop codon.
     * The container passed as input is not modified, all sites are copied.
     *
     * @param sites The container to analyse.
     * @param gCode the genetic code to use to determine stop codons.
     * @return A pointer toward a new SiteContainer.
     */
    static SiteContainer* removeStopCodonSites(const SiteContainer& sites, const GeneticCode& gCode) throw (AlphabetException);

    /**
     * @brief Create a new container with a specified set of sites.
     *
     * A new VectorSiteContainer is created with specified sites.
     * The destruction of the container is up to the user.
     * Sites are specified by their indice, beginning at 0.
     * No position verification is performed, based on the assumption that
     * the container passed as an argument is a correct one.
     * Redundant selection is not checked, so be careful with what you're doing!
     *
     * @param sequences The container from wich sequences are to be taken.
     * @param selection The positions of all sites to retrieve.
     * @return A new container with all selected sites.
     */
    static SiteContainer* getSelectedSites(const SiteContainer& sequences, const SiteSelection& selection);

    /**
     * @brief create the consensus sequence of the alignment.
     *
     * In case of ambiguity (for instance a AATT site), one state will be chosen arbitrarily.
     *
     * @param sc a site container
     * @param name the name of the sequence object that will be created. 
     * @param ignoreGap Tell if gap must be counted or not. If not (true option), only fully gapped sites will result in a gap in the consensus sequence. 
     * @param resolveUnknown Tell is unknnown characters must resolved. In a DNA sequence for instance, N will be counted as A=1/4, T=1/4, G=1/4 and C=1/4. Otherwise it will be counted as N=1.
     * If this option is set to true, a consensus sequence will never contain an unknown character.
     * @return A new Sequence object with the consensus sequence.
     */
    static Sequence* getConsensus(const SiteContainer& sc, const std::string& name = "consensus", bool ignoreGap = true, bool resolveUnknown = false);
    
    /**
     * @brief Change all gaps to unknown state in a container, according to its alphabet.
     *
     * For DNA alphabets, this change all '-' to 'N'.
     *
     * @param sites The container to be modified.
     */
    static void changeGapsToUnknownCharacters(SiteContainer& sites);

    /**
     * @brief Change all unresolved characters to gaps in a container, according to its alphabet.
     *
     * For DNA alphabets, this change all 'N', 'M', 'R', etc.  to '-'.
     *
     * @param sites The container to be modified.
     */
    static void changeUnresolvedCharactersToGaps(SiteContainer& sites);

    /**
     * @brief Resolve a container with "." notations.
     *
     * @code
     * ATGCCGTTGG
     * .C...A..C.
     * ..A....C..
     * @endcode
     * will results in
     * @code
     * ATGCCGTTGG
     * ACCCCATTCG
     * ATACCGTCGG
     * @endcode
     * for instance.
     * The first sequence is here called the "reference" sequence.
     * It need not be the first in the container.
     * The alphabet of the input alignment must be an instance of the DefaultAlphabet class, the only one which support dot characters.
     * A new alignment is created and returned, with the specified alphabet.
     *
     * If several sequences that may be considered as reference are found, the first one is used.
     * 
     * @param dottedAln The input alignment.
     * @param resolvedAlphabet The alphabet of the output alignment.
     * @return A pointer toward a dynamically created SiteContainer with the specified alphabet (can be a DefaultAlphabet).
     * @throw AlphabetException If the alphabet of the input alignment is not of class DefaultAlphabet, or if one character does not match with the output alphabet.
     * @throw Exception If no reference sequence was found, or if the input alignment contains no sequence.
     */
    static SiteContainer* resolveDottedAlignment(const SiteContainer& dottedAln, const Alphabet* resolvedAlphabet) throw (AlphabetException, Exception);

    /**
     * @name Sequences coordinates.
     *
     * @see SequenceWalker For an alternative approach.
     * @{
     */

    /**
     * @brief Get the index of each sequence position in an aligned sequence.
     *
     * If the sequence contains no gap, the translated and the original positions are the same.
     * Position numbers start at 1.
     *
     * @param seq The sequence to translate.
     * @return A map with original sequence positions as keys, and translated positions as values.
     */
    static std::map<size_t, size_t> getSequencePositions(const Sequence& seq);

    /**
     * @brief Get the index of each alignment position in an aligned sequence.
     *
     * If the sequence contains no gap, the translated and the original positions are the same.
     * Position numbers start at 1.
     *
     * @param seq The sequence to translate.
     * @return A map with original alignement positions as keys, and translated positions as values.
     */
    static std::map<size_t, size_t> getAlignmentPositions(const Sequence& seq);

    /**
     * @brief Fill a numeric matrix with the size of the alignment, containing the each sequence position.
     *
     * Positions start at 1, gaps have "position" 0.
     *
     * @param sites The input alignment.
     * @param positions A matrix object which is going to be resized and filled with the corresponding positions.
     * @author Julien Dutheil
     */
    static void getSequencePositions(const SiteContainer& sites, Matrix<size_t>& positions);
    /** @} */

    /**
     * @brief Translate alignement positions from an aligned sequence to the same sequence in a different alignment.
     *
     * Takes each position (starting at 1) in sequence 1, and look for the corresponding position in sequence 2.
     * The two sequences must be the same, excepted for the gaps.
     * If no sequence contains gaps, or if the gaps are at the same place in both sequences, the translated postion will be the same as the original positions.
     * 
     * @param seq1 The sequence to translate.
     * @param seq2 The reference sequence.
     * @return A map with original alignement positions as keys, and translated positions as values.
     * @throw AlphabetMismatchException If the sequences do not share the same alphabet.
     * @throw Exception If the sequence do not match.
     */
    static std::map<size_t, size_t> translateAlignment(const Sequence& seq1, const Sequence& seq2) throw (AlphabetMismatchException, Exception);

    /**
     * @brief Translate sequence positions from a sequence to another in the same alignment.
     *
     * Takes each position (starting at 1) in sequence 1, and look for the corresponding position in sequence 2 at the same site.
     * If no corresponding position is available (i.e. if there is a gap in sequence 2 at the corresponding position), 0 is returned.
     *
     * @param sequences The alignment to use.
     * @param i1 The index of the sequence to translate.
     * @param i2 The index of the reference sequence.
     * @return A map with original sequence positions as keys, and translated positions as values.
     */
    static std::map<size_t, size_t> translateSequence(const SiteContainer& sequences, size_t i1, size_t i2);

    /**
     * @brief Align two sequences using the Needleman-Wunsch dynamic algorithm.
     *
     * If the input sequences contain gaps, they will be ignored.
     *
     * @see BLOSUM50, DefaultNucleotideScore for score matrices.
     *
     * @param seq1 The first sequence.
     * @param seq2 The second sequence.
     * @param s The score matrix to use.
     * @param gap Gap penalty.
     * @return A new SiteContainer instance.
     * @throw AlphabetMismatchException If the sequences and the score matrix do not share the same alphabet.
     */
    static AlignedSequenceContainer* alignNW(const Sequence& seq1, const Sequence& seq2, const AlphabetIndex2& s, double gap) throw (AlphabetMismatchException);

    /**
     * @brief Align two sequences using the Needleman-Wunsch dynamic algorithm.
     *
     * If the input sequences contain gaps, they will be ignored.
     *
     * @see BLOSUM50, DefaultNucleotideScore for score matrices.
     *
     * @param seq1 The first sequence.
     * @param seq2 The second sequence.
     * @param s The score matrix to use.
     * @param opening Gap opening penalty.
     * @param extending Gap extending penalty.
     * @return A new SiteContainer instance.
     * @throw AlphabetMismatchException If the sequences and the score matrix do not share the same alphabet.
     */
    static AlignedSequenceContainer* alignNW(const Sequence& seq1, const Sequence& seq2, const AlphabetIndex2& s, double opening, double extending) throw (AlphabetMismatchException);

    /**
     * @brief Sample sites in an alignment.
     *
     * Original site positions will be kept. The resulting container will hence probably have duplicated
     * positions. You may wish to call the reindexSites() method on the returned container.
     *
     * Note: This method will be optimal with a container with vertical storage like VectorSiteContainer.
     *
     * @param sites An input alignment to sample.
     * @param nbSites The size of the resulting container.
     * @param index [out] If non-null the underlying vector will be appended with the original site indices.
     * @return A sampled alignment with nbSites sites taken from the input one.
     */
    static VectorSiteContainer* sampleSites(const SiteContainer& sites, size_t nbSites, std::vector<size_t>* index = 0);

    /**
     * @brief Bootstrap sites in an alignment.
     *
     * Original site positions will be kept. The resulting container will hence probably have duplicated
     * positions. You may wish to call the reindexSites() method on the returned container.
     *
     * Note: This method will be optimal with a container with vertical storage like VectorSiteContainer.
     *
     * @param sites An input alignment to sample.
     * @return A sampled alignment with the same number of sites than the input one.
     */
    static VectorSiteContainer* bootstrapSites(const SiteContainer& sites);

    /**
     * @brief Compute the similarity/distance score between two aligned sequences.
     *
     * The similarity measures are computed as the proportion of identical match.
     * The distance between the two sequences is defined as 1 - similarity.
     * This function can be used with any type of alphabet.
     *
     * @param seq1 The first sequence.
     * @param seq2 The second sequence.
     * @param dist Shall we return a distance instead of similarity?
     * @param gapOption How to deal with gaps:
     * - SIMILARITY_ALL: all positions are used.
     * - SIMILARITY_NODOUBLEGAP: ignore all positions with a gap in the two sequences.
     * - SIMILARITY_NOGAP: ignore all positions with a gap in at least one of the two sequences.
     * @param unresolvedAsGap Tell if unresolved characters must be considered as gaps when counting.
     * If set to yes, the gap option will also apply to unresolved characters.
     * @return The proportion of matches between the two sequences.
     * @throw SequenceNotAlignedException If the two sequences do not have the same length.
     * @throw AlphabetMismatchException If the two sequences do not share the same alphabet type.
     * @throw Exception If an invalid gapOption is passed.
     */
    static double computeSimilarity(const Sequence& seq1, const Sequence& seq2, bool dist = false, const std::string& gapOption = SIMILARITY_NODOUBLEGAP, bool unresolvedAsGap = true) throw (SequenceNotAlignedException, AlphabetMismatchException, Exception);

    /**
     * @brief Compute the similarity matrix of an alignment.
     *
     * The similarity measures are computed as the proportion of identical match.
     * The distance between the two sequences is defined as 1 - similarity.
     * This function can be used with any type of alphabet.
     * Several options concerning gaps and unresolved characters are proposed:
     * - SIMILARITY_ALL: all positions are used.
     * - SIMILARITY_NOFULLGAP: ignore positions with a gap in all the sequences in the alignment.
     * - SIMILARITY_NODOUBLEGAP: ignore all positions with a gap in the two sequences for each pair.
     * - SIMILARITY_NOGAP: ignore all positions with a gap in at least one of the two sequences for each pair.
     *
     *
     * @see computeSimilarityMatrix
     *
     * @param sites The input alignment.
     * @param dist Shall we return a distance instead of similarity?
     * @param gapOption How to deal with gaps.
     * @param unresolvedAsGap Tell if unresolved characters must be considered as gaps when counting.
     * If set to yes, the gap option will also apply to unresolved characters.
     * @return All pairwise similarity measures.
     */
    static DistanceMatrix* computeSimilarityMatrix(const SiteContainer& sites, bool dist = false, const std::string& gapOption = SIMILARITY_NOFULLGAP, bool unresolvedAsGap = true);

    static const std::string SIMILARITY_ALL;
    static const std::string SIMILARITY_NOFULLGAP;
    static const std::string SIMILARITY_NODOUBLEGAP;
    static const std::string SIMILARITY_NOGAP;

    /**
     * @brief Add the content of a site container to an exhisting one.
     *
     * The input containers are supposed to have unique sequence names.
     * If it is not the case, several things can happen:
     * - If the two containers have exactly the same names in the same order, then the content of the second one will be added as is to the first one.
     * - If the second container does not have exactly the same sequences names or in a different order, then a reordered selection of the second contianer is created first,
     *   and in that case, only the first sequence with a given name will be used and duplicated.
     * In any case, note that the second container should always contains all the sequence names from the first one,
     * otherwise an exception will be thrown.
     *
     * @author Julien Dutheil
     *
     * @param seqCont1 First container.
     * @param seqCont2 Second container. This container must contain sequences with the same names as in seqcont1.
     * Additional sequences will be ignored.
     * @param leavePositionAsIs Tell is site position should be unchanged. Otherwise (the default) is to add the size of container 1 to the positions in container 2.
     * @throw AlphabetMismatchException If the alphabet in the 2 containers do not match.
     * @throw Exception If sequence names do not match.
     */
    static void merge(SiteContainer& seqCont1, const SiteContainer& seqCont2, bool leavePositionAsIs = false) throw (AlphabetMismatchException, Exception);

    /**
     * @brief Compare an alignment to a reference alignment, and compute the column scores.
     *
     * Calculations are made according to formula for the "CS" score in Thompson et al 1999, Nucleic Acids Research (1999):27(13);2682–2690.
     *
     * @param positions1 Alignment index for the test alignment.
     * @param positions2 Alignment index for the reference alignment.
     * @param na         The score to use if the tested column is full of gap.
     * @return A vector of score, as 0 or 1.
     * @see getSequencePositions for creating the alignment indexes.
     * @warning The indexes for the two alignments must have the sequences in the exact same order!
     * @author Julien Dutheil
     */
    static std::vector<int> getColumnScores(const Matrix<size_t>& positions1, const Matrix<size_t>& positions2, int na = 0);
   
    /**
     * @brief Compare an alignment to a reference alignment, and compute the sum-of-pairs scores.
     *
     * Calculations are made according to formula for the "SPS" score in Thompson et al 1999, Nucleic Acids Research (1999):27(13);2682–2690.
     *
     * @param positions1 Alignment index for the test alignment.
     * @param positions2 Alignment index for the reference alignment.
     * @param na         The score to use if the tested column is not testable, that is not containing at least to residues.
     * @return A vector of score, between 0 and 1 (+ na value).
     * @see getSequencePositions for creating the alignment indexes.
     * @warning The indexes for the two alignments must have the sequences in the exact same order!
     * @author Julien Dutheil
     */
    static std::vector<double> getSumOfPairsScores(const Matrix<size_t>& positions1, const Matrix<size_t>& positions2, double na = 0);
  };

} //end of namespace bpp.

#endif	//_SITECONTAINERTOOLS_H_

