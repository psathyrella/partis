//
// File: SequenceTools.h
// Authors: Guillaume Deuchst
//          Julien Dutheil
//          Sylvain Gaillard
// Created on: Tue Aug 21 2003
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

#ifndef _SEQUENCETOOLS_H_
#define _SEQUENCETOOLS_H_

#include "Alphabet/Alphabet.h"
#include "Alphabet/DNA.h"
#include "Alphabet/RNA.h"
#include "Alphabet/RNY.h"
#include "GeneticCode/GeneticCode.h"
#include "Sequence.h"
#include "SymbolListTools.h"
#include "NucleicAcidsReplication.h"
#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Stat/StatTest.h>

// From the STL:
#include <string>
#include <map>
#include <vector>
#include <algorithm>

namespace bpp
{
/**
 * @brief Bowker's homogeneity test results class.
 */
class BowkerTest :
  public StatTest
{
private:
  double pvalue_;
  double stat_;

public:
  BowkerTest() : pvalue_(1.),
    stat_(0.) {}

  virtual ~BowkerTest() {}

  BowkerTest* clone() const { return new BowkerTest(*this); }

public:
  std::string getName() const { return "Bowker's test for homogeneity."; }
  double getStatistic() const { return stat_; }
  double getPValue() const { return pvalue_; }

  void setStatistic(double stat) { stat_ = stat; }
  void setPValue(double pvalue) { pvalue_ = pvalue; }
};

/**
 * @brief SequenceTools static class
 *
 * Implement methods to manipulate sequences
 */
class SequenceTools :
  public SymbolListTools
{
private:
  static DNA _DNA;
  static RNA _RNA;
  static RNY _RNY;
  static NucleicAcidsReplication _DNARep;
  static NucleicAcidsReplication _RNARep;
  static NucleicAcidsReplication _transc;

public:
  SequenceTools() {}
  virtual ~SequenceTools() {}

public:
  /**
   * @brief Get a sub-sequence.
   *
   * @param sequence The sequence to trunc.
   * @param begin The first position of the subsequence.
   * @param end   The last position of the subsequence.
   * @return A new sequence object with the given subsequence.
   * @throw IndexOutOfBoundsException, Exception In case of bad indices.
   */
  static Sequence* subseq(const Sequence& sequence, size_t begin, size_t end) throw (IndexOutOfBoundsException, Exception);

  /**
   * @brief Concatenate two sequences.
   *
   * Sequences must have the same name and alphabets.
   * Only first sequence's commentaries are kept.
   *
   * @param seq1 The first sequence.
   * @param seq2 The second sequence.
   * @return A new sequence object with the concatenation of the two sequences.
   * @throw AlphabetMismatchException If the two alphabets do not match.
   * @throw Exception If the sequence names do not match.
   */
  static Sequence* concatenate(const Sequence& seq1, const Sequence& seq2)
  throw (AlphabetMismatchException, Exception);

  /**
   * @brief Complement the nucleotide sequence itself
   *
   * @param seq The sequence to be complemented.
   * @return A ref toward the complemented sequence.
   * @throw AlphabetException if the sequence is not a nucleotide sequence.
   * @author Sylvain Gaillard
   */
  static Sequence& complement(Sequence& seq) throw (AlphabetException);

  /**
   * @brief Get the complementary sequence of a nucleotide sequence.
   *
   * @see DNAReplication
   * @return A new sequence object with the complementary sequence.
   * @param sequence The sequence to complement.
   * @throw AlphabetException If the sequence is not a nucleotide sequence.
   */
  static Sequence* getComplement(const Sequence& sequence) throw (AlphabetException);

  /**
   * @brief Get the transcription sequence of a DNA sequence.
   *
   * Translate DNA sequence into RNA sequence.
   *
   * @see DNAReplication
   * @return sequence A new sequence object with the transcription sequence.
   * @param sequence The sequence to transcript.
   * @throw AlphabetException If the sequence is not a DNA sequence.
   */
  static Sequence* transcript(const Sequence& sequence) throw (AlphabetException);

  /**
   * @brief Get the reverse-transcription sequence of a RNA sequence.
   *
   * Translate RNA sequence into DNA sequence.
   *
   * @see DNAReplication
   * @return sequence A new sequence object with the reverse-transcription sequence.
   * @param sequence The sequence to reverse-transcript.
   * @throw AlphabetException If the sequence is not a RNA sequence.
   */
  static Sequence* reverseTranscript(const Sequence& sequence) throw (AlphabetException);

  /**
   * @brief Inverse a sequence from 5'->3' to 3'->5' and vice-versa.
   *
   * ABCDEF becomes FEDCBA, and the sense attribute is changed (may be
   * inhibited).
   *
   * @param seq The sequence to inverse.
   * @return A ref toward the sequence.
   * @author Sylvain Gaillard
   */
  static Sequence& invert(Sequence& seq);

  /**
   * @brief Inverse a sequence from 5'->3' to 3'->5' and vice-versa.
   *
   * ABCDEF becomes FEDCBA, and the sense attribute is changed (may be
   * inhibited).
   *
   * @param sequence The sequence to inverse.
   * @return A new sequence object containing the inverted sequence.
   * @author Sylvain Gaillard
   */
  static Sequence* getInvert(const Sequence& sequence);

  /**
   * @brief Inverse and complement a sequence.
   *
   * This methode is more accurate than calling invert and complement
   * separatly.
   *
   * @param seq The sequence to inverse and complement.
   * @return A ref toward the sequence.
   * @author Sylvain Gaillard
   */
  static Sequence& invertComplement(Sequence& seq);

  /**
   * @return The identity percent of 2 sequence.
   * One match is counted if the two sequences have identical states.
   * @param seq1 The first sequence.
   * @param seq2 The second sequence.
   * @param ignoreGaps If true, only positions without gaps will be used for the counting.
   * @throw AlphabetMismatchException If the two sequences do not have the same alphabet.
   * @throw SequenceNotAlignedException If the two sequences do not have the same length.
   */
  static double getPercentIdentity(const Sequence& seq1, const Sequence& seq2, bool ignoreGaps = false) throw (AlphabetMismatchException, SequenceNotAlignedException);

  /**
   * @return The number of sites in the sequences, <i>i.e.</i> all positions without gaps.
   *
   * @param seq The sequence to analyse.
   */
  static size_t getNumberOfSites(const Sequence& seq);

  /**
   * @return The number of complete sites in the sequences, <i>i.e.</i> all positions without gaps and unresolved states (generic characters).
   *
   * @param seq The sequence to analyse.
   */
  static size_t getNumberOfCompleteSites(const Sequence& seq);

  /**
   * @return The number of unresolved sites in the sequence.
   *
   * @param seq The sequence to analyse.
   *
   * @author Sylvain Gaillard
   */
  static size_t getNumberOfUnresolvedSites(const Sequence& seq);

  /**
   * @brief Remove gaps from a sequence.
   *
   * The deleteElement method of the Sequence object will be used where appropriate.
   * @param seq The sequence to analyse.
   */
  static void removeGaps(Sequence& seq);

  /**
   * @brief Get a copy of the sequence without gaps.
   *
   * A whole new sequence will be created by adding all non-gap positions.
   * The original sequence will be cloned to serve as a template.
   *
   * @param seq The sequence to analyse.
   * @return A new sequence object without gaps.
   */
  static Sequence* getSequenceWithoutGaps(const Sequence& seq);

  /**
   * @brief Remove stops from a codon sequence.
   *
   * The deleteElement method of the Sequence object will be used where appropriate.
   * @param seq The sequence to analyse.
   * @param gCode The genetic code according to which stop codons are specified.
   * @throw Exception if the input sequence does not have a codon alphabet.
   */
  static void removeStops(Sequence& seq, const GeneticCode& gCode) throw (Exception);

  /**
   * @brief Get a copy of the codon sequence without stops.
   *
   * A whole new sequence will be created by adding all non-stop positions.
   * The original sequence will be cloned to serve as a template.
   *
   * @param seq The sequence to analyse.
   * @param gCode The genetic code according to which stop codons are specified.
   * @return A new sequence object without stops.
   * @throw Exception if the input sequence does not have a codon alphabet.
   */
  static Sequence* getSequenceWithoutStops(const Sequence& seq, const GeneticCode& gCode) throw (Exception);

  /**
   * @brief Replace stop codons by gaps.
   *
   * The setElement method of the Sequence object will be used where appropriate.
   * @param seq The sequence to analyse.
   * @param gCode The genetic code according to which stop codons are specified.
   * @throw Exception if the input sequence does not have a codon alphabet.
   */
  static void replaceStopsWithGaps(Sequence& seq, const GeneticCode& gCode) throw (Exception);

  /**
   * @brief Bowker's test for homogeneity.
   *
   * Computes the contingency table of occurrence of all pairs of states and test its symmetry using Bowker's (1948) test.
   *
   * Reference:<br>
   * @code
   * Ababneh F. Bioinformatics 2006 22(10) 1225-1231
   * @endcode
   *
   * @param seq1 The first sequence.
   * @param seq2 The second sequence.
   * @return A BowkerTest object with the computed statistic and p-value (computed from a chi square distribution).
   * @throw SequenceNotAlignedException If the two sequences do not have the same length.
   */
  static BowkerTest* bowkerTest(const Sequence& seq1, const Sequence& seq2) throw (SequenceNotAlignedException);

  /**
   * @brief Get all putatives haplotypes from an heterozygous sequence.
   *
   * @param seq The sequence to resolve
   * @param hap The vector to fill with the new sequences
   * @param level The maximum number of states that a generic char must code
   * (if this number is higher than level, the state will not be resolved).
   * For instance if level = 3 and Alphabet is DNA, all generic char will be
   * resolved but N.
   *
   * @author Sylvain Gaillard
   */
  static void getPutativeHaplotypes(const Sequence& seq, std::vector<Sequence*>& hap, unsigned int level = 2);

  /**
   * @brief Combine two sequences.
   *
   * @author Sylvain Gaillard
   */

  static Sequence* combineSequences(const Sequence& s1, const Sequence& s2) throw (AlphabetMismatchException);

  /**
   * @brief Subtract haplotype from an heterozygous sequence.
   *
   * Subtract an haplotype (i.e. a fully resolved sequence) from an heterozygous
   * sequence to get the other haplotype. The new haplotype could be an unresolved
   * sequence if unresolved characters in the sequence code for more than 2 states.
   *
   * For example:<br>
   * @code
   * >heterozygous sequence
   * ATTCGGGKWTATRYRM
   * >haplotype
   * ATTCGGGTATATGCAA
   * >subtracted haplotype
   * ATTCGGGGTTATATGC
   * @endcode
   *
   * @param s The heterozygous sequence.
   * @param h The haplotype to subtract.
   * @param name The name of the new computed haplotype.
   * @param level The number of states from which the site is set to fully unresolved.
   * @throw SequenceNotAlignedException if s and h don't have the same size.
   *
   * @author Sylvain Gaillard
   */
  static Sequence* subtractHaplotype(const Sequence& s, const Sequence& h, std::string name = "", unsigned int level = 1) throw (SequenceNotAlignedException);

  /**
   * @brief Get the RNY decomposition of a DNA sequence; with a given
   * phase between 1 and 3, it gives the decomposition in this phase;
   * in phase 1, the first triplet is centered on the first character.
   * Without a phase the function gives the alternative succession in
   * phases 1, 2 and 3.
   *
   * @return sequence A new sequence object with the transcription sequence.
   * @param sequence The sequence to transcript.
   * @param ph The phase to use (1,2 or 3).
   * @throw AlphabetException If the sequence is not a DNA sequence.
   *
   * @author Laurent Guéguen
   */
  static Sequence* RNYslice(const Sequence& sequence, int ph) throw (AlphabetException);
  static Sequence* RNYslice(const Sequence& sequence) throw (AlphabetException);

  /**
   * @brief Extract CDS part from a codon sequence. Optionally check for intiator and stop codons, or both.
   *
   * @param sequence The sequence to be reduced to CDS part.
   * @param gCode The genetic code according to which start and stop codons are specified.
   * @param checkInit If true, then everything before the initiator codon will be removed, together with the initiator codon if includeInit is false.
   * @param checkStop If true, then everything after the first stop codon will be removed, together with the stop codon if includeStop is false.
   * @param includeInit Tell if initiator codon should be kept or removed. No effect if checkInit is false.
   * @param includeStop Tell if stop codon should be kept or removed. No effect if checkStop is false.
   */
  static void getCDS(Sequence& sequence, const GeneticCode& gCode, bool checkInit, bool checkStop, bool includeInit = true, bool includeStop = true);

  /**
   * @brief Find the position of a motif in a sequence
   *
   * @param seq The reference sequence
   * @param motif The motif to find
   * @param strict If true (default) find exactly the motif
   *               If false find compatible match
   * @return The position of the first occurence of the motif or the seq
   * length.
   */
  static size_t findFirstOf(const Sequence& seq, const Sequence& motif, bool strict = true);

  /**
   * @brief Get a random sequence of given size and alphabet, with all state with equal probability.
   *
   * @param alphabet The alphabet to use.
   * @param length The length of the sequence to generate.
   * @return A pointer toward a new Sequence object.
   */
  static Sequence* getRandomSequence(const Alphabet* alphabet, size_t length);
};
} // end of namespace bpp.

#endif // _SEQUENCETOOLS_H_

