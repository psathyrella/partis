//
// File:       SequenceWithQualityTools.h
// Authors:    Vincent Cahais
//             Sylvain Gaillard
// Created on: 16 Apr 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (Apr 16, 2010)

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

#ifndef SEQUENCEWITHQUALITYTOOLS_H_
#define SEQUENCEWITHQUALITYTOOLS_H_

#include "SequenceTools.h"
#include "SequenceWithQuality.h"

namespace bpp {
  /**
   * @brief The SequenceWithQualityTools static class
   *
   * Implement methods to manipulate SequencesWithQuality
   *
   * @todo
   * - does complement, transcript and reverseTranscript have a meaning on
   *   SequenceWithQuality as the quality is not transferable
   *
   * @author Cahais Vincent
   */

  class SequenceWithQualityTools {
    private:
      static DNA _DNA;
      static RNA _RNA;
      static NucleicAcidsReplication _DNARep;
      static NucleicAcidsReplication _RNARep;
      static NucleicAcidsReplication _transc;

    public:

      /**
       * @brief Get a sub-sequence.
       *
       * @param sequence The sequence to trunc.
       * @param begin The first position of the subsequence.
       * @param end   The last position of the subsequence.
       * @return A new SequenceWithQuality object with the given subsequence.
       * @throw IndexOutOfBoundsException, Exception In case of bad indices.
       */
      static SequenceWithQuality* subseq(
          const SequenceWithQuality& sequence,
          unsigned int begin,
          unsigned int end
          ) throw (IndexOutOfBoundsException, Exception) ;

      /**
       * @brief Concatenate two sequences.
       *
       * Sequences must have the same name and alphabets.
       * Only first sequence's commentaries are kept.
       *
       * @param seqwq1 The first SequenceWithQuality.
       * @param seqwq2 The second SequenceWithQuality.
       * @return A new SequenceWithQuality object with the concatenation of the
       * two sequences.
       * @throw AlphabetMismatchException If the two alphabets do not match.
       * @throw Exception If the sequence names do not match.
       */
      static SequenceWithQuality* concatenate(
          const SequenceWithQuality& seqwq1,
          const SequenceWithQuality& seqwq2
          ) throw (AlphabetMismatchException, Exception) ;

      /**
       * @brief Get the complementary sequence of a nucleotide sequence.
       *
       * @see DNAReplication
       * @return sequence A new SequenceWithQuality object with the
       * complementary sequence.
       * @param sequence The sequence to complement.
       * @throw AlphabetException If the sequence is not a nucleotide sequence.
       */
      static SequenceWithQuality* complement(
          const SequenceWithQuality& sequence
          ) throw (AlphabetException);

      /**
       * @brief Get the transcription sequence of a DNA sequence.
       *
       * @see DNAReplication
       * @return sequence A new SequenceWithQuality object with the
       * transcription sequence.
       * @param sequence The sequence to transcript.
       * @throw AlphabetException If the sequence is not a DNA sequence.
       */
      static SequenceWithQuality* transcript(
          const SequenceWithQuality& sequence
          ) throw (AlphabetException);

      /**
       * @brief Get the reverse-transcription sequence of a RNA sequence.
       *
       * @see DNAReplication
       * @return sequence A new SequenceWithQuality object with the reverse-
       * transcription sequence.
       * @param sequence The SequenceWithQuality to reverse-transcript.
       * @throw AlphabetException If the sequence is not a RNA sequence.
       */

      static SequenceWithQuality* reverseTranscript(
          const SequenceWithQuality& sequence
          ) throw (AlphabetException);
      /**
       * @brief Inverse a sequence from 5'->3' to 3'->5' and vice-versa.
       *
       * ABCDEF becomes FEDCBA, and the sense attribute is changed (may be
       * inhibited).
       *
       * @return A new SequenceWithQuality object containing the inverted
       * sequence.
       * @param sequence The SequenceWithQuality to inverse.
       */
      static SequenceWithQuality* invert(
          const SequenceWithQuality& sequence
          );

      /**
       * @brief Remove gaps from a SequenceWithQuality.
       *
       * @param seq The sequence to analyse.
       * @return A new SequenceWithQuality object without gaps.
       */
      static SequenceWithQuality* removeGaps(const SequenceWithQuality& seq);

      /**
       * @brief Trim the left part of the sequence according to quality
       *
       * @param seq The sequence to analyse.
       * @return The modified sequence.
       */
      static SequenceWithQuality& trimLeft(SequenceWithQuality& seq);

  };
}

#endif /* SEQUENCEWITHQUALITYTOOLS_H_ */
