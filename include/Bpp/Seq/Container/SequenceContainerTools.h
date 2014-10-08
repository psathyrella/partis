//
// File: SequenceContainerTools.h
// Created by: Julien Dutheil
//             Sylvain Gaillard
// Created on: Sat Oct  4 09:18:34 2003
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

#ifndef _SEQUENCECONTAINERTOOLS_H_
#define _SEQUENCECONTAINERTOOLS_H_

// From the STL:
#include <string>
#include <vector>
#include <map>
#include <memory>

#include "SequenceContainer.h"
#include "OrderedSequenceContainer.h"

namespace bpp
{

typedef std::vector<size_t> SequenceSelection;

/**
 * @brief Utilitary methods dealing with sequence containers.
 */
class SequenceContainerTools
{

  public:
    SequenceContainerTools() {}
     virtual ~SequenceContainerTools() {}
  
  public:
    /**
     * @brief Create a container with @f$n@f$ void sequences.
     *
     * A new VectorSequenceContainer with the specified alphabet is created.
     * The destruction of this new container is up to the user.
     * Sequences have name "0", "1",... "n-1" and no content and comments.
     *
     * @param alphabet The alphabet to use in the container.
     * @param size     The number of sequences in the container.
     * @return A pointer toward a newly created container.
     */
    static SequenceContainer* createContainerOfSpecifiedSize(const Alphabet* alphabet, size_t size);

    /**
     * @brief Create a container with specified names.
     *
     * A new VectorSequenceContainer with the specified alphabet is created.
     * The destruction of this new container is up to the user.
     * Sequences have the specified names and no content and comments.
     *
     * @param alphabet The alphabet to use in the container.
     * @param seqNames The names of the sequences.
     * @return A pointer toward a newly created container.
     * @throw Exception If two sequence names are not unique.
     */
    static SequenceContainer* createContainerWithSequenceNames(
      const Alphabet* alphabet,
      const std::vector<std::string>& seqNames)
      throw (Exception);
 
    /**
     * @brief Generic function which creates a new container from another one,
     * by specifying the class of sequence to be stored.
     *
     * Compared to several copy constructors, this function allows to change the class of
     * the inner sequence class used for storing sequences.
     * The function used the addSequence method, so that it can also be used to
     * concatenate containers.
     *
     * @param input The container to copy.
     * @param output The container where new sequences will be appended.
     */
    template<class ContFrom, class ContTo, class Seq>
    static void convertContainer(const ContFrom& input, ContTo& output) {
      for (size_t i = 0; i < input.getNumberOfSequences(); ++i) {
        std::auto_ptr<Seq> seq(new Seq(input.getSequence(i)));
        output.addSequence(*seq);
      }
    }

    /**
     * @brief Add a specified set of sequences from a container to another.
     *
     * Sequences are specified by their position, beginning at 0.
     * Name verification will be performed, only if the output container is not empty,
     * based on the assumption that the container passed as argument is a correct one.
     * Redundant selection is not checked, so be careful with what you're doing!
     *
     * @author Julien Dutheil
     *
     * @param sequences The container from wich sequences are to be taken.
     * @param selection The positions of all sequences to retrieve.
     * @param outputCont A container where the selection should be added.
     * @throw Exception In case of bad sequence name, alphabet mismatch, etc.
     */
    static void getSelectedSequences(const OrderedSequenceContainer& sequences, const SequenceSelection& selection, SequenceContainer& outputCont) throw (Exception);

    /**
     * @brief Add a specified set of sequences from a container to another.
     *
     * Sequences are specified by their names.
     * Name verification will be performed, only if the output container is not empty,
     * based on the assumption that the container passed as argument is a correct one.
     * Redundant selection is not checked, so be careful with what you're doing!
     *
     * @author Julien Dutheil
     *
     * @param sequences The container from wich sequences are to be taken.
     * @param selection The names of all sequences to retrieve.
     * @param outputCont A container where the selection should be added.
     * @param strict If yes, trying to select a sequence that is not present
     * will raise an exception. If no, only available sequence will be added.
     * @throw Exception In case of bad sequence name, alphabet mismatch, etc.
     */
    static void getSelectedSequences(const SequenceContainer& sequences, const std::vector<std::string>& selection, SequenceContainer& outputCont, bool strict = true) throw (Exception);

    /**
     * @brief Remove all sequences that are not in a given selection from a given container.
     *
     * A new VectorSequenceContainer is created with specified sequences.
     * The destruction of the container is up to the user.
     * Sequences are specified by their position, beginning at 0.
     * Redundant selection is not checked, so be careful with what you're doing!
     *
     * @param sequences The container from wich sequences are to be taken.
     * @param selection The positions of all sequences to retrieve.
     * @return A new container with all selected sequences.
     */
    static void keepOnlySelectedSequences(OrderedSequenceContainer& sequences, const SequenceSelection& selection);
    
    /**
     * @brief Check if all sequences in a SequenceContainer have the same length.
     *
     * @param sequences The container to check.
     * @return True is all sequence have the same length.
     */
    static bool sequencesHaveTheSameLength(const SequenceContainer& sequences);
  
    /**
     * @brief Compute base counts
     *
     * Example of usage: getting the GC count from a sequence container.
     * <code>
     *   map<int, int> counts;
     *  SequenceContainerTools::getCounts(myContainer, counts); //My container is previously defined.
     *   int GCcontent = counts[1] + counts[2] ;
     * </code>
     *
     * States are stored as their int code.
     */

  static void getCounts(const SequenceContainer& sequences, std::map<int, int>&);

  /**
   * @brief Compute base frequencies.
   *
   * Example of usage: getting the GC content from a sequence container.
   * <code>
   *  map<int, double> freqs;
   *  SequenceContainerTools::getFrequencies(myContainer, freqs); //My container is previously defined.
   *   double GCcontent = (freqs[1] + freqs[2]) / (freqs[0] + freqs[1] + freqs[2] + freqs[3]);
   * </code>
   *
   * States are stored as their int code.
   */
  
  static void  getFrequencies(const SequenceContainer& sequences, std::map<int, double>& f, double pseudoCount = 0);
  
    /**
     * @brief Append all the sequences of a SequenceContainer to the end of another.
     *
     * @param seqCont1 The SequenceContainer in which the sequences will be added.
     * @param seqCont2 The SequenceContainer from which the sequences are taken.
     * @param checkNames Tell if the sequence names should be check for unicity.
     */
    static void append(SequenceContainer& seqCont1, const SequenceContainer& seqCont2, bool checkNames = true)
    throw (Exception)
    {
      std::vector<std::string> seqNames = seqCont2.getSequencesNames();
      for (size_t i = 0; i < seqNames.size(); i++)
        seqCont1.addSequence(seqCont2.getSequence(seqNames[i]), checkNames);
    }
    /**
     * @brief Append all the sequences of a SequenceContainer to the end of another, OrderedSequenceContainer implementation.
     *
     * @param seqCont1 The SequenceContainer in which the sequences will be added.
     * @param seqCont2 The SequenceContainer from which the sequences are taken.
     * @param checkNames Tell if the sequence names should be check for unicity.
     */
    static void append(SequenceContainer& seqCont1, const OrderedSequenceContainer& seqCont2, bool checkNames=true)
    throw (Exception)
    {
      for (size_t i = 0; i < seqCont2.getNumberOfSequences(); i++)
        seqCont1.addSequence(seqCont2.getSequence(i), checkNames);
    }
    
    /**
     * @brief Concatenate the sequences from two containers.
     *
     * This method will not check the original sequence names for unicity. If sequences do not have a unique name,
     * then the resulting merged container will contain the first sequence with the given duplicated name.
     *
     * @author Julien Dutheil
     *
     * @param seqCont1 First container.
     * @param seqCont2 Second container. This container must contain sequences with the same names as in seqcont1.
     * Additional sequences will be ignored.
     * @param outputCont Output sequence container to which concatenated sequences will be added.
     * @throw AlphabetMismatchException If the alphabet in the 3 containers do not match.
     */
    static void merge(const SequenceContainer& seqCont1, const SequenceContainer& seqCont2, SequenceContainer& outputCont)
    throw (Exception)
    {
      if (seqCont1.getAlphabet()->getAlphabetType() != seqCont2.getAlphabet()->getAlphabetType())
        throw AlphabetMismatchException("SequenceContainerTools::merge.", seqCont1.getAlphabet(), seqCont2.getAlphabet());

      std::vector<std::string> seqNames = seqCont1.getSequencesNames();
      for (size_t i = 0; i < seqNames.size(); i++)
      {
        BasicSequence tmp = seqCont1.getSequence(seqNames[i]);
        tmp.append(seqCont2.getContent(seqNames[i]));
        outputCont.addSequence(tmp, false);
      }
    }

    /**
     * @brief Convert a SequenceContainer with a new alphabet.
     *
     * This method assume that the original container has proper sequence names.
     * Names will be checked only if the output container is not empty.
     * @param seqCont The container to convert.
     * @param outputCont A container (most likely empty) with an alphabet into which the container will be converted.
     */
    static void convertAlphabet(const SequenceContainer& seqCont, SequenceContainer& outputCont)
    throw (Exception)
    {  
      std::vector<std::string> seqNames = seqCont.getSequencesNames();
      bool checkNames = outputCont.getNumberOfSequences() > 0;
      for (size_t i = 0; i < seqNames.size(); i++)
      {
        BasicSequence seq(seqNames[i], seqCont.toString(seqNames[i]), outputCont.getAlphabet());
        outputCont.addSequence(seq, checkNames);
      }
    }

    /**
     * @brief Extract a certain position (1, 2 or 3) from a container of codon sequences and returns the resulting nucleotide container.
     *
     * @param sequences The input sequence container, with codon alphabet.
     * @param pos       The codon position to retrieve.
     * @return          A SequenceContainer with a nucleotide alphabet.
     * @throw AlphabetException If input sequences are not registered with a codon alphabet.
     */
    static SequenceContainer* getCodonPosition(const SequenceContainer& sequences, size_t pos) throw (AlphabetException);

};

} //end of namespace bpp.

#endif  //_SEQUENCECONTAINERTOOLS_H_

