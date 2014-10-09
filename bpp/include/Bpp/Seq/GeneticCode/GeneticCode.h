//
// File: GeneticCode.h
// Created by: Julien Dutheil
// Created on: Mon Oct 13 15:37:25 2003
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

#ifndef _GENETICCODE_H_
#define _GENETICCODE_H_

#include "../Transliterator.h"
#include "../Alphabet/CodonAlphabet.h"
#include "../Alphabet/ProteicAlphabet.h"
#include <Bpp/Exceptions.h>

namespace bpp
{

  /**
   * @brief Exception thrown when a stop codon is found.
   */
  class StopCodonException:
    public Exception
  {
  private:
    std::string codon_;
			
  public:
    // Class constructor
    StopCodonException(const std::string& text, const std::string& codon);
	
    // Class destructor
    virtual ~StopCodonException() throw () {}
		
  public:
    virtual const std::string& getCodon() const { return codon_; }
  };

  /**
   * @brief Partial implementation of the Transliterator interface for genetic code object.
   *
   * A genetic code object if a translator from a codon alphabet to a proteic alphabet.
   * Depending on the codon alphabet used, several genetic code can be implemented.
   *
   * @see CodonAlphabet, ProteicAlphabet
   */
  class GeneticCode:
    public AbstractTransliterator,
    public virtual Clonable
  {
  protected:
    CodonAlphabet codonAlphabet_;
    ProteicAlphabet proteicAlphabet_;
    std::map<int, int> tlnTable_;
	
  public:
    GeneticCode(const NucleicAlphabet* alphabet):
      AbstractTransliterator(),
      codonAlphabet_(alphabet),
      proteicAlphabet_(),
      tlnTable_()
    {}

    virtual ~GeneticCode() {}
	
    virtual GeneticCode* clone() const = 0;

  public:
    /**
     * @name Methods form the Transliterator interface.
     *
     * @{
     */
    const CodonAlphabet* getSourceAlphabet() const { return &codonAlphabet_; }
    const ProteicAlphabet* getTargetAlphabet() const { return &proteicAlphabet_; }
    virtual int translate(int state) const throw (BadIntException, Exception);		
    virtual std::string translate(const std::string& state) const throw (BadCharException, Exception);
    virtual Sequence* translate(const Sequence& sequence) const throw (Exception)
    {
      return AbstractTransliterator::translate(sequence);	
    }
    /** @} */
		
  public:
    /**
     * @name Specific methods.
     *
     * @{
     */
    
    /**
     * @return The number of stop codons.
     */
    virtual size_t getNumberOfStopCodons() const = 0;

    /**
     * @return A vector will all int codes for stop codons.
     */
    virtual std::vector<int> getStopCodonsAsInt() const = 0;

    /**
     * @return A vector will all char codes for stop codons.
     */
    virtual std::vector<std::string> getStopCodonsAsChar() const = 0;

    /**
     * @brief Tells is a particular codon is a stop codon.
     *
     * @param state The numeric code for the state to test.
     * @return True if the state corresponds to a stop codon.
     * @throw BadIntException if the state is not supported by the alphabet.
     */
    virtual bool isStop(int state) const throw (BadIntException) = 0;

    /**
     * @brief Tells is a particular codon is a stop codon.
     *
     * @param state The character code for the state to test.
     * @return True if the state corresponds to a stop codon.
     * @throw BadCharException if the state is not supported by the alphabet.
     */
    virtual bool isStop(const std::string& state) const throw (BadCharException) = 0;
    
    /**
     * @brief Tells is a particular codon is a start codon.
     *
     * @param state The numeric code for the state to test.
     * @return True if the state corresponds to a start codon.
     * @throw BadIntException if the state is not supported by the alphabet.
     */
    virtual bool isStart(int state) const throw (BadIntException) {
      //Test:
      codonAlphabet_.intToChar(state); //throw exception if invalid state!
      return (state == 14);
    }

    /**
     * @brief Tells is a particular codon is a start codon.
     *
     * @param state The character code for the state to test.
     * @return True if the state corresponds to a start codon.
     * @throw BadCharException if the state is not supported by the alphabet.
     */
    virtual bool isStart(const std::string& state) const throw (BadCharException) {
      return isStart(codonAlphabet_.charToInt(state));
    }
 
    /**
     * @brief Tells is a particular codon is an alternative start codon.
     *
     * @param state The numeric code for the state to test.
     * @return True if the state corresponds to an alternative start codon.
     * @throw BadIntException if the state is not supported by the alphabet.
     */
    virtual bool isAltStart(int state) const throw (BadIntException) = 0;

    /**
     * @brief Tells is a particular codon is an alternative start codon.
     *
     * @param state The character code for the state to test.
     * @return True if the state corresponds to an alternative start codon.
     * @throw BadCharException if the state is not supported by the alphabet.
     */
    virtual bool isAltStart(const std::string& state) const throw (BadCharException) = 0;
    
    /**
     * @brief Tell if two codons are synonymous, that is, if they encode the same amino-acid.
     *
     * @param i The numeric code for the first codon.
     * @param j The numeric code for the second codon.
     * @return True if the two codons are synonymous.
     * @throw BadIntException if at least one of the states is not supported by the alphabet.
     */
    bool areSynonymous(int i, int j) const throw (BadIntException)
    {
      return (translate(i) == translate(j));
    }

    /**
     * @brief Tell if two codons are synonymous, that is, if they encode the same amino-acid.
     *
     * @param i The character code for the first codon.
     * @param j The character code for the second codon.
     * @return True if the two codons are synonymous.
     * @throw BadCharException if at least one of the states is not supported by the alphabet.
     */
    bool areSynonymous(const std::string & i, const std::string & j) const throw (BadCharException)
    {
      return (translate(i) == translate(j));
    }

    std::vector<int> getSynonymous(int aminoacid) const throw (BadIntException);
    
    std::vector<std::string> getSynonymous(const std::string & aminoacid) const throw (BadCharException);

    /**
     * @return True if the specified codon is fourfold degenerated
     * (that is, if a mutation in the fourth position does not change the aminoacid).
     * @author Benoit Nabholz, Annabelle Haudry
     * @param codon The codon to test.
     */
    bool isFourFoldDegenerated(int codon) const;

    /**
     * @brief Get the subsequence corresponding to the coding part of a given sequence.
     *
     * If lookForInitCodon if set to 'true', the subsequence will start at the first AUG motif,
     * otherwise the subsequence will start at the begining of the sequence.
     * The subsequence ends at the first stop codon (excluded) found, or the end of the sequence.
     *
     * The sequence may have a nucleotide or codon alphabet.
     * The subsequence has the same alphabet, name and comments of the input sequence.
     * In case of nucleotide sequence and if the lookForInitCodon option is checked, the phase
     * will be determined from the sequence.
     *
     * @param sequence The sequence to parse.
     * @param lookForInitCodon Tell if the AUG codon must be found.
     * @param includeInitCodon (if lookForInitCodon is true) tell if the init codon must be included in the subsequence.
     * @return A nucleotide/codon subsequence.
     */
    Sequence* getCodingSequence(const Sequence& sequence, bool lookForInitCodon = false, bool includeInitCodon = false) const throw (Exception);
    /** @} */
  };

} //end of namespace bpp.

#endif	//_GENETICCODE_H_

