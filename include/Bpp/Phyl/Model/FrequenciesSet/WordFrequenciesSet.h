//
// File: WordFrequenciesSet.h
// Created by: Laurent Gueguen
// Created on: lundi 2 avril 2012, à 13h 59
//

/*
   Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _WORDFREQUENCIESSET_H_
#define _WORDFREQUENCIESSET_H_

#include <Bpp/Seq/Alphabet/WordAlphabet.h>
#include "FrequenciesSet.h"

namespace bpp
{

  /*********************************************************************/
/****   Frequencies Set in Words *****/
/*********************************************************************/


/**
 * @brief Frequencies in words computed from the  frequencies on
 * letters. The parameters are the parameters of the Frequencies on
 * letters.
 * The WordFrequenciesSet owns the FrequenciesSet* it is built on.
 * Interface class.
 * @author Laurent Guéguen
 */

class WordFrequenciesSet :
  public virtual FrequenciesSet
{
protected:
  
  virtual size_t getSizeFromVector(const std::vector<FrequenciesSet*>& freqVector) = 0;
  
public:
#ifndef NO_VIRTUAL_COV
  WordFrequenciesSet* clone() const = 0;

  const WordAlphabet* getAlphabet() const = 0;
#endif

  /**
   *@ brief Returns the n-th FrequenciesSet&
   **/

  virtual const FrequenciesSet& getFrequenciesSetForLetter(size_t i) const = 0;

  /**
   *@ brief Returns the length of the words
   **/

  virtual size_t getLength() const = 0;
};


class AbstractWordFrequenciesSet :
  public virtual WordFrequenciesSet,
  public AbstractFrequenciesSet
{
protected:
  size_t getSizeFromVector(const std::vector<FrequenciesSet*>& freqVector);
  
public:
  AbstractWordFrequenciesSet(size_t size, const Alphabet* palph, const std::string& prefix = "", const std::string& name="");

#ifndef NO_VIRTUAL_COV
  AbstractWordFrequenciesSet*
#else
  Clonable*
#endif
  clone() const = 0;

  AbstractWordFrequenciesSet(const AbstractWordFrequenciesSet& af) :
    AbstractFrequenciesSet(af) {}

  AbstractWordFrequenciesSet & operator=(const AbstractWordFrequenciesSet& af)
  {
    AbstractFrequenciesSet::operator=(af);
    return *this;
  }

#ifndef NO_VIRTUAL_COV
  const WordAlphabet* getAlphabet() const
  {
    return dynamic_cast<const WordAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif

  virtual ~AbstractWordFrequenciesSet();
  
  /**
   *@ brief Return the length of the words
   **/
  
  size_t getLength() const;
};


/**
 * @brief the Frequencies in words are the product of Independent Frequencies in letters
 * @author Laurent Guéguen
 */

class WordFromIndependentFrequenciesSet :
    public AbstractWordFrequenciesSet
{
protected:
  std::vector<FrequenciesSet*> vFreq_;
  std::vector<std::string> vNestedPrefix_;

public:
  /**
   * @brief Constructor from a WordAlphabet* and a vector of different FrequenciesSet*.
   * Throws an Exception if their lengths do not match.
   */
  WordFromIndependentFrequenciesSet(const WordAlphabet* pWA, const std::vector<FrequenciesSet*>& freqVector, const std::string& prefix = "", const std::string& name="WordFromIndependent");

  WordFromIndependentFrequenciesSet(const WordFromIndependentFrequenciesSet& iwfs);

  ~WordFromIndependentFrequenciesSet();

  WordFromIndependentFrequenciesSet& operator=(const WordFromIndependentFrequenciesSet& iwfs);

  WordFromIndependentFrequenciesSet* clone() const { return new WordFromIndependentFrequenciesSet(*this); }

public:
  void fireParameterChanged(const ParameterList& pl);

  virtual void updateFrequencies();

  /**
   *@ brief Independent letter frequencies from given word frequencies.
   * The frequencies of a letter at a position is the sum of the
   *    frequencies of the words that have this letter at this
   *    position.
   */
  virtual void setFrequencies(const std::vector<double>& frequencies);

  /**
   *@ brief Return the n-th FrequenciesSet&
   **/
  const FrequenciesSet& getFrequenciesSetForLetter(size_t i) const { return *vFreq_[i]; }

  /**
   *@ brief Return the length of the words
   **/

  virtual size_t getLength() const;

  void setNamespace(const std::string& prefix);

  std::string getDescription() const;
};

class WordFromUniqueFrequenciesSet :
  public AbstractWordFrequenciesSet
{
protected:
  FrequenciesSet* pFreq_;
  std::string NestedPrefix_;
  size_t length_;

public:
  /**
   * @brief Constructor from a WordAlphabet* and a FrequenciesSet*
   *  repeated as many times as the length of the words.
   */
  WordFromUniqueFrequenciesSet(const WordAlphabet* pWA, FrequenciesSet* pabsfreq, const std::string& prefix = "", const std::string& name = "WordFromUnique");

  WordFromUniqueFrequenciesSet(const WordFromUniqueFrequenciesSet& iwfs);

  WordFromUniqueFrequenciesSet& operator=(const WordFromUniqueFrequenciesSet& iwfs);

  ~WordFromUniqueFrequenciesSet();

  WordFromUniqueFrequenciesSet* clone() const { return new WordFromUniqueFrequenciesSet(*this); }

public:
  virtual void fireParameterChanged(const ParameterList& pl);

  /**
   *@ brief letter frequencies from given word frequencies. The
   * frequencies of a letter at a position is the sum of the
   * frequencies of the words that have this letter at this position.
   * The frequencies of each letter is the average of the frequencies
   * of that letter at all positions.
   */
  virtual void setFrequencies(const std::vector<double>& frequencies);

  virtual void updateFrequencies();

  /**
   *@ brief Return the n-th FrequenciesSet&
   **/
  const FrequenciesSet& getFrequenciesSetForLetter(size_t i) const { return *pFreq_; }

  size_t getLength() const { return length_; }

  void setNamespace(const std::string& prefix);

  std::string getDescription() const;
};

} // end of namespace bpp.

#endif // _WORDFREQUENCIESSET_H_


