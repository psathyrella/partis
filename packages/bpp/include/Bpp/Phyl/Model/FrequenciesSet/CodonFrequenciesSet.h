//
// File: CodonFrequenciesSet.h
// Created by: laurent Gueguen
// Created on: lundi 2 avril 2012, à 14h 03
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

#ifndef _CODONFREQUENCIESSET_H_
#define _CODONFREQUENCIESSET_H_

#include "WordFrequenciesSet.h"
#include "FrequenciesSet.h"
#include "ProteinFrequenciesSet.h"

#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Numeric/Prob/Simplex.h>

namespace bpp
{

  /**
   * @brief Parametrize a set of state frequencies for codons.
   */
  class CodonFrequenciesSet :
    public virtual FrequenciesSet
  {
  public:

#ifndef NO_VIRTUAL_COV
    CodonFrequenciesSet* clone() const = 0;
    
    const CodonAlphabet* getAlphabet() const = 0;
#endif
    
  public:
    /**
     * @return The associated genetic code.
     */
    virtual const GeneticCode* getGeneticCode() const = 0;

    /**
     * @brief A helper function that provide frequencies set for codon models
     * according to PAML option.
     *
     * @param option A code describing the option, one of F61, F1X4 or F3X4.
     * @param gCode The genetic code to use. The underlying codon alphabet object will be passed to the FrequenciesSet instance.
     * @param mgmtStopFreq the optional way the frequencies assigned
     * to the stop codons are redistributed to the other codons, with
     * F1X4 and F3X4 options. The available values are:
     *  - uniform : each stop frequency is distributed evenly
     *  - linear : each stop frequency is distributed to the neighbour
     *     codons (ie 1 substitution away), in proportion to each
     *     target codon frequency.
     *  - quadratic (default): each stop frequency is distributed to the
     *     neighbour codons (ie 1 substitution away), in proportion to
     *     the square of each target codon frequency.
     *
     */
    static FrequenciesSet* getFrequenciesSetForCodons(short option, const GeneticCode* gCode, const std::string& mgmtStopFreq = "quadratic");
    
    static const short F0;
    static const short F1X4;
    static const short F3X4;
    static const short F61;
    
  };


  /**
   * @brief A generic FrequenciesSet for Full Codon alphabets.
   *
   * It is very similar to FullFrequencySet, but only the non-stop codon
   *   frequencies are parameterized.
   */
  class FullCodonFrequenciesSet :
    public virtual CodonFrequenciesSet,
    public AbstractFrequenciesSet
  {
  protected:
    const GeneticCode* pgc_;

  public:

    /**
     * @brief Construction with uniform frequencies on the letters of
     * the alphabet. The stop codon frequencies are null.
     */
    FullCodonFrequenciesSet(const GeneticCode* gCode, bool allowNullFreqs = false, const std::string& name = "Full");
    FullCodonFrequenciesSet(const GeneticCode* gCode, const std::vector<double>& initFreqs, bool allowNullFreqs = false, const std::string& name = "Full");
    
    FullCodonFrequenciesSet(const FullCodonFrequenciesSet& fcfs);
    FullCodonFrequenciesSet& operator=(const FullCodonFrequenciesSet& fcfs);

#ifndef NO_VIRTUAL_COV
    FullCodonFrequenciesSet*
#else
    Clonable*
#endif
    clone() const { return new FullCodonFrequenciesSet(*this); }

  public:
    const GeneticCode* getGeneticCode() const { return pgc_; }

    /**
     * @brief the given frequencies are normalized such that the sum of
     * the frequencies on the non-stop codons equals 1.
     *
     */
    void setFrequencies(const std::vector<double>& frequencies);

#ifndef NO_VIRTUAL_COV
    const CodonAlphabet* getAlphabet() const
    {
      return dynamic_cast<const CodonAlphabet*>(AbstractFrequenciesSet::getAlphabet());
    }
#endif

  protected:
    void fireParameterChanged(const ParameterList& parameters);
  };


  /**
   * @brief FrequenciesSet useful for homogeneous and stationary models, codon implementation
   *
   * This set contains no parameter.
   */
  class FixedCodonFrequenciesSet :
    public virtual CodonFrequenciesSet,
    public AbstractFrequenciesSet
  {
  protected:
    const GeneticCode* pgc_;

  public:
    FixedCodonFrequenciesSet(const GeneticCode* gCode, const std::vector<double>& initFreqs, const std::string& name = "Fixed");

    /**
     * @brief Construction with uniform frequencies on the letters of
     * the alphabet. The stop codon frequencies are null.
     */
    FixedCodonFrequenciesSet(const GeneticCode* gCode, const std::string& name = "Fixed");

    FixedCodonFrequenciesSet(const FixedCodonFrequenciesSet& fcfs):
      AbstractFrequenciesSet(fcfs),
      pgc_(fcfs.pgc_)
    {}

    FixedCodonFrequenciesSet& operator=(const FixedCodonFrequenciesSet& fcfs)
    {
      AbstractFrequenciesSet::operator=(fcfs);
      pgc_ = fcfs.pgc_;
      return *this;
    }

#ifndef NO_VIRTUAL_COV
    FixedCodonFrequenciesSet*
#else
    Clonable*
#endif
    clone() const { return new FixedCodonFrequenciesSet(*this); }

  public:
    const GeneticCode* getGeneticCode() const { return pgc_; }

#ifndef NO_VIRTUAL_COV
    const CodonAlphabet* getAlphabet() const
    {
      return dynamic_cast<const CodonAlphabet*>(AbstractFrequenciesSet::getAlphabet());
    }
#endif
    /**
     * @brief the given frequencies are normalized such thaat the sum of
     * the frequencies on the non-stop codons equals 1.
     *
     */
    void setFrequencies(const std::vector<double>& frequencies);

  protected:
    void fireParameterChanged(const ParameterList& parameters) {}
  };

  /**
   * @brief FrequenciesSet integrating ProteinFrequenciesSet inside
   * CodonFrequenciesSet. In this case, FrequencieSet defined inside
   * each amino acid is parametrized as a FullFrequenciesSet. Hence
   * there are 61-20=41 parameters in addition of the parameters of the
   * ProteinFrequenciesSet.
   *
   * The parametrization depends on the method used.
   * Default method is 1 (ie global ratio).
   *
   * @see Simplex
   *
   */

  class FullPerAACodonFrequenciesSet :
    public virtual CodonFrequenciesSet,
    public AbstractFrequenciesSet
  {
  
  private:
    const GeneticCode* pgc_;
    std::auto_ptr<ProteinFrequenciesSet> ppfs_;

    /*
     *@ brief vector of the simplexes, one for each AA
     *
     */
  
    std::vector<Simplex> vS_;

    void updateFrequencies();
  
  public:
  
    /**
     * @brief Create a new FullPerAACodonFrequenciesSet object.
     *
     * @param gencode The genetic code to use.
     * @param ppfs The protein frequencies to use. The codon
     * frequencies set will own the instance of the protein
     * frequencies set.
     * @param method the method used for parametrization.
     */
    FullPerAACodonFrequenciesSet(const GeneticCode* gencode, ProteinFrequenciesSet* ppfs, unsigned short method = 1);

    /**
     * @brief Construction with fixed uniform frequencies on the amino acids.
     * The stop codon frequencies are null.
     * @param gencode The genetic code to use.
     * @param method the method used for parametrization.
     */

    FullPerAACodonFrequenciesSet(const GeneticCode* gencode, unsigned short method = 1);

    FullPerAACodonFrequenciesSet(const FullPerAACodonFrequenciesSet& ffs);

    FullPerAACodonFrequenciesSet& operator=(const FullPerAACodonFrequenciesSet& ffs);

    virtual ~FullPerAACodonFrequenciesSet() {}
 
#ifndef NO_VIRTUAL_COV
    FullPerAACodonFrequenciesSet*
#else
    Clonable*
#endif
    clone() const { return new FullPerAACodonFrequenciesSet(*this); }

  public:
#ifndef NO_VIRTUAL_COV
    const CodonAlphabet* getAlphabet() const
    {
      return dynamic_cast<const CodonAlphabet*>(AbstractFrequenciesSet::getAlphabet());
    }
#endif
    
    const GeneticCode* getGeneticCode() const { return pgc_; }
    
    /**
     * @brief the given frequencies are normalized such thaat the sum of
     * the frequencies on the non-stop codons equals 1.
     *
     */
    void setFrequencies(const std::vector<double>& frequencies);

    void setNamespace(const std::string& prefix);

    const ProteinFrequenciesSet* getProteinFrequenciesSet() const
    {
      return ppfs_.get();
    }

    unsigned short getMethod() const
    {
      return (vS_.size()?vS_[0].getMethod():1);
    }
    
  protected:
    void fireParameterChanged(const ParameterList& parameters);
  };


  /**
   * @brief the Frequencies in codons are the product of Independent
   * Frequencies in letters with the frequencies of stop codons set to
   * zero.
   *
   * @author Laurent Guéguen
   */
  class CodonFromIndependentFrequenciesSet :
    public virtual CodonFrequenciesSet,
    public WordFromIndependentFrequenciesSet
  {
  private:

    // a map associating stop codons numbers with numbers of neighbour non-stop codons
    std::map<int, Vint> mStopNeigh_;

    unsigned short mgmtStopFreq_;

    const GeneticCode* pgc_;
    
  public:
    /**
     * @brief Constructor from a CodonAlphabet* and a vector of different FrequenciesSet*.
     * Throws an Exception if their lengths do not match.
     *
     * @param gCode a pointer to the genetic code to use.
     * @param freqvector a vector of pointers to the phase specific FrequenciesSets
     * @param name the optional name of the FrequenciesSet (default codon)
     * @param mgmtStopFreq the optional way the frequencies assigned to the
     * stop codons are redistributed to the other codons. The
     * available values are:
     *  - uniform : each stop frequency is distributed evenly
     *  - linear : each stop frequency is distributed to the neighbour
     *     codons (ie 1 substitution away), in proportion to each
     *     target codon frequency.
     *  - quadratic (default): each stop frequency is distributed to the
     *     neighbour codons (ie 1 substitution away), in proportion to
     *     the square of each target codon frequency.
     *
     */
    CodonFromIndependentFrequenciesSet(const GeneticCode* gCode, const std::vector<FrequenciesSet*>& freqvector, const std::string& name = "Codon", const std::string& mgmtStopFreq = "quadratic");
  
    CodonFromIndependentFrequenciesSet(const CodonFromIndependentFrequenciesSet& iwfs);

    virtual ~CodonFromIndependentFrequenciesSet(){};
  
    CodonFromIndependentFrequenciesSet& operator=(const CodonFromIndependentFrequenciesSet& iwfs);
  
    CodonFromIndependentFrequenciesSet* clone() const { return new CodonFromIndependentFrequenciesSet(*this); }

    const CodonAlphabet* getAlphabet() const;
    
    const GeneticCode* getGeneticCode() const { return pgc_; }
  
    /**
     * @ brief Update the frequencies given the parameters.
     */
    void updateFrequencies();
  };


  /**
   * @brief the Frequencies in codons are the product of the frequencies
   * for a unique FrequenciesSet in letters, with the frequencies of
   * stop codons set to zero.
   *
   * @author Laurent Guéguen
   */

  class CodonFromUniqueFrequenciesSet :
    public virtual CodonFrequenciesSet,
    public WordFromUniqueFrequenciesSet
  {
  private:

    // a map associating stop codons numbers with numbers of neighbour non-stop codons
    std::map<int, Vint> mStopNeigh_;

    unsigned short mgmtStopFreq_;
   
    const GeneticCode* pgc_;

  public:
    /**
     * @brief Constructor from a CodonAlphabet* and a FrequenciesSet*
     *  repeated three times.
     *
     * @param gCode a pointer to a genetic code.
     * @param pfreq a pointer to the nucleotidic FrequenciesSet
     * @param name the optional name of the FrequenciesSet (default codon)
     * @param mgmtStopFreq the optional way the frequencies assigned to the
     * stop codons are redistributed to the other codons. The
     * available values are:
     *  - uniform : each stop frequency is distributed evenly
     *  - linear : each stop frequency is distributed to the neighbour
     *      codons (ie 1 substitution away), in proportion to each
     *      target codon frequency.
     *  - quadratic (default): each stop frequency is distributed to the
     *      neighbour codons (ie 1 substitution away), in proportion to
     *      the square of each target codon frequency.
     */
    CodonFromUniqueFrequenciesSet(
        const GeneticCode* gCode,
        FrequenciesSet* pfreq,
        const std::string& name = "Codon",
        const std::string& mgmtStopFreq = "quadratic");
  
    CodonFromUniqueFrequenciesSet(const CodonFromUniqueFrequenciesSet& iwfs);
  
    virtual ~CodonFromUniqueFrequenciesSet() {};
  
    CodonFromUniqueFrequenciesSet& operator=(const CodonFromUniqueFrequenciesSet& iwfs);
  
    CodonFromUniqueFrequenciesSet* clone() const { return new CodonFromUniqueFrequenciesSet(*this); }
  
    const CodonAlphabet* getAlphabet() const;

    const GeneticCode* getGeneticCode() const { return pgc_; }
    
    /**
     * @brief Update the frequencies given the parameters.
     *
     */
    void updateFrequencies();
  };


} // end of namespace bpp.

#endif // _CODONFREQUENCIESSET_H_


