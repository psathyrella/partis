//
// File: Bpp0FrequenciesSetFormat.h
// Created by: Laurent Guéguen
// Created on: lundi 9 juillet 2012, à 12h 57
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

#ifndef _BPPOFREQUENCIESSETFORMAT_H_
#define _BPPOFREQUENCIESSETFORMAT_H_

#include "IoFrequenciesSetFactory.h"

//From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{
/**
 * @brief Frequencies set I/O in BppO format.
 *
 * Allow to create a new frequencies set object according to model description syntax
 * (see the Bio++ Progam Suite manual for a detailed description of this syntax).
 */
class BppOFrequenciesSetFormat :
  public virtual IFrequenciesSet,
  public virtual OFrequenciesSet
{
public:
  static unsigned char DNA;
  static unsigned char RNA;
  static unsigned char NUCLEOTIDE;
  static unsigned char PROTEIN;
  static unsigned char CODON;
  static unsigned char WORD;
  static unsigned char ALL;

private:
  unsigned char alphabetCode_;
  bool verbose_;
  std::map<std::string, std::string> unparsedArguments_;
  const GeneticCode* geneticCode_;

public:
  BppOFrequenciesSetFormat(unsigned char alphabetCode, bool verbose):
    alphabetCode_(alphabetCode),
    verbose_(verbose),
    unparsedArguments_(),
    geneticCode_(0)
  {}

  BppOFrequenciesSetFormat(const BppOFrequenciesSetFormat& format):
    alphabetCode_(format.alphabetCode_),
    verbose_(format.verbose_),
    unparsedArguments_(format.unparsedArguments_),
    geneticCode_(format.geneticCode_)
  {}

  BppOFrequenciesSetFormat& operator=(const BppOFrequenciesSetFormat& format)
  {
    alphabetCode_      = format.alphabetCode_;
    verbose_           = format.verbose_;
    unparsedArguments_ = format.unparsedArguments_;
    geneticCode_       = format.geneticCode_;
    return *this;
  }

  virtual ~BppOFrequenciesSetFormat() {}

public:
  const std::string getFormatName() const { return "BppO"; }

  const std::string getFormatDescription() const { return "Bpp Options format."; }

  /**
   * @brief Set the genetic code to use in case a codon frequencies set should be built.
   *
   * @param gCode The genetic code to use.
   */
  void setGeneticCode(const GeneticCode* gCode) {
    geneticCode_ = gCode;
  }

  FrequenciesSet* read(
      const Alphabet* alphabet,
      const std::string& freqDescription,
      const SiteContainer* data,
      bool parseArguments = true);

  const std::map<std::string, std::string>& getUnparsedArguments() const { return unparsedArguments_; }

  void write(
      const FrequenciesSet* pfreqset,
      OutputStream& out,
      std::vector<std::string>& writtenNames) const;

private:
  void initialize_(FrequenciesSet& freqSet, const SiteContainer* data);
};

} // end of namespace bpp.

#endif // _BPPOFREQUENCIESSETFORMAT_H_

