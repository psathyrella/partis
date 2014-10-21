//
// File: Bpp0SubstitutionModelFormat.h
// Created by: Laurent Guéguen
// Created on: mercredi 4 juillet 2012, à 13h 26
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

#ifndef _BPPOSUBSTITUTIONMODELFORMAT_H_
#define _BPPOSUBSTITUTIONMODELFORMAT_H_

#include "IoSubstitutionModelFactory.h"
#include "../Model/MixedSubstitutionModel.h"

// From bpp-seq
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{
/**
 * @brief Substitution model I/O in BppO format.
 *
 * Creates a new substitution model object according to model description syntax
 * (see the Bio++ Progam Suite manual for a detailed description of this syntax).
 *
 */

class BppOSubstitutionModelFormat :
  public ISubstitutionModel,
  public OSubstitutionModel
{
public:
  static unsigned char DNA;
  static unsigned char RNA;
  static unsigned char NUCLEOTIDE;
  static unsigned char PROTEIN;
  static unsigned char CODON;
  static unsigned char WORD;
  static unsigned char BINARY;
  static unsigned char ALL;

private:
  unsigned char alphabetCode_;
  bool allowCovarions_;
  bool allowMixed_;
  bool allowGaps_;
  bool verbose_;
  std::map<std::string, std::string> unparsedArguments_;
  const GeneticCode* geneticCode_;

public:
  /**
   * @brief Create a new BppOSubstitutionModelFormat object.
   *
   * @param alphabetCode     Bit saying which alphabets are allowed in the model specification.
   * @param allowCovarions   Tell is a covarion model can be returned.
   * @param allowMixed       Tell is a mixture model can be returned.
   * @param allowGaps        Tell is a gap model can be returned.
   * @param verbose          Tell if the construction is verbose.
   */
  BppOSubstitutionModelFormat(unsigned char alphabetCode, bool allowCovarions, bool allowMixed, bool allowGaps, bool verbose):
    alphabetCode_(alphabetCode),
    allowCovarions_(allowCovarions),
    allowMixed_(allowMixed),
    allowGaps_(allowGaps),
    verbose_(verbose),
    unparsedArguments_(),
    geneticCode_(0)
  {}

  BppOSubstitutionModelFormat(const BppOSubstitutionModelFormat& format):
    alphabetCode_(format.alphabetCode_),
    allowCovarions_(format.allowCovarions_),
    allowMixed_(format.allowMixed_),
    allowGaps_(format.allowGaps_),
    verbose_(format.verbose_),
    unparsedArguments_(format.unparsedArguments_),
    geneticCode_(format.geneticCode_)
  {}

  BppOSubstitutionModelFormat& operator=(const BppOSubstitutionModelFormat& format)
  {
    alphabetCode_      = format.alphabetCode_;
    allowCovarions_    = format.allowCovarions_;
    allowMixed_        = format.allowMixed_;
    allowGaps_         = format.allowGaps_;
    verbose_           = format.verbose_;
    unparsedArguments_ = format.unparsedArguments_;
    geneticCode_       = format.geneticCode_;
    return *this;
  }

  virtual ~BppOSubstitutionModelFormat() {}

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

  SubstitutionModel* read(const Alphabet* alphabet, const std::string& modelDescription, const SiteContainer* data = 0, bool parseArguments = true);

  const std::map<std::string, std::string>& getUnparsedArguments() const { return unparsedArguments_; }

  /**
   * @brief Write a substitution model to a stream.
   *
   * @param model A substitution model object;
   * @param out The output stream;
   * @param globalAliases parameters linked to global alias. The
   * output will be "name=alias_name";
   * @param writtenNames is the vector of the written
   * parameters so far [in, out];
   * @throw Exception If an error occured.
   */
  void write(const SubstitutionModel& model,
             OutputStream& out,
             std::map<std::string, std::string>& globalAliases,
             std::vector<std::string>& writtenNames) const;

  void setVerbose(bool verbose) { verbose_=verbose;}
  
private:
  MixedSubstitutionModel* readMixed_(const Alphabet* alphabet, const std::string& modelDescription, const SiteContainer* data);

  SubstitutionModel* readWord_(const Alphabet* alphabet, const std::string& modelDescription, const SiteContainer* data);

  /**
   * @brief Set parameter initial values of a given model according to options.
   *
   * Parameters actually depends on the model passed as argument.
   * See getSubstitutionModel for more information.
   *
   * This function is mainly for internal usage, you're probably looking for the getSubstitutionModel or getSubstitutionModelSet function.
   *
   * @param model                   The model to set.
   * @param data   A pointer toward the SiteContainer for which the substitution model is designed.
   *               The alphabet associated to the data must be of the same type as the one specified for the model.
   *               May be equal to NULL, but in this case use_observed_freq option will be unavailable.
   * @throw Exception if an error occured.
   */
  void initialize_(SubstitutionModel& model, const SiteContainer* data) throw (Exception);

  void writeMixed_(const MixedSubstitutionModel& model,
             OutputStream& out,
             std::map<std::string, std::string>& globalAliases,
             std::vector<std::string>& writtenNames) const;
};

} // end of namespace bpp.

#endif // _BPPOSUBSTITUTIONMODELFORMAT_H_

