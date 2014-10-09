//
// File: Bpp0DiscreteDistributionFormat.h
// Created by: Laurent Guéguen
// Created on: lundi 3 septembre 2012, à 14h 42
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

#ifndef _BPPODISCRETEDISTRIBUTIONFORMAT_H_
#define _BPPODISCRETEDISTRIBUTIONFORMAT_H_

#include "IoDiscreteDistributionFactory.h"

namespace bpp
{
/**
 * @brief Discrete Distribution I/O in BppO format.
 *
 * Creates a new discrete distribution object according to
 * distribution description syntax (see the Bio++ Progam Suite
 * manual for a detailed description of this syntax).
 *
 */
class BppODiscreteDistributionFormat :
  public virtual IDiscreteDistribution,
  public virtual ODiscreteDistribution
{
protected:
  bool verbose_;
  std::map<std::string, std::string> unparsedArguments_;

public:
  BppODiscreteDistributionFormat(bool verbose = true) : verbose_(verbose), unparsedArguments_() {}
  virtual ~BppODiscreteDistributionFormat() {}

public:
  const std::string getFormatName() const { return "BppO"; }

  const std::string getFormatDescription() const { return "Bpp Options format."; }

  DiscreteDistribution* read(const std::string& distDescription, bool parseArguments = true);

  const std::map<std::string, std::string>& getUnparsedArguments() const { return unparsedArguments_; }

  void write(const DiscreteDistribution& dist,
             OutputStream& out,
             std::map<std::string, std::string>& globalAliases,
             std::vector<std::string>& writtenNames) const;

protected:
  /**
   * @brief Set parameter initial values of a given distribution according to options.
   *
   * Parameters actually depends on the distribution passed as argument.
   *
   * @param rDist The distribution to set up.
   * @throw Exception if an error occured.
   */
  void initialize_(DiscreteDistribution& rDist) throw (Exception);
};

} // end of namespace bpp.

#endif // _BPPODISCRETEDISTRIBUTIONFORMAT_H_

