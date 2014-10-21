//
// File: BppORateDistributionFormat.h
// Created by: Laurent Guéguen and Julien Dutheil
// Created on: Fri 16 november 2012, at 13:44
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

#ifndef _BPPORATEDISTRIBUTIONFORMAT_H_
#define _BPPORATEDISTRIBUTIONFORMAT_H_

#include <Bpp/Io/BppODiscreteDistributionFormat.h>

namespace bpp
{
/**
 * @brief Rate Distribution I/O in BppO format.
 *
 * Creates a new discrete distribution object according to
 * distribution description syntax (see the Bio++ Progam Suite
 * manual for a detailed description of this syntax).
 *
 * Rate distributions are normalized and have a mean of 1, so that branch lengths are measured in mean number of substitutions per site.
 *
 * @see BppODiscreteDistribtution for a more generic parser.
 *
 */

class BppORateDistributionFormat:
  public BppODiscreteDistributionFormat
{
private:
  bool allowConstant_;

public:
  /**
   * @brief Build a new BppORateDistributionFormat object.
   *
   * @param allowConstant Is contant distribution allowed.
   */
  BppORateDistributionFormat(bool allowConstant):
    BppODiscreteDistributionFormat(),
    allowConstant_(allowConstant)
  {}

  virtual ~BppORateDistributionFormat() {}

public:

  DiscreteDistribution* read(const std::string& distDescription, bool parseArguments);

  void write(const DiscreteDistribution& dist,
             OutputStream& out,
             std::map<std::string, std::string>& globalAliases,
             std::vector<std::string>& writtenNames) const;
};

} // end of namespace bpp.

#endif // _BPPORATEDISTRIBUTIONFORMAT_H_

