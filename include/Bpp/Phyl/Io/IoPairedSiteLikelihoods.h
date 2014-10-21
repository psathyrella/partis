//
// File: PairedSiteLikelihoods.h
// Created by: Nicolas Rochette
// Created on: January 10, 2011
//

/*
   Copyright or © or Copr. Bio++ Developent Team, (November 16, 2004)

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

#ifndef _IOPAIREDSITELIKELIHOODS_H_
#define _IOPAIREDSITELIKELIHOODS_H_

// From the STL:
#include <vector>
#include <string>
#include <sstream>

// From Bio++
#include <Bpp/Io/IoFormat.h>
#include "../Likelihood/PairedSiteLikelihoods.h"

namespace bpp
{
/**
 * @brief Base class for I/O on paired-site likelihoods.
 *
 * @author Nicolas Rochette
 */
class IOPairedSiteLikelihoods : public virtual IOFormat
{};


/**
 * @brief This class provides I/O for the Tree-Puzzle/RAxML (phylip-like) paired-site likelihoods format.
 *
 * @author Nicolas Rochette
 */
class IOTreepuzzlePairedSiteLikelihoods : public virtual IOPairedSiteLikelihoods
{
public:
  /**
   * @brief Read paired-site likelihoods from a Treepuzzle/RAxML-formatted stream.
   *
   * @throw Exception If the format is not recognized.
   */
  static PairedSiteLikelihoods read(std::istream& is) throw (Exception);

  /**
   * @brief Read paired-site likelihoods from a Treepuzzle/RAxML-formatted file.
   *
   * @throw Exception If the format is not recognized.
   */
  static PairedSiteLikelihoods read(const std::string& path) throw (Exception);

  /**
   * @brief Write paired-site likelihoods to a stream.
   *
   * @param psl The PairedSiteLikelihoods object to write.
   * @param os The output stream.
   * @param delim The delimiter between model names and likelihoods. The defaut is a tab but two spaces might be used.
   */
  static void write(const PairedSiteLikelihoods& psl, std::ostream& os, const std::string& delim = "\t");

  /**
   * @brief Write paired-site likelihoods to a file.
   *
   * @param psl The PairedSiteLikelihoods object to write.
   * @param path The path of the output file.
   * @param delim The delimiter between model names and likelihoods. (The defaut is a tab but two spaces might be used.)
   */
  static void write(const PairedSiteLikelihoods& psl, const std::string& path, const std::string& delim = "\t");
};


/**
 * @brief This class provides input for the Phyml paired-site likelihoods format.
 *
 * @author Nicolas Rochette
 */
class IOPhymlPairedSiteLikelihoods : public virtual IOPairedSiteLikelihoods
{
public:
  /**
   * @brief Read paired-site likelihoods from a Phyml-formatted stream.
   * @throw Exception If the format is not recognized.
   */
  static std::vector<double> read(std::istream& is) throw (Exception);

  /**
   * @brief Read Phyml paired-site likelihoods from a Phyml-formatted file.
   * @throw Exception If the format is not recognized.
   */
  static std::vector<double> read(const std::string& path) throw (Exception);
};
} // namespace bpp
#endif
