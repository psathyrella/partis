//
// File CodonSiteTools.h
// Author : Sylvain Glémin
// Last modification : November 2005
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

#ifndef _CODONSITETOOLS_H_
#define _CODONSITETOOLS_H_

#include "SymbolListTools.h"
#include "Site.h"
#include "Alphabet/CodonAlphabet.h"
#include "GeneticCode/GeneticCode.h"

#include <Bpp/Exceptions.h>

// From the STL:
#include <map>

namespace bpp
{

/**
 * @brief Utilitary functions for codon sites.
 */

class CodonSiteTools:
  public SymbolListTools
{
	public:
		CodonSiteTools() {}
		virtual ~CodonSiteTools() {}

	public:

		/**
		 * @brief Method to know if a codon site contains gap(s) or stop codons.
		 *
		 * @param site a Site
     * @param gCode The genetic code according to which stop codons are specified.
		 * @throw AlphabetException If the alphabet associated to the site is not a codon alphabet.
		 */
		static bool hasGapOrStop(const Site& site, const GeneticCode& gCode) throw (AlphabetException);

		/**
		 * @brief Method to know if a codon site contains stop codon or not
		 *
		 * @param site a Site
     * @param gCode The genetic code according to which stop codons are specified.
		 * @throw AlphabetException If the alphabet associated to the site is not a codon alphabet.
		 */
		static bool hasStop(const Site& site, const GeneticCode& gCode) throw (AlphabetException);

		/**
		 * @brief Method to know if a polymorphic codon site is polymorphic at only one site
		 *
		 * @param site a Site
		 * @throw AlphabetException  If the alphabet associated to the site is not a codon alphabet.
		 * @throw EmptySiteException If the site has size 0.
		 */
		static bool isMonoSitePolymorphic(const Site& site) throw (AlphabetException, EmptySiteException);

		/**
		 * @brief Method to know if polymorphism at a codon site is synonymous
		 *
		 * @param site a Site
		 * @param gCode a GeneticCode
		 * @throw AlphabetException          If the alphabet associated to the site is not a codon alphabet.
		 * @throw AlphabetMismatchException  If the codon alphabet of the site do not match the codon alphabet of the genetic code.
		 * @throw EmptySiteException         If the site has size 0.
		 */
    static bool isSynonymousPolymorphic(const Site& site, const GeneticCode& gCode)
			  throw (AlphabetException, AlphabetMismatchException, EmptySiteException);

		/**
		 * @brief generate a codon site without rare variants
		 *
		 * Rare variants are replaced by the most frequent allele.
		 * This method is used to exclude rare variants in some analyses as in McDonald-Kreitman Test
		 * (McDonald & Kreitman, 1991, Nature 351 pp652-654).
		 * For an application, see for example (Fay et al. 2001, Genetics 158 pp 1227-1234).
		 *
		 * @param site a Site
     * @param gCode The genetic code according to which stop codons are specified.
		 * @param freqmin a double, allele in frequency stricly lower than freqmin are replaced
		 * @throw AlphabetException  If the alphabet associated to the site is not a codon alphabet.
		 * @throw EmptySiteException If the site has size 0.
		 *
		 */
	  static Site* generateCodonSiteWithoutRareVariant(const Site& site, const GeneticCode& gCode, double freqmin)
			  throw(AlphabetException, EmptySiteException);


		/**
		 * @brief Compute the number of differences between two codons
		 *
		 * @param i a int
		 * @param j a int
		 * @param ca a CodonAlphabet
		 */
	  static size_t numberOfDifferences(int i, int j, const CodonAlphabet& ca);

		/**
		 * @brief Compute the number of synonymous differences between two codons
		 *
		 * For complex codon:
		 * If minchange = false (default option) the different paths are equally weighted.
		 * If minchange = true the path with the minimum number of non-synonymous change is chosen.
		 * Paths included stop codons are excluded.
		 * @param i a int
		 * @param j a int
		 * @param gCode a GeneticCode
		 * @param minchange a boolean set by default to false
		 */
   	static double numberOfSynonymousDifferences(int i, int j, const GeneticCode& gCode, bool minchange=false);

		/**
		 * @brief Compute the synonymous pi per codon site
		 *
		 * The following formula is used:
		 * @f[
		 * pi = frac{n}{n-1}\sum_{i,j}x_{i}x_{j}P_{ij}
		 * @f]
		 * where n is the number of sequence, \f$x_i\f$ and \f$x_j\f$ the frequencies of each codon type occuring at the site
		 * \f$P_{i,j}\f$ the number of synonymous difference between these codons.
		 * Be careful: here, pi is not normalized by the number of synonymous sites.
		 *
		 * If minchange = false (default option) the different paths are equally weighted.
		 * If minchange = true the path with the minimum number of non-synonymous change is chosen.
		 *
		 * @param site a Site
		 * @param gCode a GeneticCode
		 * @param minchange a boolean set by default to false
		 * @throw AlphabetException          If the alphabet associated to the site is not a codon alphabet.
		 * @throw AlphabetMismatchException  If the codon alphabet of the site do not match the codon alphabet of the genetic code.
		 * @throw EmptySiteException         If the site has size 0.
		 */
		static double piSynonymous(const Site& site, const GeneticCode& gCode, bool minchange = false)
			  throw (AlphabetException, AlphabetMismatchException, EmptySiteException);

		/**
		 * @brief Compute the non-synonymous pi per codon site
		 *
		 * The following formula is used:
		 * @f[
		 * pi = frac{n}{n-1}\sum_{i,j}x_{i}x_{j}P_{ij}
		 * @f]
		 * where n is the number of sequence, \f$x_i\f$ and \f$x_j\f$ the frequencies of each codon type occuring at the site
		 * \f$P_{i,j}\f$ the number of nonsynonymous difference between these codons.
		 * Be careful: here, pi is not normalized by the number of non-synonymous sites.
		 * If minchange = false (default option) the different paths are equally weighted.
		 * If minchange = true the path with the minimum number of non-synonymous change is chosen.
		 *
		 * @param site a Site
		 * @param gCode a GeneticCode
		 * @param minchange a boolean set by default to false
		 * @throw AlphabetException         If the alphabet associated to the site is not a codon alphabet.
		 * @throw AlphabetMismatchException If the codon alphabet of the site do not match the codon alphabet of the genetic code.
		 * @throw EmptySiteException        If the site has size 0.
		 */
		static double piNonSynonymous(const Site& site, const GeneticCode& gCode, bool minchange = false)
			  throw (AlphabetException, AlphabetMismatchException, EmptySiteException);


		/**
		 * @brief Return the number of synonymous positions of a codon
		 *
		 * A site is consider as x% synonymous if x% of the possible mutations are synonymous
		 * Transition/transversion ratio can be taken into account (use the variable ratio)
		 *
     * Unresolved codons and stop codon will return a value of 0.
     *
		 * @param i a int
		 * @param gCode a GeneticCode
		 * @param ratio a double set by default to 1
		 */
		static double numberOfSynonymousPositions(int i, const GeneticCode& gCode,  double ratio = 1.0) throw (Exception);

		/**
		 * @brief Return the mean number of synonymous positions per codon site
		 *
		 * A site is consider as x% synonymous if x% of the possible mutations are synonymous
		 * Transition/transversion ratio can be taken into account (use the variable ratio)
		 * The mean is computed over the VectorSite.
		 *
     * Unresolved and stop codons are counted as 0.
     *
		 * @param site a Site
		 * @param gCode a GeneticCode
		 * @param ratio a double 			Set by default to 1
		 * @throw AlphabetException         If the alphabet associated to the site is not a codon alphabet.
		 * @throw AlphabetMismatchException If the codon alphabet of the site do not match the codon alphabet of the genetic code.
		 * @throw EmptySiteException        If the site has size 0.
		 */
		static double meanNumberOfSynonymousPositions(const Site& site, const GeneticCode& gCode, double ratio = 1)
			  throw (AlphabetException, AlphabetMismatchException, EmptySiteException);

		/**
		 * @brief Return the number of subsitutions per codon site
		 *
		 * No recombination is assumed, that is in complex codon homoplasy is assumed.
		 * Example:
     * @code
		 * ATT
		 * ATT
		 * ATT
		 * ATC
		 * ATC
		 * AGT
		 * AGT
		 * AGC
     * @endcode
		 * Here, 3 substitutions are counted. Assuming that the last codon (AGC) is a recombinant between ATC and AGT
		 * would have lead to counting only 2 subsitutions.
		 *
		 * Rare variants (<= freqmin) can be excluded.
		 *
		 * @param site a Site
		 * @param gCode a GeneticCode
		 * @param freqmin a double 	 To exclude snp in frequency strictly lower than freqmin (by default freqmin = 0)
		 * @throw AlphabetException  If the alphabet associated to the site is not a codon alphabet.
		 * @throw EmptySiteException If the site has size 0.
		 */
		static size_t numberOfSubsitutions(const Site& site, const GeneticCode& gCode, double freqmin = 0.)
			  throw(AlphabetException, EmptySiteException);

		/**
		 * @brief Return the number of Non Synonymous subsitutions per codon site.
		 *
		 * It is assumed that the path linking amino acids only involved one substitution by step.
		 *
		 * Rare variants (<= freqmin) can be excluded.
		 * In case of complex codon, the path that gives the minimum number of non-synonymous changes
		 * is chosen. The argument minchange=true is sent to numberOfSynonymousDifferences used in this method.
		 * Otherwise, a non-integer number could be return.
		 *
		 * @param site a Site
		 * @param gCode a GeneticCode
		 * @param freqmin a double 			To exclude snp in frequency strictly lower than freqmin (by default freqmin = 0).
		 * @throw AlphabetException         If the alphabet associated to the site is not a codon alphabet.
		 * @throw AlphabetMismatchException If the codon alphabet of the site do not match the codon alphabet of the genetic code.
		 * @throw EmptySiteException        If the site has size 0.
		 */
		static size_t numberOfNonSynonymousSubstitutions(const Site& site, const GeneticCode& gCode, double freqmin = 0.)
			  throw (AlphabetException, AlphabetMismatchException, EmptySiteException);

		/**
		 * @brief Return a vector with the number of fixed synonymous and non-synonymous differences per codon site
		 *
		 * Compute the number of synonymous and non-synonymous differences between
		 * the concensus codon of SiteIn (i) and SiteOut (j), which are fixed within each alignement.
		 * Example:
     * @code
		 * SiteIn
		 *	ATT
		 *	ATT
		 *	ATC
     * @endcode
     * @code
		 * SiteOut
		 *	CTA
		 *	CTA
		 *	CTA
     * @endcode
		 * Here, the first position is non-synonymous different and fixed,
		 * the third position is synonymous different but not fixed (polymorphic in SiteIn).
		 * The return vector is thus [0,1].
		 * In case of complex codon, the path that gives the minimum number of non-synonymous changes
		 * is chosen. The argument minchange=true is sent to numberOfSynonymousDifferences used in this method.
		 * Otherwise, a non-integer number could be return.
		 *
		 * Rare variants (<= freqmin) can be excluded.
		 *
		 * @param siteIn a Site
		 * @param siteOut a Site
		 * @param i an integer
		 * @param j an integer
		 * @param gCode a GeneticCode
		 * @throw AlphabetException         If the alphabet associated to one of the sites is not a codon alphabet.
		 * @throw AlphabetMismatchException If the codon alphabet each the site do not match the codon alphabet of the genetic code.
		 * @throw EmptySiteException        If one of the sites has size 0.
		 */
		static std::vector<size_t> fixedDifferences(const Site& siteIn, const Site& siteOut, int i, int j, const GeneticCode& gCode)
			  throw (AlphabetException, AlphabetMismatchException, EmptySiteException);
    
    /**
     * @return True if all sequences have a fourfold degenerated codon in the site
     * (that is, if a mutation in the fourth position does not change the aminoacid).
     * @author Benoit Nabholz, Annabelle Haudry
     * @param site The site to analyze.
     * @param gCode The genetic code to use.
     */
    static bool isFourFoldDegenerated(const Site& site, const GeneticCode& gCode);

};

} //end of namespace bpp.

#endif	//_CONDONSITETOOLS_H_

