//
// File: SequenceSimulationTools.h
// Created by: Julien Dutheil
// Created on: Wed Aug  24 16:25 2005
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

#ifndef _SEQUENCESIMULATIONTOOLS_H_
#define _SEQUENCESIMULATIONTOOLS_H_

#include "SiteSimulator.h"

//From Seqlib:
#include <Bpp/Seq/Container/SiteContainer.h>

//From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief Tools for sites and sequences simulation.
 */
class SequenceSimulationTools
{
	public:
		SequenceSimulationTools() {}
		~SequenceSimulationTools() {}

	public:
		
		/**
		 * @brief Simulate a set of sites knowing their rate.
		 *
		 * This method is rather slow.
		 * consider using a discrete rate distribution and a SequenceSimulator,
		 * which is realy faster.
		 * This method should be used only for continuous rate distribution, or
		 * as estimated from posterior rates for instance.
		 *
		 * @see SequenceSimulator
		 * @param simulator A SiteSimulator object to use to simulate sites.
		 * @param rates     the rates to use, one for each site to simulate.
		 * @return          A container with all simulated sites.
		 */
		static SiteContainer* simulateSites(const SiteSimulator& simulator, const std::vector<double>& rates);

		/**
		 * @brief Simulate a set of sites knowing their rate and ancestral state.
		 *
		 * This method is rather slow.
		 * consider using a discrete rate distribution and a SequenceSimulator,
		 * which is realy faster.
		 * This method should be used only for continuous rate distribution, or
		 * as estimated from posterior rates for instance.
		 *
		 * @see SequenceSimulator
		 * @param simulator A SiteSimulator object to use to simulate sites.
		 * @param rates     the rates to use, one for each site to simulate.
		 * @param states    the ancestral states to use, one for each site to simulate.
		 * @return          A container with all simulated sites.
		 */
		static SiteContainer* simulateSites(const SiteSimulator& simulator, const std::vector<double>& rates, const std::vector<int>& states)
			throw (Exception);

    /**
		 * @brief Simulate a set of sites knowing ancestral state.
		 *
		 * @see SequenceSimulator
		 * @param simulator A SiteSimulator object to use to simulate sites.
		 * @param states    the ancestral states to use, one for each site to simulate.
		 * @return          A container with all simulated sites.
		 */
		static SiteContainer* simulateSites(const SiteSimulator& simulator, const std::vector<int>& states)
			throw (Exception);

};

} //end of namespace bpp.

#endif //_SEQUENCESIMULATIONTOOLS_H_

