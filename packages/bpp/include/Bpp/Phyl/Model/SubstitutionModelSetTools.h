//
// File: SubstitutionModelSetTools.h
// Created by: Julien Dutheil
// Created on: tue Sep 17 16:57 2007
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _SUBSTITUTIONMODELSETTOOLS_H_
#define _SUBSTITUTIONMODELSETTOOLS_H_

#include "SubstitutionModelSet.h"
#include "../Tree.h"

//From the STL:
#include <vector>
#include <map>

namespace bpp
{

/**
 * @brief Tools for automatically creating SubstitutionModelSet objects.
 */
class SubstitutionModelSetTools
{

  public:

    /**
     * @brief Create a SubstitutionModelSet object, corresponding to the homogeneous case.
     *
     * This class is mainly for testing purpose.
     *
     * @param model     The model to use.
     * @param rootFreqs A FrequenciesSet object to parametrize root frequencies.
     * @param tree      The tree to use for the construction of the set.
     */
    static SubstitutionModelSet* createHomogeneousModelSet(
        SubstitutionModel* model,
        FrequenciesSet* rootFreqs,
        const Tree* tree
      ) throw (AlphabetException, Exception);

    /**
     * @brief Create a SubstitutionModelSet object, with one model per branch.
     *
     * All branches share the same type of model, but allow one set of parameters per branch.
     * This is also possible to specify some parameters to be common to all branches.
     *
     * @param model                The model to use.
     * @param rootFreqs            A FrequenciesSet object to parametrize root frequencies.
     * @param tree                 The tree to use for the construction of the set.
     * @param globalParameterNames Common parameters for all branches.
     * All other parameters will be considered distinct for all branches.
     */
    static SubstitutionModelSet* createNonHomogeneousModelSet(
        SubstitutionModel* model,
        FrequenciesSet* rootFreqs,
        const Tree* tree,
        const std::vector<std::string>& globalParameterNames
      ) throw (AlphabetException, Exception);

};

} //end of namespace bpp.

#endif //_SUBSTITUTIONMODELSETTOOLS_H_

