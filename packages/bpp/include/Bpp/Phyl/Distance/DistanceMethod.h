//
// File: DistanceMethod.h
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Wed jun 22 10:00 2005
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

#ifndef _DISTANCEMETHOD_H_
#define _DISTANCEMETHOD_H_

//From bpp-core:
#include <Bpp/Clonable.h>

//From the STL:
#include <string>

namespace bpp
{

class DistanceMatrix;
class Node;
class Tree;

/**
 * @brief General interface for distance-based phylogenetic reconstruction methods.
 */
class DistanceMethod:
  public virtual Clonable
{
	public:
		DistanceMethod() {}
		virtual ~DistanceMethod() {}

	public:
    /**
     * @brief Set the distance matrix to use.
     *
     * @param matrix The matrix to use.
     */
		virtual void setDistanceMatrix(const DistanceMatrix& matrix) = 0;
		
    /**
     * @brief Perform the clustering.
     */
    virtual void computeTree() = 0;
    
    /**
     * @return The computed tree.
     */
		virtual Tree* getTree() const = 0;

    /**
     * @return The name of the distance method.
     */
    virtual std::string getName() const = 0;

    /**
     * @param yn Enable/Disable verbose mode.
     */
    virtual void setVerbose(bool yn) = 0;

    /**
     * @return True if verbose mode is enabled.
     */
    virtual bool isVerbose() const = 0;
};

/**
 * @brief Interface for agglomerative distance methods.
 *
 * This interface does not contain any specific method and
 * is there only for "ontology" purposes. Specific methods
 * might be added later.
 */
class AgglomerativeDistanceMethod:
  public DistanceMethod
{
	public:
		AgglomerativeDistanceMethod() {}
		virtual ~AgglomerativeDistanceMethod() {}

};

} //end of namespace bpp.

#endif //_DISTANCEMETHOD_H_

