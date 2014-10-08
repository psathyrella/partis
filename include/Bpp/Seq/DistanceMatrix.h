//
// File: DistanceMatrix.h
// Created on: Wed jun 08 10:39 2005
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

#ifndef _DISTANCEMATRIX_H_
#define _DISTANCEMATRIX_H_

// From the STL:
#include <vector>
#include <string>
#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/VectorExceptions.h> //DimensionException
#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{

/**
 * @brief A Matrix class to store phylogenetic distances.
 */
class DistanceMatrix:
  public virtual RowMatrix<double>
{

	private:
    std::vector<std::string> names_;

	public:

    /**
     * @brief Build a new distance matrix with specified names.
     *
     * The dimension of the matrix will be equal to the number of names
     *
     * @param names The names to use.
     */
		DistanceMatrix(const std::vector<std::string>& names):
      RowMatrix<double>(names.size(), names.size()), names_(names)
		{
			reset();
		}

		/**
     * @brief Build a new distance matrix with specified size.
     *
     * Row names will be named 'Taxon 0', 'Taxon 1', and so on.
     *
     * @param n The size of the matrix.
     */
    DistanceMatrix(size_t n):
      RowMatrix<double>(n, n), names_(n)
		{
      resize(n);
		}

		virtual ~DistanceMatrix() {}

		DistanceMatrix(const DistanceMatrix& dist): RowMatrix<double>(dist), names_(dist.names_)	{}

		DistanceMatrix& operator=(const DistanceMatrix& dist)
		{
			size_t n = dist.size();
			resize(n);
			for(size_t i = 0; i < n; ++i)
      {
				for(size_t j = 0; j < n; ++j)
        {
					operator()(i, j) = dist(i, j);
				}
			}
			names_ = dist.names_;
			return *this;
		}
		
	public:

    /**
     * @brief Reset the distance matrix: all distances are set to 0.
     */
		void reset()
		{
			size_t n = size();
			for (size_t i = 0; i < n; i++)
      {
				for (size_t j = 0; j < n; j++)
        {
					operator()(i, j) = 0;
				}
			}
		}
		
    /**
     * @return The dimension of the matrix.
     */
		size_t size() const { return names_.size(); }

    /**
     * @return The names associated to the matrix.
     */
    const std::vector<std::string>& getNames() const { return names_; }

    /**
     * @return The ith name.
     * @param i Name index.
     * @throw IndexOutOfBoundsException If i is not a valid index.
     */
    const std::string& getName(size_t i) const throw (IndexOutOfBoundsException)
    { 
      if (i >= size()) throw IndexOutOfBoundsException("DistanceMatrix::getName. Invalid indice.", i, 0, size());
      return names_[i];
    }
    
    /**
     * @brief Set the ith name.
     * 
     * @param i Name index.
     * @param name The new name.
     * @throw IndexOutOfBoundsException If i is not a valid index.
     */
		void setName(size_t i, const std::string& name) throw (IndexOutOfBoundsException)
		{
			if (i >= size()) throw IndexOutOfBoundsException("DistanceMatrix::setName. Invalid indice.", i, 0, size());
			names_[i] = name;
		}

    /**
     * @brief Set the names associated to the matrix.
     * 
     * @param names Matrix names.
     * @throw DimensionException If 'names' have not the same size as the matrix.
     */
		void setNames(const std::vector<std::string>& names) throw (DimensionException)
		{
			if (names.size() != names_.size()) throw DimensionException("DistanceMatrix::setNames. Invalid number of names.", names.size(), names_.size());
			names_ = names;
		}

    /**
     * @brief Get the index of a given name.
     *
     * @param name The name to look for.
     * @return The position of the name.
     * @throw Exception If no names are attached to this matrix, or if the name was not found.
     */
    size_t getNameIndex(const std::string& name) const throw (Exception);

    /**
     * @brief Change the dimension of the matrix.
     *
     * @param n the new dimension of the matrix.
     */
    void resize(size_t n) {
      RowMatrix<double>::resize(n, n);
      names_.resize(n);
			for (size_t i = 0; i < n; ++i)
        names_[i] = "Taxon " + TextTools::toString(i);
      reset();
    }

    /**
     * @brief Access by name.
     *
     * @param iName Name 1 (row)
     * @param jName Name 2 (column)
     * @return A reference toward the specified distance.
     * @throw Exception if the matrix has no name of if one of the name do not match existing names.
     */
    virtual const double& operator()(const std::string& iName, const std::string& jName) const throw (Exception)
    {
      size_t i = getNameIndex(iName);
      size_t j = getNameIndex(jName);
      return operator()(i,j);
    }

    /**
     * @brief Access by name.
     *
     * @param iName Name 1 (row)
     * @param jName Name 2 (column)
     * @return A reference toward the specified distance.
     * @throw Exception if the matrix has no name of if one of the name do not match existing names.
     */
    virtual double& operator()(const std::string& iName, const std::string& jName) throw (Exception)
    {
      size_t i = getNameIndex(iName);
      size_t j = getNameIndex(jName);
      return operator()(i,j);
    }

    virtual const double& operator()(size_t i, size_t j) const
    {
      return RowMatrix<double>::operator()(i, j);
    }
    virtual double& operator()(size_t i, size_t j)
    {
      return RowMatrix<double>::operator()(i, j);
    }
};

} //end of namespace bpp.

#endif //_DISTANCEMATRIX_H_

