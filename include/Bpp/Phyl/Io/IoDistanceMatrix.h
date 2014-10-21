//
// File: IODistanceMatrix.h
// Created by: Julien Dutheil
// Created on: Wed Jun 08 15:43 2005
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

#ifndef _IODISTANCEMATRIX_H_
#define _IODISTANCEMATRIX_H_

#include <Bpp/Io/IoFormat.h>

// From the STL:
#include <iostream>
#include <fstream>

namespace bpp
{

class DistanceMatrix;

/**
 * @brief General interface for distance matrix I/O.
 */
class IODistanceMatrix:
  public virtual IOFormat
{
	public:
		IODistanceMatrix() {}
		virtual ~IODistanceMatrix() {}

	public:
		virtual const std::string getDataType() const { return "Distance matrix"; }
};

/**
 * @brief General interface for distance matrix readers.
 */
class IDistanceMatrix:
  public virtual IODistanceMatrix
{
	public:
		IDistanceMatrix() {}
		virtual ~IDistanceMatrix() {}

	public:
    /**
     * @brief Read a distance matrix from a file.
     *
     * @param path The file path.
     * @return A new distance matrix object.
     * @throw Exception If an error occured.
     */
		virtual DistanceMatrix* read(const std::string& path) const throw (Exception) = 0;
    /**
     * @brief Read a distance matrix from a stream.
     *
     * @param in The input stream.
     * @return A new distance matrix object.
     * @throw Exception If an error occured.
     */
		virtual DistanceMatrix* read(std::istream& in) const throw (Exception) = 0;
};

/**
 * @brief General interface for distance matrix writers.
 */
class ODistanceMatrix:
  public virtual IODistanceMatrix
{
	public:
		ODistanceMatrix() {}
		virtual ~ODistanceMatrix() {}

	public:
    /**
     * @brief Write a distance matrix to a file.
     *
     * @param dist A distance matrix object.
     * @param path The file path.
     * @param overwrite Tell if existing file must be overwritten.
     * Otherwise append to the file.
     * @throw Exception If an error occured.
     */
		virtual void write(const DistanceMatrix& dist, const std::string& path, bool overwrite) const throw (Exception) = 0;
    /**
     * @brief Write a distance matrix to a stream.
     *
     * @param dist A distance matrix object.
     * @param out The output stream.
     * @throw Exception If an error occured.
     */
		virtual void write(const DistanceMatrix& dist, std::ostream& out) const throw (Exception) = 0;
};

/**
 * @brief Partial implementation of the IDistanceMatrix interface.
 */
class AbstractIDistanceMatrix:
  public virtual IDistanceMatrix
{
	public:
		AbstractIDistanceMatrix() {}
		virtual ~AbstractIDistanceMatrix() {}

	public:
		virtual DistanceMatrix* read(const std::string& path) const throw (Exception)
		{
      std::ifstream input(path.c_str(), std::ios::in);
			DistanceMatrix* mat = read(input);
			input.close();
			return mat;
		}
		virtual DistanceMatrix* read(std::istream& in) const throw (Exception) = 0;
};

/**
 * @brief Partial implementation of the ODistanceMatrix interface.
 */
class AbstractODistanceMatrix:
  public virtual ODistanceMatrix
{
	public:
		AbstractODistanceMatrix() {}
		virtual ~AbstractODistanceMatrix() {}

	public:
		virtual void write(const DistanceMatrix& dist, const std::string& path, bool overwrite) const throw (Exception)
		{
			// Open file in specified mode
      std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out|std::ios::app));
			write(dist, output);
			output.close();
		}
		virtual void write(const DistanceMatrix& dist, std::ostream& out) const throw (Exception) = 0;
};

} //end of namespace bpp.

#endif //_IODISTANCEMATRIX_H_

