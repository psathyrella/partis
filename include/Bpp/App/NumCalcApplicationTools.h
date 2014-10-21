//
// File: NumCalcApplicationTools.h
// Created by: Sylvain Gaillard
// Created on: Tue Jan 14:58:50 CET 2009
//

/*
Copyright or © or Copr. CNRS, (January 13, 2009)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

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

#ifndef _NUMCALCAPPLICATIONTOOLS_H_
#define _NUMCALCAPPLICATIONTOOLS_H_

#include "../Text/StringTokenizer.h"
#include "ApplicationTools.h"
#include "../Numeric/VectorTools.h"
#include "../Numeric/Function/FunctionTools.h"

namespace bpp
{

class NumCalcApplicationTools
{
  public:
    NumCalcApplicationTools();
    virtual ~NumCalcApplicationTools();

  public:
    /**
     * @brief Build a vector of integers as described by a string
     *
     * Build a vector of integers following a description like:
     * "2, 5, 7-10, 4" => [2, 5, 7, 8, 9, 10, 4]
     *
     * @author Sylvain Gaillard
     * @param s The string to parse.
     * @param delim Delimiter between elements.
     * @param seqdelim Delimiter between min and max for a sequence.
     * @return A vector containing the integers
     */
    static std::vector<int> seqFromString(const std::string& s, const std::string& delim = ",", const std::string& seqdelim = "-");

    /**
     * @brief Build a vector of double from a structured text description.
     *
     * The syntax may be one of the following:
     * - Specified values: 1.23, 2.34, 3.45, 4.56
     * - Sequence macro: seq(from=1.23,to=2.45,step=0.1)
     *   or              seq(from=1.23,to=2.45,size=5)
     *   The meaning of these to form is equivalent as the R function:
     *   The first one start from 1.23 and increment 0.1 until it reaches
     *   the 2.45 value, wheras the seocnd one will compute 3 values at
     *   equal distance from 1.23 and 2.45. the 'from' and 'to' values are
     *   included, except for the first syntax when the interval is not an
     *   exact multiple of the 'step' argument.
     *
     * @author Julien Dutheil
     * @param desc The string to parse.
     * @return A vector containing the corresponding values as double.
     * @throw Exception If the syntax describing the set is not correct.
     */
    static std::vector<double> getVector(const std::string& desc) throw (Exception);

    /**
     * @brief Returns the value of the Parameter of the given name
     *  if it exists; otherwise returns the default value.
     *
     * @author Laurent Gueguen
     * @param pl A parameter list to look in.
     * @param name A string name
     * @param x A double value
     */
    static double getDefaultValue(const ParameterList& pl, const std::string& name, double x);

    /**
     * @brief Design a parameter grid from input options.
     *
     * Example:
     * @code
     * grid.number_of_parameters=3
     * grid.parameter1.name=x
     * grid.parameter1.values=0.1,0.2,0.3,0.4,0.5
     * grid.parameter2.name=y
     * grid.parameter2.values=seq(from=0.1,to=0.5,step=0.1)
     * grid.parameter3.name=z
     * grid.parameter3.values=seq(from=0.1,to=0.5,size=5)
     * @endcode
     *
     * @param params           The attribute map where options may be found.
     * @param suffix           A suffix to be applied to the parameter name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param warn             Tell if a warning must be sent in case the parameter is not found.
     * @return a parameter grid object.
     */
    static ParameterGrid* getParameterGrid(
        std::map<std::string, std::string>& params,
        const std::string& suffix = "",
        bool suffixIsOptional = true,
        bool warn = true) throw (Exception);

};

} //End of namespace bpp.

#endif  //_NUMCALCAPPLICATIONTOOLS_H_
