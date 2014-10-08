//
// File ISequenceStream.h
// Author: Sylvain Gaillard
// Created: 18/08/2009
//

/*
Copyright or © or Copr. Bio++ Development Team, (August 18, 2009)

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

#ifndef _ISEQUENCESTREAM_H_
#define _ISEQUENCESTREAM_H_

#include "IoSequenceStream.h"
#include "../Sequence.h"
#include "../Alphabet/Alphabet.h"
#include <Bpp/Exceptions.h>

namespace bpp
{

/**
 * @brief The ISequenceStream interface.
 *
 * Interface for streaming sequences input.
 *
 * @author Sylvain Gaillard
 */
class ISequenceStream: public virtual IOSequenceStream
{
	public:
		ISequenceStream() {}
		virtual ~ISequenceStream() {}

	public:
    /**
     * @brief Read sequence from stream.
     *
     * Read one sequence from a stream.
     *
     * @param input The stream to read.
     * @param seq The sequence to fill.
     * @return true if a sequence was read or false if not.
     * @throw Exception IOExecption and Sequence related Exceptions.
     */
    virtual bool nextSequence(std::istream& input, Sequence& seq) const throw (Exception) = 0;

};

} //end of namespace bpp.

#endif	// _ISEQUENCESTREAM_H_

