//
// File OSequenceStream.h
// Author: Sylvain Gaillard
// Created: 19/08/2009
//

/*
Copyright or © or Copr. CNRS, (August 19, 2009)

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

#ifndef _OSEQUENCESTREAM_H_
#define _OSEQUENCESTREAM_H_

#include "IoSequenceStream.h"
#include "../Sequence.h"
#include "../Alphabet/Alphabet.h"
#include <Bpp/Exceptions.h>

namespace bpp
{

/**
 * @brief The OSequenceStream interface.
 *
 * Interface for streaming sequences output.
 *
 * @author Sylvain Gaillard
 */
class OSequenceStream: public virtual IOSequenceStream
{
	public:
		OSequenceStream() {}
		virtual ~OSequenceStream() {}

	public:
    /**
     * @brief Read sequence from stream.
     *
     * Read one sequence from a stream.
     *
     * @param output The stream where write.
     * @param seq The sequence to write.
     * @throw Exception IOExecption.
     */
    virtual void writeSequence(std::ostream& output, const Sequence& seq) const throw (Exception) = 0;

};

} //end of namespace bpp.

#endif	// _ISEQUENCESTREAM_H_

