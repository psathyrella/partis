//
// File: DefaultColorSet.h
// Created by: Julien Dutheil
// Created on: Mon Apr 14 2008
//

/*
Copyright or © or Copr. CNRS, (November 17, 2008)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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

#ifndef _DEFAULTCOLORSET_H_
#define _DEFAULTCOLORSET_H_

#include "ColorSet.h"
#include "ColorTools.h"

namespace bpp
{

/**
 * @brief Default color definitions.
 */
class DefaultColorSet:
  public AbstractColorSet
{
  public:
    DefaultColorSet()
    {
      colors_["black"]   = ColorTools::BLACK;
      colors_["red"]     = ColorTools::RED;
      colors_["green"]   = ColorTools::GREEN;
      colors_["blue"]    = ColorTools::BLUE;
      colors_["yellow"]  = ColorTools::YELLOW;
      colors_["magenta"] = ColorTools::MAGENTA;
      colors_["cyan"]    = ColorTools::CYAN;
      colors_["white"]   = ColorTools::WHITE;
    }
    virtual ~DefaultColorSet() {}

};

} // end of namespace bpp;

#endif //_DEFAULTCOLORSET_H_

