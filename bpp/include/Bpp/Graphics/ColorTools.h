//
// File: ColorTools.h
// Created by: Julien Dutheil
// Created on: Thu Mar 16 2006
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#ifndef _COLORTOOLS_H_
#define _COLORTOOLS_H_

#include "../Text/TextTools.h"
#include "../Exceptions.h"
#include "RgbColor.h"

namespace bpp
{

/**
 * @brief Provide tools to deal with color objects.
 */
class ColorTools
{
  public:
    ColorTools() {}
    virtual ~ColorTools() {}

  public:
    /**
     * @brief Create a set of colors according to a gradient defined by two extrema.
     *
     * @param n Number of colors to output.
     * @param low First color in gradient.
     * @param high Last color in gradient.
     * @return A set of ordered colors defining a gradient.
     */
    static std::vector<RGBColor> gradient(unsigned int n, const RGBColor & low, const RGBColor & high);
 
    /**
     * @brief Create a set of colors according to a gradient defined by two extrema and a midpoint.
     *
     * @param n Number of colors to output.
     * @param low First color in gradient.
     * @param mid Midpoint color.
     * @param high Last color in gradient.
     * @return A set of ordered colors defining a gradient.
     */
    static std::vector<RGBColor> gradient(unsigned int n, const RGBColor & low, const RGBColor & mid, const RGBColor & high);
    
    /**
     * @return A gray color.
     * @param level Gray intensity ([0,1]).
     */
    static RGBColor gray(double level)
    { 
      unsigned int i = (unsigned int)round(255*level);
      return RGBColor(i, i, i);
    }

    /**
     * @brief Get a RGBColor from a cyan-magenta-yellow-key description.
     *
     * The following formula are used for the transformation:
     * @f[
     * \begin{array}{rcl}
     * r &=& 255 * (1 - c)(1 - k)\\
     * g &=& 255 * (1 - m)(1 - k)\\
     * b &=& 255 * (1 - y)(1 - k)
     * \end{array}
     * @f]
     *
     * @param c Cyan proportion.
     * @param m Magenta proportion.
     * @param y Yellow proportion.
     * @param k Black proportion.
     * @return A RGBColor object.
     */
    static RGBColor cmyk2rgb(double c, double m, double y, double k)
    {
      unsigned int r = static_cast<unsigned int>(round(255 * (1. - c) * (1. - k)));
      unsigned int g = static_cast<unsigned int>(round(255 * (1. - m) * (1. - k)));
      unsigned int b = static_cast<unsigned int>(round(255 * (1. - y) * (1. - k)));
      return RGBColor(r, g, b);
    }
 
  public:
    static const RGBColor RED;
    static const RGBColor GREEN;
    static const RGBColor BLUE;
    static const RGBColor BLACK;
    static const RGBColor WHITE;
    static const RGBColor YELLOW;
    static const RGBColor CYAN;
    static const RGBColor MAGENTA;
    static const RGBColor ORANGE;
    
};

} // end of namespace bpp.

#endif //_COLORTOOLS_H_

