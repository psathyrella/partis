//
// File: SvgGraphicDevice.h
// Created by: Julien Dutheil
// Created on: Mon Mar 10 2008
//

/*
Copyright or © or Copr. CNRS, (November 16, 2006)

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

#ifndef _SVGGRAPHICDEVICE_H_
#define _SVGGRAPHICDEVICE_H_

#include "../AbstractGraphicDevice.h"
#include "../ColorTools.h"

// From the STL:
#include <map>

namespace bpp
{

/**
 * @brief SVG plotting format.
 */
class SvgGraphicDevice:
  public AbstractGraphicDevice
{
  private:
    std::ostream& out_;
    std::map<int, std::vector<std::string>, std::greater<int> > layers_; //Layer display as in xfig
    bool inkscapeEnabled_;
    double minx_, maxx_, miny_, maxy_;
    std::map<short int, std::string> fontStyles_;
    std::map<short int, std::string> fontWeights_;

  public:
    SvgGraphicDevice(std::ostream& out, bool inkscapeEnabled = false):
      out_(out),
      layers_(),
      inkscapeEnabled_(inkscapeEnabled),
      minx_(0), maxx_(0), miny_(0), maxy_(0),
      fontStyles_(), fontWeights_()
    {
      fontStyles_[Font::STYLE_NORMAL] = "";
      fontStyles_[Font::STYLE_ITALIC] = "italic";
      fontWeights_[Font::WEIGHT_NORMAL] = "";
      fontWeights_[Font::WEIGHT_BOLD] = "bold";
    }

    virtual ~SvgGraphicDevice() {}

  public:
    void begin();
    void end();

    void drawLine(double x1, double y1, double x2, double y2);
    void drawRect(double x, double y, double width, double height, short fill = FILL_EMPTY);
    void drawCircle(double x, double y, double radius, short fill = FILL_EMPTY);
    void drawText(double x, double y, const std::string& text, short hpos = TEXT_HORIZONTAL_LEFT, short vpos = TEXT_VERTICAL_BOTTOM, double angle = 0) throw (UnvalidFlagException);
    void comment(const std::string& text)
    {
      layers_[getCurrentLayer()].push_back("<!-- " + text + " -->");
    }

  public:
    static std::string colorToText(const RGBColor& color)
    {
      return "rgb(" + TextTools::toString(color[0]) + "," + TextTools::toString(color[1]) + "," + TextTools::toString(color[2]) + ")";
    }

};

} // end of namespace bpp.

#endif //_SVGGRAPHICDEVICE_H_


