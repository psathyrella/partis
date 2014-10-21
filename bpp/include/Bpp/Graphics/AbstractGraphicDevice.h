//
// File: AbstractGraphicDevice.h
// Created by: Julien Dutheil
// Created on: Fri Jul 24 2009
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

#ifndef _ABSTRACTGRAPHICDEVICE_H_
#define _ABSTRACTGRAPHICDEVICE_H_

#include "GraphicDevice.h"

namespace bpp
{

/**
 * @brief Partial implementation of the GraphicDevice interface.
 *
 * Implement this interface to support new formats.
 */
class AbstractGraphicDevice:
  public virtual GraphicDevice
{
  private:
    double xUnit_;
    double yUnit_;
    RGBColor fgColor_;
    RGBColor bgColor_;
    Font font_;
    unsigned int pointSize_;
    short lineType_;
    int currentLayer_;

  public:
    AbstractGraphicDevice(): xUnit_(1.), yUnit_(1.),
        fgColor_(0, 0, 0), bgColor_(0, 0, 0), font_(), pointSize_(1), lineType_(LINE_SOLID), currentLayer_(-1) 
    {}

    virtual ~AbstractGraphicDevice() {}

  public:
    void setXUnit(double xu) { xUnit_ = xu; }
    void setYUnit(double yu) { yUnit_ = yu; }
    double getXUnit() const { return xUnit_; }
    double getYUnit() const { return yUnit_; }

    void setCurrentForegroundColor(const RGBColor& color) { fgColor_ = color; }
    void setCurrentBackgroundColor(const RGBColor& color) { bgColor_ = color; }
    void setCurrentFont(const Font& font) { font_ = font; }
    void setCurrentPointSize(unsigned int size) { pointSize_ = size; }
    void setCurrentLineType(short type) throw (Exception)
    { 
      if       (type == LINE_SOLID) lineType_ = type;
      else if (type == LINE_DASHED) lineType_ = type;
      else if (type == LINE_DOTTED) lineType_ = type;
      else throw Exception("AbstractGraphicDevice::setCurrentLineType. Unknown line type: " + TextTools::toString(type));
    }
    void setCurrentLayer(int layerIndex) { currentLayer_ = layerIndex; }

    const RGBColor& getCurrentForegroundColor() const { return fgColor_; }
    const RGBColor& getCurrentBackgroundColor() const { return bgColor_; }
    const Font& getCurrentFont() const { return font_; }
    unsigned int getCurrentPointSize() const { return pointSize_; }
    short getCurrentLineType() const { return lineType_; }
    int getCurrentLayer() const { return currentLayer_; }
 

  protected:
    double x_(double x) const { return x * xUnit_; }
    double y_(double y) const { return y * yUnit_; }

    double revx_(double x) const { return x / xUnit_; }
    double revy_(double y) const { return y / yUnit_; }

};

} //end of namespace bpp;

#endif //_ABSTRACTGRAPHICDEVICE_H_

