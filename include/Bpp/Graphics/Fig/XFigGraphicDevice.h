//
// File: XFigGraphicDevice.h
// Created by: Julien Dutheil
// Created on: Mon Mar 03 2008
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#ifndef _XFIGGRAPHICDEVICE_H_
#define _XFIGGRAPHICDEVICE_H_

#include "../AbstractGraphicDevice.h"
#include "../ColorManager.h"
#include "XFigLaTeXFontManager.h"
#include "XFigPostscriptFontManager.h"

// From the STL:
#include <map>
#include <iostream>

namespace bpp
{

/**
 * @brief XFig plotting format.
 */
class XFigGraphicDevice:
  public AbstractGraphicDevice
{
  private:
    std::ostream& out_;
    std::vector<std::string> content_;
    XFigColorManager colorManager_;
    XFigLaTeXFontManager latexFontManager_;
    XFigPostscriptFontManager postscriptFontManager_;
    unsigned int fgColorCode_;
    unsigned int bgColorCode_;
    int fontCode_; 
    unsigned int fontSize_;
    unsigned int fontFlag_;
    unsigned int lineTypeCode_;

  public:
    XFigGraphicDevice(std::ostream& out):
      out_(out),
      content_(),
      colorManager_(),
      latexFontManager_(),
      postscriptFontManager_(),
      fgColorCode_(0),
      bgColorCode_(0),
      fontCode_(-1),
      fontSize_(12),
      fontFlag_(FONTFLAG_POSTSCRIPT),
      lineTypeCode_(LINE_SOLID)
    {
      setCurrentLayer(0);
    }

    virtual ~XFigGraphicDevice() {}

  public:
    void begin();
    void end();

    void setCurrentForegroundColor(const RGBColor& color);
    void setCurrentBackgroundColor(const RGBColor& color);
    void setCurrentFont(const Font& font);
    void setCurrentLineType(short type) throw (Exception)
    { 
      if(type == LINE_SOLID) lineTypeCode_ = 0;
      else if(type == LINE_DASHED) lineTypeCode_ = 1;
      else if(type == LINE_DOTTED) lineTypeCode_ = 2;
      else throw Exception("XFigGraphicDevice::setCurrentLineType. Unknown line type: " + TextTools::toString(type));
      AbstractGraphicDevice::setCurrentLineType(type);
    }
    void drawLine(double x1, double y1, double x2, double y2);
    void drawRect(double x, double y, double width, double height, short fill = FILL_EMPTY);
    void drawCircle(double x, double y, double radius, short fill = FILL_EMPTY);
    void drawText(double x, double y, const std::string& text, short hpos = TEXT_HORIZONTAL_LEFT, short vpos = TEXT_VERTICAL_BOTTOM, double angle = 0) throw (UnvalidFlagException);
    void comment(const std::string& text)
    {
      content_.push_back("#" + text);
    }

    //Specific:
    void setFontFlag(unsigned int flag) { fontFlag_ = flag; }

  protected:
    int getFillCode(short fill);

  public:
    static const unsigned int FONTFLAG_LATEX;
    static const unsigned int FONTFLAG_POSTSCRIPT;

};

} // end of namespace bpp;

#endif //_XFIGGRAPHICDEVICE_H_


