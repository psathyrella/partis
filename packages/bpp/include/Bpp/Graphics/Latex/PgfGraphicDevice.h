//
// File: PgfGraphicDevice.h
// Created by: Julien Dutheil
// Created on: Thu Jun 19 2008
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

#ifndef _PGFGRAPHICDEVICE_H_
#define _PGFGRAPHICDEVICE_H_

#include "../AbstractGraphicDevice.h"
#include "../ColorTools.h"

// From the STL:
#include <map>
#include <vector>
#include <string>

namespace bpp
{

/**
 * @brief LaTeX Portable Graphic Format (Pgf) plotting format.
 */
class PgfGraphicDevice:
  public AbstractGraphicDevice
{
  private:
    std::ostream& out_;
    std::string fgColorStr_;
    std::string bgColorStr_;
    std::vector<std::string> content_;
    std::vector<int> layers_;
    std::map<const RGBColor, std::string> colorIndex_;
    unsigned int colorCount_;
    bool useLayers_;
    bool contentStarted_;
    mutable std::map<short int, std::string> fontShapes_;
    mutable std::map<short int, std::string> fontSeries_;

  public:
    /**
     * @brief Build a new Pgf device object.
     *
     * Coordinates in Pgf are the same as in LaTeX, so it can be any of cm, mm, in, pt, px, etc.
     * For compatibility with other devices, the constructor takes as input the scale of the drawing, as cm per points.
     * All coordinates and widths will be multiplied by the factor in the output file.
     *
     * @param out The output stream.
     * @param unit The unit length.
     */
    PgfGraphicDevice(std::ostream& out, double unit);

    virtual ~PgfGraphicDevice() {}

  public:
    void begin();
    void end();

    void setCurrentForegroundColor(const RGBColor& color);
    void setCurrentBackgroundColor(const RGBColor& color);
    void setCurrentFont(const Font& font);
    void setCurrentPointSize(unsigned int size);
    void setCurrentLineType(short type) throw (Exception);
    void setCurrentLayer(int layerIndex);
   
    void drawLine(double x1, double y1, double x2, double y2);
    void drawRect(double x, double y, double width, double height, short fill = FILL_EMPTY);
    void drawCircle(double x, double y, double radius, short fill = FILL_EMPTY);
    void drawText(double x, double y, const std::string & text, short hpos = TEXT_HORIZONTAL_LEFT, short vpos = TEXT_VERTICAL_BOTTOM, double angle = 0) throw (UnvalidFlagException);
    void comment(const std::string& text)
    {
      content_.push_back("%" + text);
    }

  public:
    const std::string& getFontShape(const Font& font);

};

} // end of namespace bpp.

#endif //_PGFGRAPHICDEVICE_H_


