//
// File: FontManager.h
// Created by: Julien Dutheil
// Created on: Sat Mar 08 2008
//

/*
Copyright or © or Copr. Bio++ Developement Team, (November 17, 2004)

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

#ifndef _FONTMANAGER_H_
#define _FONTMANAGER_H_

#include "Font.h"

// From the STL:
#include <vector>

#include "../../Text/TextTools.h"

namespace bpp
{

/**
 * @brief Associate special fonts to a code.
 *
 * Instances of this interface are used in some vector format.
 */
template<class CodeType>
class FontManager
{
public:
  FontManager() {}
  virtual ~FontManager() {}

public:

  /**
   * @param font The font to look for.
   * @return The code associated to a given font.
   */
  virtual CodeType getCode(const Font& font) const throw (Exception) = 0;

  /**
   * @param code The code to look for.
   * @return The font associated to a given code.
   * @throw exception if the code is not valid.
   */
  virtual const Font& getFont(CodeType& code) const throw (Exception) = 0;

  /**
   * @return All valid codes.
   */
  virtual std::vector<CodeType> getCodes() const = 0;

  /**
   * @return All available fonts.
   */
  virtual std::vector<Font> getFonts() const = 0;

  /**
   * @return The total number of fonts available.
   */
  virtual size_t getNumberOfFonts() const = 0;
};



template<class CodeType>
class AbstractFontManager :
  public virtual FontManager<CodeType>
{
private:
  std::vector<Font> fonts_;
  std::vector<CodeType> codes_;

public:
  AbstractFontManager() : fonts_(), codes_() {}

public:
  CodeType getCode(const Font& font) const throw (Exception)
  {
    for (unsigned int i = 0; i < fonts_.size(); i++)
    {
      if (fonts_[i] == font)
      {
        return codes_[i];
      }
    }
    throw Exception("AbstractFontManager::getCode. Unknown font: " + font.toString());
  }

  const Font& getFont(int &code) const throw (Exception)
  {
    for (unsigned int i = 0; i < codes_.size(); i++)
    {
      if (codes_[i] == code)
      {
        return fonts_[i];
      }
    }
    throw Exception("AbstractFontManager::getColor. No font associated to this code: " + TextTools::toString(code));
  }
  std::vector<CodeType> getCodes() const { return codes_; }
  std::vector<Font> getFonts() const { return fonts_; }
  size_t getNumberOfFonts() const { return fonts_.size(); }

protected:
  void registerFont_(const Font& font, int code)
  {
    codes_.push_back(code);
    fonts_.push_back(font);
  }

};

} // end of namespace.

#endif //_FONTMANAGER_H_

