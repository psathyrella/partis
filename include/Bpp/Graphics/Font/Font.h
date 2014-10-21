//
// File: Font.h
// Created by: Julien Dutheil
// Created on: Mon Mar 03 2008
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

#ifndef _FONT_H_
#define _FONT_H_

#include "../../Clonable.h"
#include "../../Text/TextTools.h"

//From the STL:
#include <map>
#include <string>

namespace bpp
{

/**
 * @brief Data structure for fonts.
 */
class Font:
  public virtual Clonable
{
  private:
    std::string family_;
    short int style_;
    short int weight_;
    unsigned int size_;
    mutable std::map<short int, std::string> styleDesc_;
    mutable std::map<short int, std::string> weightDesc_;

  public:
    Font(const std::string& family = "Default", short int style = STYLE_NORMAL, short int weight = WEIGHT_NORMAL, unsigned int size = 12):
      family_(family), style_(style), weight_(weight), size_(size), styleDesc_(), weightDesc_()
    {
      init_();
    }

    virtual ~Font() {}

#ifdef NO_VIRTUAL_COV
    Clonable*
#else
    Font*
#endif
    clone() const { return new Font(*this); }

  private:
    void init_();

  public:
    bool operator==(const Font& font) const
    {
      return family_ == font.family_
          &&  style_ == font.style_
          && weight_ == font.weight_;
    }

    /**
     * @return The family component of this font.
     */
    const std::string& getFamily() const { return family_; }

    /**
     * @return The style component of this font.
     */
    const short int getStyle() const { return style_; }

    /**
     * @brief Alias function for getStyle.
     * @return The shape component of this font.
     */
    const short int getShape() const { return style_; }

    /**
     * @return The weight component of this font.
     */
    const short int getWeight() const { return weight_; }

    /**
     * @brief Alias function for getWeight
     * @return The series component of this font.
     */
    const short int getSeries() const { return weight_; }

    /**
     * @return The size component of this font.
     */
    const unsigned int& getSize() const { return size_; }

    /**
     * @param family The family component of this font.
     */
    void setFamily(const std::string& family) { family_ = family; }

    /**
     * @param style The style component of this font.
     */
    void setStyle(short int style) { style_ = style; }

    /**
     * @brief Alias function for setStyle.
     * @param shape The shape component of this font.
     */
    void setShape(short int shape) { style_ = shape; }


    /**
     * @param weight The weight component of this font.
     */
    void setWeight(short int weight) { weight_ = weight; }

    /**
     * @brief Alias function for setWeight.
     * @param series The series component of this font.
     */
    void setSeries(short int series) { weight_ = series; }

    /**
     * @param size The size component of this font.
     */
    void setSize(unsigned int size) { size_ = size; }

    /**
     * @return A text description of this font (family type size).
     */
    std::string toString() const
    {
      return family_ + "/" + styleDesc_[style_] + "/" + weightDesc_[weight_] + "/" + TextTools::toString(size_);
    }

  public:
    static const short int STYLE_NORMAL;
    static const short int STYLE_ITALIC;
    
    static const short int WEIGHT_NORMAL;
    static const short int WEIGHT_BOLD;
};

} // end of namespace bpp.

#endif //_FONT_H_


