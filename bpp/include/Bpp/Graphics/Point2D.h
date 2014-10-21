//
// File Point2D.h (from file Coord.h)
// Author : Sylvain Gaillard
//          Julien Dutheil
//

/*
   Copyright or © or Copr. CNRS, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for population genetics analysis.

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

#ifndef _POINT2D_H_
#define _POINT2D_H_

#include "../Clonable.h"

namespace bpp
{

  /**
   * @brief The Point2D class.
   *
   * This is a simple class designed to store the coordinates of a point.
   * The type of the two coordinates is defined as a template.
   *
   * @author Sylvain Gaillard, Julien Dutheil
   */
  template<class T> class Point2D:
    public virtual Clonable
  {
    private:
      T x_;
      T y_;

    public: // Constructors and destructor :

      /**
       * @brief Build a new Point2D from two values.
       * 
       * The two values are set to 0 if no parametre is given to the constructor.
       *
       * @param x The longitude or abscissa.
       * @param y The latitude or ordinate.
       */
      Point2D<T>(const T x = 0, const T y = 0): x_(x), y_(y) {}

      /**
       * @brief Destroy the Point2D object.
       */
      virtual ~Point2D() {}

    public: // Methodes

      /**
       * @brief Implement the Clonable interface.
       */
      Point2D<T>* clone() const { return new Point2D(*this); }

      /**
       * @brief Set the two values.
       */
      void setCoord(const T x, const T y);

      /**
       * @brief Set only the longitude.
       */
      void setX(const T x) { x_ = x; }

      /**
       * @brief Set only the latitude.
       */
      void setY(const T y) { y_ = y; }

      /**
       * @brief Get the longitude.
       */
      const T& getX() const { return x_; }

      /**
       * @brief Get the latitude.
       */
      const T& getY() const { return y_; }

      /**
       * @brief Compares two Point2D objets.
       *
       * Return true if the coordinates of the 2 Point2D are equals.
       */
      bool hasSameCoordsAs(const Point2D<T>& coord) const
      {
        return x_ == coord.x_ && y_ == coord.y_;
      }

      /**
       * @brief The == operator.
       *
       * Return true if the coordinates of the 2 Point2Ds are equals.
       * Does the same as the asSameCoords() methode.
       */
      virtual bool operator== (const Point2D<T>& coord) const
      {
        return hasSameCoordsAs(coord);
      }

      /**
       * @brief The != operator.
       */
      virtual bool operator!= (const Point2D<T>& coord) const
      {
        return !hasSameCoordsAs(coord);
      }

  };

} //end of namespace bpp;

#endif // _POINT2D_H_

