//
// File Point2DTools.h (from CoordsTools.h)
// Author : Sylvain Gaillard
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

#ifndef _POINT2DTOOLS_H_
#define _POINT2DTOOLS_H_

// From local
#include "Point2D.h"

namespace bpp
{

  /**
   * @brief Some functions to deal with Point2D.
   *
   * @author Sylvain Gaillard
   */

  class CoordsTools
  {
    public:
      /**
       * @brief Get the distance between two Coord objects.
       *
       * @param co1 A Point2D object.
       * @param co2 An other Point2D object.
       * @return the distance between the 2 points as a double
       */
      template<class T> static double getDistanceBetween(const Point2D<T>& co1, const Point2D<T>& co2)
      {
        T base, height;
        base = co1.getX() - co2.getX();
        height = co1.getY() - co2.getY();
        base = base * base;
        height = height * height;
        return sqrt((double)base + (double)height);
      }

  };

} //end of namespace bpp;

#endif // _POINT2DTOOLS_H_

