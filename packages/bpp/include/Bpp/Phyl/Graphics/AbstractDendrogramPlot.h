//
// File: AbstractDendrogramPlot.h
// Created by: Julien Dutheil
// Created on: Fri Jul 17 11:23 2009
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

This software is a computer program whose purpose is to provide
graphic components to develop bioinformatics applications.

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

#ifndef _ABSTRACTDENDROGRAMPLOT_H_
#define _ABSTRACTDENDROGRAMPLOT_H_

#include "AbstractTreeDrawing.h"

namespace bpp
{

/**
 * @brief Basic implementation of dendrogram plots.
 *
 * Dendrograms are oriented plots, with all the leaves on one side of the plot, and the root node at the opposite side.
 * This implementation offers to option for ploting form left to right or right to left. This will affect the direction
 * of plot annotations. The drawing can always be transformed using the regular translation/rotation operation on the
 * GraphicDevice.
 */
class AbstractDendrogramPlot:
  public AbstractTreeDrawing
{
  private:
    short horOrientation_;
    short verOrientation_;

  public:
    AbstractDendrogramPlot():
      AbstractTreeDrawing(), horOrientation_(ORIENTATION_LEFT_TO_RIGHT), verOrientation_(ORIENTATION_TOP_TO_BOTTOM)
    {}

  public:
    void setHorizontalOrientation(short orientation) { horOrientation_ = orientation; }
    void setVerticalOrientation(short orientation) { verOrientation_ = orientation; }

    short getHorizontalOrientation() const { return horOrientation_; }
    short getVerticalOrientation() const { return verOrientation_; }

    void plot(GraphicDevice& gDevice) const throw (Exception);

  protected:
    virtual void drawDendrogram_(GraphicDevice& gDevice) const throw (Exception) = 0;
   
  public:
    static short ORIENTATION_LEFT_TO_RIGHT;
    static short ORIENTATION_RIGHT_TO_LEFT;
    static short ORIENTATION_TOP_TO_BOTTOM;
    static short ORIENTATION_BOTTOM_TO_TOP;

};

} //end of namespace bpp;

#endif //_ABSTRACTDENDROGRAMPLOT_H_


