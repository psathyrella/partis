//
// File: PhylogramPlot.h
// Created by: Julien Dutheil
// Created on: Tue Oct 3 20:52 2006
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _PHYLOGRAMPLOT_H_
#define _PHYLOGRAMPLOT_H_

#include "AbstractDendrogramPlot.h"

namespace bpp
{

class PhylogramDrawBranchEvent :
  public DrawIBranchEvent
{
  private:
    double orientation_;

  public:
    PhylogramDrawBranchEvent(const TreeDrawing* source, GraphicDevice* gd, const INode* node, const Cursor& cursor, short orientation);

  public:
    Cursor getBranchCursor(double position) const;
    
};


  
/**
 * @brief Phylogram representation of trees.
 *
 * This representation is for trees with branch lengths.
 */
class PhylogramPlot:
  public AbstractDendrogramPlot
{
  private:
    double totalDepth_;
    double numberOfLeaves_;

  public:
    PhylogramPlot():
      AbstractDendrogramPlot(), totalDepth_(0), numberOfLeaves_(0)
    {}
    
    virtual ~PhylogramPlot() {}

    PhylogramPlot* clone() const { return new PhylogramPlot(*this); }

  public:
    std::string getName() const { return "Phylogram"; }
 
    void setTree(const Tree* tree = 0);
    
    double getWidth() const { return totalDepth_; }
    double getHeight() const { return numberOfLeaves_; }

    void treeHasChanged()
    {
      if (hasTree())
      {
        getTree_()->setVoidBranchLengths(0.);
        totalDepth_ = TreeTemplateTools::getHeight(*getTree_()->getRootNode());
        numberOfLeaves_ = static_cast<double>(getTree_()->getNumberOfLeaves());
      }
    }
 
  private:
    void drawDendrogram_(GraphicDevice& gDevice) const throw (Exception);
 
    void recursivePlot_(GraphicDevice& gDevice, INode& node, double x, double& y, double hDirection, double vDirection, unsigned int* tipCounter) const;

};

} //end of namespace bpp.

#endif //_PHYLOGRAMPLOT_H_

