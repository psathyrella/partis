//
// File: TreeDrawingListener.h
// Created by: Julien Dutheil
// Created on: Tue May 18 10:33 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (2010)

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

#ifndef _TREEDRAWINGLISTENER_H_
#define _TREEDRAWINGLISTENER_H_

#include "TreeDrawing.h"

#include <Bpp/Clonable.h>

namespace bpp {

/**
 * @brief Interface allowing to capture drawing events.
 *
 * Implementing this interface allows you to easily and efficiently tune a plot,
 * and/or add elements to it.
 */
class TreeDrawingListener :
  public virtual Clonable
{
public:
#ifndef NO_VIRTUAL_COV
  TreeDrawingListener*
#else
  Clonable*
#endif
  clone() const = 0;

  virtual void beforeDrawTree(const DrawTreeEvent& event) = 0;
  virtual void afterDrawTree(const DrawTreeEvent& event) = 0;
  virtual void beforeDrawNode(const DrawNodeEvent& event) = 0;
  virtual void afterDrawNode(const DrawNodeEvent& event) = 0;
  virtual void beforeDrawBranch(const DrawBranchEvent& event) = 0;
  virtual void afterDrawBranch(const DrawBranchEvent& event) = 0;

  /**
   * @brief Tells if the listener is autonomous. If so, it
   * will never be hard-copied or deleted.
   */
  virtual bool isAutonomous() const = 0;

  virtual bool isEnabled() const = 0;
  virtual void enable(bool tf) = 0;
};


/**
 * @brief An empty implementation of the TreeDrawingListener interface.
 */
class TreeDrawingListenerAdapter :
  public virtual TreeDrawingListener
{
private:
  bool autonomous_;
  bool enabled_;

public:
  TreeDrawingListenerAdapter(bool autonomous) : autonomous_(autonomous), enabled_(true) {}

public:
  void beforeDrawTree(const DrawTreeEvent& event) {}
  void afterDrawTree(const DrawTreeEvent& event) {}
  void beforeDrawNode(const DrawNodeEvent& event) {}
  void afterDrawNode(const DrawNodeEvent& event) {}
  void beforeDrawBranch(const DrawBranchEvent& event) {}
  void afterDrawBranch(const DrawBranchEvent& event) {}
  bool isAutonomous() const { return autonomous_; }
  bool isEnabled() const { return enabled_; }
  void enable(bool tf) { enabled_ = tf; }
};



/**
 * @brief A TreeDrawingListener implementation that writes nodes id.
 */
class NodesIdTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
private:
  const TreeDrawingSettings* settings_;

public:
  NodesIdTreeDrawingListener(const TreeDrawingSettings* settings, bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous),
    settings_(settings)
  {}

  NodesIdTreeDrawingListener(const NodesIdTreeDrawingListener& lntdl) :
    TreeDrawingListenerAdapter(lntdl),
    settings_(lntdl.settings_) {}
  
  NodesIdTreeDrawingListener& operator=(const NodesIdTreeDrawingListener& lntdl)
  {
    TreeDrawingListenerAdapter::operator=(lntdl);
    settings_ = lntdl.settings_;
    return *this;
  }

  NodesIdTreeDrawingListener* clone() const { return new NodesIdTreeDrawingListener(*this); }

public :    
  void afterDrawNode(const DrawNodeEvent& event);

};


/**
 * @brief A TreeDrawingListener implementation that write leaf names.
 */
class LeafNamesTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
private:
  const TreeDrawingSettings* settings_;

public:
  LeafNamesTreeDrawingListener(const TreeDrawingSettings* settings, bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous),
    settings_(settings)
  {}

  LeafNamesTreeDrawingListener(const LeafNamesTreeDrawingListener& lntdl) :
    TreeDrawingListenerAdapter(lntdl),
    settings_(lntdl.settings_)
  {}
  
  LeafNamesTreeDrawingListener& operator=(const LeafNamesTreeDrawingListener& lntdl)
  {
    TreeDrawingListenerAdapter::operator=(lntdl);
    settings_   = lntdl.settings_;
    return *this;
  }

  LeafNamesTreeDrawingListener* clone() const { return new LeafNamesTreeDrawingListener(*this); }

public :    
  void afterDrawNode(const DrawNodeEvent& event);

};


/**
 * @brief A TreeDrawingListener implementation that write the branch lengths of inner nodes.
 *
 * Collapsed nodes are not labelled.
 *
 * This listener works with TreeDrawing classes, but is more efficient when used with a class that fires DrawINodeEvent events.
 */
class BranchLengthsTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
private:
  const TreeDrawingSettings* settings_;

public:
  BranchLengthsTreeDrawingListener(const TreeDrawingSettings* settings, bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous),
    settings_(settings)
  {}

  BranchLengthsTreeDrawingListener(const BranchLengthsTreeDrawingListener& lntdl) :
    TreeDrawingListenerAdapter(lntdl),
    settings_(lntdl.settings_)
  {}
  
  BranchLengthsTreeDrawingListener& operator=(const BranchLengthsTreeDrawingListener& lntdl)
  {
    TreeDrawingListenerAdapter::operator=(lntdl);
    settings_   = lntdl.settings_;
    return *this;
  }

  BranchLengthsTreeDrawingListener* clone() const { return new BranchLengthsTreeDrawingListener(*this); }

public :    
  void afterDrawBranch(const DrawBranchEvent& event);

};


/**
 * @brief A TreeDrawingListener implementation that write the bootstrap values of inner nodes.
 *
 * Collapsed nodes are not labelled.
 *
 * This listener works with TreeDrawing classes, but is more efficient when used with a class that fires DrawINodeEvent events.
 */
class BootstrapValuesTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
private:
  const TreeDrawingSettings* settings_;

public:
  BootstrapValuesTreeDrawingListener(const TreeDrawingSettings* settings, bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous),
    settings_(settings)
  {}

  BootstrapValuesTreeDrawingListener(const BootstrapValuesTreeDrawingListener& lntdl) :
    TreeDrawingListenerAdapter(lntdl),
    settings_(lntdl.settings_) {}
  
  BootstrapValuesTreeDrawingListener& operator=(const BootstrapValuesTreeDrawingListener& lntdl)
  {
    TreeDrawingListenerAdapter::operator=(lntdl);
    settings_   = lntdl.settings_;
    return *this;
  }

  BootstrapValuesTreeDrawingListener* clone() const { return new BootstrapValuesTreeDrawingListener(*this); }

public :    
  void afterDrawBranch(const DrawBranchEvent& event);

};


/**
 * @brief A TreeDrawingListener implementation that write the names of inner nodes.
 *
 * Collapsed nodes are not labelled.
 *
 * This listener works with TreeDrawing classes, but is more efficient when used with a class that fires DrawINodeEvent events.
 */
class LabelInnerNodesTreeDrawingListener :
  public TreeDrawingListenerAdapter
{

public:
  LabelInnerNodesTreeDrawingListener(bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous) {}

  LabelInnerNodesTreeDrawingListener* clone() const { return new LabelInnerNodesTreeDrawingListener(*this); }

public :    
  void afterDrawNode(const DrawNodeEvent& event);

};



/**
 * @brief A TreeDrawingListener implementation that label the collapsed nodes.
 *
 * This listener works with TreeDrawing classes, but is more efficient when used with a class that fires DrawINodeEvent events.
 */
class LabelCollapsedNodesTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
public:
  LabelCollapsedNodesTreeDrawingListener(bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous) {}

  LabelCollapsedNodesTreeDrawingListener* clone() const { return new LabelCollapsedNodesTreeDrawingListener(*this); }

public :    
  void afterDrawNode(const DrawNodeEvent& event);

};

} //end of namespace bpp

#endif //_TREEDRAWINGLISTENER_H_

