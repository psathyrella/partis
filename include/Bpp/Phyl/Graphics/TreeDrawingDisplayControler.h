//
// File: TreeDrawingDisplayControler.h
// Created by: Julien Dutheil
// Created on: Tue May 18 12:37 2010
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

#ifndef _TREEDRAWINGDISPLAYCONTROLER_H_
#define _TREEDRAWINGDISPLAYCONTROLER_H_

#include "TreeDrawingListener.h"

//From the STL:
#include <string>
#include <vector>
#include <map>
#include <algorithm>

namespace bpp {

/**
 * @brief Easy tune of tree drawings display.
 *
 * This class maintains a set of autonomous TreeDrawing listeners that
 * are used for annotating a tree drawing.
 *
 * @author Julien Dutheil
 */
class TreeDrawingDisplayControler
{
  private:
    std::map<std::string, TreeDrawingListener*> listeners_;
    std::vector<TreeDrawing*> registeredTreeDrawings_;

  public:
    TreeDrawingDisplayControler() :
      listeners_(), registeredTreeDrawings_()
    {}

  private:
    TreeDrawingDisplayControler(const TreeDrawingDisplayControler& tddc) :
      listeners_(), registeredTreeDrawings_(tddc.registeredTreeDrawings_)
    {
      for (std::map<std::string, TreeDrawingListener*>::const_iterator it = tddc.listeners_.begin();
          it != tddc.listeners_.end(); ++it)
      {
        listeners_[it->first] = dynamic_cast<TreeDrawingListener*>(it->second->clone());
      }
    }
    TreeDrawingDisplayControler& operator=(const TreeDrawingDisplayControler& tddc)
    {
      listeners_.clear();
      registeredTreeDrawings_ = tddc.registeredTreeDrawings_;
      for (std::map<std::string, TreeDrawingListener*>::const_iterator it = tddc.listeners_.begin();
          it != tddc.listeners_.end(); ++it)
      {
        listeners_[it->first] = dynamic_cast<TreeDrawingListener*>(it->second->clone());
      }
      return *this;
    }

  public:
    virtual ~TreeDrawingDisplayControler();

  public:
    /**
     * @brief Add a listener to the controler. The controler then owns the object, and will
     * copy or delete it when needed.
     */
    void addListener(const std::string& propertyName, TreeDrawingListener* listener) throw (Exception);

    bool hasListenerFor(const std::string& propertyName) const
    {
      return listeners_.find(propertyName) != listeners_.end();
    }

    void enableListener(const std::string& propertyName, bool tf) throw (Exception)
    {
      if (!hasListenerFor(propertyName))
        throw Exception("TreeDrawingDisplayControler::enableListener. No listener is registered for property " + propertyName + ".");
      listeners_[propertyName]->enable(tf);
    }

    bool isListenerEnabled(const std::string& propertyName) const throw (Exception)
    {
      if (!hasListenerFor(propertyName))
        throw Exception("TreeDrawingDisplayControler::enableListener. No listener is registered for property " + propertyName + ".");
      return listeners_.find(propertyName)->second->isEnabled();
    }

    void registerTreeDrawing(TreeDrawing* td) throw (Exception)
    {
      if (std::find(registeredTreeDrawings_.begin(), registeredTreeDrawings_.end(), td) != registeredTreeDrawings_.end())
        throw Exception("TreeDrawingDisplayControler::registerTreeDrawing. TreeDrawing is already associated to this controler.");
      for (std::map<std::string, TreeDrawingListener*>::iterator it = listeners_.begin();
          it != listeners_.end(); ++it)
        td->addTreeDrawingListener(it->second);
      registeredTreeDrawings_.push_back(td);
    }

};



/**
 * @brief Easy tune of tree drawings display, a basic implementation:
 *
 * This class maintains several "standard" drawing listener for:
 * - Plotting node id,
 * - Plotting leaves names,
 * - Plotting branch lengths,
 * - Plotting plotting bootstrap values.
 *
 * This controler takes as an argument a TreeDrawingSettings object that is used by
 * all listeners that require one.
 */
class BasicTreeDrawingDisplayControler :
  public TreeDrawingDisplayControler
{
  public:
    static const std::string PROPERTY_NODE_IDS;
    static const std::string PROPERTY_LEAF_NAMES;
    static const std::string PROPERTY_BRANCH_LENGTHS;
    static const std::string PROPERTY_BOOTSTRAP_VALUES;

  private:
    const TreeDrawingSettings* settings_;

  public:
    BasicTreeDrawingDisplayControler(const TreeDrawingSettings* settings) throw (NullPointerException) :
      settings_(settings)
    {
      if (!settings)
        throw NullPointerException("BasicTreeDrawingDisplayControler::constructor. Trying to use NULL settings.");
      addListener(PROPERTY_NODE_IDS        , reinterpret_cast<TreeDrawingListener*>(new NodesIdTreeDrawingListener        (settings_, true)));
      addListener(PROPERTY_LEAF_NAMES      , reinterpret_cast<TreeDrawingListener*>(new LeafNamesTreeDrawingListener      (settings_, true)));
      addListener(PROPERTY_BRANCH_LENGTHS  , reinterpret_cast<TreeDrawingListener*>(new BranchLengthsTreeDrawingListener  (settings_, true)));
      addListener(PROPERTY_BOOTSTRAP_VALUES, reinterpret_cast<TreeDrawingListener*>(new BootstrapValuesTreeDrawingListener(settings_, true)));
    }

  private:
    BasicTreeDrawingDisplayControler(const BasicTreeDrawingDisplayControler&) : settings_(0) {}
    BasicTreeDrawingDisplayControler& operator=(const BasicTreeDrawingDisplayControler&) { return *this; }

};

} //end of namespace bpp.

#endif //_TREEDRAWINGDISPLAYCONTROLER_H_

