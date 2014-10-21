//
// File: AbstractParameterAliasable.h
// Created by: Julien Dutheil
// Created on: Thu May 14 17:08 2009
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 19, 2004)

  This software is a computer program whose purpose is to provide classes
  for numerical calculus.

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

#ifndef _ABSTRACTPARAMETERALIASABLE_H_
#define _ABSTRACTPARAMETERALIASABLE_H_

#include "AbstractParametrizable.h"
#include "ParameterAliasable.h"

//From the STL:
#include <map>

namespace bpp
{

  /**
   * @brief Inner listener class used by AbstractParameterAliasable.
   */
  class AliasParameterListener:
    public ParameterListener
  {
  private:
    std::string id_;
    size_t alias_;
    ParameterList *pl_;
    std::string name_;
    std::string from_;

  public:
    AliasParameterListener(const std::string& id, size_t alias, ParameterList* pl, const std::string& from):
      id_(id),
      alias_(alias),
      pl_(pl),
      name_(),
      from_(from)
    {
      //This allow us to check if the parameter position have changed at some point...
      name_ = (*pl_)[alias].getName();
    }

    AliasParameterListener(const AliasParameterListener& apl):
    id_(apl.id_),
    alias_(apl.alias_),
    pl_(apl.pl_),
    name_(apl.name_),
    from_(apl.from_)
    {}

    AliasParameterListener& operator=(const AliasParameterListener& apl)
    {
      id_    = apl.id_;
      alias_ = apl.alias_;
      pl_    = apl.pl_;
      name_  = apl.name_;
      from_  = apl.from_;
      return *this;
    }

    AliasParameterListener* clone() const { return new AliasParameterListener(*this); }

  public:
    const std::string& getId() const { return id_; }

    const std::string& getFrom() const { return from_; }
    
    void setParameterList(ParameterList* pl) { pl_ = pl; }

    void parameterNameChanged(ParameterEvent& event) throw (Exception) {}
    
    void parameterValueChanged(ParameterEvent& event) throw (Exception)
    {
      Parameter* p = &(*pl_)[alias_];
      if (p->getName() != name_)
        throw Exception("AbstractParameterAliasable::AliasParameterListener::parameterValueChanged. Error, aliased parameter have change, maybe because it was renamed, or a parameter was removed?");
      p->setValue(event.getParameter()->getValue());
    }

    const std::string& getName() const { return name_; }

    void rename(const std::string& name) { name_ = name; }

    const std::string& getAlias() const { return (*pl_)[alias_].getName(); }
      
  };

  /**
   * @brief A partial implementation of the Parametrizable interface.
   *
   * Parameters are stored in a protected ParameterList object.
   *
   * The abstract fireParameterChanged() method is provided so that the derived class
   * know when a parameter has changed, and can be updated.
   * All methods call the corresponding method in ParameterList and then call the
   * fireParameterChanged() method.
   */
  class AbstractParameterAliasable:
    public AbstractParametrizable,
    public virtual ParameterAliasable
  {
  private:

    mutable ParameterList independentParameters_;

    /**
     * Contains all parameter listeners for maintening alias relationships.
     * The registry will be updated appropriately upon cloning and deleting.
     */
    std::map<std::string, AliasParameterListener *> aliasListenersRegister_;
  
  public:
    AbstractParameterAliasable(const std::string& prefix) :
      AbstractParametrizable(prefix),
      independentParameters_(),
      aliasListenersRegister_()
    {}

    AbstractParameterAliasable(const AbstractParameterAliasable& ap);
    
    AbstractParameterAliasable& operator=(const AbstractParameterAliasable& ap);

    virtual ~AbstractParameterAliasable();

  public:
    void setNamespace(const std::string& prefix);
 
    const ParameterList& getIndependentParameters() const { return independentParameters_; }
    
    size_t getNumberOfIndependentParameters() const { return independentParameters_.size(); }

    void aliasParameters(const std::string& p1, const std::string& p2) throw (ParameterNotFoundException, Exception);

    void unaliasParameters(const std::string& p1, const std::string& p2) throw (ParameterNotFoundException, Exception);

    /**
     * @brief alias the parameters following the links described in a
     * map, and update the object accordingly. Cycles in aliasing are
     * detected and forbidden.
     *
     * @param unparsedParams the map of the links : <A,B> matches for A->B aliasing.
     * @param verbose verbosity
     *
     **/
    
    void aliasParameters(std::map<std::string, std::string>& unparsedParams, bool verbose);

    /**
     * @brief Return the list of the names of the parameters that are
     * aliased (directly or not) to one of the parameters of the list.
     *
     */
    
    ParameterList getAliasedParameters(const ParameterList& pl) const;
    
    
    /**
     * @return The list of names of the parameters that are aliased with a given parameter.
     * The implementation is recursive, which means that in the case of A->B->C, getalias(C) will return both A and B.
     * @param name The name of the parameter to look for.
     */
    
    std::vector<std::string> getAlias(const std::string& name) const;

    /**
     * @return The name of the parameter from which a given parameter is aliased.
     * @param name The name of the parameter to look for.
     */

    std::string getFrom(const std::string& name) const;
      
    void fireParameterChanged(const ParameterList& parameters)
    {
      independentParameters_.matchParametersValues(getParameters());
    }

  protected:
    void addParameter_(Parameter* parameter)
    {
      AbstractParametrizable::addParameter_(parameter);
      independentParameters_.addParameter(parameter->clone());
    }

    void addParameters_(const ParameterList& parameters)
    {
      AbstractParametrizable::addParameters_(parameters);
      independentParameters_.addParameters(parameters);
    }

    void deleteParameter_(size_t index) throw (IndexOutOfBoundsException)
    {
      std::string name = getParameter_(index).getName();
      AbstractParametrizable::deleteParameter_(index);
      if (independentParameters_.hasParameter(name))
        independentParameters_.deleteParameter(name);
    }

    void deleteParameter_(std::string& name) throw (IndexOutOfBoundsException)
    {
      AbstractParametrizable::deleteParameter_(name);
      if (independentParameters_.hasParameter(name))
        independentParameters_.deleteParameter(name);
    }

    void deleteParameters_(const std::vector<std::string>& names)
    {
      std::string x;
      for (size_t i=0; i<names.size(); i++)
      {
        x=names[i];
        deleteParameter_(x);
      }
    }

    void resetParameters_()
    {
      AbstractParametrizable::resetParameters_();
      independentParameters_.reset();
    }

  };

} //end of namespace bpp.

#endif //_ABSTRACTPARAMETERALIASABLE_H_

