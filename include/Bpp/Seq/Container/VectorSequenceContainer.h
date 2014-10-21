//
// File VectorSequenceContainer.h
// Created by: Guillaume Deuchst
//             Julien Dutheil
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#ifndef _VECTORSEQUENCECONTAINER_H_
#define _VECTORSEQUENCECONTAINER_H_

#include "../Alphabet/Alphabet.h"
#include "../Sequence.h"
#include "AbstractSequenceContainer.h"
#include <Bpp/Exceptions.h>

// From the STL:
#include <algorithm>
#include <vector>

namespace bpp
{

/**
 * @brief The VectorSequenceContainer class.
 *
 * This is the simplest implementation of the OrderedSequenceContainer interface.
 * Sequences are stored in a std::vector of pointers.
 * The container is responsible for the creation and destruction of the sequence
 * objects it contains.
 */
class VectorSequenceContainer:
  public AbstractSequenceContainer
{
  private:

    /**
     * @brief A std::vector of pointers toward the sequences stored in the container.
     */
    mutable std::vector<Sequence*> sequences_;
        
  public:
    
    /**
     * @brief Build a new container from a std::vector of pointers toward sequence objects.
     *
     * The addSequence() method is called uppon each Sequence object, hence each sequence is
     * <i>copied</i> into the container.
     *
     * @param vs    The std::vector of pointers toward sequence objects.
     * @param alpha The alphabet to all sequences.
     * @throw AlphabetMismatchException if one sequence does not match the specified alphabet.
     */
    VectorSequenceContainer(
      const std::vector<const Sequence*>& vs, const Alphabet* alpha)
      throw (AlphabetMismatchException);
  
    /**
     * @brief Build an empty container that will contain sequences of a particular alphabet.
     *
     * @param alpha The alphabet of the container.
     */
    VectorSequenceContainer(const Alphabet* alpha): AbstractSequenceContainer(alpha), sequences_() {}
    
    /**
     * @name Copy contructors:
     *
     * @{
     */
    
    /**
     * @brief Copy from a VectorSequenceContainer.
     *
     * @param vsc The VectorSequenceContainer to copy into this container.
     */
    VectorSequenceContainer(const VectorSequenceContainer& vsc);
    
    /**
     * @brief Copy from an OrderedSequenceContainer.
     *
     * @param osc The OrderedSequenceContainer to copy into this container.
     */
    VectorSequenceContainer(const OrderedSequenceContainer& osc);

    /**
     * @brief Copy from a SequenceContainer.
     *
     * @param osc The SequenceContainer to copy into this container.
     */
    VectorSequenceContainer(const SequenceContainer& osc);

    /** @} */

    /**
     * @brief Assign from a VectorSequenceContainer.
     *
     * @param vsc The VectorSequenceContainer to copy into this container.
     */
    VectorSequenceContainer& operator=(const VectorSequenceContainer& vsc);

    /**
     * @brief Copy from an OrderedSequenceContainer.
     *
     * @param osc The OrderedSequenceContainer to copy into this container.
     */
    VectorSequenceContainer& operator=(const OrderedSequenceContainer& osc);
  
    /**
     * @brief Copy from a SequenceContainer.
     *
     * @param osc The SequenceContainer to copy into this container.
     */
    VectorSequenceContainer& operator=(const SequenceContainer& osc);

    /**
     * @brief Container destructor: delete all sequences in the container.
     */
    virtual ~VectorSequenceContainer() { clear(); }

  public:
    
    /**
     * @name The Clonable interface.
     *
     * @{
     */
    Clonable* clone() const { return new VectorSequenceContainer(*this); }
    /** @} */

    /**
     * @name The SequenceContainer interface.
     *
     * @{
     */
    bool hasSequence(const std::string& name) const;
  
    const Sequence& getSequence(const std::string& name) const throw (SequenceNotFoundException);

    void setSequence(const std::string& name, const Sequence& sequence, bool checkName = true) throw (Exception)
    {
      setSequence(getSequencePosition(name), sequence, checkName);
    }

    Sequence* removeSequence(const std::string& name) throw (SequenceNotFoundException)
    {
      return removeSequence(getSequencePosition(name));
    }

    void deleteSequence(const std::string& name) throw (SequenceNotFoundException)
    {
      deleteSequence(getSequencePosition(name));
    }
    
    size_t getNumberOfSequences() const { return sequences_.size(); }
    
    std::vector<std::string> getSequencesNames() const;
    void setSequencesNames(const std::vector<std::string>& names, bool checkNames = true) throw (Exception);
    void clear();
    VectorSequenceContainer * createEmptyContainer() const;
    
    int& valueAt(const std::string& sequenceName, size_t elementIndex) throw (SequenceNotFoundException, IndexOutOfBoundsException)
    {
      return getSequence_(sequenceName)[elementIndex];
    }

    const int& valueAt(const std::string& sequenceName, size_t elementIndex) const throw (SequenceNotFoundException, IndexOutOfBoundsException)
    {
      return getSequence(sequenceName)[elementIndex]; 
    }

    int& operator()(const std::string& sequenceName, size_t elementIndex)
    {
      return getSequence_(sequenceName)[elementIndex];
    }

    const int& operator()(const std::string& sequenceName, size_t elementIndex) const
    {
      return getSequence(sequenceName)[elementIndex]; 
    }
    
    int& valueAt(size_t sequenceIndex, size_t elementIndex) throw (IndexOutOfBoundsException)
    {
      return getSequence_(sequenceIndex)[elementIndex];
    }

    const int& valueAt(size_t sequenceIndex, size_t elementIndex) const throw (IndexOutOfBoundsException)
    {
      return getSequence(sequenceIndex)[elementIndex];
    }

    int& operator()(size_t sequenceIndex, size_t elementIndex)
    {
      return getSequence_(sequenceIndex)[elementIndex];
    }
    const int & operator()(size_t sequenceIndex, size_t elementIndex) const
    {
      return getSequence(sequenceIndex)[elementIndex];
    } 
    /** @} */


    /**
     * @name The OrderedSequenceContainer interface.
     *
     * @{
     */
    void setComments(const std::string & name, const Comments& comments) throw (SequenceNotFoundException)
    {
      AbstractSequenceContainer::setComments(name, comments);
    }

    void setComments(size_t sequenceIndex, const Comments& comments) throw (IndexOutOfBoundsException);
    size_t getSequencePosition(const std::string& name) const throw (SequenceNotFoundException);
    const Sequence& getSequence(size_t sequenceIndex) const throw (IndexOutOfBoundsException);
    void  setSequence(size_t sequenceIndex, const Sequence& sequence, bool checkName = true) throw (Exception);
    Sequence* removeSequence(size_t sequenceIndex) throw (IndexOutOfBoundsException);
    void deleteSequence(size_t sequenceIndex) throw (IndexOutOfBoundsException);
    /** @} */
    
    /**
     * @name Add sequence to this container.
     *
     * @{
     */
     
    /**
     * @brief Add a sequence at the end of the container.
     *
     * The sequence is copied into the container.
     * If checkNames is set to true, the method check if the name of the
     * sequence is already used in the container, and sends an exception if it
     * is the case. Otherwise, do not check the name: the method is hence faster,
     * but use it at your own risks!
     *
     * @param sequence The sequence to add.
     * @param checkName Tell if the method must check the name of the sequence
     * before adding it.
     * @throw Exception If the sequence couldn't be added to the container.
     */
    virtual void addSequence(const Sequence& sequence, bool checkName = true) throw (Exception);

    /**
     * @brief Add a sequence to the container at a particular position.
     *
     * The sequence is copied into the container.
     * If checkName is set to true, the method check if the name of the
     * sequence is already used in the container, and sends an exception if it
     * is the case. Otherwise, do not check the name: the method is hence faster,
     * but use it at your own risks!
     *
     * @param sequence The sequence to add.
     * @param sequenceIndex The position where to insert the new sequence.
     * All the following sequences will be pushed.
     * @param checkName Tell if the method must check the name of the sequence
     * before adding it.
     * @throw Exception If the sequence couldn't be added to the container.
     */
    virtual void addSequence(const Sequence& sequence, size_t sequenceIndex, bool checkName = true) throw (Exception);

  protected:

    /**
     * @name AbstractSequenceContainer methods.
     *
     * @{
     */
    Sequence& getSequence_(size_t i) throw (IndexOutOfBoundsException);
    Sequence& getSequence_(const std::string& name) throw (SequenceNotFoundException);
    /** @} */
};

} //end of namespace bpp.

#endif // _VECTORSEQUENCECONTAINER_H_

