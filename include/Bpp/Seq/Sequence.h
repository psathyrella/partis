//
// File: Sequence.h
// Created by: Guillaume Deuchst
//             Julien Dutheil
// Created on: Tue Aug 21 2003
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

#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include "SymbolList.h"
#include "SequenceExceptions.h"

// From the STL:
#include <string>
#include <vector>

namespace bpp
{

/**
 * @brief Declaration of Comments type.
 *
 * Comments are defined as a std::vector of std::strings to allow the later creation of a
 * full Comments class.
 */
typedef std::vector<std::string> Comments;

/**
 * @brief The sequence interface. 
 *
 * This is a general purpose container, containing an ordered list of states.
 * The states that allowed to be present in the sequence are defined
 * by an alphabet object.
 *
 * Sequence objets also contain a name attribute and potentially several comment lines.
 * A sequence object is also event-driven, allowing easy extension.
 *
 * @see Alphabet
 */
class Sequence:
  public virtual SymbolList
{
  public:
    virtual ~Sequence() {}

  public:
  
#ifndef NO_VIRTUAL_COV
    Sequence* clone() const = 0;
#endif
    
    /**
     * @name Setting/getting the name of the sequence.
     *
     * @{
     */
     
    /**
     * @brief Get the name of this sequence.
     *
     * @return This sequence name.
     */
    virtual const std::string& getName() const = 0;
    
    /**
     * @brief Set the name of this sequence.
     *
     * @param name The new name of the sequence.
     */
    virtual void setName(const std::string& name) = 0;    
    /** @} */
    
    /**
     * @name Setting/getting the comments associated to the sequence.
     *
     * @{
     */
     
    /**
     * @brief Get the comments associated to this sequence.
     *
     * @return The comments of the sequence.
     */
    virtual const Comments& getComments() const = 0;
    
    /**
     * @brief Set the comments associated to this sequence.
     *
     * @param comments The new comments of the sequence.
     */
    virtual void setComments(const Comments& comments) = 0;
    
    /** @} */
    
    /**
     * @name Adjusting the size of the sequence.
     *
     * @{
     */
     
    /**
     * @brief Set the whole content of the sequence.
     *
     * @param sequence The new content of the sequence.
     * @see The Sequence constructor for information about the way sequences are internaly stored.
     */
    virtual void setContent(const std::string& sequence) throw (BadCharException) = 0;
    virtual void setContent(const std::vector<int>& list) throw (BadIntException) = 0;
    virtual void setContent(const std::vector<std::string>& list) throw (BadCharException) = 0;

    /**
     * @brief Set up the size of a sequence from the right side.
     *
     * All new characters are filled with gaps.
     * If the specified size is < to the sequence size, the sequence will be truncated.
     *
     * @param newSize The new size of the sequence.
     */
    virtual void setToSizeR(size_t newSize) = 0;
    
    /**
     * @brief Set up the size of a sequence from the left side.
     *
     * All new characters are filled with gaps.
     * If the specified size is < to the sequence size, the sequence will be truncated.
     *
     * @param newSize The new size of the sequence.
     */
    virtual void setToSizeL(size_t newSize) = 0;

    /**
     * @brief Append the specified content to the sequence.
     *
     * @param content The content to append to the sequence.
     * @throw BadIntException If the content does not match the current alphabet.
     */
    virtual void append(const std::vector<int>& content) throw (BadIntException) = 0;

    /**
     * @brief Append the specified content to the sequence.
     *
     * @param content The content to append to the sequence.
     * @throw BadCharException If the content does not match the current alphabet.
     */
    virtual void append(const std::vector<std::string>& content) throw (BadCharException) = 0;

    /**
     * @brief Append the specified content to the sequence.
     *
     * @param content The content to append to the sequence.
     * @throw BadCharException If the content does not match the current alphabet.
     */
    virtual void append(const std::string& content) throw (BadCharException) = 0;
    /** @} */

};


/**
 * @brief A basic implementation of the Sequence interface. 
 *
 * This is a general purpose container, containing an ordered list of states.
 * The states that allowed to be present in the sequence are defined
 * by an alphabet object, which is passed to the sequence constructor by a pointer.
 *
 * For programming convenience, the states are stored as integers, but the translation toward
 * and from a char description is easily performed with the Alphabet classes.
 *
 * Sequence objets also contain a name attribute and potentially several comment lines.
 *
 * @see Alphabet
 */
class BasicSequence :
  public Sequence,
  public BasicSymbolList
{
  private:

    /**
     * @brief The sequence name.
     */
    std::string name_;

    /**
     * @brief The sequence comments.
     */
    Comments comments_;

  public:

    /**
     * @brief Empty constructor: build a void Sequence with just an Alphabet
     *
     * You can use it safely for all type of Alphabet in order to build an
     * empty Sequence i.e. without name nor sequence data.
     *
     * @param alpha    A pointer toward the Alphabet to be used with this Sequence.
     */
    BasicSequence(const Alphabet* alpha);

    /**
     * @brief Direct constructor: build a Sequence object from a std::string
     * You can use it safely for RNA, DNA and protein sequences.
     *
     * It can be used with codon sequences too, the std::string will be cut into
     * parts of size 3. But for more complicated alphabets, you should use one
     * complete constructors.
     *
     * @param name     The sequence name.
     * @param sequence The whole sequence to be parsed as a std::string.
     * @param alpha    A pointer toward the alphabet to be used with this sequence.
     */
    BasicSequence(const std::string& name, const std::string& sequence, const Alphabet* alpha) throw (BadCharException);
  
    /**
     * @brief Direct constructor: build a Sequence object from a std::string.
     * 
     * You can use it safely for RNA, DNA and protein sequences.
     *
     * It can be used with codon sequences too, the std::string will be cut into
     * tokens of size 3. But for more complicated alphabets, you should use one
     * complete constructors.
     *
     * @param name     The sequence name.
     * @param sequence The whole sequence to be parsed as a std::string.
     * @param comments Comments to add to the sequence.
     * @param alpha    A pointer toward the alphabet to be used with this sequence.
     */
    BasicSequence(const std::string& name, const std::string& sequence, const Comments& comments, const Alphabet* alpha) throw (BadCharException);
  
    /**
     * @brief General purpose constructor, can be used with any alphabet.
     *
     * You should note that the sequence is stored as a std::vector of int.
     * Hence each std::string in the std::vector will be translated using the alphabet object.
     *
     * @param name     The sequence name.
     * @param sequence The sequence content.
     * @param alpha    A pointer toward the alphabet to be used with this sequence.
     */
    BasicSequence(const std::string& name, const std::vector<std::string>& sequence, const Alphabet* alpha) throw (BadCharException);
    
    /**
     * @brief General purpose constructor, can be used with any alphabet.
     *
     * You should note that the sequence is stored as a std::vector of int.
     * Hence each std::string in the std::vector will be translated using the alphabet object.
     *
     * @param name     The sequence name.
     * @param sequence The sequence content.
     * @param comments Comments to add to the sequence.
     * @param alpha    A pointer toward the alphabet to be used with this sequence.
     */
    BasicSequence(const std::string& name, const std::vector<std::string>& sequence, const Comments& comments, const Alphabet* alpha) throw (BadCharException);
  
    /**
     * @brief General purpose constructor, can be used with any alphabet.
     *
     * @param name     The sequence name.
     * @param sequence The sequence content.
     * @param alpha    A pointer toward the alphabet to be used with this sequence.
     */
    BasicSequence(const std::string& name, const std::vector<int>& sequence, const Alphabet* alpha) throw (BadIntException);
    
    /**
     * @brief General purpose constructor, can be used with any alphabet.
     *
     * @param name     The sequence name.
     * @param sequence The sequence content.
     * @param comments Comments to add to the sequence.
     * @param alpha    A pointer toward the alphabet to be used with this sequence.
     */
    BasicSequence(const std::string& name, const std::vector<int>& sequence, const Comments& comments, const Alphabet* alpha) throw (BadIntException);

    /**
     * @brief The Sequence generic copy constructor. This does not perform a hard copy of the alphabet object.
     */
    BasicSequence(const Sequence& s);
   
    /**
     * @brief The Sequence copy constructor. This does not perform a hard copy of the alphabet object.
     */
    BasicSequence(const BasicSequence& s);
 
    /**
     * @brief The Sequence generic assignment operator. This does not perform a hard copy of the alphabet object.
     *
     * @return A ref toward the assigned Sequence.
     */
    BasicSequence& operator=(const Sequence& s);
   
    /**
     * @brief The Sequence assignment operator. This does not perform a hard copy of the alphabet object.
     *
     * @return A ref toward the assigned Sequence.
     */
    BasicSequence& operator=(const BasicSequence& s);

    virtual ~BasicSequence() {}

  public:
  
    /**
     * @name The Clonable interface
     *
     * @{
     */
    BasicSequence* clone() const { return new BasicSequence(*this); }
    /** @} */
        
    
    /**
     * @name Setting/getting the name of the sequence.
     *
     * @{
     */
     
    /**
     * @brief Get the name of this sequence.
     *
     * @return This sequence name.
     */
    const std::string& getName() const { return name_; }
    
    /**
     * @brief Set the name of this sequence.
     *
     * @param name The new name of the sequence.
     */
    void setName(const std::string& name) { name_ = name; }
    
    /** @} */
    
    /**
     * @name Setting/getting the comments associated to the sequence.
     *
     * @{
     */
     
    /**
     * @brief Get the comments associated to this sequence.
     *
     * @return The comments of the sequence.
     */
    const Comments& getComments() const { return comments_; }
    
    /**
     * @brief Set the comments associated to this sequence.
     *
     * @param comments The new comments of the sequence.
     */
    void setComments(const Comments& comments) { comments_ = comments; }
    
    /** @} */
    
    
    /**
     * @name Adjusting the size of the sequence.
     *
     * @{
     */
     
    /**
     * @brief Set the whole content of the sequence.
     *
     * @param sequence The new content of the sequence.
     * @see The Sequence constructor for information about the way sequences are internaly stored.
     */
    virtual void setContent(const std::string& sequence) throw (BadCharException);
    void setContent(const std::vector<int>& list) throw (BadIntException)
    {
      BasicSymbolList::setContent(list);
    }
    void setContent(const std::vector<std::string>& list) throw (BadCharException)
    {
      BasicSymbolList::setContent(list);
    }


    /**
     * @brief Set up the size of a sequence from the right side.
     *
     * All new characters are filled with gaps.
     * If the specified size is < to the sequence size, the sequence will be truncated.
     *
     * @param newSize The new size of the sequence.
     */
    virtual void setToSizeR(size_t newSize);
    
    /**
     * @brief Set up the size of a sequence from the left side.
     *
     * All new characters are filled with gaps.
     * If the specified size is < to the sequence size, the sequence will be truncated.
     *
     * @param newSize The new size of the sequence.
     */
    virtual void setToSizeL(size_t newSize);

    /**
     * @brief Append the specified content to the sequence.
     *
     * @param content The content to append to the sequence.
     * @throw BadIntException If the content does not match the current alphabet.
     */
    virtual void append(const std::vector<int>& content) throw (BadIntException);

    /**
     * @brief Append the specified content to the sequence.
     *
     * @param content The content to append to the sequence.
     * @throw BadCharException If the content does not match the current alphabet.
     */
    virtual void append(const std::vector<std::string>& content) throw (BadCharException);

    /**
     * @brief Append the specified content to the sequence.
     *
     * @param content The content to append to the sequence.
     * @throw BadCharException If the content does not match the current alphabet.
     */
    virtual void append(const std::string& content) throw (BadCharException);

    /** @} */

};

} //end of namespace bpp.

#endif // _SEQUENCE_H_

