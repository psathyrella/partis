//
// File: UserAlphabetIndex1.h
// Created by: Laurent Guéguen
// Created on: vendredi 29 mars 2013, à 13h 05
//

/*
   Copyright or © or Copr. CNRS, (November 17, 2004)

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

#ifndef _USERALPHABETINDEX1_H_
#define _USERALPHABETINDEX1_H_

#include "AlphabetIndex1.h"

// From the STL:
#include <vector>

namespace bpp
{
  /**
   * @brief Alphabet index given by user.
   */

  class UserAlphabetIndex1 :
  public virtual AlphabetIndex1
{
private:
  const Alphabet* alph_;
  std::vector<double> index_;
  
public:
  UserAlphabetIndex1(const Alphabet* alph) :
  alph_(alph),
  index_(alph->getSize(),0)
  {    
  }

  UserAlphabetIndex1(const UserAlphabetIndex1& uAlph) :
    alph_(uAlph.alph_),//->clone()),
    index_(uAlph.index_)
  {    
  }

  UserAlphabetIndex1& operator=(const UserAlphabetIndex1& uAlph)
  {
    alph_=uAlph.alph_;//->clone();
    index_=uAlph.index_;

    return *this;
  }    
  
  virtual ~UserAlphabetIndex1() {}

  UserAlphabetIndex1* clone() const { return new UserAlphabetIndex1(*this); }

public:
  double getIndex(int state) const
  {
    if (state < 0 || state >= (int)alph_->getSize())
      throw BadIntException(state, "UserAlphabetIndex1::getIndex(). Invalid state.", alph_);
    return index_[state];
  }

  void setIndex(int state, double val) 
  {
    if (state < 0 || state >= (int)alph_->getSize())
      throw BadIntException(state, "UserAlphabetIndex1::getIndex(). Invalid state.", alph_);
    index_[state]=val;
  }

  double getIndex(const std::string& state) const
  {
    return index_[alph_->charToInt(state)];
  }

  void setIndex(const std::string& state, double val) 
  {
    index_[alph_->charToInt(state)]=val;
  }

  std::vector<double>* getIndexVector() const { return new std::vector<double>(index_); }

  const Alphabet* getAlphabet() const { return alph_; }
};
} // end of namespace bpp.

#endif // _USERALPHABETINDEX1_H_

