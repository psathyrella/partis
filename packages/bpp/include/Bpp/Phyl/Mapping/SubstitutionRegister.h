//
// File: SubstitutionRegister.h
// Created by: Julien Dutheil
// Created on: Mon Dec 6 16:32 2010
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

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

#ifndef _SUBSTITUTIONREGISTER_H_
#define _SUBSTITUTIONREGISTER_H_

// From bpp-core:
#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief The SubstitutionRegister interface.
 *
 * Substitution registers are simple classes that define categories of substitutions, and assign them an index.
 *
 * @author Julien Dutheil
 */
class SubstitutionRegister :
  public virtual Clonable
{
public:
  SubstitutionRegister() {}
  virtual ~SubstitutionRegister() {}

#ifndef NO_VIRTUAL_COV
  virtual SubstitutionRegister* clone() const = 0;
#endif

public:
  /**
   * @return The alphabet associated to this instance.
   */
  virtual const Alphabet* getAlphabet() const = 0;

  /**
   * @return The number of substitution types supported by this class.
   */
  virtual size_t getNumberOfSubstitutionTypes() const = 0;

  /**
   * @brief Get the substitution type far a given pair of states.
   *
   * @param fromState Initial state (should be a state supported by the specified alphabet).
   * @param toState   Final state (should be a state supported by the specified alphabet).
   * @return The index of the corresponding substitution type, ranging from 0 to 'getNumberOfSubstitutionTypes' + 1,
   * as non-substitution (that is when fromState == toState) will always return 0.
   */
  virtual size_t getType(int fromState, int toState) const = 0;

  /**
   * @brief Get the name of a given substitution type.
   *
   * This method is only used for user-friendlyness purposes, not computational goal.
   * I can therefore be left unimplemented in some cases.
   *
   * @param type Index of the substitution (should be an size_t contained in the register).
   * @return A string describing the substitution type.
   */
  virtual std::string getTypeName(size_t type) const = 0;
};

class AbstractSubstitutionRegister :
  public virtual SubstitutionRegister
{
protected:
  const Alphabet* alphabet_;

public:
  AbstractSubstitutionRegister(const Alphabet* alphabet) :
    alphabet_(alphabet)
  {}

  AbstractSubstitutionRegister(const AbstractSubstitutionRegister& asr) :
    alphabet_(asr.alphabet_)
  {}

  AbstractSubstitutionRegister& operator=(const AbstractSubstitutionRegister& asr)
  {
    alphabet_ = asr.alphabet_;
    return *this;
  }

  virtual ~AbstractSubstitutionRegister() {}

public:
  const Alphabet* getAlphabet() const { return alphabet_; }
};

/**
 * @brief Gather states into defined categories, and count the changes between categories.
 *
 * Optionally allows for within categories substitutions.
 */
class CategorySubstitutionRegister :
  public AbstractSubstitutionRegister
{
protected:
  bool within_;
  size_t nbCategories_;
  mutable std::map<int, size_t> categories_;
  std::vector<std::string> categoryNames_;
  std::vector< std::vector<size_t> > index_;
  std::vector< std::vector<size_t> > revIndex_;

public:
  /**
   * @brief Build a new substitution register with categories. This class is mean to be inherited.
   *
   * @param alphabet The input alphabet.
   * @param within Specifies if within categories substitutions should be counted as well.
   */
  CategorySubstitutionRegister(const Alphabet* alphabet, bool within = false) :
    AbstractSubstitutionRegister(alphabet),
    within_(within),
    nbCategories_(0),
    categories_(),
    categoryNames_(),
    index_(),
    revIndex_()
  {}

protected:
  template<class T>
  void setCategories(const std::map<int, T>& categories)
  {
    // First index categories:
    nbCategories_ = 0;
    std::map<T, size_t> cats;
    for (typename std::map<int, T>::const_iterator it = categories.begin(); it != categories.end(); ++it)
    {
      if (cats.find(it->second) == cats.end())
      {
        ++nbCategories_;
        cats[it->second] = nbCategories_;
      }
    }

    // Now creates categories:
    categories_.clear();
    categoryNames_.resize(nbCategories_);
    std::vector<int> types = alphabet_->getSupportedInts();
    for (size_t i = 0; i < types.size(); ++i)
    {
      typename std::map<int, T>::const_iterator it = categories.find(types[i]);
      if (it != categories.end())
      {
        categories_[types[i]] = cats[it->second];
        categoryNames_[cats[it->second] - 1] += alphabet_->intToChar(types[i]);
      }
      else
      {
        categories_[types[i]] = 0;
      }
    }

    size_t count = 1;
    index_.resize(nbCategories_);
    for (size_t i = 0; i < index_.size(); ++i)
    {
      index_[i].resize(nbCategories_);
      for (size_t j = 0; j < index_.size(); ++j)
      {
        if (j != i)
        {
          index_[i][j] = count++;
          std::vector<size_t> pos(2);
          pos[0] = i;
          pos[1] = j;
          revIndex_.push_back(pos);
        }
      }
    }
    if (within_)
    {
      for (size_t i = 0; i < index_.size(); ++i)
      {
        index_[i][i] = count++;
        std::vector<size_t> pos(2);
        pos[0] = i;
        pos[1] = i;
        revIndex_.push_back(pos);
      }
    }
  }

public:
  virtual size_t getCategory(int state) const
  {
    if (!alphabet_->isIntInAlphabet(state))
      throw Exception("CategorySubstitutionRegister::getCategory(). State is not supported by alphabet.");
    return categories_[state];
  }

  virtual size_t getCategoryFrom(size_t type) const
  {
    if (type <= nbCategories_ * (nbCategories_ - 1))
    {
      return revIndex_[type - 1][0] + 1;
    }
    else
    {
      if (within_)
        return revIndex_[type - 1][0] + 1;
      else
        throw Exception("CategorySubstitutionRegister::getCategoryFrom. Bad substitution type.");
    }
  }

  virtual size_t getCategoryTo(size_t type) const
  {
    if (type <= nbCategories_ * (nbCategories_ - 1))
    {
      return revIndex_[type - 1][1] + 1;
    }
    else
    {
      if (within_)
        return revIndex_[type - 1][1] + 1;
      else
        throw Exception("CategorySubstitutionRegister::getCategoryTo. Bad substitution type.");
    }
  }

  virtual std::string getCategoryName(size_t category) const
  {
    return categoryNames_[category - 1];
  }

  virtual bool allowWithin() const { return within_; }

  size_t getNumberOfCategories() const { return nbCategories_; }

  size_t getNumberOfSubstitutionTypes() const { return nbCategories_ * (nbCategories_ - 1) + (within_ ? nbCategories_ : 0); }

  virtual size_t getType(int fromState, int toState) const
  {
    size_t fromCat = categories_[fromState];
    size_t toCat   = categories_[toState];
    if (fromCat > 0 && toCat > 0)
      return index_[fromCat - 1][toCat - 1];
    else
      return 0;
  }

  std::string getTypeName(size_t type) const
  {
    return getCategoryName(getCategoryFrom(type)) + "->" +  getCategoryName(getCategoryTo(type));
  }
};


/**
 * @brief Count all substitutions.
 *
 * This register has only 1 substitution type, mapped as:
 * - 0 not a substitution
 * - 1 a substitution
 */
class TotalSubstitutionRegister :
  public AbstractSubstitutionRegister
{
public:
  TotalSubstitutionRegister(const Alphabet* alphabet) :
    AbstractSubstitutionRegister(alphabet)
  {}

  TotalSubstitutionRegister* clone() const { return new TotalSubstitutionRegister(*this); }

public:
  size_t getNumberOfSubstitutionTypes() const { return 1; }

  size_t getType(int fromState, int toState) const
  {
    return fromState == toState ? 0 : 1;
  }

  std::string getTypeName(size_t type) const
  {
    if (type == 0)
    {
      return "no substitution";
    }
    else if (type == 1)
    {
      return "substitution";
    }
    else
    {
      throw Exception("TotalSubstitutionRegister::getTypeName. Bad substitution type.");
    }
  }
};

 /**
   * @brief Completion of a given substitution register to consider
   * all substitutions. The new substitutions are considered in an
   * additional type.
   *
   */
  class CompleteSubstitutionRegister :
    public AbstractSubstitutionRegister
  {
  private:
    const SubstitutionRegister* preg_;

    bool isRegComplete_;
    
  public:
    CompleteSubstitutionRegister(const SubstitutionRegister& reg) :
      AbstractSubstitutionRegister(reg.getAlphabet()),
      preg_(reg.clone()), isRegComplete_(true)
    {
      size_t size=reg.getAlphabet()->getSize();
      for (int i=0; i< (int)size; i++)
        for (int j=0; j< (int)size; j++)
          if ((i!=j) && reg.getType(i,j)==0){
            isRegComplete_=false;
            return;
          }
    }

    CompleteSubstitutionRegister* clone() const { return new CompleteSubstitutionRegister(*this); }

    CompleteSubstitutionRegister(const CompleteSubstitutionRegister& csr) :
      AbstractSubstitutionRegister(csr),
      preg_(csr.preg_->clone()),
      isRegComplete_(csr.isRegComplete_)
    {}

    CompleteSubstitutionRegister& operator=(const CompleteSubstitutionRegister& csr)
    {
      AbstractSubstitutionRegister::operator=(csr);
      preg_=csr.preg_->clone();
      isRegComplete_=csr.isRegComplete_;
      return *this;
    }
    
    ~CompleteSubstitutionRegister() {
      if (preg_)
        delete preg_;
      preg_=0;
    }
    
  public:
    size_t getNumberOfSubstitutionTypes() const {
      return preg_->getNumberOfSubstitutionTypes()+(isRegComplete_?0:1);
    }

    size_t getType(int fromState, int toState) const
    {
      size_t t=preg_->getType(fromState,toState);
      if (t==0)
        return getNumberOfSubstitutionTypes();
      else
        return t;
    }

    std::string getTypeName(size_t type) const
    {
      try {
        return preg_->getTypeName(type);
      }
      catch (Exception& e) {
        if (type == getNumberOfSubstitutionTypes())
          return "Completion substitution";
        else
          throw Exception("CompleteSubstitutionRegister::getTypeName. Bad substitution type.");
      }
    }
    
  };

/**
 * @brief Distinguishes all types of substitutions.
 *
 * This register has all n * (n-1) substitution type, where n is the size of the alphabet, mapped as:
 * - 0 not a substitution
 * - x in [1, n(n-1)] a substitution
 */
class ComprehensiveSubstitutionRegister :
  public CategorySubstitutionRegister
{
public:
  ComprehensiveSubstitutionRegister(const Alphabet* alphabet, bool within = false) :
    CategorySubstitutionRegister(alphabet, within)
  {
    std::map<int, int> categories;
    for (int i = 0; i < static_cast<int>(alphabet->getSize()); ++i)
    {
      categories[i] = i;
    }
    setCategories<int>(categories);
  }

  ComprehensiveSubstitutionRegister* clone() const { return new ComprehensiveSubstitutionRegister(*this); }
};

/**
 * @brief Sets a Register based on a matrix of integers. If M is the
 *  matrix, M[i,j] is the number of the substitution type from i to j,
 *  or 0 if there is no substitution type from i to j.
 *
 * @author Laurent Guéguen
 */
class GeneralSubstitutionRegister :
  public AbstractSubstitutionRegister
{
protected:
  /**
   * @brief The size of the matrix, i.e. the number of states.
   */
  size_t size_;

  /**
   * @brief The matrix of the substitution register.
   */
  
  RowMatrix<size_t> matrix_;

  /**
   * @brief The map from substitution types to the map of from states
   * to the vector of target states.
   *
   * This is the reverse information of matrix_
   *
   */
  
  std::map<size_t, std::map<size_t, std::vector<size_t> > > types_;
  
public:
  GeneralSubstitutionRegister(const Alphabet* alphabet) :
    AbstractSubstitutionRegister(alphabet),
    size_(alphabet->getSize()),
    matrix_(size_,size_),
    types_()
  {}

  GeneralSubstitutionRegister(const Alphabet* alphabet, const RowMatrix<size_t>& matrix) :
    AbstractSubstitutionRegister(alphabet),
    size_(alphabet->getSize()),
    matrix_(matrix),
    types_()
  {
    if (matrix_.getNumberOfRows()!=size_)
      throw DimensionException("GeneralSubstitutionRegister", size_, matrix_.getNumberOfRows());
    if (matrix_.getNumberOfColumns()!=size_)
        throw DimensionException("GeneralSubstitutionRegister", size_, matrix_.getNumberOfColumns());
    updateTypes_();
  }

  GeneralSubstitutionRegister(const GeneralSubstitutionRegister& gsr) :
    AbstractSubstitutionRegister(gsr),
    size_(gsr.size_),
    matrix_(gsr.matrix_),
    types_(gsr.types_)
  {}

  GeneralSubstitutionRegister& operator=(const GeneralSubstitutionRegister& gsr)
  {
    AbstractSubstitutionRegister::operator=(gsr);
    size_=gsr.size_;
    matrix_=gsr.matrix_;
    types_=gsr.types_;
    return *this;
  }

  GeneralSubstitutionRegister* clone() const { return new GeneralSubstitutionRegister(*this); }

  virtual ~GeneralSubstitutionRegister() {}

  size_t getType(int i,int j) const
  {
    return matrix_(i,j);
  }

  size_t getNumberOfSubstitutionTypes() const
  {
    return types_.size();
  }

  
  /**
   *@brief names of the types are their number.
   *
   */
  
  std::string getTypeName(size_t type) const
  {
    if (types_.find(type)!=types_.end())
      return TextTools::toString(type);

    throw Exception("Bad type number " + TextTools::toString(type) + " in GeneralSubstitutionRegister::getTypeName.");
  }

private:
  void updateTypes_();
};
  
/**
 * @brief Distinguishes AT<->GC from GC<->AT.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a substitution
 * - 1 a AT->GC substitution
 * - 2 a GC->AT substitution
 */
class GCSubstitutionRegister :
  public CategorySubstitutionRegister
{
public:
  GCSubstitutionRegister(const NucleicAlphabet* alphabet, bool within = false) :
    CategorySubstitutionRegister(alphabet, within)
  {
    std::map<int, int> categories;
    categories[0] = 1;
    categories[1] = 2;
    categories[2] = 2;
    categories[3] = 1;
    setCategories<int>(categories);
  }

  GCSubstitutionRegister* clone() const { return new GCSubstitutionRegister(*this); }
};

/**
 * @brief Distinguishes transitions from transversions.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a substitution
 * - 1 a transition
 * - 2 a transversion
 */
class TsTvSubstitutionRegister :
  public AbstractSubstitutionRegister
{
public:
  TsTvSubstitutionRegister(const NucleicAlphabet* alphabet) :
    AbstractSubstitutionRegister(alphabet)
  {}

  TsTvSubstitutionRegister* clone() const { return new TsTvSubstitutionRegister(*this); }

public:
  size_t getNumberOfSubstitutionTypes() const { return 2; }

  size_t getType(int fromState, int toState) const
  {
    if (fromState == toState)
      return 0;  // nothing happens
    if ((fromState == 0 && toState == 2)
        || (fromState == 2 && toState == 0)
        || (fromState == 1 && toState == 3)
        || (fromState == 3 && toState == 1))
      return 1;  // This is a transition
    return 2; // This is a transversion
  }

  std::string getTypeName(size_t type) const
  {
    if (type == 0)
    {
      return "no substitution";
    }
    else if (type == 1)
    {
      return "transition";
    }
    else if (type == 2)
    {
      return "transversion";
    }
    else
    {
      throw Exception("TsTvSubstitutionRegister::getTypeName. Bad substitution type.");
    }
  }
};

/**
 * @brief Distinguishes synonymous from non-synonymous substitutions.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a substitution
 * - 1 a synonymous substitution
 * - 2 a non-synonymous substitution
 */
class DnDsSubstitutionRegister :
  public AbstractSubstitutionRegister
{
private:
  const GeneticCode* code_;
  bool countMultiple_;

public:
  DnDsSubstitutionRegister(const GeneticCode* gc, bool countMultiple = false) :
    AbstractSubstitutionRegister(gc->getSourceAlphabet()),
    code_(gc),
    countMultiple_(countMultiple)
  {}

  DnDsSubstitutionRegister(const DnDsSubstitutionRegister& reg) :
    AbstractSubstitutionRegister(reg),
    code_(reg.code_),
    countMultiple_(reg.countMultiple_)
  {}

  DnDsSubstitutionRegister& operator=(const DnDsSubstitutionRegister& reg)
  {
    AbstractSubstitutionRegister::operator=(reg);
    code_ = reg.code_;
    countMultiple_ = reg.countMultiple_;
    return *this;
  }

  DnDsSubstitutionRegister* clone() const { return new DnDsSubstitutionRegister(*this); }

public:
  size_t getNumberOfSubstitutionTypes() const { return 2; }

  size_t getType(int fromState, int toState) const
  {
    const CodonAlphabet* cAlpha = dynamic_cast<const CodonAlphabet*>(alphabet_);
    if (code_->isStop(fromState) || code_->isStop(toState))
      return 0;
    if (fromState == toState)
      return 0;  // nothing happens
    if (!countMultiple_)
    {
      size_t countPos = 0;
      if (cAlpha->getFirstPosition(fromState) != cAlpha->getFirstPosition(toState))
        countPos++;
      if (cAlpha->getSecondPosition(fromState) != cAlpha->getSecondPosition(toState))
        countPos++;
      if (cAlpha->getThirdPosition(fromState) != cAlpha->getThirdPosition(toState))
        countPos++;
      if (countPos > 1)
        return 0;
    }
    return code_->areSynonymous(fromState, toState) ? 1 : 2;
  }

  std::string getTypeName (size_t type) const
  {
    if (type == 0)
    {
      return "no substitution";
    }
    else if (type == 1)
    {
      return "synonymous";
    }
    else if (type == 2)
    {
      return "non synonymous";
    }
    else
    {
      throw Exception("DnDsSubstitutionRegister::getTypeName. Bad substitution type.");
    }
  }
};

/**
 * @brief Distinguishes AT->GC vs GC->AT inside synonymous
 * substitutions on third codon position.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a counted substitution
 * - 1 a AT->GC synonymous substitution
 * - 2 a GC->AT synonymous substitution
 *
 * Multiple substitutions are forbidden.
 *
 */

class GCSynonymousSubstitutionRegister :
  public CategorySubstitutionRegister
{
private:
  const GeneticCode* code_;

public:
  GCSynonymousSubstitutionRegister(const GeneticCode* gc, bool within = false) :
    CategorySubstitutionRegister(gc->getSourceAlphabet(), within),
    code_(gc)
  {
    const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(gc->getSourceAlphabet());

    std::map<int, int> categories;
    for (int i = 0; i < static_cast<int>(pCA->getSize()); i++)
    {
      int n = pCA->getThirdPosition(i);
      switch (n)
      {
      case 0:
      case 3:
        categories[i] = 1;
        break;
      case 1:
      case 2:
        categories[i] = 2;
        break;
      }
    }
    setCategories<int>(categories);
  }

  GCSynonymousSubstitutionRegister(const GCSynonymousSubstitutionRegister& reg) :
    CategorySubstitutionRegister(reg),
    code_(reg.code_)
  {}

  GCSynonymousSubstitutionRegister& operator=(const GCSynonymousSubstitutionRegister& reg)
  {
    CategorySubstitutionRegister::operator=(reg);
    code_ = reg.code_;
    return *this;
  }

  GCSynonymousSubstitutionRegister* clone() const { return new GCSynonymousSubstitutionRegister(*this); }

public:
  size_t getNumberOfSubstitutionTypes() const { return 2; }

  size_t getType(int fromState, int toState) const
  {
    const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(code_->getSourceAlphabet());
    if (code_->isStop(fromState) || code_->isStop(toState) || !code_->areSynonymous(fromState, toState))
      return 0;

    // only substitutions between 3rd positions

    if ((pCA->getFirstPosition(fromState) != pCA->getFirstPosition(toState)) ||
        (pCA->getSecondPosition(fromState) != pCA->getSecondPosition(toState)))
      return 0;

    size_t fromCat = categories_[fromState];
    size_t toCat   = categories_[toState];

    if (fromCat > 0 && toCat > 0)
      return index_[fromCat - 1][toCat - 1];
    else
      return 0;
  }

  std::string getTypeName (size_t type) const
  {
    if (type == 0)
    {
      return "no AT<->GC substitution or non-synonymous substitution";
    }
    else if (type == 1)
    {
      return "AT->GC synonymous";
    }
    else if (type == 2)
    {
      return "GC->AT synonymous";
    }
    else
    {
      throw Exception("GCSynonymousSubstitutionRegister::getTypeName. Bad substitution type.");
    }
  }
};
} // end of namespace bpp.

#endif // _SUBSTITUTIONREGISTER_H_

