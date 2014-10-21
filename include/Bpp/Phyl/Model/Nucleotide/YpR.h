//
// File: YpR.h
// Created by: Laurent Gueguen
// Created on: Wed Aug 3 2007
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _YpR_H_
#define _YpR_H_

#include "../AbstractSubstitutionModel.h"

#include <Bpp/Seq/Alphabet/RNY.h>


// From Utils:
#include <Bpp/Exceptions.h>

using namespace std;

namespace bpp
{
/**
 * @brief YpR  model.
 * @author Laurent Guéguen
 *
 * Model YpR, on RNY triplets, with independent positions.
 *
 * This model is made from neighbourhood parameters and
 * a nucleotidic probabilistic model such that the R/Y
 * condition is respected : with letters A, T, C, G
 * @f[
 * M=\begin{pmatrix}
 * . & \beta_T & \beta_C & \alpha_G \\
 * \beta_A & . & \alpha_C & \beta_G \\
 * \beta_A & \alpha_T & . & \beta_G \\
 * \alpha_A & \beta_T & \beta_C & .
 * \end{pmatrix}
 * @f]
 *
 * From this model, the rates are a multiplication of the rates for the
 * letters of the triplet.
 * For the first letter, on alphabet R, C, T:
 * @f[
 * M_1=\begin{pmatrix}
 * . & \beta_C & \beta_T \\
 * \beta_A + \beta_G & . & \alpha_T \\
 * \beta_A + \beta_G & \alpha_C & .
 * \end{pmatrix}
 * @f]
 * For the second letter, on alphabet A, G, C, T:
 * @f[
 * M_2=\begin{pmatrix}
 * . & \alpha_G & \beta_C & \beta_T \\
 * \alpha_A & . & \beta_C & \beta_T \\
 * \beta_A & \beta_G & . & \alpha_T \\
 * \beta_A & \beta_G & \alpha_C & .
 * \end{pmatrix}
 * @f]
 * For the third letter, on alphabet A, G, Y:
 * @f[
 * M_3 = \begin{pmatrix}
 * . & \alpha_G & \beta_C + \beta_T\\
 * \alpha_A & . & \beta_C + \beta_T\\
 * \beta_A & \beta_G & .
 * \end{pmatrix}.
 * @f]
 * And the model is the "union" of the 3 matrices, since in the
 * model  only the mutations concerning only one position  are
 * not null.
 * In addition to this, neighbour dependency parameters
 * (in inherited models) give the extra mutation rates for transitions
 * inside YpR dinucleotides.
 * For example from CpG to CpA and TpG, in relative proportion
 * to the single transition rate:
 *  CGx -> CAx   += rcGA*G->A
 *  CGx -> TGx   += rCgT*C->T
 *  xCG -> xCA   += rcGA*G->A
 *  xCG -> xTG   += rCgT*C->T
 *
 * Other YpR neighbour dependency rates are:
 *   TG -> CG, TG -> TA, CA -> TA, CA -> CG, TA -> CA, TA -> TC
 *
 * The generator in normalized such that there is on average
 *  a substitution per site per unit of time ON THE CENTRAL POSITION
 *  of the triplet.
 * @see AbstractSubstitutionModel
 *
 */

class YpR :
  public AbstractSubstitutionModel
{
protected:
  SubstitutionModel*  _pmodel;

  // Check that the model is good for YpR
  void check_model(SubstitutionModel* const) const
  throw (Exception);

  std::string _nestedPrefix;

protected:
  /**
   * @brief Build a new YpR substitution model, with no dependency
   *   parameters
   */

  YpR(const RNY*, SubstitutionModel* const, const std::string& prefix);

  YpR(const YpR&, const std::string& prefix);

  YpR(const YpR& ypr);

  YpR& operator=(const YpR& ypr)
  {
    AbstractParameterAliasable::operator=(ypr);
    AbstractSubstitutionModel::operator=(ypr);
    _nestedPrefix = ypr._nestedPrefix;
    _pmodel = ypr._pmodel->clone();
    return *this;
  }

public:
  virtual ~YpR()
  {
    if (_pmodel)
      delete _pmodel;
    _pmodel = 0;
  }

protected:
  void updateMatrices(double, double, double, double,
                      double, double, double, double);

  string getNestedPrefix() const
  {
    return _nestedPrefix;
  }

public:
  //  virtual std::string getName() const;

  const SubstitutionModel* getNestedModel() const {return _pmodel;}
  
  size_t getNumberOfStates() const { return 36; }

  virtual void updateMatrices();

  virtual void setNamespace(const std::string&);

  void fireParameterChanged(const ParameterList& parameters)
  {
    AbstractSubstitutionModel::fireParameterChanged(parameters);
   _pmodel->matchParametersValues(parameters);
   updateMatrices();
  }
};
}


// //////////////////////////////////////
// //////// YpR_symetrical

namespace bpp
{
/**
 * @brief symetrical YpR  model.
 *
 * Model YpR, on RNY triplets with symetrical dependency parameters
 *
 * The neighbour dependency parameters are noted as: XyZ for
 *  XpY -> ZpY substitution. Since the process are symetrical,
 *  only these are necessary.
 * @see YpR
 *
 */

class YpR_Sym :
  public YpR
{
public:
  /**
   * @brief Build a new YpR_Sym substitution model.
   * @param CgT, TgC, CaT, TaC neighbour dependency parameters
   * @param alph RNY alphabet
   * @param pm Substitution model.
   */

  YpR_Sym(const RNY* alph,
          SubstitutionModel* pm,
          double CgT = 0., double TgC = 0.,
          double CaT = 0., double TaC = 0.);

  YpR_Sym(const YpR_Sym&);

  virtual ~YpR_Sym() {}

  YpR_Sym* clone() const { return new YpR_Sym(*this); }

  std::string getName() const;

  void updateMatrices();
};
}

// //////////////////////////////////////
// //////// YpR_general

namespace bpp
{
/**
 * @brief General YpR  model.
 *
 * Model YpR, on RNY triplets with general dependency parameters
 *
 * The neighbour dependency parameters are noted as: XyZ for
 *  XpY -> ZpY substitution and xYZ for XpY -> XpZ substitution.
 *
 * @see YpR
 *
 */

class YpR_Gen :
  public YpR
{
public:
  /**
   * @brief Build a new YpR_Gen substitution model.
   * @param CgT, cGA, TgC, tGA, CaT, cAG, TaC, tAG neighbour
   * dependency parameters
   * @param alph RNY alphabet
   * @param pm Substitution model.
   */

  YpR_Gen(const RNY* alph,
          SubstitutionModel* pm,
          double CgT = 0., double cGA = 0.,
          double TgC = 0., double tGA = 0.,
          double CaT = 0., double cAG = 0.,
          double TaC = 0., double tAG = 0.);

  YpR_Gen(const YpR_Gen&);

  virtual ~YpR_Gen() {}

  YpR_Gen* clone() const { return new YpR_Gen(*this); }

  std::string getName() const;

  void updateMatrices();
};
}


#endif // _YpR_H_


