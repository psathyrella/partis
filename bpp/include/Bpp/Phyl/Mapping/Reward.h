//
// File: Reward.h
// Created by: Laurent Guéguen
// Created on: mercredi 27 mars 2013, à 09h 58
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

#ifndef _REWARD_H_
#define _REWARD_H_

#include "../Model/SubstitutionModel.h"

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex1.h>

//From the STL:
#include <vector>

namespace bpp
{

  /**
   * @brief The Reward interface.
   *
   * Provide a method to compute the reward of a real-valued function
   * @f$f@f$ of the states on a branch. Namely, on a branch of length
   * @f$t@f$, with initial state @f$x@f$ and final state @f$y@f$, if
   * @f$t_s@f$ is the time spent in each state @f$s@f$, the reward of
   * @f$f@f$ is @f$\sum_{s} f(s).E(t_s) @f$.
   * 
   * @author Laurent Guéguen
   *
   * See:
   * Minin, V.N. and Suchard, M.A., 
   * Fast, accurate and simulation-free stochastic mapping
   * Philosophical Transactions of the Royal Society B 2008 363:3985-95.
   */
  class Reward:
    public virtual Clonable
  {
  public:
    Reward() {}
    virtual ~Reward() {}
    virtual Reward* clone() const = 0;
	
  public:
    /**
     * @return Tell if an alphabet index  has been attached to this class.
     */
    virtual bool hasAlphabetIndex() const = 0;

    /**
     * @return The AlphabetIndex1 object associated to this instance.
     * The alphabet index contains the value associated to each state.
     */
    virtual const AlphabetIndex1* getAlphabetIndex() const = 0;

    /**
     * @return The AlphabetIndex1 object associated to this instance.
     * The alphabet index contains the value associated to each state.
     */
    
    virtual AlphabetIndex1* getAlphabetIndex() = 0;

    /**
     * @param alphind The new AlphabetIndex1 object to be associated to this instance.
     */
    
    virtual void setAlphabetIndex(AlphabetIndex1* alphind) = 0;

    /**
     * @brief Short cut function, equivalent to getSubstitutionRegister()->getAlphabet().
     *
     * @return The alphabet associated to this substitution count.
     */
    virtual const Alphabet* getAlphabet() const { return getAlphabetIndex()->getAlphabet(); }

    /**
     * @brief Short cut function, equivalent to getSubstitutionRegister()->getAlphabet()->getSize().
     *
     * @return The number of states in the model/alphabet.
     */
    virtual size_t getNumberOfStates() const { return getAlphabet()->getSize(); }


    /**
     * @brief Get the reward of susbstitutions on a branch, given the initial and final states, and the branch length.
     *
     * @param initialState The initial state.
     * @param finalState   The final state.
     * @param length       The length of the branch.
     * @return The reward of the function on a branch of specified length and
     * according to initial and final states.
     */
    virtual double getReward(int initialState, int finalState, double length) const = 0;
		
    /**
     * @brief Get the rewards on a branch, for each initial and final
     * states, and given the branch length.
     *
     * @param length       The length of the branch.
     * @return A matrix with all rewards for each initial and final states.
     */
    virtual Matrix<double>* getAllRewards(double length) const = 0;

    /**
     * @brief Set the substitution model associated with this reward, if relevant.
     *
     * @param model The substitution model to use with this reward.
     */
    virtual void setSubstitutionModel(const SubstitutionModel* model) = 0;
  };

  /**
   * @brief Basic implementation of the the Reward interface.
   *
   * This partial implementation deals with the AlphabetIndex1 gestion, by maintaining a pointer.
   */
  class AbstractReward:
    public virtual Reward
  {
  protected:
    AlphabetIndex1* alphIndex_;

  public:
    AbstractReward(AlphabetIndex1* alphIndex):
      alphIndex_(alphIndex)
    {}

    AbstractReward(const AbstractReward& ar):
      alphIndex_(dynamic_cast<AlphabetIndex1*>(ar.alphIndex_->clone()))
    {}

    AbstractReward& operator=(const AbstractReward& ar) {
      if (alphIndex_)
        delete alphIndex_;
      alphIndex_ = dynamic_cast<AlphabetIndex1*>(ar.alphIndex_->clone());
      return *this;
    }

    ~AbstractReward() {
      if (alphIndex_)
        delete alphIndex_;
    }

  public:
    bool hasAlphabetIndex() const { return (alphIndex_ != 0); }

    /*
     *@brief attribution of an AlphabetIndex1
     *
     *@param alphIndex pointer to a AlphabetIndex1
     *
     */
    
    void setAlphabetIndex(AlphabetIndex1* alphIndex) {
      if (alphIndex_) delete alphIndex_;
      alphIndex_ = alphIndex;
      alphabetIndexHasChanged();
    }

    const AlphabetIndex1* getAlphabetIndex() const { return alphIndex_; }
    
    AlphabetIndex1* getAlphabetIndex() { return alphIndex_; }

  protected:
    virtual void alphabetIndexHasChanged() = 0;
  };

} //end of namespace bpp.

#endif //_REWARD_H_

