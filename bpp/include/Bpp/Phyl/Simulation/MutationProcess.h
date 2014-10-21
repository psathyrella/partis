//
// File: MutationProcess.h
// Created by: Julien Dutheil
// Created on: Wed Mar 12 16:11:44 2003
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
 
#ifndef _MUTATIONPROCESS_H_
#define _MUTATIONPROCESS_H_

#include "../Model/SubstitutionModel.h"
#include "../Mapping/SubstitutionRegister.h"

#include <Bpp/Numeric/VectorTools.h>

namespace bpp
{

/**
 * @brief This class is used by MutationProcess to store detailed results of simulations.
 *
 * @author Julien Dutheil
 */
class MutationPath
{
	private:

    const Alphabet* alphabet_;

		/**
		 * @brief The states taken, without intiial state.
		 */
    std::vector<int> states_;

		/**
		 * @brief Times between states.
		 * The first element in array is the time between the initial state and the first state in states_.
		 */
    std::vector<double> times_;
		
		/**
		 * @brief The initial state.
		 */
		int initialState_;

		/**
		 * @brief Total time of evolution.
		 * Typically, this is a branch length.
		 */
		double totalTime_;

	public:

		/**
		 * @brief Builds a new MutationPath object with initial state 'initialState' and total time 'time'.
		 *
     * @param alphabet     The alphabet associated to the states in this path.
		 * @param initialState The initial state.
		 * @param time         The total time of evolution.
		 */
		MutationPath(const Alphabet* alphabet, int initialState, double time) :
      alphabet_(alphabet), states_(), times_(), initialState_(initialState), totalTime_(time) {};
		
    MutationPath(const MutationPath& path) :
      alphabet_(path.alphabet_), states_(path.states_), times_(path.times_), initialState_(path.initialState_), totalTime_(path.totalTime_) {};

    MutationPath& operator=(const MutationPath& path) {
      alphabet_     = path.alphabet_;
      states_       = path.states_;
      times_        = path.times_;
      initialState_ = path.initialState_;
      totalTime_    = path.totalTime_;
      return *this;
    }

		virtual ~MutationPath() {};

	public:
	
    /**
     * @return A pointer toward the alphabet associated to this path.
     */
    const Alphabet* getAlphabet() const { return alphabet_; }
		
    /**
		 * @brief Add a new mutation event.
		 *
		 * @param state The new state after mutation event.
		 * @param time  The time between this mutation and previous mutation (or initial state).
		 */
		void addEvent(int state, double time) {
			states_.push_back(state);
			times_.push_back(time);
		}

		/**
		 * @brief Retrieve the initial state.
		 *
		 * @return The initial state of this path.
		 */
		int getInitialState() const { return initialState_; }

		/**
		 * @brief Retrieve the total time of evolution.
		 *
		 * @return The total time of evolution.
		 */
		double getTotalTime() const { return totalTime_; }
		
		/**
		 * @brief Retrieve the number of substitution events.
		 *
		 * @return The number of substitution events, i.e. the number of states (without initial state).
		 */
		size_t getNumberOfEvents() const { return states_.size(); }

    /**
     * @brief Retrieve the number of substitution events per type of substitution.
     *
     * @param counts A matrix with the same size as the alphabet. The substitution counts will be incremented according to the mutation path, which allows to efficiently sum various mutation paths with a look.
     */
    template<class Scalar>
    void getEventCounts(Matrix<Scalar>& counts) const {
      if (counts.getNumberOfRows()    != alphabet_->getSize()
       || counts.getNumberOfColumns() != alphabet_->getSize())
        throw Exception("MutationPath::getEventCounts. Incorrect input matrix, does not match alphabet size.");
      int currentState = initialState_;
      for (size_t i = 0; i < states_.size(); ++i) {
        int newState = states_[i];
        counts(currentState, newState)++;
        currentState = newState;
      }
    }

    /**
     * @brief Retrieve the number of substitution events per type of substitution, defined by a SubstitutionRegister object.
     *
     * @param counts A vector with the appropriate size, as defined by SubstitutionRegister::getNumberOfSubstitutionTypes(). The substitution counts will be incremented according to the mutation path, which allows to efficiently sum various mutation paths with a look.
     * @param reg The substitution register to use to categorize substitutions.
     */
    template<class Scalar>
    void getEventCounts(std::vector<Scalar>& counts, const SubstitutionRegister& reg) const {
      if (counts.size() != reg.getNumberOfSubstitutionTypes())
        throw Exception("MutationPath::getEventCounts. Incorrect input vector, does not match alphabet size.");
      int currentState = initialState_;
      for (size_t i = 0; i < states_.size(); ++i) {
        int newState = states_[i];
        size_t type = reg.getType(currentState, newState);
        if (type > 0) counts[type - 1]++;
        currentState = newState;
      }
    }

		/**
		 * @brief Retrieve the final state of this path.
		 *
		 * @return The initial state if no mutation occured, otherwise sends the state after last mutation event.
		 */
		int getFinalState() const {
			if (states_.size() == 0) return initialState_;
			else return states_[states_.size() - 1];
		}
};

/**
 * @brief Interface for simulations.
 *
 * A mutation process defines the rules for mutations to occure.
 * The MutationProcess interface provides two methods, one for mutating a character in
 * state i in another character, another for achieving this task n times.
 */
class MutationProcess
{
	public:
		MutationProcess() {};
		virtual ~MutationProcess() {};
	
	public:
		
    /**
     * @brief Mutate a character in state i.
		 *
		 * @param state The current state of the character.
     */
    virtual int mutate(int state) const = 0;

    /**
     * @brief Mutate a character in state i n times.
     * 
		 * @param state The current state of the character.
		 * @param n The number of mutations to perform.
     */
    virtual int mutate(int state, unsigned int n) const = 0;
	
		/**
		 * @brief Get the time before next mutation event.
		 *
		 * @param state The actual state of the chain;
		 * @return A random time before next mutation event.
		 */
		virtual double getTimeBeforeNextMutationEvent(int state) const = 0;
		
		/**
		 * @brief Simulation a character evolution during a specified time
		 * according to the given substitution model and send the final state.
		 *
		 * @param initialState The state before beginning evolution.
		 * @param time         The time during which evolution must occure.
		 * @return The resulting state after evolution is completed.
		 */
		virtual int evolve(int initialState, double time) const = 0;
	
		/**
		 * @brief Simulation a character evolution during a specified time
		 * according to the given substitution model and send the total path
		 * with all intermediate states and times between mutation events.
		 *
		 * @param initialState The state before beginning evolution.
		 * @param time         The time during which evolution must occure.
		 * @return The resulting mutation path.
		 */
		virtual MutationPath detailedEvolve(int initialState, double time) const = 0;

		/**
		 * @brief Get the substitution model associated to the mutation process.
		 *
		 * @return The SubstitutionModel associated to this instance.
		 */
		virtual const SubstitutionModel* getSubstitutionModel() const = 0;
};

/**
 * @brief Partial implmentation of the MutationProcess interface.
 *
 * This class provides an implementation of the MutationProcess interface.
 * It assumes that there are size_ states allowed for the character of interest,
 * and that the distribution of probabilities are in repartition_.
 * As a matter of facts, probabilities must be cumulative, so that repartition_
 * contains values of the repartition function.
 * The mutate function hence draws a random number between 0 and 1 and gives the
 * corresponding character using the bijection of the repartition function.
 *
 * All derived classes must initialize the repartition_ and size_ fields.
 */
class AbstractMutationProcess :
  public virtual MutationProcess
{
	protected:
		
		/**
		 * @brief The substitution model to use:
		 */
		const SubstitutionModel* model_;
	
		/**
		 * @brief The number of states allowed for the character to mutate.
		 */
    size_t size_;
	
		/**
		 * @brief The repartition function for states probabilities.
		 *
		 * repartition_[i][j] = probability that, being in state i at time t,
		 * we'll be in state <= j at time t+1.
		 */
    VVdouble repartition_;
	
	public:
		AbstractMutationProcess(const SubstitutionModel* model) :
      model_(model), size_(), repartition_()
    {}

    AbstractMutationProcess(const AbstractMutationProcess& amp) :
      model_(amp.model_), size_(amp.size_), repartition_(amp.repartition_)
    {}

    AbstractMutationProcess& operator=(const AbstractMutationProcess& amp)
    {
      model_       = amp.model_;
      size_        = amp.size_;
      repartition_ = amp.repartition_;
      return *this;
    }

		virtual ~AbstractMutationProcess() {}
	
	public:
    int mutate(int state) const;
    int mutate(int state, unsigned int n) const;
		double getTimeBeforeNextMutationEvent(int state) const;
		int evolve(int initialState, double time) const;
		MutationPath detailedEvolve(int initialState, double time) const;
		const SubstitutionModel* getSubstitutionModel() const { return model_; }
};

/**
 * @brief Generally used mutation process model.
 *
 * This builds a MutationProcess according to a given SubstitutionModel.
 * The underlying mutation process is the following:
 * <ol>
 * <li>Draw a random time @f$ t @f$ from an exponential law with parameter
 * @f$ - \lambda_i @f$,</li>
 * <li> Mutate the initial state. The probability of mutating state @f$i@f$ 
 * to state @f$j@f$ is:
 * @f[ \frac{Q_{i,j}}{\sum_k Q_{i,k}}. @f]</li>
 * </ol>
 */
class SimpleMutationProcess : public AbstractMutationProcess
{
	public: // Constructor and destructor:
		
		/**
		 * @brief Build a new SimpleMutationProcess object.
		 *
		 * @param model The substitution model to use.
		 */
  	SimpleMutationProcess(const SubstitutionModel* model);
	
		virtual ~SimpleMutationProcess();

    /**
     * @brief Method redefinition for better performance.
     *
		 * @param initialState The state before beginning evolution.
		 * @param time         The time during which evolution must occure.
		 * @return The resulting state after evolution is completed.
		 */
		int evolve(int initialState, double time) const;
};

/**
 * @brief This class is mainly for testing purpose.
 * It allow "self" mutation of the kind i->i;
 */
class SelfMutationProcess : public AbstractMutationProcess
{
  	public:
  		SelfMutationProcess(int alphabetSize);
	
			virtual ~SelfMutationProcess();
};

} //end of namespace bpp.

#endif	//_MUTATIONPROCESS_H_

