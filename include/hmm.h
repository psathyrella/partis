#ifndef HAM_HMM_H
#define HAM_HMM_H

#include "state.h"
#include "yaml-cpp/yaml.h"

using namespace std;
namespace ham {

class model {
public:
  model();
  void parse(string);
  void finalize();

  inline string& name() { return name_; }
  inline size_t n_states() { return states_.size(); }
  // inline string& getStateName(size_t iter){ return states_[iter]->getName(); }
  // inline string& getStateLabel(size_t iter){ return states_[iter]->getLabel(); }
  State* state(const string&);
  inline State* state(size_t iter) { return states_[iter]; }
  inline bitset<STATE_MAX>* initial_to_states() { return &(initial_->to); }  //!Get vector of states that the initial state transitions to
  inline State* init_state() { return initial_; }  //!Get pointer to the initial state
  inline State* end_state() { return ending_; }  //!Get pointer to the ending state
  inline double overall_gene_prob() { return overall_gene_prob_; }

  void add_state(State*);
  inline void set_init(State* st) { initial_=st; }
  inline void set_end(State* st) { ending_=st; }
  bool checkTopology();
		
private:
  double overall_gene_prob_;  // prob to select this gene

  //!Flag set to tell whether the transitions bitsets have been set foreach
  //!state.  Model is also checked for correct order of states
  bool finalized;   
		
  string name_;	//! Model Name
		
  Tracks tracks_;

  vector<State*> states_; //!  All the states contained in the model

  map<string,State*> stateByName; //Ptr to state stored by State name;
		
  State* initial_; //!Initial state q0
  State* ending_;	//!Ending state
		
  void _addStateToFromTransition(State*);
  void _checkTopology(State* st, vector<uint16_t>& visited); //!Checks to see that all states are connected and there
};

}
#endif
